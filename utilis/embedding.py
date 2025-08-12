# utils/embedding.py
import os
import pickle
import re
from typing import Dict, List, Optional, Union
import numpy as np
import torch
from transformers import AutoTokenizer, AutoModel
from sklearn.metrics.pairwise import cosine_similarity
import config

# Global models and tokenizers to avoid reloading
sci_tok, sci_bert = None, None
bio_tok, bio_model = None, None
current_device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"[INFO] Using device: {current_device}")

def load_bert_models():
    """Loads SciBERT and BioBERT models into global scope."""
    global sci_tok, sci_bert, bio_tok, bio_model
    print("[INFO] Loading SciBERT...")
    sci_tok = AutoTokenizer.from_pretrained("allenai/scibert_scivocab_uncased")
    sci_bert = AutoModel.from_pretrained("allenai/scibert_scivocab_uncased").to(current_device).eval()

    print("[INFO] Loading BioBERT...")
    bio_tok = AutoTokenizer.from_pretrained("dmis-lab/biobert-base-cased-v1.1")
    bio_model = AutoModel.from_pretrained("dmis-lab/biobert-base-cased-v1.1").to(current_device).eval()

def _to_device(batch):
    return {k: v.to(current_device) for k, v in batch.items()}

def embed_with_bert(text: str) -> np.ndarray:
    if sci_tok is None or sci_bert is None:
        raise RuntimeError("SciBERT models are not loaded. Call load_bert_models() first.")
    toks = sci_tok(text, truncation=True, max_length=512, padding=False, return_tensors="pt")
    toks = _to_device(toks)
    with torch.no_grad():
        vec = sci_bert(**toks).last_hidden_state[:, 0, :]
    return vec.squeeze().cpu().numpy()

def get_embedding(text: str) -> np.ndarray:
    if bio_tok is None or bio_model is None:
        raise RuntimeError("BioBERT models are not loaded. Call load_bert_models() first.")
    tks = bio_tok.tokenize(text)
    max_len = bio_tok.model_max_length - 2
    chunks = [tks[i:i+max_len] for i in range(0, len(tks), max_len)]
    embs = []
    for ch in chunks:
        ids = bio_tok.build_inputs_with_special_tokens(bio_tok.convert_tokens_to_ids(ch))
        with torch.no_grad():
            out = bio_model(_to_device({"input_ids": torch.tensor([ids])}))[0][:, 0, :]
        embs.append(out.squeeze().cpu().numpy())
    return np.mean(embs, axis=0) if embs else np.zeros(bio_model.config.hidden_size)


class AttentionSeeker:
    """Loads a fine-tuned PubMedBERT and manages keyword + embedding caches."""
    def __init__(self):
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.tokenizer: Optional[AutoTokenizer] = None
        self.model: Optional[AutoModel] = None
        self.keyword_corpus = config.KEYWORD_CORPUS_PATH
        self.keyword_cache = config.KEYWORD_CACHE_PATH
        self.keywords: List[str] = []
        self.keyword_embeddings: np.ndarray = np.empty((0, 768))

        self._load_model(str(config.MODEL_DIR))
        if self.model:
            self._load_or_build_keyword_cache()
        else:
            print("[AttentionSeeker] Model failed to load, bio-aware filtering will be disabled.")

    def _load_model(self, path: str) -> None:
        try:
            self.tokenizer = AutoTokenizer.from_pretrained(path)
            self.model = AutoModel.from_pretrained(path).to(self.device).eval()
            print(f"[AttentionSeeker INFO] PubMedBERT loaded from {path}")
        except Exception as exc:
            print(f"[AttentionSeeker ERROR] Cannot load fine-tuned model at {path}\n{exc}")
            self.model = None

    def _load_plain_keywords(self) -> List[str]:
        kws: List[str] = []
        try:
            if not os.path.exists(self.keyword_corpus):
                print(f"[AttentionSeeker WARNING] Keyword file not found: {self.keyword_corpus}.")
                return kws
            with open(self.keyword_corpus, "r", encoding="utf-8") as fh:
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) >= 1 and parts[-1].strip():
                        kws.append(parts[-1].strip().lower())
        except Exception as exc:
            print(f"[AttentionSeeker ERROR] Failed to read keyword file {self.keyword_corpus}: {exc}")
        return kws

    @torch.no_grad()
    def _embed(self, texts: Union[str, List[str]]):
        if not self.model or not self.tokenizer:
            return None if isinstance(texts, str) else [None] * len(texts)
        if isinstance(texts, str):
            texts = [texts]
        if not texts:
            return []
        try:
            inputs = self.tokenizer(texts, return_tensors="pt", padding="longest", truncation=True, max_length=128).to(self.device)
            outputs = self.model(**inputs).last_hidden_state
            mask = inputs["attention_mask"].unsqueeze(-1).expand(outputs.size()).float()
            mean_embs = (outputs * mask).sum(1) / mask.sum(1).clamp(min=1e-9)
            result = [emb.cpu().numpy().flatten() for emb in mean_embs]
            return result[0] if len(result) == 1 and isinstance(texts, list) and len(texts) == 1 else result
        except Exception as exc:
            print(f"[AttentionSeeker ERROR] Batch embedding failed for texts: '{texts[0][:50]}...': {exc}")
            return [None] * len(texts)

    def _load_or_build_keyword_cache(self) -> None:
        if not os.path.exists(self.keyword_corpus):
            print(f"[AttentionSeeker WARNING] Keyword corpus file '{self.keyword_corpus}' not found.")
            return

        corpus_mtime = os.path.getmtime(self.keyword_corpus)
        cache_valid = False
        if os.path.exists(self.keyword_cache):
            try:
                with open(self.keyword_cache, "rb") as fh:
                    cache = pickle.load(fh)
                if (isinstance(cache, dict) and "mod_time" in cache and "keywords" in cache and "embeddings" in cache and cache["mod_time"] >= corpus_mtime):
                    self.keywords = cache["keywords"]
                    self.keyword_embeddings = cache["embeddings"]
                    cache_valid = True
                    print(f"[AttentionSeeker INFO] Keyword embeddings loaded from cache ({len(self.keywords)} keywords).")
            except Exception as exc:
                print(f"[AttentionSeeker WARNING] Failed to read keyword cache '{self.keyword_cache}': {exc}. Rebuilding.")

        if cache_valid:
            return

        print(f"[AttentionSeeker INFO] Building new keyword embeddings cache from '{self.keyword_corpus}'...")
        self.keywords = self._load_plain_keywords()
        if not self.keywords:
            print("[AttentionSeeker WARNING] No keywords loaded from corpus.")
            return
        
        embeds = self._embed(self.keywords)
        valid_embeds = [emb for emb in embeds if emb is not None and not np.all(emb == 0)]
        
        if valid_embeds:
            self.keyword_embeddings = np.vstack(valid_embeds)
            config.EMBEDDING_CACHE_DIR.mkdir(parents=True, exist_ok=True)
            with open(self.keyword_cache, "wb") as fh:
                pickle.dump({"keywords": self.keywords, "embeddings": self.keyword_embeddings, "mod_time": corpus_mtime}, fh)
            print(f"[AttentionSeeker INFO] Computed and cached {len(self.keywords)} keyword embeddings.")
        else:
            print("[AttentionSeeker WARNING] No valid embeddings generated.")
            self.keyword_embeddings = np.empty((0, 768))

    def calculate_similarity(self, tokens: List[str], threshold: float, token_cache: Dict[str, np.ndarray]) -> Dict[str, float]:
        # ... (code for calculate_similarity)
        pass
