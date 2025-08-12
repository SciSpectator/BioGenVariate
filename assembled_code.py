
import torch
import os
import sys
import gzip
import sqlite3
import tempfile
import uuid
import threading
import queue
import re
import math
import json
import webbrowser
from pathlib import Path
from typing import Optional, Dict, List, Set, Union
from datetime import datetime
import pandas as pd
import numpy as np
import tkinter as tk
from tkinter import messagebox, simpledialog, filedialog, ttk
from tkinter.constants import BOTH, END, LEFT, RIGHT, TOP, W, X, Y, CENTER, NORMAL, DISABLED, NSEW
from collections import Counter
from nltk.corpus import stopwords, wordnet as wn
import nltk

from transformers import AutoTokenizer, AutoModel
from sklearn.metrics.pairwise import cosine_similarity
from rake_nltk import Rake

# --- Global Setup ---
nltk.download("punkt", quiet=True)
nltk.download("stopwords", quiet=True)
nltk.download("wordnet", quiet=True)
stop_words = set(stopwords.words("english"))
current_device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"[INFO] Using device: {current_device}")

# --- BERT Model Loading ---
print("[INFO] Loading SciBERT...")
sci_tok = AutoTokenizer.from_pretrained("allenai/scibert_scivocab_uncased")
sci_bert = AutoModel.from_pretrained("allenai/scibert_scivocab_uncased").to(current_device).eval()

print("[INFO] Loading BioBERT...")
bio_tok = AutoTokenizer.from_pretrained("dmis-lab/biobert-base-cased-v1.1")
bio_model = AutoModel.from_pretrained("dmis-lab/biobert-base-cased-v1.1").to(current_device).eval()

# --- Helper Functions ---
def _to_device(batch):
    return {k: v.to(current_device) for k, v in batch.items()}

def embed_with_bert(text: str) -> np.ndarray:
    toks = sci_tok(text, truncation=True, max_length=512, padding=False, return_tensors="pt")
    toks = _to_device(toks)
    with torch.no_grad():
        vec = sci_bert(**toks).last_hidden_state[:, 0, :].squeeze().cpu().numpy()
    return vec

def get_embedding(text: str) -> np.ndarray:
    tks = bio_tok.tokenize(text)
    max_len = bio_tok.model_max_length - 2
    chunks = [tks[i:i + max_len] for i in range(0, len(tks), max_len)]
    embs = []
    for ch in chunks:
        ids = bio_tok.build_inputs_with_special_tokens(bio_tok.convert_tokens_to_ids(ch))
        with torch.no_grad():
            out = bio_model(_to_device({"input_ids": torch.tensor([ids])}))[0][:, 0, :]
        embs.append(out.squeeze().cpu().numpy())
    return np.mean(embs, axis=0) if embs else np.zeros(bio_model.config.hidden_size)

def build_final_text(row):
    return " ".join(str(row[col]) for col in ("title", "description", "characteristics") if col in row and pd.notnull(row[col]))

# ───────────────────────────────────────────────────────────────────────────────
# ## FIX: ManualLabelingDialog DEFINED BEFORE THREADS THAT USE IT ##
# ───────────────────────────────────────────────────────────────────────────────
class ManualLabelingDialog(simpledialog.Dialog):
    def __init__(self, parent, title, text_content):
        self.text_content = text_content
        self.user_choice = -1
        super().__init__(parent, title)

    def body(self, master):
        self.resizable(True, True)
        text_widget = tk.Text(master, wrap=tk.WORD, height=20, width=80)
        text_widget.pack(padx=5, pady=5, fill=tk.BOTH, expand=True)
        text_widget.insert(tk.END, self.text_content)
        text_widget.config(state=tk.DISABLED)
        return text_widget

    def buttonbox(self):
        box = tk.Frame(self)
        tk.Button(box, text="Case (1)", width=10, command=lambda: self.ok(1), bg="#2ca02c", fg="white").pack(side=tk.LEFT, padx=5, pady=5)
        tk.Button(box, text="Control (0)", width=10, command=lambda: self.ok(0), bg="#d62728", fg="white").pack(side=tk.LEFT, padx=5, pady=5)
        tk.Button(box, text="Skip (-1)", width=10, command=lambda: self.ok(-1)).pack(side=tk.LEFT, padx=5, pady=5)
        self.bind("<Escape>", lambda event: self.ok(-1))
        box.pack()

    def ok(self, choice=None):
        self.user_choice = choice
        super().ok()

# ───────────────────────────────────────────────────────────────────────────────
# Step 1: ExtractionThread Class
# ───────────────────────────────────────────────────────────────────────────────
class ExtractionThread(threading.Thread):
    SIM_THRESHOLD = 0.65
    EXCLUDED_COLUMNS = {"gsm", "contact", "supplementary_file", "data_row_count", "channel_count", "organism_ch1", "status", "series_id", "submission_date", "last_update_date", "data_processing"}

    def __init__(self, plat_filter: str, tokens: str, log_cb, on_finish=None):
        super().__init__()
        self.plat_filter, self.tokens, self.log_cb, self.on_finish = plat_filter, tokens, log_cb, on_finish
        self.final_df, self.gse_keywords = None, {}
        self._stop_event = threading.Event()

    def stop(self): self._stop_event.set()
    def log(self, msg): self.log_cb(msg)

    def run(self):
        mem_conn = None
        tmp_sql_path = None
        try:
            self.log("PROGRESS: 0")
            gz_path = "./GEOmetadb.sqlite.gz"
            if not os.path.exists(gz_path): return self.log("[STEP 1 ERROR] GEOmetadb.sqlite.gz not found.")

            self.log("[STEP 1] Loading GEOmetadb into thread-local memory...")
            with tempfile.NamedTemporaryFile(suffix=".sqlite", delete=False) as tmp:
                tmp_sql_path = tmp.name
                with gzip.open(gz_path, "rb") as gzfi: tmp.write(gzfi.read())

            disk_conn = sqlite3.connect(tmp_sql_path)
            # ## FIX: AVOID UTF-8 DECODING ERRORS ##
            disk_conn.text_factory = lambda b: b.decode('utf-8', 'replace')
            
            mem_conn = sqlite3.connect(":memory:")
            mem_conn.text_factory = disk_conn.text_factory
            disk_conn.backup(mem_conn)
            disk_conn.close()
            os.remove(tmp_sql_path)
            tmp_sql_path = None
            
            self.log("PROGRESS: 10")
            df = pd.read_sql_query("SELECT * FROM gsm", mem_conn)
            self.log(f"[STEP 1] Loaded {len(df):,} GSM rows.")
            if self._stop_event.is_set(): return

            if self.plat_filter.strip():
                wanted = {p.strip().lower() for p in self.plat_filter.split(",")}
                df = df[df["gpl"].astype(str).str.lower().isin(wanted)].copy()
            self.log(f"[STEP 1] After platform filter: {len(df):,} rows.")
            self.log("PROGRESS: 30")

            raw_toks = [t.strip().lower() for t in self.tokens.split(",") if t.strip()]
            if not raw_toks: return self.log("[STEP 1 ERROR] No filter tokens provided.")

            lex_tokens = set(raw_toks)
            try:
                for t in raw_toks:
                    for syn in wn.synsets(t):
                        lex_tokens.update(lemma.name().lower().replace("_", " ") for lemma in syn.lemmas())
            except Exception as e: self.log(f"[WARN] WordNet expansion failed: {e}")

            searchable_cols = [c for c in df.columns if c.lower() not in self.EXCLUDED_COLUMNS]
            def fast_hit(row):
                blob = " ".join(str(row[c]) for c in searchable_cols if pd.notnull(row[c])).lower()
                return any(tok in blob for tok in lex_tokens)
            
            cand = df[df.apply(fast_hit, axis=1)].copy()
            if cand.empty: return self.log("[STEP 1 ERROR] Nothing passed lexical screen.")
            self.log(f"[STEP 1] After lexical screen: {len(cand):,} rows.")
            self.log("PROGRESS: 40")

            query_emb = embed_with_bert(" ".join(raw_toks))
            def sem_hit(row):
                txt = " ".join(str(row[c]) for c in searchable_cols if pd.notnull(row[c]))
                if not txt.strip(): return False
                try:
                    sim = cosine_similarity(embed_with_bert(txt).reshape(1, -1), query_emb.reshape(1, -1))[0, 0]
                    return sim >= self.SIM_THRESHOLD
                except Exception: return False

            matched = cand[cand.apply(sem_hit, axis=1)].copy()
            if matched.empty: return self.log("[STEP 1 ERROR] No rows passed semantic check.")
            self.log(f"[STEP 1] After semantic check: {len(matched):,} rows.")
            self.log("PROGRESS: 50")

            docs = matched.groupby("series_id").apply(lambda g: "\n".join(g.apply(build_final_text, axis=1))).to_dict()
            try:
                rake = Rake()
                for i, (gse, doc) in enumerate(docs.items()):
                    if self._stop_event.is_set(): return
                    rake.extract_keywords_from_text(doc)
                    self.gse_keywords[gse] = list(set(k for k in rake.get_ranked_phrases()[:20] if not re.search(r"\d", k)))
            except NameError:
                self.log("[WARN] rake_nltk not found. Skipping keyword extraction.")
                self.gse_keywords = {gse: [] for gse in docs.keys()}
            self.log("PROGRESS: 70")
            
            gse_ids_sql = "','".join(g.replace("'", "''") for g in docs.keys())
            self.final_df = pd.read_sql_query(f"SELECT * FROM gsm WHERE series_id IN ('{gse_ids_sql}')", mem_conn)
            self.final_df["Token_Match"] = self.final_df["gsm"].isin(matched["gsm"])
            self.final_df["Final_Text"] = self.final_df.apply(build_final_text, axis=1)
            
            self.log(f"[STEP 1] Final GSM set: {len(self.final_df):,} rows.")
            self.log("PROGRESS: 80")

            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            tok_tag = re.sub(r'[^\w\-]+', '_', self.tokens or "tokens")[:50]
            plat_tag = re.sub(r'[^\w\-]+', '_', self.plat_filter or "ALL")[:50]
            out_dir = Path("NEW_RESULTS_ROOT") / f"{tok_tag}_{plat_tag}_{timestamp}"
            out_dir.mkdir(parents=True, exist_ok=True)
            self.final_df.to_csv(out_dir / "step1_extracted_gsm.csv", index=False)
            self.log(f"[STEP 1] CSV saved -> {out_dir}")
            
            self.log("PROGRESS: 100")
            self.log("[STEP 1 DONE] Extraction complete.")
        except Exception as e:
            self.log(f"[STEP 1 THREAD EXCEPTION] {type(e).__name__}: {e}")
            import traceback; self.log(traceback.format_exc())
        finally:
            if mem_conn: mem_conn.close()
            if tmp_sql_path: os.remove(tmp_sql_path)
            if self.on_finish: app.after(0, self.on_finish)

# ───────────────────────────────────────────────────────────────────────────────
# Step 2: LabelingThread Class
# ───────────────────────────────────────────────────────────────────────────────
class LabelingThread(threading.Thread):
    def __init__(self, df_sub: pd.DataFrame, folder_name: str, log_cb, mode: str, case_toks: str = "", ctrl_toks: str = ""):
        super().__init__()
        self.df_sub, self.folder_name, self.log_cb, self.mode = df_sub.copy(), folder_name, log_cb, mode
        self.case_toks, self.ctrl_toks, self.user_labels = case_toks, ctrl_toks, []
        self._stop_event = threading.Event()

    def stop(self): self._stop_event.set()
    def log(self, msg): self.log_cb(msg)

    def run(self):
        try:
            if self.df_sub.empty: return self.log("[STEP 2 ERROR] No GSM rows to label.")
            self.log("PROGRESS: 0")
            if self.mode == 'automatic': self.run_automatic_mode()
            elif self.mode == 'manual': self.run_manual_mode()
        except Exception as e:
            self.log(f"[STEP 2 EXCEPTION] {type(e).__name__}: {e}")
            import traceback; self.log(traceback.format_exc())

    def run_automatic_mode(self):
        self.log("[STEP 2] Starting Automatic Labeling...")
        case_emb = get_embedding(" ".join(t.strip() for t in self.case_toks.split(",")))
        ctrl_emb = get_embedding(" ".join(t.strip() for t in self.ctrl_toks.split(",")))
        self.log("PROGRESS: 20")
        def auto_classify(text):
            if not isinstance(text, str) or not text.strip(): return 0.5
            emb = get_embedding(text)
            sc = cosine_similarity(emb.reshape(1, -1), case_emb.reshape(1, -1))[0, 0]
            sb = cosine_similarity(emb.reshape(1, -1), ctrl_emb.reshape(1, -1))[0, 0]
            return 1 if sc > sb else 0
        self.df_sub["Predicted_Label"] = self.df_sub["Final_Text"].apply(auto_classify)
        self.log("PROGRESS: 60")
        self.save_results()

    def run_manual_mode(self):
        self.log("[STEP 2] Starting Manual Labeling...")
        total_samples = len(self.df_sub)
        result_queue = queue.Queue()
        for index, row in self.df_sub.iterrows():
            if self._stop_event.is_set(): return self.log("[INFO] Manual labeling stopped.")
            self.log(f"PROGRESS: {(index + 1) / total_samples * 80:.1f}")
            def ask_label(r, i):
                dialog = ManualLabelingDialog(app, f"Label Sample {i + 1}/{total_samples} (GSM: {r['gsm']})", r['Final_Text'])
                result_queue.put(dialog.user_choice)
            app.after(0, lambda r=row, i=index: ask_label(r, i))
            label = result_queue.get()
            if label is None: return self.log("[STEP 2] Labeling aborted.")
            self.user_labels.append(label)
        self.df_sub["User_Label"] = self.user_labels
        self.log("PROGRESS: 90")
        self.save_results()

    def save_results(self):
        mode_tag = "manual_labeled" if self.mode == "manual" else "auto_labeled"
        tok_tag = re.sub(r'[^\w\-]+', '_', app.filter_entry.get() or "tokens")[:50]
        plat_tag = re.sub(r'[^\w\-]+', '_', app.platform_entry.get() or "ALL")[:50]
        step1_base_dir = Path("NEW_RESULTS_ROOT")
        matching_dirs = sorted([d for d in step1_base_dir.iterdir() if d.is_dir() and d.name.startswith(f"{tok_tag}_{plat_tag}")], reverse=True)
        out_dir = matching_dirs[0] if matching_dirs else step1_base_dir / f"{tok_tag}_{plat_tag}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        out_dir.mkdir(exist_ok=True)
        out_path = out_dir / f"step2_{self.folder_name}_{mode_tag}.csv"
        try:
            self.df_sub.to_csv(out_path, index=False)
            self.log(f"[STEP 2] Saved -> {out_path}")
        except Exception as e: self.log(f"[STEP 2 ERROR] Could not write file:\n{e}")
        self.log("PROGRESS: 100")
        self.log(f"[STEP 2 DONE] {self.mode.capitalize()} labeling complete.")

# ───────────────────────────────────────────────────────────────────────────────
# Ontology related helpers
CACHE_DIR = Path("./.ontology_cache")
CACHE_DIR.mkdir(exist_ok=True)

ONTOLOGIES = {
    "uberon.obo": "https://purl.obolibrary.org/obo/uberon.obo",
    "doid.obo":    "https://purl.obolibrary.org/obo/doid.obo",
}
SPECIES_TERMS = {"homo sapiens", "mus musculus", "rattus norvegicus", "homo", "sapiens", "human", "mouse", "rat"}

def cached_path(fname: str, url: str) -> Path:
    path = CACHE_DIR / fname
    if not path.exists():
        print(f"[ONTOLOGY] Downloading {fname} ...")
        tmp = path.with_suffix(".tmp")
        try:
            urllib.request.urlretrieve(url, tmp)
            tmp.rename(path)
            print(f"[ONTOLOGY] Saved -> {path.name} ({path.stat().st_size / 1e6:.1f} MB)")
        except Exception as e:
            print(f"[ONTOLOGY ERROR] Failed to download {url} to {path}: {e}")
            if tmp.exists():
                tmp.unlink() # Clean up temp file
            raise
    return path

def smart_open(path: Path) -> IO[str]:
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "rt", encoding="utf-8", errors="ignore")

def parse_obo(path: Path) -> List[Tuple[str, List[str]]]:
    terms: List[Tuple[str, List[str]]] = []
    try:
        with smart_open(path) as fh:
            in_term = False
            label, syns, obsolete = None, [], False
            for line in fh:
                if line.startswith("[Term]"):
                    if in_term and label and not obsolete:
                        terms.append((label, syns))
                    label, syns, obsolete = None, [], False
                    in_term = True
                    continue
                if not in_term:
                    continue
                if line.startswith("name:"):
                    label = line.partition("name:")[2].strip()
                elif line.startswith("synonym:"):
                    syns.append(line.split('"', 2)[1])
                elif line.startswith("is_obsolete: true"):
                    obsolete = True
            if in_term and label and not obsolete:
                terms.append((label, syns))
    except Exception as e:
        print(f"[ONTOLOGY ERROR] Failed to parse OBO file {path}: {e}")
        return [] # Return empty on error
    return terms

def load_ontology_sets() -> Tuple[Set[str], Set[str], Dict, Dict]:
    """
    Loads Uberon and DOID ontologies. 
    Returns sets of canonical terms and pre-processed dictionaries for fast lookup.
    """
    print("[ONTOLOGY] Loading and pre-processing ontologies for categorization...")
    uberon_map, disease_map = defaultdict(list), defaultdict(list)
    uberon_canons, disease_canons = set(), set()
    
    try:
        uberon_terms = parse_obo(cached_path("uberon.obo", ONTOLOGIES["uberon.obo"]))
        for canon, synonyms in uberon_terms:
            uberon_canons.add(canon.lower())
            all_terms = [canon] + synonyms
            for term in all_terms:
                for word in re.split(r'\W+', term.lower()):
                    if word and len(word) > 2:
                        uberon_map[word].append(canon)

        disease_terms = parse_obo(cached_path("doid.obo", ONTOLOGIES["doid.obo"]))
        for canon, synonyms in disease_terms:
            disease_canons.add(canon.lower())
            all_terms = [canon] + synonyms
            for term in all_terms:
                for word in re.split(r'\W+', term.lower()):
                    if word and len(word) > 2:
                        disease_map[word].append(canon)

        # Make mappings unique
        for word in uberon_map: uberon_map[word] = sorted(list(set(uberon_map[word])), key=len)
        for word in disease_map: disease_map[word] = sorted(list(set(disease_map[word])), key=len)
        
        print(f"[ONTOLOGY] Loaded {len(uberon_canons)} Uberon terms and {len(disease_canons)} Disease terms.")
    except Exception as e:
        print(f"[ONTOLOGY ERROR] Could not load all ontology files: {e}. Categorization may be incomplete.")
    
    return uberon_canons, disease_canons, uberon_map, disease_map

def find_canonical_term(token: str, disease_map: Dict, uberon_map: Dict) -> str:
    """Finds the best matching canonical term for a given token."""
    token_words = {word for word in re.split(r'\W+', token.lower()) if len(word) > 2}
    
    # Prioritize disease matches
    for word in token_words:
        if word in disease_map:
            return disease_map[word][0] # Return the shortest canonical term
            
    # Then check tissue matches
    for word in token_words:
        if word in uberon_map:
            return uberon_map[word][0] # Return the shortest canonical term
            
    return token # Fallback to the token itself if no match is found

def get_canon_color(canon_name: str, disease_canons: Set[str], uberon_canons: Set[str]) -> str:
    """Determines the color for a canonical token based on its ontology."""
    canon_lower = canon_name.lower()
    if canon_lower in SPECIES_TERMS:
        return 'grey' # Changed to grey for species
    if canon_lower in disease_canons:
        return 'red'
    elif canon_lower in uberon_canons:
        return 'green'
    else:
        return 'black'

# Repaired `CompareDistributionsWindow`
class CompareDistributionsWindow(tk.Toplevel):
    def __init__(self, parent, app_ref):
        super().__init__(parent)
        self.parent = parent
        self.app_ref = app_ref
        self.title("Compare Distributions")
        self.geometry("1200x850")
        self.user_defined_groups = {}

        self.fig_density, self.ax_density, self.canvas_density, self.toolbar_density = None, None, None, None
        self.fig_hist, self.ax_hist, self.canvas_hist, self.toolbar_hist = None, None, None, None

        try:
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
            self.FigureCanvasTkAgg = FigureCanvasTkAgg
            self.NavigationToolbar2Tk = NavigationToolbar2Tk
            self.modules_loaded = True
        except ImportError as e:
            self.app_ref.enqueue_log(f"[ERROR] Matplotlib TkAgg backend not available. {e}")
            self.modules_loaded = False

        self._setup_ui()
        self.protocol("WM_DELETE_WINDOW", self._on_closing)

    def _setup_ui(self):
        main_frame = ttk.Frame(self, padding="10")
        main_frame.pack(fill=BOTH, expand=True)
        main_frame.grid_columnconfigure(0, weight=2, minsize=400) 
        main_frame.grid_columnconfigure(1, weight=5) 
        main_frame.grid_rowconfigure(0, weight=1)

        control_pane = ttk.Frame(main_frame, padding=5)
        control_pane.grid(row=0, column=0, sticky=NSEW)
        control_pane.grid_rowconfigure(5, weight=1)
        control_pane.grid_columnconfigure(0, weight=1)

        data_loading_frame = ttk.LabelFrame(control_pane, text="1. Load Case-Control File(s)")
        data_loading_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        data_loading_frame.grid_columnconfigure(0, weight=1)
        
        load_button_frame = ttk.Frame(data_loading_frame)
        load_button_frame.pack(fill=X, padx=5, pady=5)
        ttk.Button(load_button_frame, text="Load File(s)", command=self._load_case_control_file).pack(side=LEFT, fill=X, expand=True)
        ttk.Button(load_button_frame, text="Clear All", command=self._clear_user_data).pack(side=LEFT, padx=5)

        ttk.Label(control_pane, text="2. Select Groups to Compare", font=("Helvetica", 10, "bold")).grid(row=1, column=0, sticky=W, padx=5, pady=(10, 2))
        self.loaded_files_listbox = tk.Listbox(control_pane, height=6, selectmode=tk.EXTENDED)
        self.loaded_files_listbox.grid(row=2, column=0, sticky="ew", padx=5, pady=2)

        mode_frame = ttk.LabelFrame(control_pane, text="3. Choose Comparison Mode")
        mode_frame.grid(row=3, column=0, sticky="ew", padx=5, pady=10)
        self.comparison_mode = tk.StringVar(value="groups_only")
        ttk.Radiobutton(mode_frame, text="Compare Selected Case/Control Groups Only", variable=self.comparison_mode, value="groups_only").pack(anchor=W, padx=10)
        ttk.Radiobutton(mode_frame, text="Compare Groups vs. Gene(s) on Platform(s)", variable=self.comparison_mode, value="vs_gene").pack(anchor=W, padx=10)
        ttk.Radiobutton(mode_frame, text="Compare Groups vs. Entire Platform(s) (All Genes)", variable=self.comparison_mode, value="vs_platform").pack(anchor=W, padx=10)

        spec_frame = ttk.LabelFrame(control_pane, text="4. Specify Gene(s) and/or Platform(s)")
        spec_frame.grid(row=4, column=0, sticky="ew", padx=5, pady=5)

        ttk.Label(spec_frame, text="Gene Symbol(s) (Optional for Loaded Groups):").pack(anchor=W, padx=5, pady=(5,0))
        self.gene_entry = ttk.Entry(spec_frame)
        self.gene_entry.pack(fill=X, padx=5, pady=2)
        
        ttk.Label(spec_frame, text="Comparison Platform(s) (for modes 2 & 3):").pack(anchor=W, padx=5, pady=(5,0))
        self.platform_frame = ttk.Frame(spec_frame)
        self.platform_frame.pack(fill=X, padx=5, pady=2)
        self.platform_vars = {}
        if self.app_ref.gpl_datasets:
            for p in sorted(self.app_ref.gpl_datasets.keys()):
                var = tk.BooleanVar(value=False)
                ttk.Checkbutton(self.platform_frame, text=p, variable=var).pack(side=LEFT)
                self.platform_vars[p] = var

        ttk.Button(control_pane, text="Plot & Analyze Distributions", command=self._plot_and_analyze).grid(row=5, column=0, sticky="ew", padx=5, pady=15, ipady=5)
        self.status_label = ttk.Label(control_pane, text="Ready.", foreground="blue")
        self.status_label.grid(row=6, column=0, sticky=W, padx=5, pady=5)

        output_pane = ttk.Frame(main_frame)
        output_pane.grid(row=0, column=1, sticky=NSEW, padx=5)
        output_pane.grid_rowconfigure(0, weight=3); output_pane.grid_rowconfigure(1, weight=2)
        output_pane.grid_columnconfigure(0, weight=1); output_pane.grid_columnconfigure(1, weight=1)

        self.density_plot_frame = ttk.LabelFrame(output_pane, text="Overlaid Density Plots")
        self.density_plot_frame.grid(row=0, column=0, sticky=NSEW, padx=2, pady=5)
        self.hist_plot_frame = ttk.LabelFrame(output_pane, text="Overlaid Histograms & Classification")
        self.hist_plot_frame.grid(row=0, column=1, sticky=NSEW, padx=2, pady=5)
        self.stats_frame = ttk.LabelFrame(output_pane, text="Statistical Test Results (Pairwise Wilcoxon rank-sum)")
        self.stats_frame.grid(row=1, column=0, columnspan=2, sticky=NSEW, padx=2, pady=5)

    def _clear_user_data(self):
        self.user_defined_groups.clear()
        self.loaded_files_listbox.delete(0, END)
        self.status_label.config(text="Cleared all loaded files.")

    def _on_closing(self):
        self._clear_figure()
        self.destroy()

    def _clear_figure(self):
        for fig, canvas, toolbar in [(self.fig_density, self.canvas_density, self.toolbar_density), (self.fig_hist, self.canvas_hist, self.toolbar_hist)]:
            if fig: plt.close(fig)
            if canvas: canvas.get_tk_widget().destroy()
            if toolbar: toolbar.destroy()
        self.fig_density, self.canvas_density, self.toolbar_density = None, None, None
        self.fig_hist, self.canvas_hist, self.toolbar_hist = None, None, None
        
    def _clear_stats_table(self):
        for widget in self.stats_frame.winfo_children(): widget.destroy()

    def _ask_for_platform(self, platforms: List[str]) -> Optional[str]:
        dialog = tk.Toplevel(self)
        dialog.title("Select Platform")
        dialog.transient(self); dialog.grab_set()
        ttk.Label(dialog, text="Select the platform for the loaded samples:", wraplength=300).pack(pady=10)
        var = tk.StringVar(value=platforms[0])
        for p in platforms:
            ttk.Radiobutton(dialog, text=p, variable=var, value=p).pack(anchor=W, padx=30)
        result = {"platform": None}
        def on_ok():
            result["platform"] = var.get()
            dialog.destroy()
        ttk.Button(dialog, text="OK", command=on_ok).pack(pady=10, ipadx=20)
        self.wait_window(dialog)
        return result["platform"]

    def _load_case_control_file(self):
        filepaths = filedialog.askopenfilenames(
            title="Select Case-Control CSV or gzipped CSV file(s)",
            filetypes=[("CSV files", "*.csv"), ("Gzipped CSV files", "*.csv.gz")]
        )
        if not filepaths: return

        for filepath in filepaths:
            try:
                df_annot = pd.read_csv(filepath, compression='gzip' if filepath.endswith('.gz') else None, low_memory=False)
                if df_annot.empty: continue

                columns = df_annot.columns.tolist()
                formatted_columns = '\n'.join(columns)
                gsm_col = simpledialog.askstring("GSM Column", f"For file '{Path(filepath).name}', enter column name for GSM IDs.\n\nAvailable:\n{formatted_columns}", parent=self)
                if not gsm_col or gsm_col not in columns: continue

                cc_col = simpledialog.askstring("Case/Control Column", f"For file '{Path(filepath).name}', enter column for labels (0=Control, 1=Case).\n\nAvailable:\n{formatted_columns}", parent=self)
                if not cc_col or cc_col not in columns: continue

                loaded_platforms = list(self.app_ref.gpl_datasets.keys())
                if not loaded_platforms:
                    messagebox.showerror("Error", "No platforms are loaded in the main app.", parent=self)
                    return
                
                platform = self._ask_for_platform(loaded_platforms)
                if not platform: continue

                df_annot.dropna(subset=[gsm_col, cc_col], inplace=True)
                df_annot[gsm_col] = df_annot[gsm_col].astype(str).str.upper()
                df_annot[cc_col] = pd.to_numeric(df_annot[cc_col], errors='coerce')
                
                if not set(df_annot[cc_col].unique()).issubset({0.0, 1.0}):
                    messagebox.showerror("Data Error", f"File '{Path(filepath).name}': The case/control column '{cc_col}' must contain only 0s and 1s.", parent=self)
                    continue

                case_gsms = df_annot[df_annot[cc_col] == 1.0][gsm_col].tolist()
                control_gsms = df_annot[df_annot[cc_col] == 0.0][gsm_col].tolist()
                
                filename = Path(filepath).name
                if case_gsms:
                    label = f"{filename} (Case)"
                    self.user_defined_groups[label] = {"platform": platform, "gsms": case_gsms}
                    self.loaded_files_listbox.insert(END, label)
                if control_gsms:
                    label = f"{filename} (Control)"
                    self.user_defined_groups[label] = {"platform": platform, "gsms": control_gsms}
                    self.loaded_files_listbox.insert(END, label)
            except Exception as e:
                messagebox.showerror("Error", f"Failed to process file {Path(filepath).name}:\n{e}", parent=self)

    def _plot_and_analyze(self):
        if not self.modules_loaded: return
        self._clear_figure()
        self._clear_stats_table()
        
        data_to_plot = {}
        
        selected_indices = self.loaded_files_listbox.curselection()
        selected_groups = [self.loaded_files_listbox.get(i) for i in selected_indices]
        mode = self.comparison_mode.get()
        gene_input_str = self.gene_entry.get().strip()
        genes_to_plot = [g.strip().upper() for g in gene_input_str.split(',') if g.strip()]
        selected_platforms = [p for p, var in self.platform_vars.items() if var.get()]

        # --- Main Logic Branching ---

        # --- Use Case 1: User has selected loaded Case/Control groups ---
        if selected_groups:
            if mode == "vs_gene" and not genes_to_plot:
                messagebox.showerror("Input Required", "Please enter gene(s) to compare your groups against platforms.", parent=self)
                return
            if mode != "groups_only" and not selected_platforms:
                messagebox.showerror("Input Required", f"Please select at least one Comparison Platform for the chosen mode.", parent=self)
                return

            # 1a. Add data for selected user-defined groups
            for label in selected_groups:
                group_info = self.user_defined_groups.get(label)
                if not group_info: continue
                platform, gsms = group_info["platform"], group_info["gsms"]
                df_plat = self.app_ref.gpl_datasets.get(platform)
                gmap_plat = self.app_ref.gpl_gene_mappings.get(platform, {})
                if df_plat is None or not gmap_plat: continue
                
                if genes_to_plot:
                    for gene in genes_to_plot:
                        gene_col = gmap_plat.get(gene)
                        if gene_col and gene_col in df_plat.columns:
                            subset_df = df_plat[df_plat['GSM'].isin(gsms)]
                            expr = pd.to_numeric(subset_df[gene_col], errors='coerce').dropna()
                            if not expr.empty: data_to_plot[f"{gene} - {label}"] = expr
                else:
                    all_gene_cols = [col for col in gmap_plat.values() if col in df_plat.columns]
                    subset_df = df_plat[df_plat['GSM'].isin(gsms)]
                    all_expr = pd.concat([pd.to_numeric(subset_df[col], errors='coerce').dropna() for col in all_gene_cols], ignore_index=True)
                    if not all_expr.empty: data_to_plot[label] = all_expr
            
            # 1b. Add data for comparison platforms based on mode
            if mode == "vs_gene":
                for platform in selected_platforms:
                    df, gmap = self.app_ref.gpl_datasets.get(platform), self.app_ref.gpl_gene_mappings.get(platform, {})
                    if df is None or not gmap: continue
                    for gene in genes_to_plot:
                        gene_col = gmap.get(gene)
                        if gene_col in df.columns:
                            expr = pd.to_numeric(df[gene_col], errors='coerce').dropna()
                            if not expr.empty: data_to_plot[f"{gene} ({platform})"] = expr
            elif mode == "vs_platform":
                for platform in selected_platforms:
                    df, gmap = self.app_ref.gpl_datasets.get(platform), self.app_ref.gpl_gene_mappings.get(platform, {})
                    if df is None or not gmap: continue
                    all_expr = pd.concat([pd.to_numeric(df[col], errors='coerce').dropna() for col in gmap.values() if col in df.columns], ignore_index=True)
                    if not all_expr.empty: data_to_plot[f"All Genes ({platform})"] = all_expr

        # --- Use Case 2: No groups selected, direct platform/gene comparison ---
        else:
            if not genes_to_plot or not selected_platforms:
                self.status_label.config(text="No valid data to plot. Select genes and platforms.")
                return

            for platform in selected_platforms:
                df, gmap = self.app_ref.gpl_datasets.get(platform), self.app_ref.gpl_gene_mappings.get(platform, {})
                if df is None or not gmap: continue
                for gene in genes_to_plot:
                    gene_col = gmap.get(gene)
                    if gene_col in df.columns:
                        expr = pd.to_numeric(df[gene_col], errors='coerce').dropna()
                        if not expr.empty: 
                            data_to_plot[f"{gene} ({platform})"] = expr
        
        # --- Final validation and plotting ---
        if not data_to_plot:
            self.status_label.config(text="No valid data to plot based on selections.")
            return

        self.status_label.config(text=f"Plotting {len(data_to_plot)} distributions...")
        
        self.fig_density, self.ax_density = plt.subplots(figsize=(7, 5))
        for label, data in data_to_plot.items():
            sns.kdeplot(data, ax=self.ax_density, label=label, fill=True, alpha=0.4)
        self.ax_density.set_title("Overlaid Density", fontsize=12)
        self.ax_density.legend(fontsize='x-small'); self.fig_density.tight_layout()
        color_map = {text.get_text(): line.get_color() for text, line in zip(self.ax_density.get_legend().get_texts(), self.ax_density.get_lines())}
        
        self.canvas_density = self.FigureCanvasTkAgg(self.fig_density, master=self.density_plot_frame)
        self.canvas_density.draw(); self.canvas_density.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
        self.toolbar_density = self.NavigationToolbar2Tk(self.canvas_density, self.density_plot_frame); self.toolbar_density.update()
        
        self.fig_hist, self.ax_hist = plt.subplots(figsize=(7, 5))
        if data_to_plot:
            all_values = pd.concat(data_to_plot.values())
            bins = np.linspace(all_values.min(), all_values.max(), 50)
            classifications = {label: analyze_gene_distribution(series) for label, series in data_to_plot.items()}
            for label, data in data_to_plot.items():
                self.ax_hist.hist(data, bins=bins, color=color_map.get(label), label=f"{label} ({classifications.get(label, 'N/A')})", alpha=0.6)
        
        self.ax_hist.set_title("Overlaid Histograms", fontsize=12); self.ax_hist.set_xlabel("Expression"); self.ax_hist.set_ylabel("Frequency")
        self.ax_hist.legend(fontsize='x-small'); self.fig_hist.tight_layout()
        
        self.canvas_hist = self.FigureCanvasTkAgg(self.fig_hist, master=self.hist_plot_frame)
        self.canvas_hist.draw(); self.canvas_hist.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
        self.toolbar_hist = self.NavigationToolbar2Tk(self.canvas_hist, self.hist_plot_frame); self.toolbar_hist.update()

        self._run_stats_analysis(data_to_plot)
        self.status_label.config(text="Analysis complete.")

    def _run_stats_analysis(self, data_for_stats):
        if len(data_for_stats) < 2: return
        self._clear_stats_table()
        
        columns = ["Distribution 1", "Distribution 2", "Wilcoxon Score (U)"]
        stats_tree = ttk.Treeview(self.stats_frame, columns=columns, show="headings")
        stats_tree.pack(side=LEFT, fill=BOTH, expand=True, padx=5, pady=5)
        
        vsb = ttk.Scrollbar(self.stats_frame, orient="vertical", command=stats_tree.yview)
        vsb.pack(side=RIGHT, fill=Y)
        stats_tree.configure(yscrollcommand=vsb.set)

        stats_tree.column("Distribution 1", width=300, anchor=W)
        stats_tree.column("Distribution 2", width=300, anchor=W)
        stats_tree.column("Wilcoxon Score (U)", width=150, anchor=CENTER)
        for col in columns:
            stats_tree.heading(col, text=col)

        for combo in itertools.combinations(data_for_stats.keys(), 2):
            label1, label2 = combo
            data1, data2 = data_for_stats[label1], data_for_stats[label2]
            try:
                if len(data1) > 1 and len(data2) > 1:
                    stat, _ = ranksums(data1, data2)
                    row_data = [label1, label2, f"{stat:.2f}"]
                else:
                    row_data = [label1, label2, "N/A"]
            except ValueError:
                row_data = [label1, label2, "Error"]
            stats_tree.insert("", "end", values=row_data)

# Main GUI class
class GeoWorkflowGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("GEO Samples Labeling Workflow")
        self.geometry("1024x960")

        # --- PATHS ---
        self.data_dir = "/home/mwinn99/GEOMETADB_TOKENS_PLATFORMS"
        self.model_dir = "/home/mwinn99/code/notebooks/MateuszSz/GPL570/finetuned_pubmedbert"
        self.keywords_dir = "/home/mwinn99/code/notebooks/MateuszSz/GPL570"
        self.keyword_cache_file = os.path.join(self.data_dir, "embeddings_cache", "biomed_keyword_embeddings.pkl")

        # --- Instance Variables ---
        self.log_queue = queue.Queue()
        self.gds_conn = None
        self.step1_results_df = None
        self.step2_data_df = None
        self.step1_gse_keywords = {}
        self.gse_to_keep_for_step2 = []
        self.current_extraction_thread = None
        self.current_labeling_thread = None
        self.gpl_datasets = {}
        self.gpl_gene_mappings = {}
        self.is_closing = False
        
        # --- UI SETUP ---
        self._load_geometadb_connection()

        self.model_dir_valid = os.path.exists(self.model_dir)
        self.keywords_dir_valid = os.path.exists(self.keywords_dir) and os.path.exists(os.path.join(self.keywords_dir, "biomed_corpus.txt"))
        
        # ... (The rest of your __init__ continues from here)
        
        # Frame: Load GPL Datasets
        df_frame = tk.LabelFrame(self, text="Load GPL Datasets", padx=10, pady=10)
        df_frame.pack(fill=tk.X, padx=5, pady=5)
        df_frame_row1 = tk.Frame(df_frame); df_frame_row1.pack(fill=tk.X)
        tk.Button(df_frame_row1, text="Load GPL570", command=self.load_gpl570_dataset, width=25).pack(side=tk.LEFT, padx=5, pady=5)
        tk.Button(df_frame_row1, text="Load GPL96", command=self.load_gpl96_dataset, width=25).pack(side=tk.LEFT, padx=5, pady=5)
        # ... (and so on for all your GPL load buttons)

        # Step 1: GSE Extraction
        step1 = tk.LabelFrame(self, text="Step 1: GSE Extraction", padx=10, pady=10)
        step1.pack(fill=tk.X, padx=5, pady=5)
        tk.Label(step1, text="Platform Filter (e.g., GPL570,GPL96):").grid(row=0, column=0, sticky=tk.W)
        self.platform_entry = tk.Entry(step1, width=60); self.platform_entry.grid(row=0, column=1, padx=3, pady=2, sticky=tk.EW)
        tk.Label(step1, text="Filtering Tokens (comma-separated):").grid(row=1, column=0, sticky=tk.W)
        self.filter_entry = tk.Entry(step1, width=60); self.filter_entry.grid(row=1, column=1, padx=3, pady=2, sticky=tk.EW)
        step1.grid_columnconfigure(1, weight=1)
        tk.Button(step1, text="Run GSE Extraction", command=self.start_extraction, bg="green", fg="white").grid(row=2, column=0, columnspan=2, pady=8)

        # Step 1.5: GSE Selection (initially hidden)
        self.gse_frame = tk.LabelFrame(self, text="Step 1.5: GSEs Selected for Step 2", padx=10, pady=10)
        self.gse_listbox = tk.Listbox(self.gse_frame, selectmode=tk.MULTIPLE, width=80, height=4)
        self.gse_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb_gse = tk.Scrollbar(self.gse_frame, command=self.gse_listbox.yview); sb_gse.pack(side=tk.RIGHT, fill=tk.Y)
        self.gse_listbox.config(yscrollcommand=sb_gse.set)
        
        # Step 2: Case-Control Labeling
        lab = tk.LabelFrame(self, text="Step 2: Case-Control Labeling", padx=10, pady=10)
        lab.pack(fill=tk.X, padx=5, pady=5)
        self.step2_status_label = tk.Label(lab, text="Status: Waiting for data from Step 1 or an external file.", fg="blue")
        self.step2_status_label.pack(pady=2)
        bf = tk.Frame(lab); bf.pack(pady=5)
        tk.Button(bf, text="Load External File", width=20, command=self.load_external_file_for_step2, bg="teal", fg="white").pack(side=tk.LEFT, padx=10)
        tk.Button(bf, text="Automatic Labeling", width=20, command=self.run_auto_labeling, bg="blue", fg="white").pack(side=tk.LEFT, padx=10)
        tk.Button(bf, text="Manual Labeling", width=20, command=self.run_manual_labeling, bg="darkgreen", fg="white").pack(side=tk.LEFT, padx=10)

        # Interactive Console & Log
        cons = tk.LabelFrame(self, text="Interactive Console", padx=10, pady=10)
        cons.pack(fill=tk.X, padx=5, pady=5)
        tk.Label(cons, text="Enter command:").pack(side=tk.LEFT)
        self.console_entry = tk.Entry(cons, width=70); self.console_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        tk.Button(cons, text="Submit", command=self.process_console_command).pack(side=tk.LEFT, padx=5)
        self.progressbar = ttk.Progressbar(self, orient="horizontal", mode="determinate", length=600)
        self.progressbar.pack(pady=10, fill=tk.X, padx=5)
        log_frame = tk.LabelFrame(self, text="LOG OUTPUT", padx=10, pady=10)
        log_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        self.log_text = tk.Text(log_frame, wrap=tk.WORD, height=15); self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        sb_log = ttk.Scrollbar(log_frame, command=self.log_text.yview); sb_log.pack(side=tk.RIGHT, fill=tk.Y)
        self.log_text.config(yscrollcommand=sb_log.set)
        
        # Final setup
        if BERTOPIC_AVAILABLE:
            self.enqueue_log("Loading BERTopic model...")
            try:
                self.bertopic_model = BERTopic(embedding_model=SentenceTransformer("all-MiniLM-L6-v2", device=current_device), verbose=False)
                self.enqueue_log("BERTopic model loaded successfully.")
            except Exception as e: self.enqueue_log(f"[CRITICAL ERROR] Failed to load BERTopic: {e}.")
        
        self.uberon_canons, self.disease_canons, self.uberon_map, self.disease_map = load_ontology_sets()
        if self.model_dir_valid and self.keywords_dir_valid:
            self._load_attention_seeker_async()
        
        # Final Setup Calls
        self.after(100, self.process_log_queue)
        self.protocol("WM_DELETE_WINDOW", self.on_closing)

    def _load_geometadb_connection(self):
        """Loads the GEOmetadb.sqlite.gz into a persistent in-memory connection."""
        gz_path = "./GEOmetadb.sqlite.gz"
        if not os.path.exists(gz_path):
            messagebox.showerror("Error", f"GEOmetadb.sqlite.gz not found.")
            self.after(10, self.quit)
            return
        tmp_sql_path = None
        try:
            with tempfile.NamedTemporaryFile(suffix=".sqlite", delete=False) as tmp:
                tmp_sql_path = tmp.name
                with gzip.open(gz_path, "rb") as gzfi:
                    tmp.write(gzfi.read())
            
            disk_conn = sqlite3.connect(tmp_sql_path)
            # --- FIX FOR UTF-8 DECODING ERROR ---
            disk_conn.text_factory = lambda b: b.decode('utf-8', 'replace')
            # --- END FIX ---
            
            self.gds_conn = sqlite3.connect(":memory:")
            self.gds_conn.text_factory = disk_conn.text_factory # Inherit the safe text factory
            disk_conn.backup(self.gds_conn)
            disk_conn.close()
            os.remove(tmp_sql_path)
            print("[INFO] GEOmetadb loaded into memory.")
        except Exception as e:
            messagebox.showerror("DB Error", f"Could not load GEOmetadb into memory: {e}.")
            if self.gds_conn: self.gds_conn.close()
            self.gds_conn = None
            if tmp_sql_path and os.path.exists(tmp_sql_path):
                os.remove(tmp_sql_path)
            self.after(10, self.quit)

    def on_closing(self):
        self.is_closing = True
        if self.current_extraction_thread and self.current_extraction_thread.is_alive():
            if messagebox.askyesno("Exit", "Extraction in progress. Stop and exit?"):
                self.current_extraction_thread.stop()
            else: return
        if self.current_labeling_thread and self.current_labeling_thread.is_alive():
            if messagebox.askyesno("Exit", "Labeling in progress. Stop and exit?"):
                self.current_labeling_thread.stop()
            else: return
        if self.gds_conn: self.gds_conn.close()
        self.destroy()
        

    def enqueue_log(self, msg):
        self.log_queue.put(msg)

    def process_log_queue(self):
        try:
            while True:
                msg = self.log_queue.get_nowait()
                if msg.startswith("PROGRESS:"):
                    try:
                        self.progressbar["value"] = float(msg.split(":",1)[1])
                    except: pass
                else:
                    self.log_text.insert(tk.END, msg + "\n")
                    self.log_text.see(tk.END)

                if self.current_extraction_thread and not self.current_extraction_thread.is_alive():
                    if self.current_extraction_thread.final_df is not None:
                        self.step1_results_df = self.current_extraction_thread.final_df
                    self.current_extraction_thread = None

                if self.current_labeling_thread and not self.current_labeling_thread.is_alive():
                    self.enqueue_log("[INFO] Labeling thread finished.")
                    self.current_labeling_thread = None
        except queue.Empty:
            pass
        if self.winfo_exists():
            self.after(100, self.process_log_queue)

    def _load_attention_seeker_async(self):
        """Loads the AttentionSeeker model and keyword embeddings in a separate thread."""
        self.enqueue_log("Loading AttentionSeeker model and keywords asynchronously…")
        def job():
            if not self.model_dir_valid or not self.keywords_dir_valid:
                self.after(0, lambda: self.enqueue_log("[INFO] AttentionSeeker loading skipped: prerequisite directories are invalid."))
                self.attn = None
                return

            try:
                # Instantiate AttentionSeeker with correct paths
                self.attn = AttentionSeeker(self.model_dir, self.keywords_dir, self.keyword_cache_file)
                if self.attn.model is not None and len(self.attn.keywords) > 0 and self.attn.keyword_embeddings.size > 0:
                    self.after(0, lambda: self.enqueue_log("AttentionSeeker loaded successfully. Bio-aware filtering enabled."))
                else:
                    self.after(0, lambda: self.enqueue_log("AttentionSeeker failed to fully load. Bio-aware filtering will be disabled."))
                    self.attn = None # Explicitly set to None if not fully loaded
            except FileNotFoundError as e_fnf:
                self.after(0, lambda: self.enqueue_log(f"[ERROR] File not found during AttentionSeeker load. Check paths in code. Details: {e_fnf}"))
                self.attn = None
            except Exception as e:
                self.after(0, lambda: self.enqueue_log(f"[ERROR] Error loading AttentionSeeker: {e}. Bio-aware filtering disabled."))
                self.attn = None # Ensure attn is None on failure
        threading.Thread(target=job, daemon=True).start()

        
    def start_extraction(self):
        if self.current_extraction_thread and self.current_extraction_thread.is_alive():
            messagebox.showwarning("Busy", "Extraction already running.")
            return
        pf = self.platform_entry.get()
        toks = self.filter_entry.get()
        if not toks.strip():
            messagebox.showerror("Error", "Enter at least one token.")
            return
        
        self.step1_results_df = None
        self.step2_data_df = None
        self.step1_gse_keywords = {}
        self.gse_to_keep_for_step2 = []
        self.gse_listbox.delete(0, tk.END)
        self.gse_frame.pack_forget()
        self.step2_status_label.config(text="Status: Waiting for data from Step 1...", fg="blue")
        
        self.enqueue_log("[UI] Starting Step 1…")
        self.current_extraction_thread = ExtractionThread(pf, toks, self.enqueue_log, on_finish=self.on_extraction_finish)
        self.current_extraction_thread.start()

    def on_extraction_finish(self):
        """Callback to open the review dialog once the extraction thread is done."""
        if self.current_extraction_thread:
            self.step1_gse_keywords = self.current_extraction_thread.gse_keywords
            self.open_gse_review_dialog()

    def open_gse_review_dialog(self):
        """Opens a Toplevel window to review extracted GSEs in a detailed table."""
        if not self.step1_gse_keywords or self.step1_results_df is None:
            messagebox.showinfo("Nothing to Review", "No GSEs were extracted in Step 1.")
            return

        win = tk.Toplevel(self)
        win.title("Step 1 Results: Review & Select GSEs to Keep")
        win.geometry("900x600")

        gse_summary_data = []
        for gse_id, keywords in self.step1_gse_keywords.items():
            gse_df = self.step1_results_df[self.step1_results_df['series_id'] == gse_id]
            total_samples = len(gse_df)
            matching_samples = int(gse_df['Token_Match'].sum())
            percentage = (matching_samples / total_samples * 100) if total_samples > 0 else 0
            gse_summary_data.append({
                "GSE ID": gse_id, "Matching Samples": matching_samples,
                "Total Samples": total_samples, "% Matched": f"{percentage:.1f}%",
                "Keywords": ", ".join(keywords)
            })
        
        gse_summary_data.sort(key=lambda x: x["Matching Samples"], reverse=True)

        top_frame = ttk.Frame(win, padding="5")
        top_frame.pack(fill=tk.BOTH, expand=True)
        columns = ["GSE ID", "Matching Samples", "Total Samples", "% Matched", "Keywords"]
        tree = ttk.Treeview(top_frame, columns=columns, show="headings", selectmode="extended")
        
        for col in columns:
            tree.heading(col, text=col)
            width = 400 if col == "Keywords" else 100 if col == "GSE ID" else 80
            anchor = tk.W if col in ["GSE ID", "Keywords"] else tk.CENTER
            tree.column(col, width=width, anchor=anchor, minwidth=60)

        for item in gse_summary_data:
            tree.insert("", tk.END, values=[item[col] for col in columns])

        vsb = ttk.Scrollbar(top_frame, orient="vertical", command=tree.yview)
        hsb = ttk.Scrollbar(top_frame, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        top_frame.grid_rowconfigure(0, weight=1); top_frame.grid_columnconfigure(0, weight=1)

        bottom_frame = ttk.Frame(win, padding="5")
        bottom_frame.pack(fill=tk.X)

        def save_selection_and_close():
            selected_items = tree.selection()
            if not selected_items:
                messagebox.showerror("No Selection", "Please select one or more GSEs to continue.", parent=win)
                return
            
            kept_gses = [tree.item(item, "values")[0] for item in selected_items]
            self.gse_to_keep_for_step2 = kept_gses
            self.enqueue_log(f"[UI] User selected {len(kept_gses)} GSEs to keep: {', '.join(kept_gses)}")
            
            self.gse_listbox.delete(0, tk.END)
            for gse_id in kept_gses:
                kws = self.step1_gse_keywords.get(gse_id, [])
                self.gse_listbox.insert(tk.END, f"{gse_id}: {', '.join(kws[:5])}...")
            
            self.gse_frame.pack(fill=tk.X, padx=5, pady=5)
            self.step2_status_label.config(text=f"Status: Loaded {len(kept_gses)} GSEs from Step 1.", fg="green")
            win.destroy()

        ttk.Button(bottom_frame, text="Use Selected GSEs for Step 2", command=save_selection_and_close).pack(pady=5)
        win.transient(self); win.grab_set(); self.wait_window(win)

    def load_external_file_for_step2(self):
        """Loads an external CSV file for direct use in Step 2."""
        filepath = filedialog.askopenfilename(
            title="Select CSV file for Step 2 Labeling",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        if not filepath: return

        try:
            self.step2_data_df = pd.read_csv(filepath)
            self.step1_results_df = None
            self.gse_to_keep_for_step2 = []
            self.gse_listbox.delete(0, tk.END)
            self.gse_frame.pack_forget()
            
            self.step2_status_label.config(text=f"Loaded: {os.path.basename(filepath)} ({len(self.step2_data_df)} rows)", fg="green")
            self.enqueue_log(f"[UI] Loaded external file. Ready for Step 2.")
        except Exception as e:
            messagebox.showerror("File Load Error", f"Failed to load file:\n{e}")
            self.step2_data_df = None
            self.step2_status_label.config(text="Error loading file.", fg="red")

    def _start_labeling(self, mode: str):
        """A unified function to prepare data and launch the labeling thread."""
        if self.current_labeling_thread and self.current_labeling_thread.is_alive():
            messagebox.showwarning("Busy", "A labeling process is already running.")
            return

        source_df, data_source_name = (None, "")
        if self.step1_results_df is not None and self.gse_to_keep_for_step2:
            source_df = self.step1_results_df[self.step1_results_df["series_id"].isin(self.gse_to_keep_for_step2)].copy()
            data_source_name = "Step 1 Results"
        elif self.step2_data_df is not None:
            source_df = self.step2_data_df.copy()
            data_source_name = "External File"
        else:
            messagebox.showerror("No Data", "Please run Step 1 and select GSEs, or load an external file.")
            return
            
        if source_df.empty:
            messagebox.showerror("Empty Dataset", "The selected data source is empty.")
            return

        available_cols_str = "\n".join(source_df.columns)
        
        gsm_col = simpledialog.askstring("Input Required", f"Which column contains GSM IDs?\n\nAvailable:\n{available_cols_str}", parent=self)
        if not gsm_col or gsm_col not in source_df.columns:
            messagebox.showerror("Invalid Column", f"Column '{gsm_col}' not found.", parent=self)
            return
        
        source_df.rename(columns={gsm_col: "gsm"}, inplace=True)

        if 'Final_Text' not in source_df.columns:
            text_cols_str = simpledialog.askstring("Input Required", "Which column(s) for descriptive text?\n(Separate with commas)", parent=self)
            if not text_cols_str: return
            text_cols = [c.strip() for c in text_cols_str.split(',') if c.strip() in source_df.columns]
            if not text_cols:
                messagebox.showerror("Invalid Columns", "None of the specified columns were found.", parent=self)
                return
            source_df['Final_Text'] = source_df[text_cols].astype(str).agg(' '.join, axis=1)
        
        case_toks, ctrl_toks = "", ""
        if mode == 'automatic':
            case_toks = simpledialog.askstring("Case Tokens", "Enter CASE tokens:", parent=self)
            if not case_toks: return
            ctrl_toks = simpledialog.askstring("Control Tokens", "Enter CONTROL tokens:", parent=self)
            if not ctrl_toks: return

        self.enqueue_log(f"[UI] Starting Step 2 ({mode.capitalize()}) using {data_source_name}...")
        folder_base = "external_data" if data_source_name == "External File" else "_".join(sorted(self.gse_to_keep_for_step2))
        
        self.current_labeling_thread = LabelingThread(
            df_sub=source_df, folder_name=folder_base, log_cb=self.enqueue_log,
            mode=mode, case_toks=case_toks, ctrl_toks=ctrl_toks
        )
        self.current_labeling_thread.start()

    def run_auto_labeling(self):
        """Wrapper to start automatic labeling."""
        self._start_labeling(mode='automatic')

    def run_manual_labeling(self):
        """Wrapper to start manual labeling."""
        self._start_labeling(mode='manual')

    def process_console_command(self):
        """Processes commands entered in the interactive console."""
        cmd = self.console_entry.get().strip()
        self.console_entry.delete(0, tk.END)
        self.enqueue_log(f"[Console] {cmd}")
        c = cmd.lower()
        if c == "help":
            self.enqueue_log("Commands: help, list_gses, clear_log, test_cuda")
        elif c == "list_gses":
            if self.step1_gse_keywords:
                self.enqueue_log(f"GSEs: {', '.join(self.step1_gse_keywords.keys())}")
            else:
                self.enqueue_log("No GSEs extracted.")
        elif c == "clear_log":
            self.log_text.delete(1.0, tk.END)
        elif c == "test_cuda":
            self.enqueue_log(f"Cuda available: {torch.cuda.is_available()}")
            if torch.cuda.is_available():
                self.enqueue_log(f"Device count: {torch.cuda.device_count()}")
                self.enqueue_log(f"Current device: {torch.cuda.current_device()}")
                self.enqueue_log(f"Device name: {torch.cuda.get_device_name(torch.cuda.current_device())}")
        else:
            self.enqueue_log(f"[Console] Unknown: {cmd}")

    def _load_gpl_data(self, gpl_name: str, file_path: str, metadata_path: Optional[str]=None):
        """
        Loads GPL gene expression data from a gzipped CSV file.
        Occasionally merges with metadata and performs min-max scaling for numerical columns.
        Identifies gene columns based on numerical type and common exclusions.
        """
        try:
            self.enqueue_log(f"[UI] Loading {gpl_name} from {file_path}...")
            if gpl_name in self.gpl_datasets:
                if not messagebox.askyesno("Reload?", f"{gpl_name} already loaded. Reload?"):
                    self.enqueue_log(f"[UI] {gpl_name} load cancelled.")
                    return
            data_df = pd.read_csv(file_path, compression="gzip", low_memory=False)
            self.enqueue_log(f"[DEBUG] Loaded {gpl_name} with initial columns: {list(data_df.columns)}")

            if gpl_name == "GPL570" and metadata_path:
                if not os.path.exists(metadata_path):
                    messagebox.showerror(f"{gpl_name} Error", f"Metadata not found: {metadata_path}")
                    self.enqueue_log(f"[UI] {gpl_name} metadata missing.")
                    return
                meta_df = pd.read_csv(metadata_path, compression="gzip", low_memory=False)
                num_cols_data = data_df.select_dtypes(include=np.number).columns
                if len(meta_df) == len(data_df):
                    scaler = MinMaxScaler().fit_transform(data_df[num_cols_data]) if len(num_cols_data)>0 else None
                    if scaler is not None:
                        df_scaled = pd.DataFrame(scaler, columns=num_cols_data, index=data_df.index)
                        rest = data_df.select_dtypes(exclude=np.number)
                        # Concatenate metadata, scaled numerical data, and remaining non-numerical columns
                        # Ensure no duplicate columns from `rest` if already in `meta_df`
                        data_df = pd.concat([meta_df, df_scaled, rest.loc[:, ~rest.columns.isin(meta_df.columns)]], axis=1)
                    else:
                        data_df = pd.concat([meta_df, data_df], axis=1)
                else:
                    self.enqueue_log(f"[UI] {gpl_name}: metadata row mismatch; loading only data.")
                    
            # Standardize GSM column name
            if "gsm" in data_df.columns and "GSM" not in data_df.columns:
                data_df.rename(columns={"gsm":"GSM"}, inplace=True)
            if "GSM" not in data_df.columns:
                cand = [col for col in data_df.columns if "GSM" in col.upper() or "SAMPLE" in col.upper() or col.upper()=="ID"]
                if cand:
                    data_df.rename(columns={cand[0]:"GSM"}, inplace=True)
                    self.enqueue_log(f"[UI] Renamed {cand[0]}→GSM")

            data_df["GSM"] = data_df["GSM"].astype(str).str.upper()

            self.gpl_datasets[gpl_name] = data_df
            loaded_df = data_df # Reference to the fully loaded and processed DataFrame

            gene_map = {}
            all_numeric_cols = loaded_df.select_dtypes(include=np.number).columns.tolist()

            # --- FIX APPLIED HERE: Start gene data from 8th column (index 7) ---
            if len(all_numeric_cols) >= 8: # Ensure there are at least 8 numerical columns
                # Take numerical columns starting from the 8th column (index 7) onwards
                gene_cols_to_map = all_numeric_cols[7:]  
                self.enqueue_log(f"[DEBUG] Identifying gene columns from index 7 onwards. Found {len(gene_cols_to_map)} gene columns.")
            else:
                self.enqueue_log(f"[WARN] Less than 8 numerical columns in {gpl_name}. Using all numerical columns as genes. Columns: {all_numeric_cols}")
                gene_cols_to_map = all_numeric_cols
            # ------------------------------------------------------------------

            for col in gene_cols_to_map:
                # Still filter out common non-gene identifiers, just in case they are numeric and slipped through
                # This makes it more robust if the 8-column assumption is slightly off for some GPLs.
                if not any(x in col.upper() for x in ["GSM","ID","SAMPLE","_X","_Y","COORD","UNNAMED","AVG_SIGNAL","DETECTION_PVAL"]): # Added more exclusions
                    gene_map[col.upper()] = col
                    
            if not gene_map:
                self.enqueue_log(f"[ERROR] No valid gene columns identified for {gpl_name} after filtering and starting from column 8.")


            self.gpl_gene_mappings[gpl_name] = gene_map
            self.enqueue_log(f"[UI] {gpl_name} loaded: {len(loaded_df)} rows, {len(loaded_df.columns)} cols; {len(gene_map)} genes identified.")
            messagebox.showinfo(gpl_name, f"Loaded {gpl_name}: {len(loaded_df)} rows, {len(loaded_df.columns)} cols. Identified {len(gene_map)} gene columns.")
        except FileNotFoundError:
            messagebox.showerror(f"{gpl_name} Error", f"File not found: {file_path}")
            self.enqueue_log(f"[UI] {gpl_name} file missing.")
        except Exception as e:
            messagebox.showerror(f"{gpl_name} Error", f"Error loading {gpl_name}:\n{e}")
            self.enqueue_log(f"[UI] {gpl_name} load error: {e}")
            import traceback
            self.enqueue_log(traceback.format_exc()) # Log full traceback for load errors

    def load_gpl570_dataset(self):
        meta_path = "/home/mwinn99/code/notebooks/MateuszSz/GPL570/GPL570_metadata.csv.gz"
        data_path = "/home/mwinn99/code/notebooks/MateuszSz/GPL570/gpl570_all_samples_normalized_non_scaled_with_nans.csv.gz"
        self._load_gpl_data("GPL570", data_path, metadata_path=meta_path)
    def load_gpl96_dataset(self):
        file_path = "/home/mwinn99/code/notebooks/MateuszSz/GPL96/gpl96_normalized_scaled.csv.gz"
        self._load_gpl_data("GPL96", file_path)
    def load_gpl6947_dataset(self):
        file_path = "/home/mwinn99/hf_datasets_biovdb/GPL6947/gpl6947_all_samples_normalized_scaled_with_nans.csv.gz"
        self._load_gpl_data("GPL6947", file_path)
    def load_gpl7202_dataset(self):
        file_path = "/home/mwinn99/hf_datasets_biovdb/GPL7202/gpl7202_all_samples_normalized_scaled_with_nans.csv.gz"
        self._load_gpl_data("GPL7202", file_path)
    def load_gpl6885_dataset(self):
        file_path = "/home/mwinn99/hf_datasets_biovdb/GPL6885/gpl6885_all_samples_normalized_scaled_with_nans.csv.gz"
        self._load_gpl_data("GPL6885", file_path)
    def load_gpl1261_dataset(self):
        file_path = "/home/mwinn99/hf_datasets_biovdb/GPL1261/gpl1261_all_samples_normalized_scaled_with_nans.csv.gz"
        self._load_gpl_data("GPL1261", file_path)
    def load_gpl10558_dataset(self):
        file_path = "/home/mwinn99/hf_datasets_biovdb/GPL10558/gpl10558_all_samples_normalized_scaled_with_nans.csv.gz"
        self._load_gpl_data("GPL10558", file_path)

    def open_compare_window(self):
            """
            Instantiates and opens the Compare Distributions window.
            Allows opening even if no platforms are loaded, for custom file analysis.
            """
            if not self.gpl_datasets:
                messagebox.showinfo(
                    "No Platforms Loaded", 
                    "No GPL datasets are loaded in the main app.\n\nYou can still use this feature by loading your own custom data file from within the comparison window.", 
                    parent=self
                )
            
            # This prevents opening multiple instances of the same window
            if not hasattr(self, 'compare_window') or not self.compare_window.winfo_exists():
                # Pass self as parent and as the app_ref for callbacks
                self.compare_window = CompareDistributionsWindow(self, self)
            
            # Bring the window to the front
            self.compare_window.lift()
            self.compare_window.focus_force()


    def show_gene_distribution_popup(self):
        """Opens the Gene Distribution Explorer popup window."""
        import tkinter as tk_local
        from tkinter import messagebox as mb_local

        if not self.gds_conn:
            mb_local.showerror("Database Error", "GEOmetadb connection not established.", parent=self)
            return

        self.gene_dist_popup_root = tk_local.Toplevel(self)
        self.gene_dist_popup_root.title("Gene Distribution Explorer")
        self.gene_dist_popup_root.geometry("1100x800")

        top_controls_frame = tk_local.Frame(self.gene_dist_popup_root)
        top_controls_frame.pack(fill=X, padx=5, pady=5)

        plat_frame = tk_local.Frame(top_controls_frame)
        plat_frame.pack(fill=X, pady=2)
        tk_local.Label(plat_frame, text="Platforms:").pack(side=LEFT, padx=(0, 5))
        gpls = ["GPL96","GPL570","GPL6947","GPL7202","GPL6885","GPL1261","GPL10558"]
        self.gpl_selection_vars = {p: tk_local.BooleanVar() for p in gpls}
        for p in gpls:
            cb = tk_local.Checkbutton(plat_frame, text=p, variable=self.gpl_selection_vars[p],
                                     state=NORMAL if p in self.gpl_datasets and self.gpl_datasets[p] is not None else DISABLED)
            cb.pack(side=LEFT, padx=3)

        gene_sim_frame = tk_local.Frame(top_controls_frame)
        gene_sim_frame.pack(fill=X, pady=2)

        # Updated label for multiple genes
        tk_local.Label(gene_sim_frame, text="Gene Symbol(s) (comma-separated):").pack(side=LEFT, padx=(0, 5))
        self.current_gene_entry = tk_local.Entry(gene_sim_frame, width=40) # Increased width
        self.current_gene_entry.pack(side=LEFT, padx=5)

        tk.Label(gene_sim_frame, text="Similarity Threshold (%):").pack(side=LEFT, padx=(15, 5))
        self.similarity_thr_var_popup = tk.IntVar(value=75)
        self.similarity_thr_label_popup = ttk.Label(gene_sim_frame, text="75%", width=4)
        ttk.Scale(gene_sim_frame, from_=0, to=100, variable=self.similarity_thr_var_popup,
                  length=160,
                  command=lambda v: self.similarity_thr_label_popup.config(text=f"{int(float(v))}%")
                 ).pack(side=LEFT, padx=5)
        self.similarity_thr_label_popup.pack(side=LEFT)

        btn_frame = tk_local.Frame(self.gene_dist_popup_root)
        btn_frame.pack(fill=X, pady=10)

        tk_local.Button(btn_frame,text="Plot Distributions",bg="#2ca02c",fg="white",command=self._plot_histograms).pack(side=tk_local.LEFT,padx=5)
        tk_local.Button(btn_frame,text="Analyze LEFT Tail",bg="#d62728",fg="white",command=lambda: self._analyze_tail("left")).pack(side=tk_local.LEFT,padx=5)
        tk_local.Button(btn_frame,text="Analyze RIGHT Tail",bg="#d62728",fg="white",command=lambda: self._analyze_tail("right")).pack(side=tk_local.LEFT,padx=5)
        tk_local.Button(btn_frame,text="Analyze Full Distribution",bg="DarkGoldenrod",fg="white",command=self._analyze_full_distribution).pack(side=tk_local.LEFT,padx=5)
        tk_local.Button(btn_frame,text="Run Batch Tail Analysis for ALL Genes",bg="teal",fg="white",command=self._batch_process_all_genes_for_tails).pack(side=tk_local.LEFT,padx=5)
        tk_local.Button(btn_frame, text="Batch Genes - Negative Control", bg="gray", fg="white", command=self._batch_negative_control).pack(side=tk.LEFT, padx=5)
        btn_frame.pack(fill=X, pady=10)

        self.gene_dist_out_frame = tk_local.Frame(self.gene_dist_popup_root)
        self.gene_dist_out_frame.pack(fill=BOTH, expand=True, padx=5, pady=5)

        self._current_popup_figs = {}

        def _on_close_popup_handler():
            import matplotlib.pyplot as plt_local
            for key in list(self._current_popup_figs.keys()):
                fig, canv_widget, tool = self._current_popup_figs.pop(key)
                canv_widget.destroy(); tool.destroy(); plt_local.close(fig) # Close Matplotlib figure
            self._axis_map_dist_plot = {}
            self._threshold_map_dist_plot = {}
            self._tok_df_cache = {}
            self.current_genes = [] # Reset to empty list
            self.gene_dist_popup_root.destroy()

        self.gene_dist_popup_root.protocol("WM_DELETE_WINDOW", _on_close_popup_handler)
        self.gene_dist_popup_root.transient(self); self.gene_dist_popup_root.grab_set()

    def _pack_figure_into_popup(self, fig, parent_frame, key):
        """Packs a matplotlib figure into the Tkinter popup frame."""
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        import matplotlib.pyplot as plt_local
        if key in self._current_popup_figs:
            old_fig, old_canv, old_tool = self._current_popup_figs.pop(key)
            old_canv.destroy(); old_tool.destroy();
            plt_local.close(old_fig)
        canv = FigureCanvasTkAgg(fig, master=parent_frame)
        widget = canv.get_tk_widget(); widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        tool = NavigationToolbar2Tk(canv, parent_frame); tool.update(); tool.pack(side=tk.TOP, fill=tk.X)
        canv.draw()
        self._current_popup_figs[key] = (fig, widget, tool)

    def _plot_histograms(self):
        """
        Plots gene expression histograms for selected platforms and multiple genes.
        Includes dynamic binning and individual y-axis scaling for better visibility.
        """
        import matplotlib.pyplot as plt_local

        # Close any previously opened matplotlib figures associated with this popup
        for key in list(self._current_popup_figs.keys()):
            fig, canv, tool = self._current_popup_figs.pop(key)
            canv.destroy(); tool.destroy()
            plt_local.close(fig) # Explicitly close the figure

        gene_input_str = self.current_gene_entry.get().strip()
        if not gene_input_str:
            messagebox.showerror("Input Missing", "Enter at least one gene symbol.", parent=self.gene_dist_popup_root); return
            
        self.current_genes = [g.strip().upper() for g in gene_input_str.split(',') if g.strip()]
        if not self.current_genes:
            messagebox.showerror("Input Missing", "Invalid gene symbol(s).", parent=self.gene_dist_popup_root); return

        sel_plats = [p for p,v in self.gpl_selection_vars.items() if v.get()]
        if not sel_plats:
            messagebox.showerror("Input Missing", "Select platform(s).", parent=self.gene_dist_popup_root); return

        # Determine grid size for subplots
        num_genes = len(self.current_genes)
        num_platforms = len(sel_plats)

        # Adjust figsize dynamically for better readability with more plots
        fig_width = 6 * num_platforms
        fig_height = 5.5 * num_genes 
        fig, axs = plt_local.subplots(num_genes, num_platforms, 
                                     figsize=(fig_width, fig_height), 
                                     squeeze=False) # squeeze=False ensures axs is always 2D

        self._axis_map_dist_plot = {} # Clear previous mappings for interactive plots

        valid_plots_count = 0
        for r_idx, gene in enumerate(self.current_genes):
            for c_idx, plat in enumerate(sel_plats):
                ax = axs[r_idx, c_idx] # Get the specific subplot axis

                dfg = self.gpl_datasets.get(plat)
                gmap = self.gpl_gene_mappings.get(plat,{})
                col = gmap.get(gene)

                if dfg is None or dfg.empty:
                    ax.set_title(f"{plat}\nNot loaded", color='red'); ax.axis("off"); continue
                if col is None:
                    ax.set_title(f"{plat}-{gene}\nNot found", color='red'); ax.axis("off"); continue
                
                expr = dfg[col].dropna().astype(float)
                if expr.empty:
                    ax.set_title(f"{plat}-{col}\nNo data", color='orange'); ax.axis("off"); continue
                
                # --- Dynamic Binning Calculation ---
                # Use Freedman-Diaconis rule for bin width
                q25, q75 = np.percentile(expr, [25, 75])
                iqr = q75 - q25
                bin_width = 2 * iqr / (len(expr)**(1/3))
                
                if bin_width == 0: # Handle cases with no variance, fall back to a reasonable small number of bins
                    num_bins = 10 
                else:
                    num_bins = int(np.ceil((expr.max() - expr.min()) / bin_width))
                num_bins = max(10, min(num_bins, 100)) # Ensure minimum 10, maximum 100 bins for visual clarity
                self.enqueue_log(f"[DEBUG] {plat}-{gene}: Calculated {num_bins} bins.")
                # -----------------------------------

                dist = analyze_gene_distribution(expr)
                mu, sd = expr.mean(), expr.std()
                low, high = mu - 3*sd, mu + 3*sd
                self._threshold_map_dist_plot[ax] = (low, high) # Store thresholds per axis for click analysis

                counts, bins, patches = ax.hist(expr, bins=num_bins, edgecolor="black", alpha=0.7)
                
                # --- AUTOMATIC Y-AXIS SCALING PER SUBPLOT (Enhanced) ---
                if len(counts) > 0 and max(counts) > 0:
                    ax.set_ylim(0, max(counts) * 1.1) # Scale Y-axis to fit its own data with a 10% buffer
                else:
                    ax.set_ylim(0, 1) # Fallback if no counts
                # --------------------------------------------------------

                for p in patches: p.set_facecolor("cornflowerblue")
                
                ax.set_title(f"{plat}-{col}\n({dist}, n={len(expr)})", fontsize=10)
                ax.set_xlabel("Expression"); ax.set_ylabel("Freq")
                ax.axvline(low, linestyle="--", color='red', label=f'Low: {low:.2f}'); 
                ax.axvline(high, linestyle="--", color='red', label=f'High: {high:.2f}')
                ax.legend()
                
                # Store data needed for _on_click_histogram_bin
                self._axis_map_dist_plot[ax] = (bins, patches, dfg.copy(), plat, col, gene)
                valid_plots_count += 1
        
        if valid_plots_count == 0:
            messagebox.showinfo("No Data", "Nothing to plot for the selected genes and platforms.", parent=self.gene_dist_popup_root)
            plt_local.close(fig) # Close the empty figure
            return

        fig.tight_layout() # Adjust layout to prevent overlapping titles/labels
        self._pack_figure_into_popup(fig, self.gene_dist_out_frame, "dist_plot")
        fig.canvas.mpl_connect("button_press_event", self._on_click_histogram_bin)

    def _on_click_histogram_bin(self, event):
        """
        Callback for clicking on a histogram bin to analyze samples within that bin.
        This will trigger the full analysis and saving for the *specific* gene/platform
        and bin clicked.
        """
        ax = event.inaxes
        if ax not in self._axis_map_dist_plot or event.xdata is None:
            return

        # Unpack the stored data, including the gene symbol
        bins, patches, dfg, plat, col, gene_symbol = self._axis_map_dist_plot[ax]
        x = event.xdata
        idx = None

        # Find which bin was clicked
        for i in range(len(bins) - 1):
            if bins[i] <= x < bins[i+1]:
                idx = i
                break
        if idx is None:
            return

        # Highlight the clicked bin visually
        for p in patches:
            p.set_facecolor("cornflowerblue") # Reset all to default
        patches[idx].set_facecolor("salmon") # Highlight clicked bin
        ax.figure.canvas.draw_idle() # Redraw the canvas

        low, high = bins[idx], bins[idx+1]
        
        # Filter samples for the selected bin
        sel_df_base = dfg[(dfg[col] >= low) & (dfg[col] < high)].copy()
        if sel_df_base.empty:
            messagebox.showinfo("Empty Bin", "No samples in this bin.", parent=self.gene_dist_popup_root)
            return
        
        sel_df_base["Expression"] = sel_df_base[col]
        sel_df_base["Platform"] = plat

        # --- Create a unique output directory for this specific bin analysis ---
        # The directory now includes gene symbol, platform, and bin info
        analysis_folder_name = f"{gene_symbol}_{plat}_BIN_{idx}_{low:.2f}-{high:.2f}".replace('.', 'p').replace('-', 'm') # Sanitize for path
        out_dir = os.path.join("NEW_RESULTS_ROOT", "GeneDistributionAnalysis", analysis_folder_name)
        os.makedirs(out_dir, exist_ok=True)
        base_name_for_files = f"{gene_symbol}_{plat}_BIN_{idx}"
        self.enqueue_log(f"[INFO] Analyzing bin {idx} for gene {gene_symbol} ({plat}). Results saving to: {out_dir}")

        # --- 1. Density Plot (REORDERED to be FIRST) ---
        background_expressions = dfg.loc[~dfg.index.isin(sel_df_base.index), col].dropna()
        self._display_density_plot_popup(
            selected_expressions=sel_df_base[col].dropna(),
            background_expressions=background_expressions,
            title=f"Density Plot: {gene_symbol} - Bin {idx} ({plat})",
            output_dir=out_dir,
            base_name=base_name_for_files
        )

        # Prepare DataFrame with tokens for word cloud generation
        sel_df_with_tokens = self._prepare_df_with_tokens(sel_df_base, plat)
        self.show_table_popup_only_table(sel_df_with_tokens, plat, gene_symbol, f"Bin {idx} Samples (n={len(sel_df_with_tokens)} GSMs)")

        # Get current similarity threshold from UI
        similarity_threshold = self.similarity_thr_var_popup.get() / 100.0
        attn_available = self.attn is not None and self.attn.model is not None and self.attn.keywords and self.attn.keyword_embeddings.size > 0
        
        # Get background token data
        bg_df = self._get_background_df(plat, sel_df_with_tokens["GSM"])

        # --- Generate and save word cloud images and stats ---
        self._generate_and_save_wordcloud_images(
            sel_df_with_tokens, bg_df, plat, gene_symbol, f"Bin {idx}",
            out_dir, similarity_threshold, attn_available, display_mode=True
        )
        
        # --- New: Call for cluster analysis and display after word clouds ---
        if attn_available:
            temp_sel_tokens = self.flatten_tokens(sel_df_with_tokens["tokens"])
            temp_cnt_sel, temp_N_sel = Counter(temp_sel_tokens), len(temp_sel_tokens)
            temp_cnt_bg, temp_N_bg = Counter(self.flatten_tokens(bg_df["tokens"])), len(self.flatten_tokens(bg_df["tokens"]))

            # Convert counters to pandas Series for vectorized operations
            sel_counts_series = pd.Series(temp_cnt_sel, name='sel_count')
            bg_counts_series = pd.Series(temp_cnt_bg, name='bg_count')
            all_tokens_combined = pd.concat([sel_counts_series, bg_counts_series], axis=1).fillna(0)
            all_tokens_combined.columns = ['sel_count', 'bg_count']

            total_counts = all_tokens_combined['sel_count'] + all_tokens_combined['bg_count']
            
            # Calculate Z-score: (count_in_sel - count_in_bg) / sqrt(count_in_sel + count_in_bg)
            # Handle cases where total_counts is 0 to avoid ZeroDivisionError
            z_scores = np.where(total_counts > 0,
                                 (all_tokens_combined['sel_count'] - all_tokens_combined['bg_count']) / np.sqrt(total_counts),
                                 0) # If total_counts is 0, z-score is 0 (no difference)

            # Calculate one-tailed p-value for enrichment (sel_count > bg_count)
            p_values = 1 - norm.cdf(z_scores)

            # Calculate percentages
            sel_pct = np.where(temp_N_sel > 0, all_tokens_combined['sel_count'] / temp_N_sel, 0)
            bg_pct = np.where(temp_N_bg > 0, all_tokens_combined['bg_count'] / temp_N_bg, 0)

            temp_full_stats_df = pd.DataFrame({
                'token': all_tokens_combined.index,
                'sel_count': all_tokens_combined['sel_count'].astype(int),
                'bg_count': all_tokens_combined['bg_count'].astype(int),
                'p_value': p_values,
                'sel_pct': sel_pct,
                'bg_pct': bg_pct
            })

            if not temp_full_stats_df.empty:
                temp_full_stats_df = temp_full_stats_df.sort_values(by=['p_value', 'sel_count'], ascending=[True, False])
                
            enriched_tokens_df = temp_full_stats_df[(temp_full_stats_df['p_value'] < 0.05) & (temp_full_stats_df['sel_pct'] > temp_full_stats_df['bg_pct'])].copy()

            bio_specific_tokens_df = pd.DataFrame()
            if not enriched_tokens_df.empty:
                enriched_tok_list_for_sim = enriched_tokens_df['token'].tolist()
                sim_scores_dict = self.attn.calculate_similarity(enriched_tok_list_for_sim, similarity_threshold, self.token_cache)
                bio_specific_tokens_df = enriched_tokens_df[enriched_tokens_df['token'].isin(sim_scores_dict.keys())].copy()
                if not bio_specific_tokens_df.empty:
                    bio_specific_tokens_df['sim_score'] = bio_specific_tokens_df['token'].map(sim_scores_dict)
                    bio_specific_tokens_df = bio_specific_tokens_df.sort_values(by=['sim_score', 'p_value'], ascending=[False, True])

            if not bio_specific_tokens_df.empty:
                self.enqueue_log(f"[Clustering] Found {len(bio_specific_tokens_df)} bio-specific tokens for clustering for {gene_symbol} Bin {idx} ({plat}).")
                self._perform_token_clustering_and_display(
                    bio_specific_tokens_df,
                    f"{gene_symbol} Bin {idx} ({plat})",
                    out_dir
                )
            else:
                self.enqueue_log(f"[Clustering] No bio-specific tokens found for clustering for {gene_symbol} Bin {idx} ({plat}). Skipping cluster plots.")


    def _analyze_tail(self, which_tail: str):
        """
        Analyzes samples in the left or right tails of the gene distribution for multiple genes.
        Results (density plots, word clouds, cluster analysis) are saved for each gene/platform.
        """
        if not self._axis_map_dist_plot:
            messagebox.showerror("Plot Not Generated", "Please plot distributions first to define gene/platform contexts.", parent=self.gene_dist_popup_root)
            return

        gene_input_str = self.current_gene_entry.get().strip()
        if not gene_input_str:
            messagebox.showerror("Input Missing", "Enter at least one gene symbol.", parent=self.gene_dist_popup_root); return
            
        genes_to_analyze = [g.strip().upper() for g in gene_input_str.split(',') if g.strip()]
        if not genes_to_analyze:
            messagebox.showerror("Input Missing", "Invalid gene symbol(s).", parent=self.gene_dist_popup_root); return

        sel_plats = [p for p,v in self.gpl_selection_vars.items() if v.get()]
        if not sel_plats:
            messagebox.showerror("Input Missing", "Select platform(s).", parent=self.gene_dist_popup_root); return

        similarity_threshold = self.similarity_thr_var_popup.get() / 100.0
        attn_available = self.attn is not None and self.attn.model is not None and self.attn.keywords and self.attn.keyword_embeddings.size > 0

        self.enqueue_log(f"[Tail Analysis] Starting {which_tail.upper()} tail analysis for selected genes and platforms...")
        
        # Iterate through all selected genes and platforms to perform analysis
        for gene_symbol in genes_to_analyze:
            for plat in sel_plats:
                # Find the stored data for the current gene and platform from _plot_histograms run
                current_ax_data = None
                # We need to iterate through the stored axis map to find the correct data for this gene/platform pair
                # self._axis_map_dist_plot stores: {ax_object: (bins, patches, dfg_copy, plat_name, col_name, gene_symbol)}
                for ax_obj, data_tuple in self._axis_map_dist_plot.items():
                    # Check if platform and gene symbol match
                    if data_tuple[3] == plat and data_tuple[5] == gene_symbol:
                        current_ax_data = data_tuple
                        break

                if current_ax_data is None:
                    self.enqueue_log(f"[Tail] No plot data (or gene/column) found for gene '{gene_symbol}' on platform '{plat}'. Skipping tail analysis for this combination.")
                    continue

                # Unpack the relevant data
                bins, patches, df_plat, plat_name_retrieved, gene_col, _ = current_ax_data
                
                # Retrieve the thresholds previously calculated and stored for this axis
                low_thr, high_thr = self._threshold_map_dist_plot.get(ax_obj, (None, None))
                if low_thr is None or high_thr is None:
                    self.enqueue_log(f"[Tail] Thresholds not calculated for '{gene_symbol}' on '{plat}'. Skipping tail analysis for this combination.")
                    continue

                # Re-color histograms for emphasis (only for the *current* plot within the UI)
                # This visual feedback only applies to the last plot generated by _plot_histograms
                # but the underlying logic applies to all gene/platform combinations.
                if ax_obj in self._axis_map_dist_plot: # Only if this subplot is still rendered
                    current_patches = self._axis_map_dist_plot[ax_obj][1] # Get patches reference
                    for i, p in enumerate(current_patches):
                        mid = (bins[i] + bins[i+1]) / 2
                        is_left, is_right = mid <= low_thr, mid >= high_thr
                        highlight = (is_left if which_tail == "left" else is_right)
                        p.set_facecolor("orangered" if highlight else "cornflowerblue")
                    ax_obj.figure.canvas.draw_idle()

                # Filter samples for the selected tail
                tail_mask = (df_plat[gene_col] <= low_thr) if which_tail == "left" else (df_plat[gene_col] >= high_thr)
                sel_df_base = df_plat.loc[tail_mask].copy()

                if sel_df_base.empty:
                    self.enqueue_log(f"[Tail] No {which_tail}-tail samples found for gene '{gene_symbol}' in '{plat}'. Skipping further analysis for this tail.")
                    continue

                sel_df_base["Expression"] = sel_df_base[gene_col]
                sel_df_base["Platform"] = plat

                # --- Create a unique output directory for this specific tail analysis ---
                analysis_folder_name = f"{gene_symbol}_{plat}_{which_tail.upper()}_TAIL"
                out_dir = os.path.join("NEW_RESULTS_ROOT", "GeneDistributionAnalysis", analysis_folder_name)
                os.makedirs(out_dir, exist_ok=True)
                base_name_for_files = f"{gene_symbol}_{plat}_{which_tail.upper()}_TAIL"
                self.enqueue_log(f"[INFO] Analyzing {which_tail.upper()} tail for gene {gene_symbol} ({plat}). Results saving to: {out_dir}")

                # --- 1. Density Plot (REORDERED to be FIRST) ---
                background_expressions = df_plat.loc[~df_plat.index.isin(sel_df_base.index), gene_col].dropna()
                self._display_density_plot_popup(
                    selected_expressions=sel_df_base[gene_col].dropna(),
                    background_expressions=background_expressions,
                    title=f"Density Plot: {gene_col} - {which_tail.capitalize()} Tail ({plat})",
                    output_dir=out_dir,
                    base_name=base_name_for_files
                )

                # Prepare DataFrame with tokens for word cloud generation
                sel_df_with_tokens = self._prepare_df_with_tokens(sel_df_base, plat)
                
                # Show sample table for this specific tail (optional, could be removed for batch-like behavior if too many popups)
                self.show_table_popup_only_table(sel_df_with_tokens, plat, gene_col,
                                                 f"{which_tail.capitalize()} Tail Samples (n={len(sel_df_with_tokens)} GSMs)")
                
                # Get background token data
                bg_df = self._get_background_df(plat, sel_df_with_tokens["GSM"])

                # --- Generate and save word cloud images and stats ---
                self._generate_and_save_wordcloud_images(
                    sel_df_with_tokens, bg_df, plat, gene_col, f"{which_tail}_tail",
                    out_dir, similarity_threshold, attn_available, display_mode=True
                )
                
                # --- New: Call for cluster analysis and display after word clouds ---
                if attn_available:
                    temp_sel_tokens = self.flatten_tokens(sel_df_with_tokens["tokens"])
                    temp_cnt_sel, temp_N_sel = Counter(temp_sel_tokens), len(temp_sel_tokens)
                    temp_cnt_bg, temp_N_bg = Counter(self.flatten_tokens(bg_df["tokens"])), len(self.flatten_tokens(bg_df["tokens"]))

                    # Convert counters to pandas Series for vectorized operations
                    sel_counts_series = pd.Series(temp_cnt_sel, name='sel_count')
                    bg_counts_series = pd.Series(temp_cnt_bg, name='bg_count')
                    all_tokens_combined = pd.concat([sel_counts_series, bg_counts_series], axis=1).fillna(0)
                    all_tokens_combined.columns = ['sel_count', 'bg_count']

                    total_counts = all_tokens_combined['sel_count'] + all_tokens_combined['bg_count']
                    
                    # Calculate Z-score: (count_in_sel - count_in_bg) / sqrt(count_in_sel + count_in_bg)
                    # Handle cases where total_counts is 0 to avoid ZeroDivisionError
                    z_scores = np.where(total_counts > 0,
                                         (all_tokens_combined['sel_count'] - all_tokens_combined['bg_count']) / np.sqrt(total_counts),
                                         0) # If total_counts is 0, z-score is 0 (no difference)

                    # Calculate one-tailed p-value for enrichment (sel_count > bg_count)
                    p_values = 1 - norm.cdf(z_scores)

                    # Calculate percentages
                    sel_pct = np.where(temp_N_sel > 0, all_tokens_combined['sel_count'] / temp_N_sel, 0)
                    bg_pct = np.where(temp_N_bg > 0, all_tokens_combined['bg_count'] / temp_N_bg, 0)

                    temp_full_stats_df = pd.DataFrame({
                        'token': all_tokens_combined.index,
                        'sel_count': all_tokens_combined['sel_count'].astype(int),
                        'bg_count': all_tokens_combined['bg_count'].astype(int),
                        'p_value': p_values,
                        'sel_pct': sel_pct,
                        'bg_pct': bg_pct
                    })

                    if not temp_full_stats_df.empty:
                        temp_full_stats_df = temp_full_stats_df.sort_values(by=['p_value', 'sel_count'], ascending=[True, False])
                        
                    enriched_tokens_df = temp_full_stats_df[(temp_full_stats_df['p_value'] < 0.05) & (temp_full_stats_df['sel_pct'] > temp_full_stats_df['bg_pct'])].copy()

                    bio_specific_tokens_df = pd.DataFrame()
                    if not enriched_tokens_df.empty:
                        enriched_tok_list_for_sim = enriched_tokens_df['token'].tolist()
                        sim_scores_dict = self.attn.calculate_similarity(enriched_tok_list_for_sim, similarity_threshold, self.token_cache)
                        bio_specific_tokens_df = enriched_tokens_df[enriched_tokens_df['token'].isin(sim_scores_dict.keys())].copy()
                        if not bio_specific_tokens_df.empty:
                            bio_specific_tokens_df['sim_score'] = bio_specific_tokens_df['token'].map(sim_scores_dict)
                            bio_specific_tokens_df = bio_specific_tokens_df.sort_values(by=['sim_score', 'p_value'], ascending=[False, True])

                    if not bio_specific_tokens_df.empty:
                        self.enqueue_log(f"[Clustering] Found {len(bio_specific_tokens_df)} bio-specific tokens for clustering for {gene_symbol} {which_tail.capitalize()} Tail ({plat}).")
                        self._perform_token_clustering_and_display(
                            bio_specific_tokens_df,
                            f"{gene_symbol} {which_tail.capitalize()} Tail ({plat})",
                            out_dir
                        )
                    else:
                        self.enqueue_log(f"[Clustering] No bio-specific tokens found for clustering for {gene_symbol} {which_tail.capitalize()} Tail ({plat}). Skipping cluster plots.")
                
        self.enqueue_log(f"[Tail Analysis] Completed analysis for all selected genes and platforms.")


    def _analyze_full_distribution(self):
        """
        Analyzes the full distribution of multiple genes across selected platforms.
        Results (density plots, word clouds, cluster analysis) are saved for each gene/platform.
        """
        if not self.gpl_datasets:
            messagebox.showerror("Data Missing", "No GPL datasets loaded.", parent=self.gene_dist_popup_root)
            return

        gene_input_str = self.current_gene_entry.get().strip()
        if not gene_input_str:
            messagebox.showerror("Input Missing", "Enter at least one gene symbol.", parent=self.gene_dist_popup_root)
            return
            
        genes_to_analyze = [g.strip().upper() for g in gene_input_str.split(',') if g.strip()]
        if not genes_to_analyze:
            messagebox.showerror("Input Missing", "Invalid gene symbol(s).", parent=self.gene_dist_popup_root); return

        sel_plats = [p for p,v in self.gpl_selection_vars.items() if v.get()]
        if not sel_plats:
            messagebox.showerror("Input Missing", "Select platform(s).", parent=self.gene_dist_popup_root)
            return

        similarity_threshold = self.similarity_thr_var_popup.get() / 100.0
        attn_available = self.attn is not None and self.attn.model is not None and self.attn.keywords and self.attn.keyword_embeddings.size > 0

        self.enqueue_log("[Full Distribution Analysis] Starting full distribution analysis for selected genes and platforms...")

        for gene_symbol in genes_to_analyze:
            all_selected_gene_samples_for_gene, all_selected_gsms_for_gene = [], set()

            for plat in sel_plats:
                dfg = self.gpl_datasets.get(plat)
                gmap = self.gpl_gene_mappings.get(plat, {})
                gene_col = gmap.get(gene_symbol)
                
                if dfg is None or dfg.empty or gene_col is None:
                    self.enqueue_log(f"[Full Dist] Gene '{gene_symbol}' not found or platform '{plat}' not loaded. Skipping this combination.")
                    continue

                current_plat_data = dfg.copy()
                current_plat_data["Expression"] = current_plat_data[gene_col]
                current_plat_data["Platform"] = plat

                current_plat_data_with_tokens = self._prepare_df_with_tokens(current_plat_data, plat)
                if not current_plat_data_with_tokens.empty:
                    all_selected_gene_samples_for_gene.append(current_plat_data_with_tokens)
                    all_selected_gsms_for_gene.update(current_plat_data_with_tokens["GSM"].unique())

            if not all_selected_gene_samples_for_gene:
                self.enqueue_log(f"[Full Dist] No samples found across selected platforms for gene '{gene_symbol}'. Skipping full distribution analysis for this gene.")
                continue

            combined_sel_df = pd.concat(all_selected_gene_samples_for_gene, ignore_index=True)
            
            # For full distribution analysis, background is typically all other samples
            # This logic gathers all tokens from other platforms/samples not in the current selection
            all_bg_dfs_raw = []
            for p in self.gpl_datasets.keys(): # Iterate through all currently loaded platforms
                full_plat_df = self._get_tok_df(p)
                if full_plat_df is not None:
                    # Exclude samples that are part of the current gene's "selected" set (across all its platforms)
                    bg_subset = full_plat_df[~full_plat_df["GSM"].isin(all_selected_gsms_for_gene)].copy()
                    bg_subset["tokens"] = bg_subset["tokens"].fillna("[]").apply(lambda x: json.dumps(x) if isinstance(x, list) else str(x))
                    all_bg_dfs_raw.append(bg_subset)

            combined_bg_df = pd.concat(all_bg_dfs_raw, ignore_index=True) if all_bg_dfs_raw else pd.DataFrame()


            # --- Create unique output directory for this specific full distribution analysis ---
            analysis_folder_name = f"{gene_symbol}_FULL_DISTRIBUTION"
            out_dir = os.path.join("NEW_RESULTS_ROOT", "GeneDistributionAnalysis", analysis_folder_name)
            os.makedirs(out_dir, exist_ok=True)
            base_name_for_files = f"{gene_symbol}_FULL_DISTRIBUTION"
            self.enqueue_log(f"[INFO] Analyzing full distribution for gene {gene_symbol}. Results saving to: {out_dir}")

            # --- 1. Density Plot (REORDERED to be FIRST) ---
            # For "Full Distribution", the "background" is truly 'all other samples NOT in this gene's set'.
            # This logic gathers all expression data from all other loaded genes on all loaded platforms
            all_platforms_all_expr = []
            for p, df in self.gpl_datasets.items():
                gmap = self.gpl_gene_mappings.get(p, {})
                for g_sym, g_col in gmap.items():
                    # Only include other genes
                    if g_sym != gene_symbol:
                        all_platforms_all_expr.append(df[g_col].dropna())
            
            all_bg_expr_series = pd.concat(all_platforms_all_expr, ignore_index=True) if all_platforms_all_expr else pd.Series(dtype='float64')

            self._display_density_plot_popup(
                selected_expressions=combined_sel_df["Expression"].dropna(),
                background_expressions=all_bg_expr_series,
                title=f"Density Plot: {gene_symbol} - Full Distribution",
                output_dir=out_dir,
                base_name=base_name_for_files
            )


            self.show_table_popup_only_table(combined_sel_df, "All Platforms", gene_symbol, f"Full Dist Samples (n={len(combined_sel_df)})")
            
            # Generate and save word cloud images and stats
            self._generate_and_save_wordcloud_images(
                combined_sel_df, combined_bg_df, "All Platforms", gene_symbol, "Full_Distribution",
                out_dir, similarity_threshold, attn_available, display_mode=True
            )
            
            # New: Call for cluster analysis and display after word clouds
            if attn_available:
                temp_sel_tokens = self.flatten_tokens(combined_sel_df["tokens"])
                temp_cnt_sel, temp_N_sel = Counter(temp_sel_tokens), len(temp_sel_tokens)
                temp_cnt_bg, temp_N_bg = Counter(self.flatten_tokens(combined_bg_df["tokens"])), len(self.flatten_tokens(combined_bg_df["tokens"]))

                # Convert counters to pandas Series for vectorized operations
                sel_counts_series = pd.Series(temp_cnt_sel, name='sel_count')
                bg_counts_series = pd.Series(temp_cnt_bg, name='bg_count')
                all_tokens_combined = pd.concat([sel_counts_series, bg_counts_series], axis=1).fillna(0)
                all_tokens_combined.columns = ['sel_count', 'bg_count']

                total_counts = all_tokens_combined['sel_count'] + all_tokens_combined['bg_count']
                
                # Calculate Z-score: (count_in_sel - count_in_bg) / sqrt(count_in_sel + count_in_bg)
                # Handle cases where total_counts is 0 to avoid ZeroDivisionError
                z_scores = np.where(total_counts > 0,
                                     (all_tokens_combined['sel_count'] - all_tokens_combined['bg_count']) / np.sqrt(total_counts),
                                     0) # If total_counts is 0, z-score is 0 (no difference)

                # Calculate one-tailed p-value for enrichment (sel_count > bg_count)
                p_values = 1 - norm.cdf(z_scores)

                # Calculate percentages
                sel_pct = np.where(temp_N_sel > 0, all_tokens_combined['sel_count'] / temp_N_sel, 0)
                bg_pct = np.where(temp_N_bg > 0, all_tokens_combined['bg_count'] / temp_N_bg, 0)

                temp_full_stats_df = pd.DataFrame({
                    'token': all_tokens_combined.index,
                    'sel_count': all_tokens_combined['sel_count'].astype(int),
                    'bg_count': all_tokens_combined['bg_count'].astype(int),
                    'p_value': p_values,
                    'sel_pct': sel_pct,
                    'bg_pct': bg_pct
                })

                if not temp_full_stats_df.empty:
                    temp_full_stats_df = temp_full_stats_df.sort_values(by=['p_value', 'sel_count'], ascending=[True, False])
                    
                enriched_tokens_df = temp_full_stats_df[(temp_full_stats_df['p_value'] < 0.05) & (temp_full_stats_df['sel_pct'] > temp_full_stats_df['bg_pct'])].copy()

                bio_specific_tokens_df = pd.DataFrame()
                if not enriched_tokens_df.empty:
                    enriched_tok_list_for_sim = enriched_tokens_df['token'].tolist()
                    sim_scores_dict = self.attn.calculate_similarity(enriched_tok_list_for_sim, similarity_threshold, self.token_cache)
                    bio_specific_tokens_df = enriched_tokens_df[enriched_tokens_df['token'].isin(sim_scores_dict.keys())].copy()
                    if not bio_specific_tokens_df.empty:
                        bio_specific_tokens_df['sim_score'] = bio_specific_tokens_df['token'].map(sim_scores_dict)
                        bio_specific_tokens_df = bio_specific_tokens_df.sort_values(by=['sim_score', 'p_value'], ascending=[False, True])

                if not bio_specific_tokens_df.empty:
                    self.enqueue_log(f"[Clustering] Found {len(bio_specific_tokens_df)} bio-specific tokens for clustering for {gene_symbol} Full Distribution (All Platforms).")
                    self._perform_token_clustering_and_display(
                        bio_specific_tokens_df,
                        f"{gene_symbol} Full Distribution (All Platforms)",
                        out_dir
                    )
                else:
                    self.enqueue_log(f"[Clustering] No bio-specific tokens found for clustering for {gene_symbol} Full Distribution (All Platforms). Skipping cluster plots.")

            self.enqueue_log(f"[Full Dist] Analysis for {gene_symbol} completed.")

    def _get_tok_df(self, plat: str) -> Optional[pd.DataFrame]:
        """
        Loads and caches the token DataFrame for a given platform.
        This function assumes a directory structure like:
        `/path/to/data_dir/{PLATFORM_NAME}/geometadb_tokens_{PLATFORM_NAME}.csv.gz`
        """
        if plat in self._tok_df_cache:
            return self._tok_df_cache[plat].copy() if self._tok_df_cache[plat] is not None else None

        # Corrected path construction for token files based on the specified structure
        path = os.path.join(self.data_dir, plat, f"geometadb_tokens_{plat}.csv.gz")

        self.enqueue_log(f"[DEBUG] Attempting to load token file: {path}")

        if not os.path.exists(path):
            self.enqueue_log(f"[ERROR] Token file NOT FOUND for {plat}. Expected at: {path}")
            self.enqueue_log(f"[ERROR] Word cloud generation will fail for {plat}. Please check file path and naming convention.")
            self._tok_df_cache[plat] = None
            return None

        try:
            self.enqueue_log(f"[INFO] Loading token file: {path}")
            df = pd.read_csv(path, compression="gzip", low_memory=False)
            self.enqueue_log(f"[DEBUG] Token file loaded. Columns: {list(df.columns)}")

            # Standardize GSM column name if needed
            if "gsm" in df.columns and "GSM" not in df.columns:
                df.rename(columns={"gsm":"GSM"}, inplace=True)
            elif "GSM" not in df.columns: # Added more robust check if 'gsm' isn't there either
                found_gsm_col = False
                for col_name in df.columns:
                    if "GSM" in col_name.upper() or "SAMPLE" in col_name.upper() or col_name.upper()=="ID":
                        df.rename(columns={col_name:"GSM"}, inplace=True)
                        self.enqueue_log(f"[INFO] Renamed '{col_name}' to 'GSM' in token file for {plat}.")
                        found_gsm_col = True
                        break
                if not found_gsm_col:
                    self.enqueue_log(f"[ERROR] 'GSM' column (or suitable alternative) not found in token file for {plat}. Found columns: {list(df.columns)}")
                    self._tok_df_cache[plat] = None
                    return None

            df["GSM"] = df["GSM"].astype(str).str.upper()
            
            # Ensure 'tokens' column exists, otherwise word cloud generation will fail
            if "tokens" not in df.columns:
                self.enqueue_log(f"[ERROR] 'tokens' column not found in token file for {plat}. Word cloud generation will fail.")
                self._tok_df_cache[plat] = None
                return None

            # Attempt to add series_id if not present, for consistency
            if "series_id" not in df.columns and self.gds_conn:
                ids = df["GSM"].unique().tolist()
                if ids:
                    frames = []
                    for i in range(0, len(ids), 999):
                        chunk_ids = ids[i:i+999]
                        placeholders = ",".join("?" * len(chunk_ids))
                        query = f"SELECT gsm AS GSM, series_id FROM gsm WHERE gsm IN ({placeholders})"
                        try:
                            part = pd.read_sql_query(query, self.gds_conn, params=chunk_ids)
                            frames.append(part)
                        except Exception as sqle:
                            self.enqueue_log(f"[ERROR] SQLite query for series_id failed for {plat} chunk: {sqle}")
                            break # Stop trying to fetch series_id for this plat
                    if frames:
                        map_df = pd.concat(frames, ignore_index=True)
                        map_df["GSM"] = map_df["GSM"].astype(str).str.upper()
                        df = df.merge(map_df, on="GSM", how="left")
                        self.enqueue_log(f"[INFO] Added series_id column to token data for {plat}.")
                    else:
                        self.enqueue_log(f"[WARN] No GSMs to query for series_id mapping for {plat} or query failed.")

            self.enqueue_log(f"[INFO] Successfully loaded and processed tokens for {plat}. Final DataFrame head:\n{df.head().to_string()}")
            self._tok_df_cache[plat] = df
            return df.copy()

        except Exception as e:
            self.enqueue_log(f"[CRITICAL ERROR] Failed to read or process token file {path}. Reason: {e}")
            self._tok_df_cache[plat] = None
            return None

    def _get_background_df(self, plat: str, excluded_gsms: pd.Series) -> pd.DataFrame:
        """Retrieves background token data for a given platform, excluding specified GSMs."""
        full_plat_tok_df = self._get_tok_df(plat)
        if full_plat_tok_df is not None:
            # Ensure 'tokens' column is stringified list for consistency before filtering
            bg_df = full_plat_tok_df[~full_plat_tok_df["GSM"].isin(excluded_gsms)].copy()
            # Ensure the 'tokens' column is always a string representation of a list (or empty list)
            bg_df["tokens"] = bg_df["tokens"].fillna("[]").apply(lambda x: json.dumps(x) if isinstance(x, list) else str(x))
            return bg_df
        return pd.DataFrame({"GSM": [], "tokens": []})

    def _prepare_df_with_tokens(self, df_input: pd.DataFrame, plat_name: str) -> pd.DataFrame:
        """
        Merges input DataFrame with token data from the specified platform.
        Ensures 'tokens' column is present and in string format suitable for JSON parsing.
        """
        tok_df = self._get_tok_df(plat_name)

        df_input_copy = df_input.copy()

        if "GSM" in df_input_copy.columns:
            df_input_copy["GSM"] = df_input_copy["GSM"].astype(str).str.upper()
        else:
            self.enqueue_log(f"[ERROR] Input DataFrame for _prepare_df_with_tokens is missing 'GSM' column for platform {plat_name}. Adding empty 'tokens' column.")
            df_input_copy["tokens"] = "[]" # Cannot merge without GSM, return with empty tokens
            return df_input_copy


        if tok_df is None:
            self.enqueue_log(f"[WARN] Merge skipped for {plat_name} as token data is unavailable. Adding empty 'tokens' column to output.")
            df_input_copy["tokens"] = "[]" # Add empty tokens column if token data is missing
        else:
            # Merge with tok_df, ensuring 'tokens' column exists and is string
            # Select only 'GSM' and 'tokens' from tok_df for merge to avoid column conflicts
            df_output = df_input_copy.merge(tok_df[["GSM", "tokens"]], on="GSM", how="left")
            # Ensure the 'tokens' column is always a string representation of a list (or empty list)
            df_output["tokens"] = df_output["tokens"].fillna("[]").apply(lambda x: json.dumps(x) if isinstance(x, list) else str(x))
            df_input_copy = df_output # Update the df_input_copy reference

        return df_input_copy

    def flatten_tokens(self, series) -> List[str]:
        """Flattens a Series of (potentially JSON-stringified) token lists into a single list of tokens."""
        all_tokens = []
        for item in series.dropna():
            if isinstance(item, str):
                try:
                    tokens_list = json.loads(item) # Try parsing as JSON list
                except (json.JSONDecodeError, TypeError):
                    # If JSON parsing fails, treat as a raw string and split by comma
                    tokens_list = [t.strip() for t in item.split(',') if t.strip()]
            elif isinstance(item, list):
                tokens_list = item # Already a list
            else:
                tokens_list = [] # Unexpected type

            if isinstance(tokens_list, list):
                all_tokens.extend([str(t).strip().lower() for t in tokens_list if len(str(t).strip()) > 2])
        return all_tokens


    def _generate_and_save_wordcloud_images(self, sel_df: pd.DataFrame, bg_df: pd.DataFrame,
                                             plat: str, gene: str, analysis_type: str, out_dir: str,
                                             similarity_threshold: float, attn_available: bool,
                                             display_mode: bool = False):
        """
        Generates and saves word cloud images (raw, enriched, bio-specific).
        In interactive mode, ONLY Bio-Specific word clouds are displayed (combined with stats).
        In BATCH mode, saves ONLY the bio-specific wordcloud JPG, and saves enriched
        and bio-specific token stats to separate .csv.gz files.
        """
        sel_tokens = self.flatten_tokens(sel_df["tokens"])
        bg_tokens = self.flatten_tokens(bg_df["tokens"])

        if not sel_tokens:
            self.enqueue_log(f"[WC] No tokens found in selection for {gene} {analysis_type} ({plat}). Skipping word cloud generation.")
            return

        cnt_sel, N_sel = Counter(sel_tokens), len(sel_tokens)
        cnt_bg, N_bg = Counter(bg_tokens), len(bg_tokens)

        # Pre-cache embeddings for all unique tokens if AttentionSeeker is available
        if attn_available:
            unique_toks_to_cache = set(cnt_sel.keys()) | set(cnt_bg.keys())
            # Convert to list for batch embedding
            unique_toks_to_embed_list = [tok for tok in unique_toks_to_cache if tok not in self.token_cache]
            if unique_toks_to_embed_list:
                new_embeddings = self.attn._embed(unique_toks_to_embed_list)
                if new_embeddings is None: new_embeddings = [None] * len(unique_toks_to_embed_list)
                for i, tok in enumerate(unique_toks_to_embed_list):
                    if new_embeddings[i] is not None and not np.all(new_embeddings[i] == 0):
                        self.token_cache[tok] = new_embeddings[i]
                    else:
                        self.token_cache[tok] = np.zeros(self.attn.keyword_embeddings.shape[1] if self.attn.keyword_embeddings.size > 0 else 768)


        # --- Generate data for full statistics table (VECTORIZED) ---
        sel_counts_series = pd.Series(cnt_sel, name='sel_count')
        bg_counts_series = pd.Series(cnt_bg, name='bg_count')

        all_tokens_combined = pd.concat([sel_counts_series, bg_counts_series], axis=1).fillna(0)
        all_tokens_combined.columns = ['sel_count', 'bg_count']

        total_counts = all_tokens_combined['sel_count'] + all_tokens_combined['bg_count']
        
        # Calculate Z-score: (count_in_sel - count_in_bg) / sqrt(count_in_sel + count_in_bg)
        # Handle cases where total_counts is 0 to avoid ZeroDivisionError
        z_scores = np.where(total_counts > 0,
                             (all_tokens_combined['sel_count'] - all_tokens_combined['bg_count']) / np.sqrt(total_counts),
                             0) # If total_counts is 0, z-score is 0 (no difference)

        # Calculate one-tailed p-value for enrichment (sel_count > bg_count)
        p_values = 1 - norm.cdf(z_scores)

        # Calculate percentages
        sel_pct = np.where(N_sel > 0, all_tokens_combined['sel_count'] / N_sel, 0)
        bg_pct = np.where(N_bg > 0, all_tokens_combined['bg_count'] / N_bg, 0)

        full_stats_df = pd.DataFrame({
            'token': all_tokens_combined.index,
            'sel_count': all_tokens_combined['sel_count'].astype(int),
            'bg_count': all_tokens_combined['bg_count'].astype(int),
            'p_value': p_values,
            'sel_pct': sel_pct,
            'bg_pct': bg_pct
        })
        
        if not full_stats_df.empty:
            full_stats_df = full_stats_df.sort_values(by=['p_value', 'sel_count'], ascending=[True, False])

        # --- SAVE ALL P-VALUES FILE (NEW FEATURE) ---
        if not full_stats_df.empty:
            # Sanitize gene name for the filename
            safe_gene_name = re.sub(r'[^\w\-]+', '_', gene)
            pvalues_filename = f"all_p_values_{safe_gene_name}_{plat}_{analysis_type}.csv.gz"
            pvalues_filepath = os.path.join(out_dir, pvalues_filename)
            try:
                full_stats_df.to_csv(pvalues_filepath, index=False, compression='gzip')
            except Exception as e:
                self.enqueue_log(f"[ERROR] Failed to save all-token p-values file: {e}")
        # ---------------------------------------------

        # --- Process and save CSVs (only in batch mode or when display_mode=False) ---
        os.makedirs(out_dir, exist_ok=True) # Ensure output directory exists

        # Enriched Tokens data for CSV
        self.enqueue_log(f"[WC] Generating enriched tokens for {gene} {analysis_type} ({plat})...")
        enriched_data_csv = full_stats_df[(full_stats_df['p_value'] < 0.05) & (full_stats_df['sel_pct'] > full_stats_df['bg_pct'])].copy()
        if not enriched_data_csv.empty:
            enriched_csv_path = os.path.join(out_dir, f"enriched_token_stats_{analysis_type}.csv.gz")
            enriched_data_csv.to_csv(enriched_csv_path, index=False, compression='gzip')
            self.enqueue_log(f"[WC] Enriched token statistics saved to {enriched_csv_path}")
        else:
            self.enqueue_log(f"[WC] No enriched token statistics to save for {gene} {analysis_type} ({plat}).")
            
        # Bio-Specific Tokens data for CSV and BERTopic
        self.enqueue_log(f"[WC] Generating bio-specific tokens for {gene} {analysis_type} ({plat})...")
        bio_specific_data_csv = pd.DataFrame() # Initialize
        bertopic_topics_info = [] # To store BERTopic results

        if attn_available and not enriched_data_csv.empty:
            enriched_tok_list_for_sim = enriched_data_csv['token'].tolist()
            sim_scores_dict = self.attn.calculate_similarity(enriched_tok_list_for_sim, similarity_threshold, self.token_cache)
            bio_specific_data_csv = enriched_tokens_df[enriched_tokens_df['token'].isin(sim_scores_dict.keys())].copy()
            if not bio_specific_data_csv.empty:
                bio_specific_data_csv['sim_score'] = bio_specific_data_csv['token'].map(sim_scores_dict)
                bio_specific_data_csv = bio_specific_data_csv.sort_values(by=['sim_score', 'p_value'], ascending=[False, True])
                
                # BERTopic Integration: Only if more than 3 bio-specific tokens
                if BERTOPIC_AVAILABLE and self.bertopic_model and len(bio_specific_data_csv) > 3:
                    self.enqueue_log(f"[BERTopic] Running topic modeling for {gene} {analysis_type} ({plat})...")
                    docs_for_bertopic = bio_specific_data_csv['token'].tolist()
                    try:
                        topics, probs = self.bertopic_model.fit_transform(docs_for_bertopic)
                        # Get topic information
                        freq = self.bertopic_model.get_topic_info()
                        # Filter out -1 topic (outliers)
                        freq = freq[freq.Topic != -1]
                        if not freq.empty:
                            bertopic_topics_info = []
                            for _, row in freq.head(5).iterrows(): # Get top 5 topics
                                topic_words = [word for word, _ in self.bertopic_model.get_topic(row.Topic)[:5]] # Top 5 words per topic
                                bertopic_topics_info.append(f"Topic {row.Topic} ({row.Count} tokens): {', '.join(topic_words)}")
                            self.enqueue_log(f"[BERTopic] Found {len(freq)} topics.")
                            # Add topic ID to bio_specific_data_csv
                            topic_mapping = {docs_for_bertopic[i]: topics[i] for i in range(len(docs_for_bertopic))}
                            bio_specific_data_csv['bertopic_topic'] = bio_specific_data_csv['token'].map(topic_mapping).fillna(-1).astype(int)
                        else:
                            self.enqueue_log("[BERTopic] No significant topics found.")
                    except Exception as e:
                        self.enqueue_log(f"[BERTopic ERROR] Failed to perform topic modeling: {e}")
                        bertopic_topics_info = [f"BERTopic failed: {e}"]
                elif BERTOPIC_AVAILABLE and self.bertopic_model and len(bio_specific_data_csv) <= 3:
                    self.enqueue_log(f"[BERTopic] Skipping topic modeling for {gene} {analysis_type} ({plat}): Less than 4 bio-specific tokens.")
                elif not BERTOPIC_AVAILABLE:
                    self.enqueue_log(f"[BERTopic] BERTopic library not available; skipping topic modeling for {gene} {analysis_type} ({plat}).")
                elif not self.bertopic_model:
                    self.enqueue_log(f"[BERTopic] BERTopic model not initialized; skipping topic modeling for {gene} {analysis_type} ({plat}).")


                biospecific_csv_path = os.path.join(out_dir, f"biospecific_token_stats_{analysis_type}.csv.gz")
                bio_specific_data_csv.to_csv(biospecific_csv_path, index=False, compression='gzip')
                self.enqueue_log(f"[WC] Bio-specific token statistics saved to {biospecific_csv_path}")
            else:
                self.enqueue_log(f"[WC] No bio-specific token statistics to save for {gene} {analysis_type} ({plat}).")
        elif not attn_available:
            self.enqueue_log(f"[WC] AttentionSeeker not available; skipping bio-specific token stats for {gene} {analysis_type} ({plat}).")
        else:
            self.enqueue_log(f"[WC] No enriched tokens available for bio-specific analysis for {gene} {analysis_type} ({plat}).")

        # --- Generate and save Word Cloud Images (interactive/batch) ---

        # Bio-Specific Tokens (ALWAYS generate/save; ONLY display interactively if display_mode=True)
        wc_bio_specific = None
        bio_specific_data_for_wc = pd.DataFrame() # Initialize for WC generation

        if attn_available and not full_stats_df.empty:
            self.enqueue_log(f"[WC] Attempting to generate Bio-Specific word cloud for {gene} {analysis_type} ({plat})...")
            enriched_from_full_stats = full_stats_df[(full_stats_df['p_value'] < 0.05) & (full_stats_df['sel_pct'] > full_stats_df['bg_pct'])].copy()
            
            if not enriched_from_full_stats.empty:
                enriched_tok_list = enriched_from_full_stats['token'].tolist()
                sim_scores_dict = self.attn.calculate_similarity(enriched_tok_list, similarity_threshold, self.token_cache)
                
                bio_specific_data_for_wc = enriched_from_full_stats[enriched_from_full_stats['token'].isin(sim_scores_dict.keys())].copy()
                
                if not bio_specific_data_for_wc.empty:
                    bio_specific_data_for_wc['sim_score'] = bio_specific_data_for_wc['token'].map(sim_scores_dict)
                    bio_specific_data_for_wc = bio_specific_data_for_wc.sort_values(by=['sim_score', 'p_value'], ascending=[False, True])
                    
                    bio_freq = {row['token']: row['sel_count'] for _, row in bio_specific_data_for_wc.iterrows()}
                    
                    if bio_freq:
                        wc_bio_specific = WordCloud(width=800, height=400, background_color="white", colormap="magma", random_state=42).generate_from_frequencies(bio_freq)
                        
                        # Always save bio-specific wordcloud image (for both interactive & batch)
                        save_path_bio = os.path.join(out_dir, f"wordcloud_biospecific_tokens_sim{int(similarity_threshold*100)}.jpg")
                        try:
                            wc_bio_specific.to_file(save_path_bio)
                            self.enqueue_log(f"[WC] Bio-Specific word cloud saved to {save_path_bio}")
                        except Exception as e:
                            self.enqueue_log(f"[ERROR] Failed to save bio-specific word cloud {save_path_bio}: {e}")
                            
                        if display_mode:
                            # Display combined word cloud image and its statistics table (plain text summary)
                            self._show_combined_wordcloud_and_stats_popup(
                                wc_bio_specific, # The wordcloud object
                                bio_specific_data_for_wc, # The DataFrame for statistics table
                                f"{gene} - {analysis_type} ({plat}) - Bio-Specific (≥{int(similarity_threshold*100)}%)",
                                bertopic_topics_info # Pass BERTopic topics
                            )
                    else:
                        self.enqueue_log(f"[WC] No bio-specific tokens with frequencies to generate WC for {gene} {analysis_type} ({plat}).")
                else:
                    self.enqueue_log(f"[WC] No bio-specific tokens found after semantic filtering for {gene} {analysis_type} ({plat}).")
            else:
                self.enqueue_log(f"[WC] No enriched tokens found to filter for bio-specific word cloud for {gene} {analysis_type} ({plat}).")
        elif not attn_available:
            self.enqueue_log(f"[WC] AttentionSeeker not available; skipping bio-specific word cloud for {gene} {analysis_type} ({plat}).")
        else:
            self.enqueue_log(f"[WC] No enriched tokens available for bio-specific analysis for {gene} {analysis_type} ({plat}).")

    def _show_wordcloud_image_only(self, wc_obj, title: str):
        """
        Internal helper to display a WordCloud image in a Tkinter Toplevel window.
        This is separated from the statistics table display.
        """
        win = tk.Toplevel(self)
        win.title(title)
        win.geometry("820x500") # Adjust size for image only

        fig, ax = plt.subplots(figsize=(8.2, 4.8)) # Match window geometry aspect ratio
        ax.imshow(wc_obj, interpolation="bilinear")
        ax.axis("off")
        fig.tight_layout(pad=0) # Remove padding

        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, win)
        toolbar.update(); toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.draw()
        
        # Store a reference to avoid garbage collection
        # Using a more unique key than just title in case multiple same-titled images are opened
        key = f"WC_Image_{title}_{uuid.uuid4().hex}"
        self._current_popup_figs[key] = (fig, canvas_widget, toolbar)
        
        # Ensure matplotlib figure is closed when the Tkinter window is closed
        win.protocol("WM_DELETE_WINDOW", lambda: (plt.close(fig), win.destroy(), self._current_popup_figs.pop(key, None)))
        win.transient(self)
        # Removed grab_set() here to allow multiple windows to be open
        win.lift()
        win.attributes('-topmost', True)
        win.after(50, lambda: win.attributes('-topmost', False))


    def _show_wordcloud_statistics_table_popup(self, title: str, stats_df: pd.DataFrame):
        """
        Creates a new Toplevel window to display token statistics in a ttk.Treeview.
        This is specifically for the word cloud related statistics.
        """
        if stats_df.empty:
            self.enqueue_log(f"[WC Stats Table] No statistics data for '{title}'. Skipping table display.")
            return

        win = tk.Toplevel(self)
        win.title(title)
        win.geometry("850x400") # Adjusted size for statistics table

        frame = ttk.Frame(win, padding="5")
        frame.pack(fill=BOTH, expand=True)

        # Define columns for statistics table
        cols_to_show = ["token", "sel_count", "bg_count", "p_value"]
        headers = ["Token", "Sel Count", "Bg Count", "P-Value"]

        if 'sim_score' in stats_df.columns:
            cols_to_show.append("sim_score")
            headers.append("Sim Score")
        if 'bertopic_topic' in stats_df.columns:
            cols_to_show.append("bertopic_topic")
            headers.append("BERTopic Topic")


        tree = ttk.Treeview(frame, columns=cols_to_show, show="headings")

        for i, col_name in enumerate(cols_to_show):
            tree.heading(col_name, text=headers[i]) # Use descriptive headers
            if col_name == "token":
                tree.column(col_name, width=150, minwidth=100, anchor=W)
            elif col_name in ["sel_count", "bg_count"]:
                tree.column(col_name, width=80, minwidth=50, anchor=CENTER)
            elif col_name in ["p_value", "sim_score"]:
                tree.column(col_name, width=90, minwidth=70, anchor=CENTER)
            elif col_name == "bertopic_topic":
                tree.column(col_name, width=100, minwidth=70, anchor=CENTER)
            else:
                tree.column(col_name, width=80, minwidth=50, anchor=W)

        # Populate the treeview
        for _, row in stats_df.iterrows():
            row_values = []
            for col_name in cols_to_show:
                val = row[col_name]
                if col_name == "p_value":
                    val = f"{val:.2e}"
                elif col_name == "sim_score":
                    val = f"{val:.3f}"
                elif isinstance(val, (int, float)):
                    val = f"{val:.0f}" # For counts, display as integer
                row_values.append(val)
            tree.insert("", "end", values=row_values)
            
        vsb = ttk.Scrollbar(frame, orient="vertical", command=tree.yview)
        hsb = ttk.Scrollbar(frame, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)

        vsb.pack(side=RIGHT, fill=Y)
        hsb.pack(side=BOTTOM, fill=X)
        tree.pack(side=LEFT, fill=BOTH, expand=True)

        close_button = ttk.Button(win, text="Close", command=win.destroy)
        close_button.pack(pady=5)
        
        win.transient(self)
        win.update_idletasks()
        # Removed grab_set() here
        win.lift()
        win.attributes('-topmost', True)
        win.after(50, lambda: win.attributes('-topmost', False))


    def _show_combined_wordcloud_and_stats_popup(self, wc_obj: WordCloud, stats_df: pd.DataFrame, title: str, bertopic_topics: List[str] = None):
        """
        Creates a new Toplevel window to display a WordCloud image and its statistics summary
        in plain text format within a single Matplotlib figure.
        This is used specifically for the Bio-Specific word clouds in interactive mode.
        """
        import matplotlib.pyplot as plt_local
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

        win = tk.Toplevel(self)
        win.title(title)
        win.geometry("900x800") # Adjusted size to accommodate both word cloud and text summary

        # Create a figure with two subplots: one for word cloud, one for text summary
        fig, axs = plt_local.subplots(2, 1, figsize=(9, 7.5),
                                     gridspec_kw={'height_ratios': [0.65, 0.35]}) # Adjust ratios as needed
        
        # --- Word Cloud Subplot (axs[0]) ---
        if wc_obj and wc_obj.words_:
            axs[0].imshow(wc_obj, interpolation="bilinear")
        else:
            axs[0].text(0.5, 0.5, "No Word Cloud Data", horizontalalignment='center', verticalalignment='center', transform=axs[0].transAxes, fontsize=12, color='gray')
        axs[0].set_title("Bio-Specific Word Cloud", fontsize=12)
        axs[0].axis("off")

        # --- Statistics Text Summary Subplot (axs[1]) ---
        axs[1].axis('off') # Hide axes for this subplot

        summary_text = ""
        if not stats_df.empty:
            # Prepare text data for display (top 15 tokens)
            display_df = stats_df.head(15).copy()
            
            # Create a formatted string for display
            summary_text += "Top Token Statistics:\n"
            summary_text += "{:<10} {:>5} {:>5} {:>8} {:>5}\n".format(
                "Token", "Sel", "Bg", "P-val", "Sim") # Removed Sel% and Bg%
            summary_text += "-" * 40 + "\n" # Adjust separator length based on header

            for _, row in display_df.iterrows():
                token = row['token']
                if len(token) > 9: # Truncate token for display
                    token = token[:7] + '...'
                
                p_val_str = f"{row['p_value']:.1e}"
                sim_score_str = f"{row['sim_score']:.2f}" if 'sim_score' in row else "N/A"

                summary_text += "{:<10} {:>5} {:>5} {:>8} {:>5}\n".format( # Removed Sel% and Bg%
                    token, int(row['sel_count']), int(row['bg_count']), p_val_str, sim_score_str)
        else:
            summary_text += "No Token Statistics Available.\n"

        # Add BERTopic topics if available
        if bertopic_topics:
            summary_text += "\nBERTopic Topics:\n"
            for topic_str in bertopic_topics:
                summary_text += f"- {topic_str}\n"
        elif self.bertopic_model and len(stats_df) > 3: # If BERTopic was supposed to run but found no topics
              summary_text += "\nBERTopic Topics: No significant topics found.\n"
        elif self.bertopic_model and len(stats_df) <= 3 and len(stats_df) > 0:
            summary_text += "\nBERTopic Topics: Skipped (less than 4 bio-specific tokens).\n"
        elif not self.bertopic_model and BERTOPIC_AVAILABLE:
            summary_text += "\nBERTopic Topics: Model not loaded or initialized.\n"
        elif not BERTOPIC_AVAILABLE:
            summary_text += "\nBERTopic Topics: BERTopic library not available.\n"


        axs[1].text(0.01, 0.99, summary_text, 
                                            transform=axs[1].transAxes, 
                                            fontsize=8, verticalalignment='top', 
                                            fontfamily='monospace', # Use monospace for better alignment
                                            wrap=True)

        fig.tight_layout(rect=[0, 0, 1, 0.95]) # Adjust layout to prevent title overlapping
        fig.suptitle(title, fontsize=14, y=0.98) # Main title for the entire figure
        
        # Embed the Matplotlib figure into the Tkinter window
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas_widget = canvas.get_tk_widget(); canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(canvas, win); toolbar.update(); toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.draw()
        
        # Store a reference to avoid garbage collection
        key = f"Combined_WC_Stats_{title}_{uuid.uuid4().hex}"
        self._current_popup_figs[key] = (fig, canvas_widget, toolbar)
        
        # Ensure matplotlib figure is closed when the Tkinter window is closed
        win.protocol("WM_DELETE_WINDOW", lambda: (plt_local.close(fig), win.destroy(), self._current_popup_figs.pop(key, None)))
        win.transient(self)
        win.lift()
        win.attributes('-topmost', True)
        win.after(50, lambda: win.attributes('-topmost', False))


    def _batch_process_all_genes_for_tails(self):
        """Starts batch processing of all genes for their tails across selected platforms."""
        sel_plats = [p for p,v in self.gpl_selection_vars.items() if v.get()]
        if not sel_plats:
            messagebox.showerror("Input Missing","Select platform(s).", parent=self.gene_dist_popup_root)
            return

        sim_thr = self.similarity_thr_var_popup.get() / 100.0

        if not self.attn or not self.attn.model or not self.attn.keywords:
            messagebox.showwarning("AttentionSeeker Not Ready", "Bio-aware filtering (AttentionSeeker) is not fully loaded or initialized. Bio-specific word clouds will be skipped during batch processing.", parent=self.gene_dist_popup_root)
            attn_available = False
        else:
            attn_available = True

        self.enqueue_log("[BATCH] Starting batch processing for all genes (tails) across selected platforms...")
        self.progressbar["value"] = 0
        total_genes = sum(len(self.gpl_gene_mappings.get(p,{})) for p in sel_plats)
        if total_genes == 0:
            self.enqueue_log("[BATCH] No genes to process."); return

        # Start the batch processing in a separate thread
        threading.Thread(target=self._batch_job, args=(sel_plats, total_genes, sim_thr, attn_available), daemon=True).start()

    def _batch_job(self, sel_plats: List[str], total_genes: int, sim_thr: float, attn_available: bool):
        """
        The main logic for batch processing gene tails, running in a separate thread.
        Generates and saves word clouds and statistics for each gene's tails.
        """
        processed_genes = 0
        for plat in sel_plats:
            if self.is_closing: break # Allow graceful exit if app is closing
            dfg = self.gpl_datasets.get(plat)
            gmap = self.gpl_gene_mappings.get(plat, {})
            if dfg is None or dfg.empty:
                self.enqueue_log(f"[BATCH INFO] Platform {plat} not loaded or empty. Skipping.")
                processed_genes += len(gmap) # Account for genes that would have been processed
                self.after(0, lambda p=processed_genes, tg=total_genes: self.progressbar.config(value=(p/tg)*100))
                continue

            for gene_upper, gene_col in gmap.items():
                if self.is_closing: break # Allow graceful exit if app is closing

                expr = dfg[gene_col].dropna().astype(float)
                # Skip if no data or constant values, which prevents std() from being 0
                if expr.empty or expr.std() == 0:
                    self.enqueue_log(f"[BATCH INFO] Gene {gene_upper} in {plat} has no variance or data. Skipping.")
                    processed_genes += 1
                    self.after(0, lambda p=processed_genes, tg=total_genes: self.progressbar.config(value=(p/tg)*100))
                    continue

                mu, sd = expr.mean(), expr.std()
                # Define tails as outside 3 standard deviations from the mean
                low_thr, high_thr = mu - 3*sd, mu + 3*sd

                # Ensure thresholds are within data range to avoid empty selections if data is very skewed
                min_expr, max_expr = expr.min(), expr.max()
                low_thr_adj = max(low_thr, min_expr)
                high_thr_adj = min(high_thr, max_expr)

                processed_any_tail_for_gene = False

                # Process LEFT Tail
                # Ensure at least one sample is outside the 'normal' range
                tail_mask_left = dfg[gene_col] <= low_thr_adj
                sel_df_base_left = dfg.loc[tail_mask_left].copy()
                if not sel_df_base_left.empty:
                    self.enqueue_log(f"[BATCH] Processing {plat}, gene {gene_upper}: LEFT_TAIL ({len(sel_df_base_left)} samples)")
                    sel_df_base_left["Expression"] = sel_df_base_left[gene_col]
                    sel_df_base_left["Platform"] = plat
                    sel_df_with_tokens_left = self._prepare_df_with_tokens(sel_df_base_left, plat)

                    # Check if the prepared DataFrame has actual tokens.
                    if not sel_df_with_tokens_left.empty and "tokens" in sel_df_with_tokens_left.columns and \
                       any(len(json.loads(str(x))) > 0 for x in sel_df_with_tokens_left["tokens"].dropna() if isinstance(x, str)):
                        bg_df = self._get_background_df(plat, sel_df_with_tokens_left["GSM"])
                        out_dir = os.path.join("NEW_RESULTS_ROOT", "BatchGeneAnalysis", f"{gene_upper}_{plat}_LEFT_TAIL") # FIX: Changed RESULTS_ROOT
                        os.makedirs(out_dir, exist_ok=True)
                        # Call the word cloud image generation function (not the combined one)
                        self._generate_and_save_wordcloud_images(sel_df_with_tokens_left, bg_df, plat,
                                                                 gene_upper, "LEFT_TAIL", out_dir, sim_thr, attn_available,
                                                                 display_mode=False) # Save as JPG for batch
                        processed_any_tail_for_gene = True
                    else:
                        self.enqueue_log(f"[BATCH WARN] No valid samples with tokens for {plat}, gene {gene_upper} LEFT_TAIL. Skipping word cloud generation.")


                # Process RIGHT Tail
                tail_mask_right = dfg[gene_col] >= high_thr_adj
                sel_df_base_right = dfg.loc[tail_mask_right].copy()
                if not sel_df_base_right.empty:
                    self.enqueue_log(f"[BATCH] Processing {plat}, gene {gene_upper}: RIGHT_TAIL ({len(sel_df_base_right)} samples)")
                    sel_df_base_right["Expression"] = sel_df_base_right[gene_col]
                    sel_df_base_right["Platform"] = plat
                    sel_df_with_tokens_right = self._prepare_df_with_tokens(sel_df_base_right, plat)

                    # Check if the prepared DataFrame has actual tokens.
                    if not sel_df_with_tokens_right.empty and "tokens" in sel_df_with_tokens_right.columns and \
                       any(len(json.loads(str(x))) > 0 for x in sel_df_with_tokens_right["tokens"].dropna() if isinstance(x, str)):
                        bg_df = self._get_background_df(plat, sel_df_with_tokens_right["GSM"])
                        out_dir = os.path.join("NEW_RESULTS_ROOT", "BatchGeneAnalysis", f"{gene_upper}_{plat}_RIGHT_TAIL") # FIX: Changed RESULTS_ROOT
                        os.makedirs(out_dir, exist_ok=True)
                        # Call the word cloud image generation function (not the combined one)
                        self._generate_and_save_wordcloud_images(sel_df_with_tokens_right, bg_df, plat,
                                                                 gene_upper, "RIGHT_TAIL", out_dir, sim_thr, attn_available,
                                                                 display_mode=False) # Save as JPG for batch
                        processed_any_tail_for_gene = True
                    else:
                        self.enqueue_log(f"[BATCH WARN] No valid samples with tokens for {plat}, gene {gene_upper} RIGHT_TAIL. Skipping word cloud generation.")

                if not processed_any_tail_for_gene:
                    self.enqueue_log(f"[BATCH INFO] No samples found in significant tails for {gene_upper} in {plat}.")

                processed_genes += 1
                self.after(0, lambda p=processed_genes, tg=total_genes: self.progressbar.config(value=(p/tg)*100))

        if not self.is_closing:
            self.enqueue_log("[BATCH] Batch processing finished.")
            self.after(0, lambda: self.progressbar.config(value=0)) # Reset progress bar

    def _batch_negative_control(self):
        from tkinter import messagebox as mb
        sel_plats = [p for p,v in self.gpl_selection_vars.items() if v.get()]
        if not sel_plats:
            mb.showerror("Input Missing", "Select at least one platform.", parent=self.gene_dist_popup_root)
            return
        sim_thr = self.similarity_thr_var_popup.get() / 100.0
        self.enqueue_log("[NEG_CTRL] Starting negative-control batch with shuffled tokens...")
        # Run in background
        threading.Thread(target=self._batch_neg_ctrl_job, args=(sel_plats, sim_thr), daemon=True).start()

    def _batch_neg_ctrl_job(self, sel_plats: List[str], sim_thr: float):
        """
        Worker thread: for each gene/platform, shuffle tokens across samples,
        then apply low/high tail logic and save stats+wordclouds under
        NEW_RESULTS_ROOT/BatchGeneAnalysis/NEGATIVE_CONTROL_SHUFFLED_TOKENS.
        """
        total_genes = sum(len(self.gpl_gene_mappings.get(p, {})) for p in sel_plats)
        processed = 0
        for plat in sel_plats:
            # Load expression and token data
            expr_df = self.gpl_datasets.get(plat)
            tok_df = self._get_tok_df(plat)
            gmap = self.gpl_gene_mappings.get(plat, {})
            if expr_df is None or tok_df is None:
                processed += len(gmap)
                continue

            # Shuffle only the 'tokens' column
            shuffled_tokens = tok_df['tokens'].sample(frac=1, random_state=42).reset_index(drop=True)
            tok_shuf_df = tok_df.copy()
            tok_shuf_df['tokens'] = shuffled_tokens.values

            for gene_sym, gene_col in gmap.items():
                vals = expr_df[gene_col].dropna().astype(float)
                if len(vals) < 2 or vals.std() == 0:
                    processed += 1
                    continue
                mu, sd = vals.mean(), vals.std()
                low_thr = max(vals.min(), mu - 3*sd)
                high_thr = min(vals.max(), mu + 3*sd)

                for tail_name, mask in [('LEFT', vals <= low_thr), ('RIGHT', vals >= high_thr)]:
                    sel_samples = expr_df.loc[mask.index[mask]].copy()
                    if sel_samples.empty:
                        continue
                    sel_samples['Expression'] = sel_samples[gene_col]
                    sel_samples['Platform'] = plat
                    # Merge shuffled tokens by GSM
                    merged = sel_samples[['GSM']].merge(
                        tok_shuf_df[['GSM','tokens']], on='GSM', how='left'
                    )
                    merged['Expression'] = sel_samples['Expression']
                    merged['Platform'] = plat

                    out_folder = os.path.join(
                        'NEW_RESULTS_ROOT','BatchGeneAnalysis',
                        'NEGATIVE_CONTROL_SHUFFLED_TOKENS',
                        f"{gene_sym}_{plat}_{tail_name}_TAIL"
                    )
                    os.makedirs(out_folder, exist_ok=True)
                    # Background unchanged
                    bg = self._get_background_df(plat, merged['GSM'])
                    self._generate_and_save_wordcloud_images(
                        merged, bg, plat, gene_sym,
                        f"{tail_name}_TAIL", out_folder,
                        sim_thr, attn_available=bool(self.attn),
                        display_mode=False
                    )
                processed += 1
                # update progress bar
                self.after(0, lambda p=processed, t=total_genes: self.progressbar.config(value=(p/t)*100))

        self.enqueue_log("[NEG_CTRL] Negative-control batch complete.")
        self.after(0, lambda: self.progressbar.config(value=0))

    def show_table_popup_only_table(self, df_show: pd.DataFrame, plat_name: str, gene_name: str, suffix_title: str):
        """
        Creates a new Toplevel window to display a DataFrame in a ttk.Treeview.
        This function is updated to remove the 'tokens' column and provide resizable
        panes for the sample table and the GSE distribution table.
        """
        if df_show.empty:
            self.enqueue_log(f"[INFO] No data to display in table for {gene_name} ({plat_name}).")
            return

        win_table = tk.Toplevel(self)
        win_table.title(f"{gene_name} {suffix_title} — {plat_name}")
        win_table.geometry("950x700")
        
        tk.Label(win_table, text=f"Total samples in selection: {len(df_show)}").pack(pady=4)

        # --- Data Preparation (GSE and Classification Merging) ---
        gsms_to_lookup = []
        if "GSM" in df_show.columns:
            gsms_to_lookup = df_show["GSM"].astype(str).unique().tolist()
        
        mapping_df = pd.DataFrame()
        if gsms_to_lookup and self.gds_conn:
            gsm_chunks = [gsms_to_lookup[i:i+999] for i in range(0, len(gsms_to_lookup), 999)]
            all_gse_mappings = []
            for chunk_ids in gsm_chunks:
                placeholders = ",".join("?" * len(chunk_ids))
                query = f"SELECT gsm AS GSM, series_id AS GSE FROM gsm WHERE gsm IN ({placeholders})"
                try:
                    part_map_df = pd.read_sql_query(query, self.gds_conn, params=chunk_ids)
                    all_gse_mappings.append(part_map_df)
                except Exception as e:
                    self.enqueue_log(f"[WARN] Failed to get GSE mapping from DB for a chunk: {e}")
            
            if all_gse_mappings:
                mapping_df = pd.concat(all_gse_mappings, ignore_index=True)
                if "GSM" in mapping_df.columns:
                    mapping_df["GSM"] = mapping_df["GSM"].astype(str).str.upper()
            
        if not mapping_df.empty and "GSM" in df_show.columns and "GSM" in mapping_df.columns:
            df_show["GSM_upper_for_merge"] = df_show["GSM"].astype(str).str.upper()
            mapping_df["GSM_upper_for_merge"] = mapping_df["GSM"].astype(str).str.upper()
            existing_gse_col = "GSE" if "GSE" in df_show.columns else None
            df_show = df_show.merge(mapping_df[["GSE", "GSM_upper_for_merge"]], on="GSM_upper_for_merge", how="left", suffixes=('', '_map'))
            
            if 'GSE_map' in df_show.columns:
                if existing_gse_col:
                    df_show[existing_gse_col] = df_show[existing_gse_col].fillna(df_show['GSE_map'])
                else:
                    df_show.rename(columns={'GSE_map': 'GSE'}, inplace=True)
                df_show.drop(columns=['GSE_map'], inplace=True, errors='ignore')
            
            df_show.drop(columns=["GSM_upper_for_merge"], inplace=True, errors='ignore')

        cls_cols_to_add = []
        try:
            attention_file_path = os.path.join(self.data_dir, "bio_token_outputs", "Attention_results", f"{plat_name}_attention.csv")
            if os.path.exists(attention_file_path):
                att_df = pd.read_csv(attention_file_path)
                if "GSM" in att_df.columns:
                    att_df["GSM_upper_for_merge"] = att_df["GSM"].astype(str).str.upper()
                    defined_cls_cols = ["Tissue_Classification", "Disease_Classification", "Condition_Classification"]
                    cls_cols_to_add = [c for c in defined_cls_cols if c in att_df.columns]
                    if cls_cols_to_add and "GSM" in df_show.columns:
                        df_show["GSM_upper_for_merge"] = df_show["GSM"].astype(str).str.upper()
                        df_show = df_show.merge(att_df[["GSM_upper_for_merge"] + cls_cols_to_add],
                                                 on="GSM_upper_for_merge", how="left", suffixes=('', '_att'))
                        df_show.drop(columns=["GSM_upper_for_merge"], errors='ignore', inplace=True)
        except Exception as e:
            self.enqueue_log(f"[WARN] Could not load or merge attention classification data: {e}")

        # --- Define columns for the Treeview display (MODIFIED: 'tokens' removed) ---
        cols_to_display = ["GSM"]
        if "GSE" in df_show.columns: cols_to_display.append("GSE")
        if "Platform" in df_show.columns: cols_to_display.append("Platform")
        if "Expression" in df_show.columns: cols_to_display.append("Expression")
        cols_to_display.extend(cls_cols_to_add)
        
        final_cols_to_display = [col for col in cols_to_display if col in df_show.columns]

        # --- LAYOUT MODIFICATION: Use PanedWindow for resizable sections ---
        main_pane = ttk.PanedWindow(win_table, orient=tk.VERTICAL)
        main_pane.pack(fill=BOTH, expand=True, padx=5, pady=5)

        # --- Top Pane: Samples Table ---
        tv_frame = ttk.Frame(main_pane, padding=5)
        main_pane.add(tv_frame, weight=2) # Give it a higher initial weight

        tv = ttk.Treeview(tv_frame, columns=final_cols_to_display, show="headings", selectmode="browse")
        
        for c_name in final_cols_to_display:
            tv.heading(c_name, text=c_name)
            if c_name in ["GSM", "GSE"]:
                tv.column(c_name, width=100, minwidth=80, anchor=W, stretch=tk.NO)
            elif c_name == "Expression":
                tv.column(c_name, width=80, minwidth=60, anchor=CENTER, stretch=tk.NO)
            elif c_name in ["Tissue_Classification", "Disease_Classification", "Condition_Classification"]:
                tv.column(c_name, width=200, minwidth=150, anchor=W, stretch=tk.YES)
            else:
                tv.column(c_name, width=80, minwidth=60, anchor=W, stretch=tk.NO)
        
        vsb_tv = ttk.Scrollbar(tv_frame, orient=tk.VERTICAL, command=tv.yview)
        hsb_tv = ttk.Scrollbar(tv_frame, orient=tk.HORIZONTAL, command=tv.xview)
        tv.configure(yscrollcommand=vsb_tv.set, xscrollcommand=hsb_tv.set)
        
        tv.grid(row=0, column=0, sticky="nsew")
        vsb_tv.grid(row=0, column=1, sticky="ns")
        hsb_tv.grid(row=1, column=0, sticky="ew")
        
        tv_frame.grid_rowconfigure(0, weight=1)
        tv_frame.grid_columnconfigure(0, weight=1)

        for _, row_data in df_show.iterrows():
            display_values = [row_data.get(c, "") for c in final_cols_to_display]
            for i, c_name in enumerate(final_cols_to_display):
                if c_name == "Expression" and isinstance(display_values[i], (int, float)):
                    display_values[i] = f"{display_values[i]:.4g}"
            tv.insert("", tk.END, values=display_values)

        def on_item_double_click(event):
            item_id = tv.identify_row(event.y)
            if not item_id: return
            column_clicked_id = tv.identify_column(event.x)
            column_index = int(column_clicked_id.replace('#','')) - 1
            item_values = tv.item(item_id, "values")
            gsm_col_idx = final_cols_to_display.index("GSM") if "GSM" in final_cols_to_display else -1
            gse_col_idx = final_cols_to_display.index("GSE") if "GSE" in final_cols_to_display else -1
            if column_index == gse_col_idx and gse_col_idx != -1:
                acc_val = item_values[gse_col_idx]
                if isinstance(acc_val, str) and acc_val.startswith("GSE"):
                    webbrowser.open(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc_val}")
            elif column_index == gsm_col_idx and gsm_col_idx != -1:
                acc_val = item_values[gsm_col_idx]
                if isinstance(acc_val, str) and acc_val.startswith("GSM"):
                    webbrowser.open(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc_val}")
        
        tv.bind("<Double-1>", on_item_double_click)

        # --- Bottom Pane: GSE Distribution Table ---
        gse_dist_frame = ttk.LabelFrame(main_pane, text="GSE Distribution in Selection", padding=5)
        main_pane.add(gse_dist_frame, weight=1)

        if "GSE" in df_show.columns and df_show["GSE"].notna().any():
            selected_counts_per_gse = df_show["GSE"].value_counts()
            full_token_df = self._get_tok_df(plat_name)
            gse_stats_data = []
            
            if full_token_df is not None and any(c in full_token_df.columns for c in ["series_id", "GSE"]):
                gse_col_in_tok_df = 'GSE' if 'GSE' in full_token_df.columns else 'series_id'
                for gse_id in selected_counts_per_gse.index:
                    num_selected = selected_counts_per_gse.get(gse_id, 0)
                    total_gse_samples = full_token_df[full_token_df[gse_col_in_tok_df] == gse_id].shape[0]
                    percentage = (num_selected / total_gse_samples) * 100 if total_gse_samples > 0 else 0
                    gse_stats_data.append({
                        "GSE ID": gse_id,
                        "Selected Samples %": percentage,
                        "Total Samples in GSE": total_gse_samples,
                        "Samples in Selection": num_selected
                    })
            
                gse_stats_df = pd.DataFrame(gse_stats_data).sort_values(by="Selected Samples %", ascending=False)
                
                gse_tree = ttk.Treeview(gse_dist_frame, columns=list(gse_stats_df.columns), show="headings")
                for col in gse_stats_df.columns:
                    gse_tree.heading(col, text=col)
                    gse_tree.column(col, width=150, anchor=CENTER)
                
                for _, row in gse_stats_df.iterrows():
                    gse_tree.insert("", "end", values=(
                        row["GSE ID"],
                        f"{row['Selected Samples %']:.1f}%",
                        int(row["Total Samples in GSE"]),
                        int(row["Samples in Selection"])
                    ), tags=(row["GSE ID"],))

                gse_tree_vsb = ttk.Scrollbar(gse_dist_frame, orient=tk.VERTICAL, command=gse_tree.yview)
                gse_tree_hsb = ttk.Scrollbar(gse_dist_frame, orient=tk.HORIZONTAL, command=gse_tree.xview)
                gse_tree.configure(yscrollcommand=gse_tree_vsb.set, xscrollcommand=gse_tree_hsb.set)
                
                gse_dist_frame.grid_rowconfigure(0, weight=1)
                gse_dist_frame.grid_columnconfigure(0, weight=1)
                gse_tree.grid(row=0, column=0, sticky="nsew")
                gse_tree_vsb.grid(row=0, column=1, sticky="ns")
                gse_tree_hsb.grid(row=1, column=0, sticky="ew")

                def on_gse_tree_double_click(event):
                    item_id = gse_tree.identify_row(event.y)
                    if item_id:
                        gse_id_clicked = gse_tree.item(item_id, "tags")[0]
                        if gse_id_clicked.startswith("GSE"):
                            webbrowser.open(f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id_clicked}")
                gse_tree.bind("<Double-1>", on_gse_tree_double_click)
            else:
                tk.Label(gse_dist_frame, text="Could not retrieve full GSE sample counts for distribution analysis.").pack(pady=5)
        else:
            tk.Label(gse_dist_frame, text="GSE column not available or empty in the data for distribution analysis.").pack(pady=5)

        ttk.Button(win_table, text="Close", command=win_table.destroy).pack(pady=5)
        win_table.transient(self)
        win_table.lift()
        win_table.attributes('-topmost', False)


    def _perform_token_clustering_and_display(self, bio_specific_tokens_df: pd.DataFrame, title_prefix: str, output_dir: str):
        """
        REPAIRED: This function now performs ontology-based clustering and calls the display/save functions.
        """
        if not bio_specific_tokens_df.empty:
            self.enqueue_log(f"[Clustering] Starting ontology-based clustering for {title_prefix}...")

            # Map each token to its canonical ontology term
            bio_specific_tokens_df['canonical_token'] = bio_specific_tokens_df['token'].apply(
                lambda t: find_canonical_term(t, self.disease_map, self.uberon_map)
            )

            # Group by the new canonical token to form clusters
            clusters = {}
            for canon, group in bio_specific_tokens_df.groupby('canonical_token'):
                members = []
                for _, row in group.iterrows():
                    members.append((row['token'], int(row['sel_count'])))
                clusters[canon] = members
            
            base_name = re.sub(r'[^\w\-]+', '_', title_prefix)

            # Display and SAVE canonical word cloud
            self._display_canonical_wordcloud_popup(clusters, self.disease_canons, self.uberon_canons, title_prefix, output_dir, base_name)

            # Display and SAVE graph
            self._display_graph_popup(clusters, self.disease_canons, self.uberon_canons, title_prefix, output_dir, base_name)

            # Display and SAVE frequency bar plot
            self._display_frequency_barplot_popup(clusters, self.disease_canons, self.uberon_canons, title_prefix, output_dir, base_name)

            self.enqueue_log(f"[Clustering] Cluster visualizations generated, displayed, and saved for {title_prefix}.")
        else:
            self.enqueue_log(f"[Clustering] No bio-specific tokens available for clustering for {title_prefix}. Skipping cluster plots.")

    def _display_density_plot_popup(self, selected_expressions: pd.Series, background_expressions: pd.Series, title: str, output_dir: str, base_name: str):
        """NEW: Displays a density plot and SAVES it to the specified output directory."""
        import matplotlib.pyplot as plt_local
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        import seaborn as sns

        if selected_expressions.empty:
            self.enqueue_log(f"[Density Plot] No selected expression data for '{title}'. Skipping plot.")
            return

        win = tk.Toplevel(self)
        win.title(title)
        win.geometry("800x600")

        fig, ax = plt_local.subplots(figsize=(8, 5.5))
        
        # Plot densities using seaborn for a nice look
        sns.kdeplot(selected_expressions, color='red', fill=True, alpha=0.5, label=f'Selected (n={len(selected_expressions)})', ax=ax)
        if not background_expressions.empty:
            sns.kdeplot(background_expressions, color='blue', fill=True, alpha=0.5, label=f'Background (n={len(background_expressions)})', ax=ax)

        ax.set_title(title, fontsize=12, weight='bold')
        ax.set_xlabel("Gene Expression Level", fontsize=10)
        ax.set_ylabel("Density", fontsize=10)
        ax.legend()
        ax.grid(axis='y', linestyle='--', alpha=0.7)
        fig.tight_layout()

        # --- ADD THIS BLOCK TO SAVE THE FIGURE ---
        try:
            out_file = Path(output_dir) / f"density_plot_{base_name}.png"
            fig.savefig(out_file, dpi=150, bbox_inches='tight')
            self.enqueue_log(f"[DENSITY PLOT] Plot saved to {out_file}")
        except Exception as e:
            self.enqueue_log(f"[ERROR] Failed to save density plot: {e}")
        # ------------------------------------------

        # --- Embed in Tkinter ---
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas_widget = canvas.get_tk_widget(); canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, win); toolbar.update(); toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.draw()

        key = f"Density_Plot_{title}_{uuid.uuid4().hex}"
        self._current_popup_figs[key] = (fig, canvas_widget, toolbar)
        win.protocol("WM_DELETE_WINDOW", lambda: (plt_local.close(fig), win.destroy(), self._current_popup_figs.pop(key, None)))
        win.transient(self); win.lift()
        win.attributes('-topmost', True)
        win.after(50, lambda: win.attributes('-topmost', False))


    def _display_graph_popup(self, clusters: dict, disease_canons: Set[str], uberon_canons: Set[str], title_prefix: str, output_dir: str, base_name: str):
        """REPAIRED: Displays the force-directed graph and SAVES it to the specified output directory."""
        if not NETWORKX_AVAILABLE:
            self.enqueue_log(f"[GRAPH] NetworkX not available, skipping graph display for {title_prefix}.")
            return
        import matplotlib.pyplot as plt_local
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

        if not clusters: return

        win = tk.Toplevel(self)
        win.title(f"Canonical Token Graph: {title_prefix}")
        win.geometry("800x800")

        G = nx.Graph()
        labels = {}
        canonical_nodes = list(clusters.keys())

        for canon, members in clusters.items():
            if canon not in G: G.add_node(canon); labels[canon] = canon
            for tok, cnt in members:
                if tok not in G: G.add_node(tok); labels[tok] = f"{tok} ({cnt})"
                G.add_edge(canon, tok)

        n_nodes = G.number_of_nodes()
        fig_size = max(10, n_nodes * 0.1)
        fig, ax = plt_local.subplots(figsize=(fig_size, fig_size))
        pos = nx.spring_layout(G, k=1.5/math.sqrt(n_nodes) if n_nodes > 0 else 1, iterations=50, seed=42)
        nx.draw_networkx_edges(G, pos, edge_color='gray', alpha=0.6, width=1.0, ax=ax)

        for node, (x, y) in pos.items():
            if node in canonical_nodes:
                color = get_canon_color(node, disease_canons, uberon_canons)
                circle = plt_local.Circle((x, y), radius=0.04, color=color, fill=False, lw=2, zorder=2)
                ax.add_patch(circle)
                ax.text(x, y, labels[node], color=color, fontweight='bold', ha='center', va='center', fontsize=10, zorder=3)
            else:
                ax.text(x, y, labels[node], color='black', ha='center', va='center', fontsize=8, zorder=3)

        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label='Disease', markerfacecolor='red', markersize=12),
            Line2D([0], [0], marker='o', color='w', label='Tissue', markerfacecolor='green', markersize=12),
            Line2D([0], [0], marker='o', color='w', label='Other', markerfacecolor='black', markersize=12),
            Line2D([0], [0], marker='o', color='w', label='Species', markerfacecolor='grey', markersize=12)
        ]
        ax.legend(handles=legend_elements, title="Canonical Gene Categories", loc='upper right', fontsize=10)
        ax.set_title(f"Token Clusters for {title_prefix}", fontsize=12, weight='bold')
        ax.axis("off")
        
        if pos:
            x_vals, y_vals = zip(*pos.values())
            ax.set_xlim(min(x_vals) - 0.1, max(x_vals) + 0.1)
            ax.set_ylim(min(y_vals) - 0.1, max(y_vals) + 0.1)
        fig.tight_layout()

        # REPAIRED: Save the figure to the correct directory
        out_file = Path(output_dir) / f"{base_name}_clusters_graph.png"
        fig.savefig(out_file, dpi=150, bbox_inches='tight')
        self.enqueue_log(f"[GRAPH] Graph saved to {out_file}")

        # --- Display in Tkinter ---
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas_widget = canvas.get_tk_widget(); canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, win); toolbar.update(); toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.draw()

        key = f"Graph_Plot_{title_prefix}_{uuid.uuid4().hex}"
        self._current_popup_figs[key] = (fig, canvas_widget, toolbar)
        win.protocol("WM_DELETE_WINDOW", lambda: (plt_local.close(fig), win.destroy(), self._current_popup_figs.pop(key, None)))
        win.transient(self); win.lift()
        win.attributes('-topmost', True)
        win.after(50, lambda: win.attributes('-topmost', False))

    def _display_canonical_wordcloud_popup(self, clusters: dict, disease_canons: Set[str], uberon_canons: Set[str], title_prefix: str, output_dir: str, base_name: str):
        """REPAIRED: Displays the canonical word cloud and SAVES it to the specified output directory."""
        import matplotlib.pyplot as plt_local
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

        if not clusters: return

        frequencies = {canon: sum(cnt for _, cnt in members) for canon, members in clusters.items() if sum(cnt for _, cnt in members) > 0}
        if not frequencies: return

        win = tk.Toplevel(self)
        win.title(f"Canonical Token Word Cloud: {title_prefix}")
        win.geometry("800x600")

        fig, ax = plt_local.subplots(figsize=(16, 10))
        
        def color_function(word, **kwargs):
            return get_canon_color(word, disease_canons, uberon_canons)

        wc = WordCloud(width=1600, height=1000, background_color='white', collocations=False, color_func=color_function)
        wc.generate_from_frequencies(frequencies)

        ax.imshow(wc, interpolation='bilinear'); ax.axis("off")
        ax.set_title(f"Canonical Token Word Cloud for: {title_prefix}", fontsize=24, weight='bold')

        legend_elements = [
            Line2D([0], [0], marker='s', color='w', label='Disease', markerfacecolor='red', markersize=15),
            Line2D([0], [0], marker='s', color='w', label='Tissue', markerfacecolor='green', markersize=15),
            Line2D([0], [0], marker='s', color='w', label='Other', markerfacecolor='black', markersize=15),
            Line2D([0], [0], marker='s', color='w', label='Species', markerfacecolor='grey', markersize=15)
        ]
        ax.legend(handles=legend_elements, title="Category", loc='lower left', fontsize=12, frameon=True, facecolor='white', framealpha=0.8)
        
        plt.tight_layout(pad=1)

        # REPAIRED: Save the figure to the correct directory
        out_file = Path(output_dir) / f"{base_name}_clusteredtokens_wordcloud.png"
        plt.savefig(out_file, dpi=150, bbox_inches='tight')
        self.enqueue_log(f"[WC] Canonical word cloud saved to {out_file}")

        # --- Display in Tkinter ---
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas_widget = canvas.get_tk_widget(); canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, win); toolbar.update(); toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.draw()

        key = f"Canonical_WC_{title_prefix}_{uuid.uuid4().hex}"
        self._current_popup_figs[key] = (fig, canvas_widget, toolbar)
        win.protocol("WM_DELETE_WINDOW", lambda: (plt_local.close(fig), win.destroy(), self._current_popup_figs.pop(key, None)))
        win.transient(self); win.lift()
        win.attributes('-topmost', True)
        win.after(50, lambda: win.attributes('-topmost', False))

    def _display_frequency_barplot_popup(self, clusters: dict, disease_canons: Set[str], uberon_canons: Set[str], title_prefix: str, output_dir: str, base_name: str):
        """REPAIRED: Displays the frequency bar plot and SAVES it to the specified output directory."""
        import matplotlib.pyplot as plt_local
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

        if not clusters: return

        frequencies = {canon: sum(cnt for _, cnt in members) for canon, members in clusters.items() if sum(cnt for _, cnt in members) > 0}
        if not frequencies: return

        sorted_canons = sorted(frequencies.items(), key=lambda item: item[1], reverse=False)
        
        labels = [item[0] for item in sorted_canons]
        counts = [item[1] for item in sorted_canons]
        colors = [get_canon_color(label, disease_canons, uberon_canons) for label in labels]

        win = tk.Toplevel(self)
        win.title(f"Canonical Token Frequency Plot: {title_prefix}")
        win.geometry("900x700")

        fig, ax = plt_local.subplots(figsize=(10, max(7, len(labels) * 0.3)))
        ax.barh(labels, counts, color=colors)
        
        ax.set_xlabel("Total Frequency (Sum of Member Token Counts)", fontsize=10)
        ax.set_ylabel("Clustered Tokens (based on ontology matching)", fontsize=10)
        ax.set_title(f"Canonical Token Frequencies for {title_prefix}", fontsize=12, weight='bold')
        ax.tick_params(axis='y', labelsize=8)
        ax.grid(axis='x', linestyle='--', alpha=0.7)

        legend_elements = [
            Line2D([0], [0], marker='s', color='red', label='Disease', markersize=8, linestyle='None'),
            Line2D([0], [0], marker='s', color='green', label='Tissue', markersize=8, linestyle='None'),
            Line2D([0], [0], marker='s', color='black', label='Other', markersize=8, linestyle='None'),
            Line2D([0], [0], marker='s', color='grey', label='Species', markersize=8, linestyle='None')
        ]
        ax.legend(handles=legend_elements, title="Category", loc='lower right', fontsize=8)

        fig.tight_layout()

        # REPAIRED: Save the figure to the correct directory
        out_file = Path(output_dir) / f"{base_name}_clusteredtokens_frequencyplot.png"
        plt.savefig(out_file, dpi=150, bbox_inches='tight')
        self.enqueue_log(f"[BARPLOT] Frequency bar plot saved to {out_file}")

        # --- Display in Tkinter ---
        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas_widget = canvas.get_tk_widget(); canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, win); toolbar.update(); toolbar.pack(side=tk.BOTTOM, fill=tk.X)
        canvas.draw()

        key = f"Frequency_Plot_{title_prefix}_{uuid.uuid4().hex}"
        self._current_popup_figs[key] = (fig, canvas_widget, toolbar)
        win.protocol("WM_DELETE_WINDOW", lambda: (plt_local.close(fig), win.destroy(), self._current_popup_figs.pop(key, None)))
        win.transient(self); win.lift()
        win.attributes('-topmost', True)
        win.after(50, lambda: win.attributes('-topmost', False))

if __name__ == "__main__":
    if not os.path.exists("NEW_RESULTS_ROOT"): # FIX: Changed RESULTS_ROOT
        os.makedirs("NEW_RESULTS_ROOT")

    app = GeoWorkflowGUI()
    app.mainloop()
