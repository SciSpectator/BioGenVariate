# utils/helpers.py
import re
import json
from typing import List
import pandas as pd
import nltk
from nltk.tokenize import word_tokenize
from nltk.corpus import stopwords, wordnet as wn
from nltk.stem import WordNetLemmatizer
from nltk.util import ngrams

# NLTK setup
for model in ["punkt", "stopwords", "wordnet"]:
    try:
        if model == "punkt":
             nltk.data.find(f'tokenizers/{model}')
        else:
             nltk.data.find(f'corpora/{model}')
    except nltk.downloader.DownloadError:
        nltk.download(model, quiet=True)

lemmatizer = WordNetLemmatizer()
stop_words = set(stopwords.words("english"))

def extract_search_tokens(text, keywords):
    cleaned = re.sub(r"[^A-Za-z0-9 ]", " ", str(text)).lower()
    toks = [lemmatizer.lemmatize(t) for t in word_tokenize(cleaned) if t not in stop_words]
    found = []
    for n in (2, 3):
        for ng in ngrams(toks, n):
            phrase = " ".join(ng)
            if phrase in keywords:
                found.append(phrase)
    found += [t for t in toks if t in keywords]
    for kw in keywords:
        found += re.findall(rf"\b{re.escape(kw)}\b", cleaned)
    return list(set(found))

def build_final_text(row):
    return " ".join(
        str(row[col]) for col in ("title", "description", "characteristics")
        if col in row and pd.notnull(row[col])
    )

def flatten_tokens(series) -> List[str]:
    all_tokens = []
    for item in series.dropna():
        if isinstance(item, str):
            try:
                tokens_list = json.loads(item)
            except (json.JSONDecodeError, TypeError):
                tokens_list = [t.strip() for t in item.split(',') if t.strip()]
        elif isinstance(item, list):
            tokens_list = item
        else:
            tokens_list = []

        if isinstance(tokens_list, list):
            all_tokens.extend([str(t).strip().lower() for t in tokens_list if len(str(t).strip()) > 2])
    return all_tokens
