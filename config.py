# config.py
from pathlib import Path

# --- Project Root ---
# This assumes the script is run from the project's root directory.
PROJECT_ROOT = Path(__file__).resolve().parent

# --- Core Data Paths ---
# INSTRUCTION: Replace "directory_to_this" with the actual paths to your data directories.
DATA_DIR = Path("directory_to_this/GEOMETADB_TOKENS_PLATFORMS")
MODEL_DIR = Path("directory_to_this/finetuned_pubmedbert")
KEYWORDS_DIR = Path("directory_to_this/GPL570_keywords") # Or a more general keywords directory

# --- Database Path ---
GEOMETADB_PATH = PROJECT_ROOT / "GEOmetadb.sqlite.gz"

# --- Output & Cache Directories ---
RESULTS_ROOT = PROJECT_ROOT / "NEW_RESULTS_ROOT"
CACHE_DIR = PROJECT_ROOT / ".cache"
ONTOLOGY_CACHE_DIR = CACHE_DIR / "ontology_cache"
EMBEDDING_CACHE_DIR = CACHE_DIR / "embedding_cache"

# --- Derived Paths ---
KEYWORD_CORPUS_PATH = KEYWORDS_DIR / "biomed_corpus.txt"
KEYWORD_CACHE_PATH = EMBEDDING_CACHE_DIR / "biomed_keyword_embeddings.pkl"

# --- GPL File Paths ---
GPL_PATHS = {
    "GPL570": {
        "data": DATA_DIR / "GPL570/gpl570_all_samples_normalized_non_scaled_with_nans.csv.gz",
        "metadata": DATA_DIR / "GPL570/GPL570_metadata.csv.gz"
    },
    "GPL96": {
        "data": DATA_DIR / "GPL96/gpl96_normalized_scaled.csv.gz"
    },
    "GPL6947": {
        "data": DATA_DIR / "GPL6947/gpl6947_all_samples_normalized_scaled_with_nans.csv.gz"
    },
    "GPL7202": {
        "data": DATA_DIR / "GPL7202/gpl7202_all_samples_normalized_scaled_with_nans.csv.gz"
    },
    "GPL6885": {
        "data": DATA_DIR / "GPL6885/gpl6885_all_samples_normalized_scaled_with_nans.csv.gz"
    },
    "GPL1261": {
        "data": DATA_DIR / "GPL1261/gpl1261_all_samples_normalized_scaled_with_nans.csv.gz"
    },
    "GPL10558": {
        "data": DATA_DIR / "GPL10558/gpl10558_all_samples_normalized_scaled_with_nans.csv.gz"
    }
}

# --- Analysis Settings ---
SIMILARITY_THRESHOLD = 0.65
