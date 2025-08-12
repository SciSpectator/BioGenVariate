"""
Tokenise GEOmetadb (7 GPL platforms) with BioNLP-13CG
=====================================================

This script processes sample data from specific platforms within the GEOmetadb
database. It extracts raw text from various metadata fields, tokenises this
text using a biomedical-specific NLP model, and saves the resulting tokens.

The NLP model used is `en_ner_bionlp13cg_md`, which is trained on the
BioNLP-13CG (Cancer Genetics) shared task dataset. This makes it highly
effective at recognizing biomedical entities like genes, proteins, tissues,
diseases, and cell types from scientific text.

Key Features:
- Auto-downloads the required spaCy model if it's not already installed.
- Automatically uses a GPU if available, otherwise defaults to CPU.
- For each platform, it creates a gzipped CSV file containing the sample ID (GSM)
  and a JSON list of its extracted tokens.
- Includes a resume functionality to continue processing from where it left off.
- Prints a detailed example every 200 samples for monitoring purposes.
"""

# ── imports ─────────────────────────────────────────────────────
"""
Import all necessary libraries for the script.
- os: Used for interacting with the operating system.
- sqlite3: For interacting with the GEOmetadb SQLite database.
- gzip: To handle the compressed database file.
- tempfile: To temporarily store the decompressed database on disk.
- re: For regular expression operations used in text cleaning.
- json: To format the list of tokens for storage in the CSV.
- sys: To exit the script if the database file is not found.
- textwrap: To format the raw text for clean printing in the demo output.
- warnings: To suppress specific warnings from the NLTK library.
- pathlib: For object-oriented handling of file system paths.
- pandas: For data manipulation and reading/writing CSV files.
- nltk: Used here for its list of common English "stopwords".
- spacy: The core NLP library for tokenisation and entity recognition.
- tqdm: To display progress bars for long-running loops.
"""
import os, sqlite3, gzip, tempfile, re, json, sys, textwrap, warnings
from pathlib import Path
import pandas as pd, nltk, spacy
from nltk.corpus import stopwords
from tqdm import tqdm

# --- On-demand model download ----------------------------------
def ensure_spacy_model(name: str):
    """
    Checks if a spaCy model is installed and downloads it if not.
    This avoids errors and makes the script easier to set up.

    Args:
        name (str): The name of the spaCy model to load (e.g., "en_ner_bionlp13cg_md").

    Returns:
        A loaded spaCy language model object.
    """
    import importlib
    try:
        # Try to load the model directly.
        return spacy.load(name)
    except Exception:
        # If loading fails, the model is likely not installed.
        print(f"[INFO] downloading {name} …")
        from spacy.cli import download
        download(name)
        # Invalidate caches to ensure the new model is found.
        importlib.invalidate_caches()
        # Retry loading the model.
        return spacy.load(name)

# ── config ─────────────────────────────────────────────────────
"""
Configuration settings for the script.
- GEO_DB_PATH: The file path to the compressed GEOmetadb SQLite database.
- PLATFORMS: A list of GEO platform IDs (GPLs) to be processed.
- OUTPUT_BASE: The name of the root directory where results will be saved.
- DEMO_EVERY: How often (in number of samples) to print a detailed demo output.
- SAVE_EVERY: How often (in number of samples) to save intermediate progress.
- EXCLUDE: A set of database columns to ignore when concatenating raw text.
  These columns typically contain IDs, dates, or other non-descriptive metadata.
"""
GEO_DB_PATH = "/home/mwinn99/code/notebooks/MateuszSz/GPL570/GEOmetadb.sqlite.gz"
PLATFORMS   = ["GPL96","GPL570","GPL6947","GPL7202","GPL6887","GPL1261","GPL10558"]
OUTPUT_BASE = "GEOMETADB_TOKENS_PLATFORMS"
DEMO_EVERY  = 200
SAVE_EVERY  = 1000

EXCLUDE = {
    'ID','gsm','series_id','gpl','status','submission_date',
    'last_update_date','type','contact','supplementary_file',
    'data_row_count','channel_count'
}

# ── stop-words ──────────────────────────────────────────────────
"""
Set up a list of "stopwords" – common English words (like "the", "a", "is")
that are generally uninformative for tokenisation and are therefore removed.
NLTK's stopwords list is downloaded quietly and stored in a set for fast lookups.
"""
warnings.filterwarnings("ignore", category=UserWarning, module="nltk")
nltk.download("stopwords", quiet=True)
STOP = set(stopwords.words("english"))

# ── load / fetch model ─────────────────────────────────────────
"""
Load the `en_ner_bionlp13cg_md` spaCy model. This model is pre-trained on the
BioNLP-13 Cancer Genetics (CG) task dataset, making it specialized for
recognising biomedical named entities from text, such as genes, proteins,
diseases, tissues, and cell lines. The `ensure_spacy_model` function handles
the download if necessary.

Two additional components from the `scispacy` library are added to the
NLP pipeline:
1. AbbreviationDetector: Helps link abbreviations to their definitions
   (e.g., "G-CSF" -> "granulocyte colony-stimulating factor").
2. UmlsEntityLinker: Links recognised entities to the Unified Medical
   Language System (UMLS) to resolve them to standard concepts.
"""
NLP = ensure_spacy_model("en_ner_bionlp13cg_md")

if "abbreviation_detector" not in NLP.pipe_names:
    from scispacy.abbreviation import AbbreviationDetector
    NLP.add_pipe("abbreviation_detector", last=True)
if "scispacy_linker" not in NLP.pipe_names:
    from scispacy.umls_linking import UmlsEntityLinker
    NLP.add_pipe("scispacy_linker",
                 config={"resolve_abbreviations": True,"max_entities_per_mention":1},
                 last=True)

# ── helpers ────────────────────────────────────────────────────
def clean_lower(t:str)->str:
    """
    Performs basic cleaning on a string.
    1. Removes any character that is not a word character, whitespace, or hyphen.
    2. Replaces multiple whitespace characters with a single space.
    3. Strips leading/trailing whitespace and converts to lowercase.

    Args:
        t (str): The input string.

    Returns:
        str: The cleaned and lowercased string.
    """
    t=re.sub(r"[^\w\s\-]"," ",t)
    return re.sub(r"\s+"," ",t).strip().lower()

def load_sqlite_mem(gz:str):
    """
    Decompresses a gzipped SQLite database file and loads it into memory.
    This allows for much faster database queries compared to reading from disk.

    Args:
        gz (str): Path to the gzipped SQLite file.

    Returns:
        sqlite3.Connection: A connection object to the in-memory database.
    """
    with gzip.open(gz,"rb") as f: buf=f.read()
    import tempfile, os
    # Write the decompressed content to a temporary file on disk.
    tmp=tempfile.NamedTemporaryFile(delete=False,suffix=".sqlite")
    tmp.write(buf); tmp.close()
    # Connect to the disk file, then back it up to an in-memory database.
    disk=sqlite3.connect(tmp.name); mem=sqlite3.connect(":memory:")
    mem.text_factory=lambda b:b.decode("utf-8","replace")
    disk.backup(mem); disk.close(); os.unlink(tmp.name)
    return mem

def raw_text(row):
    """
    Concatenates all relevant text fields from a pandas DataFrame row into a
    single string. It skips columns listed in the `EXCLUDE` set and ignores
    empty or null values.

    Args:
        row (pd.Series): A row from the GSM DataFrame.

    Returns:
        str: A single string containing all descriptive text for the sample.
    """
    return " ".join(str(v).strip() for c,v in row.items()
                    if c not in EXCLUDE and pd.notna(v) and str(v).strip())

def tokens(text):
    """
    The main tokenisation function. It processes a raw text string to extract
    meaningful single-word and multi-word (phrase) tokens.

    Args:
        text (str): The concatenated raw text for a sample.

    Returns:
        list[str]: A sorted list of unique, cleaned tokens.
    """
    if not text.strip(): return []
    # Process the cleaned text with the spaCy NLP pipeline.
    doc=NLP(clean_lower(text))

    # 1. Extract single-word tokens (lemmas).
    # A lemma is the base form of a word (e.g., "running" -> "run").
    # It filters for alphabetic tokens and excludes stopwords.
    singles={t.lemma_.lower() for t in doc if t.is_alpha and t.lemma_.lower() not in STOP}

    # 2. Extract multi-word tokens (named entities/chunks).
    # These are phrases identified by the model (e.g., "colon cancer", "t-cell").
    chunks=set()
    for ent in doc.ents:
        # Lemmatise and clean the words within the entity.
        lem=[t.lemma_.lower() for t in ent if t.is_alpha and t.lemma_.lower() not in STOP]
        if lem:
            ph=" ".join(lem)
            # Exclude phrases that are just numbers.
            if not re.fullmatch(r"\d+",ph): chunks.add(ph)

    # 3. Remove sub-phrases.
    # If we have both "cell" and "stem cell", we only want to keep "stem cell".
    # This loop removes any token that is a substring of another token.
    final=chunks.copy()
    for a in chunks:
        for b in chunks:
            if a!=b and a in b:
                final.discard(a)
                break

    # 4. Combine the single words and the final phrases and return a sorted list.
    return sorted(singles.union(final))

def save(rows,path):
    """
    Saves a list of processed data to a gzipped CSV file.

    Args:
        rows (list[dict]): A list of dictionaries, where each dict has 'gsm' and 'tokens'.
        path (Path): The file path to save the CSV to.
    """
    if rows:
        pd.DataFrame(rows).to_csv(path,index=False,compression="gzip")

# ── main ───────────────────────────────────────────────────────
"""
The main execution block of the script.
"""
# Exit if the database file doesn't exist.
if not Path(GEO_DB_PATH).exists():
    sys.exit(f"DB not found: {GEO_DB_PATH}")

# Load the database into memory and read the entire `gsm` table into a pandas DataFrame.
conn=load_sqlite_mem(GEO_DB_PATH)
df_all=pd.read_sql_query("SELECT * FROM gsm",conn); conn.close()
Path(OUTPUT_BASE).mkdir(exist_ok=True)

# Loop through each platform specified in the config.
for gpl in tqdm(PLATFORMS,desc="Platforms"):
    # Check if the platform exists in the database.
    if gpl not in df_all["gpl"].unique():
        tqdm.write(f"[WARN] {gpl} missing – skip"); continue

    # Set up output directory and file path.
    out_dir=Path(OUTPUT_BASE)/gpl; out_dir.mkdir(parents=True,exist_ok=True)
    out_file=out_dir/f"geometadb_tokens_{gpl}.csv.gz"

    # --- Resume Logic ---
    processed,seen=[],set()
    if out_file.exists():
        try:
            # If an output file already exists, load it.
            prev=pd.read_csv(out_file,compression="gzip")
            # Add its contents to our `processed` list and `seen` set.
            processed.extend(prev.to_dict("records"))
            seen.update(prev["gsm"].astype(str))
            tqdm.write(f"[RESUME] {gpl}: {len(seen)} done")
        except Exception:
            tqdm.write(f"[WARN] cannot resume {gpl}")

    # --- Processing Loop ---
    # Filter the DataFrame for the current platform and iterate over its rows.
    for _,row in tqdm(df_all[df_all["gpl"]==gpl].iterrows(),
                      total=(df_all["gpl"]==gpl).sum(),
                      desc=gpl,leave=False):
        gsm=str(row["gsm"])
        if gsm in seen: continue # Skip if already processed.

        # Get raw text and generate tokens.
        rt=raw_text(row); tok=tokens(rt)
        # Append the result as a dictionary. The tokens are stored as a JSON string.
        processed.append({"gsm":gsm,"tokens":json.dumps(tok,ensure_ascii=False)})
        seen.add(gsm)

        # --- Monitoring and Saving ---
        cnt=len(processed)
        # Print a demo output periodically.
        if cnt%DEMO_EVERY==0:
            tqdm.write(f"\n[{gpl}] DEMO #{cnt}  GSM {gsm}\n"
                       f"=== RAW TEXT ===\n{textwrap.fill(rt,120)}\n\n"
                       f"=== {len(tok)} TOKENS (combined) ===\n{tok}\n")

        # Save progress to the file periodically.
        if cnt%SAVE_EVERY==0:
            save(processed,out_file)

    # Final save after the platform is fully processed.
    save(processed,out_file)
    tqdm.write(f"[DONE] {gpl}: {len(processed)} GSM rows → {out_file}")

