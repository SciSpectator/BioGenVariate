import gzip
import re
import urllib.request
from collections import defaultdict
from typing import IO, Dict, List, Set, Tuple
from pathlib import Path
import config

# Ensure the cache directory exists before any operations
config.ONTOLOGY_CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Define the ontologies to be downloaded and their URLs
ONTOLOGIES = {
    "uberon.obo": "https://purl.obolibrary.org/obo/uberon.obo",
    "doid.obo": "https://purl.obolibrary.org/obo/doid.obo",
}

# Define a set of common species terms to be categorized separately
SPECIES_TERMS = {"homo sapiens", "mus musculus", "rattus norvegicus", "homo", "sapiens", "human", "mouse", "rat"}

def cached_path(fname: str, url: str) -> Path:
    """
    Checks if an ontology file exists in the cache directory.
    If not, it downloads the file from the specified URL.
    """
    path = config.ONTOLOGY_CACHE_DIR / fname
    if not path.exists():
        print(f"[ONTOLOGY] Downloading {fname}...")
        tmp_path = path.with_suffix(".tmp")
        try:
            urllib.request.urlretrieve(url, tmp_path)
            tmp_path.rename(path)
            print(f"[ONTOLOGY] Saved -> {path.name} ({path.stat().st_size / 1e6:.1f} MB)")
        except Exception as e:
            print(f"[ONTOLOGY ERROR] Failed to download {url} to {path}: {e}")
            if tmp_path.exists():
                tmp_path.unlink()  # Clean up the temporary file on failure
            raise
    return path

def smart_open(path: Path) -> IO[str]:
    """Opens a file, transparently handling gzip compression if present."""
    with open(path, "rb") as fh:
        # Read the first two bytes to check for gzip magic number
        magic_number = fh.read(2)
    if magic_number == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "rt", encoding="utf-8", errors="ignore")

def parse_obo(path: Path) -> List[Tuple[str, List[str]]]:
    """
    Parses a .obo file to extract canonical term names and their synonyms.
    Skips obsolete terms.
    """
    terms: List[Tuple[str, List[str]]] = []
    try:
        with smart_open(path) as fh:
            in_term_block = False
            current_label, current_synonyms, is_obsolete = None, [], False
            for line in fh:
                if line.startswith("[Term]"):
                    # Save the previous term if it was valid
                    if in_term_block and current_label and not is_obsolete:
                        terms.append((current_label, current_synonyms))
                    # Reset for the new term
                    current_label, current_synonyms, is_obsolete = None, [], False
                    in_term_block = True
                    continue
                
                if not in_term_block:
                    continue

                if line.startswith("name:"):
                    current_label = line.partition("name:")[2].strip()
                elif line.startswith("synonym:"):
                    # Synonyms are enclosed in double quotes
                    try:
                        synonym = line.split('"', 2)[1]
                        current_synonyms.append(synonym)
                    except IndexError:
                        continue # Skip malformed synonym lines
                elif line.startswith("is_obsolete: true"):
                    is_obsolete = True
            
            # Add the very last term in the file if it's valid
            if in_term_block and current_label and not is_obsolete:
                terms.append((current_label, current_synonyms))

    except Exception as e:
        print(f"[ONTOLOGY ERROR] Failed to parse OBO file {path}: {e}")
        return []  # Return an empty list on error
    return terms

def load_ontology_sets() -> Tuple[Set[str], Set[str], Dict, Dict]:
    """
    Loads, parses, and processes the Uberon and DOID ontologies.
    Returns sets of canonical terms and pre-processed dictionaries for fast lookup.
    """
    print("[ONTOLOGY] Loading and pre-processing ontologies for categorization...")
    uberon_map, disease_map = defaultdict(list), defaultdict(list)
    uberon_canons, disease_canons = set(), set()
    
    try:
        # Process Uberon (Tissue/Anatomy) Ontology
        uberon_terms = parse_obo(cached_path("uberon.obo", ONTOLOGIES["uberon.obo"]))
        for canon, synonyms in uberon_terms:
            uberon_canons.add(canon.lower())
            all_terms = [canon] + synonyms
            for term in all_terms:
                for word in re.split(r'\W+', term.lower()):
                    if word and len(word) > 2:
                        uberon_map[word].append(canon)

        # Process DOID (Disease) Ontology
        disease_terms = parse_obo(cached_path("doid.obo", ONTOLOGIES["doid.obo"]))
        for canon, synonyms in disease_terms:
            disease_canons.add(canon.lower())
            all_terms = [canon] + synonyms
            for term in all_terms:
                for word in re.split(r'\W+', term.lower()):
                    if word and len(word) > 2:
                        disease_map[word].append(canon)

        # Ensure canonical mappings are unique and sorted (shortest first)
        for word in uberon_map:
            uberon_map[word] = sorted(list(set(uberon_map[word])), key=len)
        for word in disease_map:
            disease_map[word] = sorted(list(set(disease_map[word])), key=len)
            
        print(f"[ONTOLOGY] Loaded {len(uberon_canons)} Uberon terms and {len(disease_canons)} Disease terms.")
    
    except Exception as e:
        print(f"[ONTOLOGY ERROR] Could not load all ontology files: {e}. Categorization may be incomplete.")
    
    return uberon_canons, disease_canons, uberon_map, disease_map

def find_canonical_term(token: str, disease_map: Dict, uberon_map: Dict) -> str:
    """Finds the best matching canonical term for a given token, prioritizing disease terms."""
    token_words = {word for word in re.split(r'\W+', token.lower()) if len(word) > 2}
    
    # Prioritize disease matches
    for word in token_words:
        if word in disease_map:
            return disease_map[word][0]  # Return the shortest canonical term
            
    # Then check tissue matches
    for word in token_words:
        if word in uberon_map:
            return uberon_map[word][0]  # Return the shortest canonical term
            
    return token  # Fallback to the token itself if no match is found

def get_canon_color(canon_name: str, disease_canons: Set[str], uberon_canons: Set[str]) -> str:
    """Determines the color for a canonical token based on its ontology."""
    canon_lower = canon_name.lower()
    if canon_lower in SPECIES_TERMS:
        return 'grey'
    if canon_lower in disease_canons:
        return 'red'
    elif canon_lower in uberon_canons:
        return 'green'
    else:
        return 'black'
