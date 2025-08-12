# utils/ontology.py
import gzip
import re
import urllib.request
from collections import defaultdict
from typing import IO, Dict, List, Set, Tuple
from pathlib import Path
import config

config.ONTOLOGY_CACHE_DIR.mkdir(parents=True, exist_ok=True)

ONTOLOGIES = {
    "uberon.obo": "https://purl.obolibrary.org/obo/uberon.obo",
    "doid.obo": "https://purl.obolibrary.org/obo/doid.obo",
}
SPECIES_TERMS = {"homo sapiens", "mus musculus", "rattus norvegicus", "human", "mouse", "rat", "homo", "sapiens"}

def cached_path(fname: str, url: str) -> Path:
    path = config.ONTOLOGY_CACHE_DIR / fname
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
                tmp.unlink()
            raise
    return path

def smart_open(path: Path) -> IO[str]:
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "rt", encoding="utf-8", errors="ignore")

def parse_obo(path: Path) -> List[Tuple[str, List[str]]]:
    # ... (code for parse_obo)
    pass

def load_ontology_sets() -> Tuple[Set[str], Set[str], Dict, Dict]:
    # ... (code for load_ontology_sets)
    pass

def find_canonical_term(token: str, disease_map: Dict, uberon_map: Dict) -> str:
    # ... (code for find_canonical_term)
    pass

def get_canon_color(canon_name: str, disease_canons: Set[str], uberon_canons: Set[str]) -> str:
    # ... (code for get_canon_color)
    pass
