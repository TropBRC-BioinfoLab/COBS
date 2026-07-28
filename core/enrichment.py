"""core.enrichment

Enrichment analysis helper functions using Enrichr (Ma'ayan Lab).

Workflow:
- addList (POST multipart) -> returns userListId
- enrich (GET) -> returns enrichment rows for a given library

This module is intentionally dependency-light: only `requests` + (optional) `pandas`.
"""

from __future__ import annotations

from typing import Any, Dict, List, Sequence
import time
import requests

try:
    import pandas as pd  # type: ignore
except Exception:  # pragma: no cover
    pd = None  # type: ignore


# -----------------------
# Enrichr (Ma'ayan Lab)
# -----------------------
ENRICHR_ADDLIST_URL = "https://maayanlab.cloud/Enrichr/addList"
ENRICHR_ENRICH_URL = "https://maayanlab.cloud/Enrichr/enrich"

# Fixed libraries requested (no UI selection)
# label -> Enrichr backgroundType
ENRICHR_FIXED_LIBS: Dict[str, str] = {
    "GO Biological Process": "GO_Biological_Process_2021",
    "GO Molecular Function": "GO_Molecular_Function_2021",
    "GO Cellular Component": "GO_Cellular_Component_2021",
    "KEGG Pathway": "KEGG_2021_Human",
}


def enrichr_add_list(genes: Sequence[str], description: str = "") -> Dict[str, Any]:
    """Upload gene list to Enrichr and return the addList response.

    Enrichr expects a newline-separated string of gene identifiers.
    """
    genes_clean = [str(g).strip() for g in genes if str(g).strip()]
    genes_str = "\n".join(genes_clean)

    # Enrichr docs often show multipart form ("files") payload.
    payload = {
        "list": (None, genes_str),
        "description": (None, description or "gene list"),
    }

    resp = requests.post(ENRICHR_ADDLIST_URL, files=payload, timeout=60)
    if not resp.ok:
        raise RuntimeError(f"Enrichr addList failed: HTTP {resp.status_code}: {resp.text[:300]}")
    return resp.json()


def enrichr_enrich(user_list_id: int, library: str) -> Dict[str, Any]:
    """Run Enrichr enrichment for a given library."""
    params = {
        "userListId": str(user_list_id),
        "backgroundType": library,
    }
    resp = requests.get(ENRICHR_ENRICH_URL, params=params, timeout=60)
    if not resp.ok:
        raise RuntimeError(f"Enrichr enrich failed: HTTP {resp.status_code}: {resp.text[:300]}")
    return resp.json()


def enrichr_to_dataframe(enrich_json: Dict[str, Any], library: str):
    """Convert Enrichr enrichment JSON to a DataFrame (if pandas is available)."""
    rows = enrich_json.get(library, []) or []
    # Enrichr returns list rows like:
    # [rank, term_name, p_value, z_score, combined_score, overlapping_genes, adjusted_p_value, ...]
    parsed: List[Dict[str, Any]] = []
    for r in rows:
        if not isinstance(r, (list, tuple)) or len(r) < 7:
            continue
        parsed.append(
            {
                "rank": r[0],
                "term": r[1],
                "p_value": r[2],
                "z_score": r[3],
                "combined_score": r[4],
                "overlap_genes": r[5],
                "adj_p_value": r[6],
            }
        )
    if pd is None:
        return parsed
    return pd.DataFrame(parsed)


def run_enrichr(
    genes: Sequence[str],
    library: str,
    description: str = "",
    sleep_s: float = 0.0,
):
    """Convenience wrapper: addList -> enrich -> DataFrame/rows."""
    add = enrichr_add_list(genes, description=description)
    user_list_id = int(add["userListId"])
    if sleep_s and float(sleep_s) > 0:
        time.sleep(float(sleep_s))
    enr = enrichr_enrich(user_list_id, library=library)
    return enrichr_to_dataframe(enr, library=library)
