import time
import io
from pathlib import Path
from typing import Iterable, Tuple, Dict, Any, List, Optional

import networkx as nx
import pandas as pd
import requests


# -------------------------
# Local / toy datasets
# -------------------------
def load_default_graph(name: str) -> nx.Graph:
    if name.startswith("Zachary Karate"):
        return nx.karate_club_graph()
    raise ValueError(f"Unknown dataset: {name}")


# -------------------------
# Helpers
# -------------------------
def _split_identifiers(text: str) -> List[str]:
    # terima koma, semicolon, whitespace, dan baris baru
    if not text:
        return []
    cleaned = (
        text.replace("\r", "\n")
            .replace(";", ",")
            .replace("\t", ",")
    )
    parts = []
    for chunk in cleaned.split("\n"):
        for p in chunk.split(","):
            p = p.strip()
            if p:
                parts.append(p)
    # deduplicate sambil menjaga urutan
    seen = set()
    out = []
    for p in parts:
        if p not in seen:
            out.append(p)
            seen.add(p)
    return out


def _string_api_url(base_url: str, output_format: str, method: str) -> str:
    base_url = base_url.rstrip("/")
    return f"{base_url}/api/{output_format}/{method}"


# -------------------------
# UniProt (REST API)
# -------------------------
def _uniprot_kb_search(
    query: str,
    size: int = 200,
    fields: str = "accession,gene_primary,protein_name",
    base_url: str = "https://rest.uniprot.org",
    timeout: int = 30,
) -> Dict[str, Any]:
    url = base_url.rstrip("/") + "/uniprotkb/search"
    params = {
        "query": query,
        "format": "json",
        "fields": fields,
        "size": int(size),
    }
    r = requests.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r.json()


def _extract_gene_primary(item: Dict[str, Any]) -> Optional[str]:
    # struktur JSON UniProt bervariasi; kita ambil yang paling umum
    genes = item.get("genes") or []
    if not genes:
        return None
    # geneName utama biasanya di index 0
    gene_name = genes[0].get("geneName", {})
    if isinstance(gene_name, dict):
        return gene_name.get("value")
    return None


def _extract_protein_name(item: Dict[str, Any]) -> Optional[str]:
    pdict = item.get("proteinDescription") or {}
    rec = pdict.get("recommendedName") or {}
    full = rec.get("fullName") or {}
    if isinstance(full, dict):
        return full.get("value")
    return None


def uniprot_genes_from_disease(
    disease: str,
    organism_id: int = 9606,
    reviewed_only: bool = True,
    size: int = 200,
) -> Tuple[List[str], List[Dict[str, Any]]]:
    """
    Mengambil daftar gen (gene symbols) dari UniProt berdasarkan string penyakit.

    Catatan:
    - UniProt punya query language yang kaya; untuk 'penyakit' kita coba memakai field disease:"...".
      Jika field tersebut tidak dikenali (mis. HTTP 400), kita fallback ke full-text query \"...\"
      agar tetap menghasilkan kandidat protein.
    """
    disease = (disease or "").strip()
    if not disease:
        raise ValueError("Nama penyakit kosong.")

    # 1) coba query berbasis field disease
    base = f'(disease:"{disease}") AND (organism_id:{int(organism_id)})'
    if reviewed_only:
        base = f"({base}) AND (reviewed:true)"

    try:
        payload = _uniprot_kb_search(base, size=size)
        used_query = base
    except requests.HTTPError:
        # 2) fallback full-text
        ft = f'("{disease}") AND (organism_id:{int(organism_id)})'
        if reviewed_only:
            ft = f"({ft}) AND (reviewed:true)"
        payload = _uniprot_kb_search(ft, size=size)
        used_query = ft

    results = payload.get("results") or []

    records: List[Dict[str, Any]] = []
    genes: List[str] = []
    seen = set()

    for it in results:
        acc = it.get("primaryAccession")
        g = _extract_gene_primary(it)
        pname = _extract_protein_name(it)
        records.append(
            {
                "accession": acc,
                "gene_primary": g,
                "protein_name": pname,
            }
        )
        if g and g not in seen:
            genes.append(g)
            seen.add(g)

    # simpan query yang dipakai untuk debugging/reproducibility (opsional)
    for r in records:
        r["uniprot_query"] = used_query

    return genes, records


# -------------------------
# STRING (REST API)
# -------------------------
def string_map_identifiers(
    identifiers: Iterable[str],
    species: int = 9606,
    base_url: str = "https://string-db.org",
    caller_identity: str = "community-detect-streamlit",
    timeout: int = 30,
) -> pd.DataFrame:
    ids = list(identifiers)
    if not ids:
        raise ValueError("Daftar identifier kosong.")

    url = _string_api_url(base_url, "tsv", "get_string_ids")

    # pisahkan dengan carriage return sesuai dokumentasi STRING
    params = {
        "identifiers": "\r".join(ids),
        "species": int(species),
        "limit": 1,
        "echo_query": 1,
        "caller_identity": caller_identity,
    }

    r = requests.post(url, data=params, timeout=timeout)
    r.raise_for_status()
    if not r.text.strip():
        return pd.DataFrame()

    return pd.read_csv(io.StringIO(r.text), sep="\t")


def string_fetch_network(
    identifiers: Iterable[str],
    species: int = 9606,
    required_score: int = 700,
    add_nodes: int = 0,
    base_url: str = "https://string-db.org",
    caller_identity: str = "community-detect-streamlit",
    timeout: int = 60,
) -> pd.DataFrame:
    ids = list(identifiers)
    if not ids:
        raise ValueError("Daftar identifier kosong.")

    url = _string_api_url(base_url, "tsv", "network")
    params = {
        "identifiers": "\r".join(ids),
        "species": int(species),
        "required_score": int(required_score),
        "add_nodes": int(add_nodes),
        "caller_identity": caller_identity,
    }

    r = requests.post(url, data=params, timeout=timeout)
    r.raise_for_status()
    if not r.text.strip():
        return pd.DataFrame()

    return pd.read_csv(io.StringIO(r.text), sep="\t")


def build_nx_graph_from_string_edges(df_edges: pd.DataFrame) -> nx.Graph:
    if df_edges is None or df_edges.empty:
        return nx.Graph()

    # kolom umum dari STRING network API
    # (tetap dibuat robust kalau format berubah sedikit)
    col_a = "preferredName_A" if "preferredName_A" in df_edges.columns else None
    col_b = "preferredName_B" if "preferredName_B" in df_edges.columns else None
    if col_a is None or col_b is None:
        # fallback ke stringId
        col_a = "stringId_A" if "stringId_A" in df_edges.columns else None
        col_b = "stringId_B" if "stringId_B" in df_edges.columns else None

    if col_a is None or col_b is None:
        raise ValueError(f"Kolom node tidak ditemukan pada output STRING: {list(df_edges.columns)}")

    score_col = "score" if "score" in df_edges.columns else None

    G = nx.Graph()
    for _, row in df_edges.iterrows():
        u = str(row[col_a])
        v = str(row[col_b])
        if u == v:
            continue
        if score_col:
            w = float(row[score_col])
            G.add_edge(u, v, weight=w, score=w)
        else:
            G.add_edge(u, v)

    return G


# -------------------------
# Public builders for app
# -------------------------
def build_graph_from_gene_list(
    genes_text: str,
    string_species: int = 9606,
    required_score: int = 700,
    add_nodes: int = 0,
    string_base_url: str = "https://string-db.org",
    caller_identity: str = "community-detect-streamlit",
) -> Tuple[nx.Graph, Dict[str, Any]]:
    genes = _split_identifiers(genes_text)
    if not genes:
        raise ValueError("Tidak ada gen yang bisa diproses.")

    # STRING docs menganjurkan menunggu ~1 detik antar call
    map_df = string_map_identifiers(
        genes,
        species=string_species,
        base_url=string_base_url,
        caller_identity=caller_identity,
    )
    time.sleep(1.0)

    # pakai STRING IDs untuk query network (lebih cepat & lebih konsisten)
    if not map_df.empty and "stringId" in map_df.columns:
        ids_for_network = map_df["stringId"].dropna().astype(str).tolist()
    else:
        ids_for_network = genes

    edges_df = string_fetch_network(
        ids_for_network,
        species=string_species,
        required_score=required_score,
        add_nodes=add_nodes,
        base_url=string_base_url,
        caller_identity=caller_identity,
    )

    G = build_nx_graph_from_string_edges(edges_df)

    meta = {
        "source": "string_gene_list",
        "input_genes": genes[:50],
        "string_mapping": map_df.head(200).to_dict(orient="records") if isinstance(map_df, pd.DataFrame) else [],
        "string_edges": edges_df.head(200).to_dict(orient="records") if isinstance(edges_df, pd.DataFrame) else [],
        "n_input_genes": len(genes),
        "n_nodes": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
    }
    return G, meta


def build_graph_from_disease(
    disease: str,
    organism_id: int = 9606,
    reviewed_only: bool = True,
    uniprot_size: int = 200,
    string_species: int = 9606,
    required_score: int = 700,
    add_nodes: int = 0,
    string_base_url: str = "https://string-db.org",
    caller_identity: str = "community-detect-streamlit",
) -> Tuple[nx.Graph, Dict[str, Any]]:
    genes, records = uniprot_genes_from_disease(
        disease=disease,
        organism_id=organism_id,
        reviewed_only=reviewed_only,
        size=uniprot_size,
    )
    if not genes:
        raise ValueError("UniProt tidak mengembalikan gene_primary. Coba mode 'Gene list → STRING'.")

    # STRING docs menganjurkan menunggu ~1 detik antar call
    map_df = string_map_identifiers(
        genes,
        species=string_species,
        base_url=string_base_url,
        caller_identity=caller_identity,
    )
    time.sleep(1.0)

    if not map_df.empty and "stringId" in map_df.columns:
        ids_for_network = map_df["stringId"].dropna().astype(str).tolist()
    else:
        ids_for_network = genes

    edges_df = string_fetch_network(
        ids_for_network,
        species=string_species,
        required_score=required_score,
        add_nodes=add_nodes,
        base_url=string_base_url,
        caller_identity=caller_identity,
    )

    G = build_nx_graph_from_string_edges(edges_df)

    meta = {
        "source": "uniprot_disease_to_string",
        "disease": disease,
        "uniprot_records": records[:200],
        "n_uniprot_records": len(records),
        "n_gene_primary": len(genes),
        "string_mapping": map_df.head(200).to_dict(orient="records") if isinstance(map_df, pd.DataFrame) else [],
        "string_edges": edges_df.head(200).to_dict(orient="records") if isinstance(edges_df, pd.DataFrame) else [],
        "n_nodes": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
    }
    return G, meta
# -------------------------
# Local disease datasets (*.txt identifiers) -> STRING
# -------------------------
def read_identifiers_from_txt_file(file_path: str) -> List[str]:
    """Baca identifier per baris / dipisah koma dari file teks."""
    p = Path(file_path)
    if not p.exists():
        raise FileNotFoundError(f"File tidak ditemukan: {file_path}")
    text = p.read_text(encoding="utf-8", errors="ignore")
    return _split_identifiers(text)


def build_graph_from_identifier_file(
    file_path: str,
    source_name: str,
    string_species: int = 9606,
    required_score: int = 700,
    add_nodes: int = 0,
    string_base_url: str = "https://string-db.org",
    caller_identity: str = "community-detect-streamlit",
) -> Tuple[nx.Graph, Dict[str, Any]]:
    """
    Bangun graf PPI dari file teks berisi identifier gen/protein (1 per baris
    atau dipisah koma), lalu dimapping ke STRING dan diambil network-nya.

    Catatan:
    - Dapat dipakai untuk UniProt accession maupun gene symbol, selama STRING
      bisa memetakan identifier tersebut untuk species yang dipilih.
    """
    identifiers = read_identifiers_from_txt_file(file_path)
    if not identifiers:
        raise ValueError("File identifier kosong / tidak ada item valid.")

    map_df = string_map_identifiers(
        identifiers,
        species=string_species,
        base_url=string_base_url,
        caller_identity=caller_identity,
    )
    time.sleep(1.0)

    # pakai STRING IDs untuk query network (lebih cepat & lebih konsisten)
    if not map_df.empty and "stringId" in map_df.columns:
        ids_for_network = map_df["stringId"].dropna().astype(str).tolist()
    else:
        ids_for_network = identifiers

    edges_df = string_fetch_network(
        ids_for_network,
        species=string_species,
        required_score=required_score,
        add_nodes=add_nodes,
        base_url=string_base_url,
        caller_identity=caller_identity,
    )

    G = build_nx_graph_from_string_edges(edges_df)

    meta = {
        "source": source_name,
        "file_path": str(file_path),
        "n_input_identifiers": len(identifiers),
        "input_identifiers": identifiers[:50],
        "string_mapping": map_df.head(200).to_dict(orient="records") if isinstance(map_df, pd.DataFrame) else [],
        "string_edges": edges_df.head(200).to_dict(orient="records") if isinstance(edges_df, pd.DataFrame) else [],
        "n_nodes": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
    }
    return G, meta


def build_graph_from_lungadenocarcinoma_file(
    file_path: str,
    string_species: int = 9606,
    required_score: int = 700,
    add_nodes: int = 0,
    string_base_url: str = "https://string-db.org",
    caller_identity: str = "community-detect-streamlit",
) -> Tuple[nx.Graph, Dict[str, Any]]:
    """
    Bangun graf PPI lung adenocarcinoma dari file teks, lalu mapping ke STRING.
    """
    return build_graph_from_identifier_file(
        file_path=file_path,
        source_name="lungadenocarcinoma_txt_to_string",
        string_species=string_species,
        required_score=required_score,
        add_nodes=add_nodes,
        string_base_url=string_base_url,
        caller_identity=caller_identity,
    )


def build_graph_from_breastcancer_file(
    file_path: str,
    string_species: int = 9606,
    required_score: int = 700,
    add_nodes: int = 0,
    string_base_url: str = "https://string-db.org",
    caller_identity: str = "community-detect-streamlit",
) -> Tuple[nx.Graph, Dict[str, Any]]:
    """
    Bangun graf PPI breast cancer dari file teks hasil filter Open Targets,
    lalu mapping ke STRING.
    """
    return build_graph_from_identifier_file(
        file_path=file_path,
        source_name="breastcancer_txt_to_string",
        string_species=string_species,
        required_score=required_score,
        add_nodes=add_nodes,
        string_base_url=string_base_url,
        caller_identity=caller_identity,
    )

def build_graph_from_breastcancerdisgenet_file(
    file_path: str,
    string_species: int = 9606,
    required_score: int = 700,
    add_nodes: int = 0,
    string_base_url: str = "https://string-db.org",
    caller_identity: str = "community-detect-streamlit",
) -> Tuple[nx.Graph, Dict[str, Any]]:
    """
    Bangun graf PPI breast cancer dari file teks hasil ekspor DisGeNET,
    lalu mapping ke STRING.
    """
    return build_graph_from_identifier_file(
        file_path=file_path,
        source_name="breastcancer_disgenet_txt_to_string",
        string_species=string_species,
        required_score=required_score,
        add_nodes=add_nodes,
        string_base_url=string_base_url,
        caller_identity=caller_identity,
    )

