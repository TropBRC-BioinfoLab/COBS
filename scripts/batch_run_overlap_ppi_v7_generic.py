#!/usr/bin/env python3
"""
Batch run algoritma deteksi komunitas overlap pada dataset PPI dari file identifier
TANPA bergantung pada core.data dan TANPA paket requests.

Yang dipakai dari project:
- core.algos.registry.list_specs
- core.metrics.postprocess_and_memberships
- core.metrics.compute_evaluations

Graf dibangun langsung dari file identifier, misalnya breastcancer.txt, menggunakan urllib + STRING API,
dengan fallback UniProt untuk accession yang tidak langsung dikenali STRING.

Contoh:
    python batch_run_overlap_ppi_v7_generic.py --input-file breastcancer.txt --dataset-label breastcancer
    python batch_run_overlap_ppi_v7_generic.py --input-file breastcancer.txt --dataset-label breastcancer --include-cobs
    python batch_run_overlap_ppi_v7_generic.py --input-file breastcancer.txt --dataset-label breastcancer --include-keys slpa,kclique,cobs --include-cobs
    python batch_run_overlap_ppi_v7_generic.py --input-file lungadenocarcinoma.txt --dataset-label lungadenocarcinoma
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import traceback
import multiprocessing as mp
from dataclasses import asdict, is_dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError

import networkx as nx
import pandas as pd


# =========================
# Import project modules (tanpa core.data)
# =========================
try:
    from core.algos.registry import list_specs
    from core.metrics import postprocess_and_memberships, compute_evaluations
except Exception as e:
    raise RuntimeError(
        "Gagal mengimpor modul project yang dibutuhkan (registry/metrics). "
        "Jalankan script ini dari root project yang memiliki paket 'core'. "
        f"Detail: {type(e).__name__}: {e}"
    ) from e


# =========================
# Konstanta API
# =========================
STRING_API_URL = "https://string-db.org/api"
UNIPROT_API_URL = "https://rest.uniprot.org/uniprotkb/search"
CALLER_IDENTITY = "batch-overlap-ppi-generic"
NETWORK_TYPE = "functional"


# =========================
# Helpers umum
# =========================
def now_jakarta_iso() -> str:
    # timestamp lokal mesin; cukup untuk logging run
    return datetime.now().astimezone().isoformat(timespec="seconds")


def safe_jsonable(obj: Any):
    if is_dataclass(obj):
        return asdict(obj)
    if isinstance(obj, dict):
        return {str(k): safe_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [safe_jsonable(v) for v in obj]
    if isinstance(obj, Path):
        return str(obj)
    try:
        json.dumps(obj)
        return obj
    except Exception:
        return repr(obj)


def parse_csv_arg(text: str | None):
    if not text:
        return []
    return [x.strip().lower() for x in text.split(",") if x.strip()]


def chunked(seq: Iterable, size: int):
    seq = list(seq)
    for i in range(0, len(seq), size):
        yield seq[i:i + size]


def read_identifiers_txt(file_path: str | Path) -> list[str]:
    p = Path(file_path)
    if not p.exists():
        raise FileNotFoundError(f"File input tidak ditemukan: {p}")
    identifiers = [line.strip() for line in p.read_text(encoding="utf-8").splitlines()]
    identifiers = [x for x in identifiers if x]
    identifiers = list(dict.fromkeys(identifiers))
    if not identifiers:
        raise RuntimeError("File input kosong / tidak berisi identifier yang valid.")
    return identifiers


# =========================
# HTTP util (stdlib only)
# =========================
def _http_post_form(url: str, data: dict, timeout: int = 180) -> str:
    body = urlencode(data).encode("utf-8")
    req = Request(url, data=body, headers={
        "User-Agent": "Python urllib",
        "Content-Type": "application/x-www-form-urlencoded",
    })
    with urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8", errors="replace")


def _http_get_json(url: str, params: dict, timeout: int = 180) -> dict:
    full_url = url + "?" + urlencode(params)
    req = Request(full_url, headers={"User-Agent": "Python urllib"})
    with urlopen(req, timeout=timeout) as resp:
        return json.loads(resp.read().decode("utf-8", errors="replace"))


def string_api_post(method: str, params: dict, output_format: str = "tsv-no-header") -> str:
    url = "/".join([STRING_API_URL.rstrip("/"), output_format, method])
    return _http_post_form(url, params, timeout=180)


# =========================
# Graph builder (standalone)
# =========================
def map_to_string_ids_batched(identifiers: list[str], species: int = 9606, batch_size: int = 250) -> pd.DataFrame:
    rows = []
    batches = list(chunked(identifiers, batch_size))
    total = len(batches)
    for bi, batch in enumerate(batches, start=1):
        print(f"[mapping {bi}/{total}] direct STRING id mapping untuk {len(batch)} identifier...")
        params = {
            "identifiers": "\r".join(batch),
            "species": int(species),
            "limit": 1,
            "echo_query": 1,
            "caller_identity": CALLER_IDENTITY,
        }
        txt = string_api_post("get_string_ids", params, output_format="tsv-no-header")
        for line in txt.strip().split("\n"):
            if not line.strip():
                continue
            parts = line.split("\t")
            # queryItem, queryIndex, stringId, ncbiTaxonId, taxonName, preferredName, annotation
            if len(parts) < 6:
                continue
            rows.append({
                "queryItem": parts[0],
                "stringId": parts[2],
                "ncbiTaxonId": parts[3],
                "preferredName": parts[5],
                "resolvedQuery": parts[0],
                "mappingSource": "direct_string",
                "mappingPriority": 1,
            })
    return pd.DataFrame(rows)


def _quote_uniprot_query_value(value: str) -> str:
    """Quote nilai query UniProt agar aman untuk simbol seperti H3-3B."""
    value = str(value).strip().replace('\\', '\\\\').replace('"', '\\"')
    return f'"{value}"'


def _looks_like_uniprot_accession(value: str) -> bool:
    """
    Heuristik accession UniProt.
    Tujuannya bukan validasi sempurna, tetapi mencegah query accession:... untuk gene symbol
    yang jelas bukan accession dan bisa memicu HTTP 400.
    """
    x = str(value).strip().upper()
    patterns = [
        r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$',
        r'^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]$',
        r'^[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]-[0-9]+$',
    ]
    return any(re.match(p, x) for p in patterns)


def _extract_uniprot_rows(data: dict) -> list[dict]:
    rows = []
    for item in data.get("results", []) or []:
        accession = item.get("primaryAccession")
        entry_name = item.get("uniProtkbId")

        gene_primary = None
        for g in item.get("genes", []) or []:
            gname = ((g or {}).get("geneName") or {}).get("value")
            if gname:
                gene_primary = gname
                break

        protein_name = None
        pdsc = item.get("proteinDescription") or {}
        rn = ((pdsc.get("recommendedName") or {}).get("fullName") or {}).get("value")
        if rn:
            protein_name = rn
        if protein_name is None:
            sn = pdsc.get("submissionNames") or []
            if sn:
                protein_name = (((sn[0] or {}).get("fullName") or {}).get("value"))

        rows.append({
            "accession": accession,
            "gene_primary": gene_primary,
            "entry_name": entry_name,
            "protein_name": protein_name,
        })
    return rows


def _uniprot_candidate_queries(identifier: str, species: int) -> list[str]:
    """
    Identifier bisa berupa accession UniProt atau gene symbol.
    Untuk dataset breast cancer, mayoritas berupa gene symbol.
    Untuk dataset lung lama, bisa berupa accession.
    """
    ident = str(identifier).strip()
    q_ident = _quote_uniprot_query_value(ident)
    queries = []

    # Untuk accession yang tampak valid, coba accession lebih dahulu.
    if _looks_like_uniprot_accession(ident):
        queries.append(f'(accession:{ident}) AND organism_id:{int(species)}')

    # Untuk gene symbol, terutama breast cancer: H3-3B, MT-CYB, dll. harus di-quote.
    queries.extend([
        f'(gene_exact:{q_ident}) AND organism_id:{int(species)}',
        f'(gene:{q_ident}) AND organism_id:{int(species)}',
        f'(id:{q_ident}) AND organism_id:{int(species)}',
    ])

    # Fallback terakhir: full text quoted, lebih longgar tetapi aman.
    queries.append(f'({q_ident}) AND organism_id:{int(species)}')
    return queries


def query_uniprot_accessions(accessions: list[str], species: int = 9606, batch_size: int = 100) -> pd.DataFrame:
    """
    Fallback UniProt yang robust.

    Versi lama mengasumsikan semua unmapped identifier adalah accession dan membuat query
    accession:<identifier>. Pada dataset breast cancer, input adalah gene symbol; sebagian
    symbol yang tidak dikenali STRING dapat membuat UniProt mengembalikan HTTP 400.

    Strategi baru:
    - jika tampak seperti accession, coba accession dulu;
    - selalu coba gene_exact/gene/id/full-text dengan quoting;
    - jika satu identifier tetap gagal, jangan hentikan proses batch; catat saja sebagai unmapped.
    """
    rows = []
    if not accessions:
        return pd.DataFrame(columns=["accession", "gene_primary", "entry_name", "protein_name"])

    # Sengaja one-by-one: jumlah unmapped biasanya kecil, dan ini menghindari seluruh batch gagal karena 1 query buruk.
    total = len(accessions)
    for i, ident in enumerate(accessions, start=1):
        print(f"[uniprot {i}/{total}] lookup fallback untuk {ident}...")
        found = False
        last_error = None

        for query in _uniprot_candidate_queries(ident, species):
            params = {
                "query": query,
                "format": "json",
                "fields": "accession,gene_primary,id,protein_name",
                "size": 5,
            }
            try:
                data = _http_get_json(UNIPROT_API_URL, params, timeout=180)
                extracted = _extract_uniprot_rows(data)
                if extracted:
                    rows.extend(extracted)
                    found = True
                    break
            except HTTPError as e:
                last_error = f"HTTPError {e.code}: {e.reason}"
                continue
            except URLError as e:
                last_error = f"URLError: {e.reason}"
                continue
            except Exception as e:
                last_error = f"{type(e).__name__}: {e}"
                continue

        if not found:
            msg = f"[uniprot] tidak ditemukan / dilewati: {ident}"
            if last_error:
                msg += f" ({last_error})"
            print(msg)

    df = pd.DataFrame(rows)
    if df.empty:
        return pd.DataFrame(columns=["accession", "gene_primary", "entry_name", "protein_name"])
    return df.drop_duplicates(subset=["accession"])

def map_unmapped_via_uniprot_gene(unmapped_accessions: list[str], species: int = 9606, batch_size: int = 250):
    uni_df = query_uniprot_accessions(unmapped_accessions, species=species)
    if uni_df.empty:
        return (
            pd.DataFrame(columns=["queryItem", "stringId", "ncbiTaxonId", "preferredName", "resolvedQuery", "mappingSource", "mappingPriority"]),
            uni_df,
        )

    gene_df = uni_df.dropna(subset=["gene_primary"]).copy()
    if gene_df.empty:
        return (
            pd.DataFrame(columns=["queryItem", "stringId", "ncbiTaxonId", "preferredName", "resolvedQuery", "mappingSource", "mappingPriority"]),
            uni_df,
        )

    gene_map = map_to_string_ids_batched(
        gene_df["gene_primary"].tolist(),
        species=species,
        batch_size=batch_size,
    )
    if gene_map.empty:
        return (
            pd.DataFrame(columns=["queryItem", "stringId", "ncbiTaxonId", "preferredName", "resolvedQuery", "mappingSource", "mappingPriority"]),
            uni_df,
        )

    gene_map = gene_map.rename(columns={"queryItem": "gene_primary"})
    fallback_df = gene_df.merge(
        gene_map[["gene_primary", "stringId", "ncbiTaxonId", "preferredName"]],
        on="gene_primary",
        how="left",
    )
    fallback_df["queryItem"] = fallback_df["accession"]
    fallback_df["resolvedQuery"] = fallback_df["gene_primary"]
    fallback_df["mappingSource"] = fallback_df["stringId"].notna().map(
        lambda ok: "uniprot_gene_fallback" if ok else "uniprot_lookup_only"
    )
    fallback_df["mappingPriority"] = fallback_df["stringId"].notna().map(lambda ok: 2 if ok else 99)

    cols = ["queryItem", "stringId", "ncbiTaxonId", "preferredName", "resolvedQuery", "mappingSource", "mappingPriority"]
    return fallback_df[cols].copy(), uni_df


def fetch_network_edges(string_ids: list[str], species: int, required_score: int) -> pd.DataFrame:
    params = {
        "identifiers": "\r".join(string_ids),
        "species": int(species),
        "required_score": int(required_score),
        "network_type": NETWORK_TYPE,
        "caller_identity": CALLER_IDENTITY,
    }
    txt = string_api_post("network", params, output_format="tsv-no-header")
    rows = []
    for line in txt.strip().split("\n"):
        if not line.strip():
            continue
        p = line.split("\t")
        if len(p) < 6:
            continue
        rows.append({
            "stringId_A": p[0],
            "stringId_B": p[1],
            "preferredName_A": p[2],
            "preferredName_B": p[3],
            "ncbiTaxonId": p[4],
            "score_raw": p[5],
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df["score_raw"] = pd.to_numeric(df["score_raw"], errors="coerce")
        df["score"] = df["score_raw"].apply(lambda x: x / 1000.0 if pd.notna(x) and x > 1.0 else x)
    return df


def fetch_partners_edges_one_by_one(query_ids: list[str], species: int, required_score: int,
                                    limit_per_protein: int = 1000, sleep_sec: float = 0.15) -> pd.DataFrame:
    all_rows = []
    total = len(query_ids)
    for i, q in enumerate(query_ids, start=1):
        print(f"[partners {i}/{total}] mengambil partner untuk {q}...")
        params = {
            "identifiers": q,
            "species": int(species),
            "required_score": int(required_score),
            "limit": int(limit_per_protein),
            "network_type": NETWORK_TYPE,
            "caller_identity": CALLER_IDENTITY,
        }
        txt = string_api_post("interaction_partners", params, output_format="tsv-no-header")
        for line in txt.strip().split("\n"):
            if not line.strip():
                continue
            p = line.split("\t")
            if len(p) < 6:
                continue
            score_raw = pd.to_numeric(p[5], errors="coerce")
            score = (score_raw / 1000.0) if (pd.notna(score_raw) and score_raw > 1.0) else score_raw
            all_rows.append({
                "stringId_A": p[0],
                "stringId_B": p[1],
                "preferredName_A": p[2],
                "preferredName_B": p[3],
                "ncbiTaxonId": p[4],
                "score": score,
            })
        if sleep_sec and sleep_sec > 0:
            time.sleep(float(sleep_sec))
    return pd.DataFrame(all_rows)


def choose_added_nodes(input_string_ids: list[str], species: int, required_score: int, add_nodes: int,
                       partner_limit_per_protein: int = 200, sleep_sec: float = 0.10) -> list[str]:
    if add_nodes <= 0:
        return []
    print(f"Menentukan {add_nodes} node tambahan dari partner STRING...")
    partner_df = fetch_partners_edges_one_by_one(
        input_string_ids,
        species=species,
        required_score=required_score,
        limit_per_protein=partner_limit_per_protein,
        sleep_sec=sleep_sec,
    )
    if partner_df.empty:
        return []

    input_set = set(input_string_ids)
    cand_rows = []
    for _, r in partner_df.iterrows():
        a = r.get("stringId_A")
        b = r.get("stringId_B")
        s = r.get("score")
        if pd.isna(s):
            continue
        if a in input_set and b not in input_set:
            cand_rows.append((b, float(s)))
        elif b in input_set and a not in input_set:
            cand_rows.append((a, float(s)))

    if not cand_rows:
        return []

    cand_df = pd.DataFrame(cand_rows, columns=["stringId", "score"])
    cand_df = cand_df.groupby("stringId", as_index=False)["score"].max().sort_values("score", ascending=False)
    added = cand_df.head(add_nodes)["stringId"].tolist()
    print(f"Node tambahan terpilih: {len(added)}")
    return added


def build_graph_from_identifier_file_local(file_path: str | Path,
                                                   string_species: int = 9606,
                                                   required_score: int = 400,
                                                   add_nodes: int = 0):
    identifiers = read_identifiers_txt(file_path)
    print(f"Identifier input unik: {len(identifiers)}")
    print(f"Contoh awal: {identifiers[:10]}")

    mapping_direct_df = map_to_string_ids_batched(identifiers, species=string_species, batch_size=250)
    if mapping_direct_df.empty:
        print("Peringatan: direct mapping STRING kosong. Lanjut fallback UniProt untuk semua identifier.")
        direct_mapped_items = set()
        unmapped_after_direct = list(identifiers)
    else:
        direct_mapped_items = set(mapping_direct_df["queryItem"].tolist())
        unmapped_after_direct = [x for x in identifiers if x not in direct_mapped_items]

    fallback_df, uniprot_fallback_df = map_unmapped_via_uniprot_gene(
        unmapped_after_direct,
        species=string_species,
        batch_size=250,
    )

    mapping_frames = []
    if not mapping_direct_df.empty:
        mapping_frames.append(mapping_direct_df)
    if not fallback_df.empty:
        mapping_frames.append(fallback_df)

    if mapping_frames:
        mapping_audit_df = pd.concat(mapping_frames, ignore_index=True, sort=False)
        mapping_audit_df = mapping_audit_df.sort_values(["queryItem", "mappingPriority", "preferredName", "stringId"])
        mapping_audit_df = mapping_audit_df.drop_duplicates(subset=["queryItem"], keep="first")
    else:
        mapping_audit_df = pd.DataFrame(columns=[
            "queryItem", "stringId", "ncbiTaxonId", "preferredName",
            "resolvedQuery", "mappingSource", "mappingPriority"
        ])

    resolved_original_df = mapping_audit_df.dropna(subset=["stringId"]).copy()
    resolved_items = set(resolved_original_df["queryItem"].tolist())
    unmapped = [x for x in identifiers if x not in resolved_items]

    if resolved_original_df.empty:
        raise RuntimeError(
            "Tidak ada identifier yang berhasil di-resolve ke STRING, bahkan setelah fallback UniProt. "
            "Cek koneksi internet, response STRING/UniProt, atau format identifier input."
        )

    mapped = resolved_original_df.sort_values(["mappingPriority", "queryItem"]).drop_duplicates(subset=["stringId"], keep="first").copy()
    collapsed_originals = len(resolved_original_df) - len(mapped)

    print(f"Mapped langsung ke STRING  : {len(direct_mapped_items)}")
    print(f"Tambahan dari fallback     : {len(resolved_items - direct_mapped_items)}")
    print(f"Resolved accession total   : {len(resolved_items)}")
    print(f"Unique STRING proteins     : {len(mapped)}")
    print(f"Collapsed (same STRING ID) : {collapsed_originals}")
    print(f"Unmapped tersisa           : {len(unmapped)}")
    if unmapped:
        print("Contoh unmapped:", unmapped[:20])

    preferred_counts = mapped["preferredName"].value_counts(dropna=False).to_dict()
    mapped["node_label"] = mapped.apply(
        lambda r: r["preferredName"] if preferred_counts.get(r["preferredName"], 0) == 1 else f"{r['preferredName']}|{r['stringId']}",
        axis=1,
    )
    string_to_node_label = dict(zip(mapped["stringId"], mapped["node_label"]))

    input_string_ids = mapped["stringId"].tolist()
    added_string_ids = choose_added_nodes(input_string_ids, string_species, required_score, add_nodes)
    all_query_string_ids = list(dict.fromkeys(input_string_ids + added_string_ids))
    all_query_string_id_set = set(all_query_string_ids)

    edges_df = pd.DataFrame()
    try:
        print("Mencoba Mode A: /network (sekali panggil)...")
        edges_df = fetch_network_edges(all_query_string_ids, species=string_species, required_score=required_score)
        print(f"Edges retrieved (raw): {len(edges_df)}")
    except Exception as e:
        print("Mode A gagal:", repr(e))

    if edges_df.empty:
        print("Pindah ke Mode B: /interaction_partners (per protein, lebih lama)...")
        edges_df = fetch_partners_edges_one_by_one(
            all_query_string_ids,
            species=string_species,
            required_score=required_score,
            limit_per_protein=max(1000, add_nodes + 50),
            sleep_sec=0.10,
        )
        print(f"Edges retrieved (raw): {len(edges_df)}")

    if edges_df.empty:
        raise RuntimeError("Tidak ada edge yang didapat dari STRING. Cek mapping / parameter / koneksi.")

    edges_df = edges_df.dropna(subset=["stringId_A", "stringId_B", "score"]).copy()
    edges_df = edges_df[
        edges_df["stringId_A"].isin(all_query_string_id_set) &
        edges_df["stringId_B"].isin(all_query_string_id_set)
    ].copy()
    edges_df = edges_df[edges_df["stringId_A"] != edges_df["stringId_B"]].copy()

    # pastikan node tambahan juga punya label aman
    if added_string_ids:
        add_map_df = map_to_string_ids_batched(added_string_ids, species=string_species, batch_size=250)
        if not add_map_df.empty:
            add_map_df = add_map_df.drop_duplicates(subset=["stringId"])
            for _, r in add_map_df.iterrows():
                sid = r["stringId"]
                if sid not in string_to_node_label:
                    pname = r.get("preferredName") or sid
                    string_to_node_label[sid] = pname

    edges_df["node_A"] = edges_df["stringId_A"].map(string_to_node_label)
    edges_df["node_B"] = edges_df["stringId_B"].map(string_to_node_label)
    edges_df = edges_df.dropna(subset=["node_A", "node_B"]).copy()

    a = edges_df["node_A"].astype(str)
    b = edges_df["node_B"].astype(str)
    edges_df["u"] = a.where(a <= b, b)
    edges_df["v"] = b.where(a <= b, a)
    edges_df = edges_df.groupby(["u", "v"], as_index=False)["score"].max()
    edges_df = edges_df.rename(columns={"u": "node_A", "v": "node_B"})

    print(f"Edges after filtering/dedup: {len(edges_df)}")

    G = nx.Graph()
    for _, r in mapped.iterrows():
        node = r["node_label"]
        G.add_node(
            node,
            stringId=r.get("stringId"),
            preferredName=r.get("preferredName"),
            original_query=r.get("queryItem"),
            mapping_source=r.get("mappingSource"),
            is_input=True,
        )

    for sid in added_string_ids:
        node = string_to_node_label.get(sid, sid)
        if node not in G:
            G.add_node(node, stringId=sid, preferredName=node, original_query=None, mapping_source="added_partner", is_input=False)

    for _, r in edges_df.iterrows():
        G.add_edge(r["node_A"], r["node_B"], weight=float(r["score"]), score=float(r["score"]))

    meta = {
        "input_file": str(file_path),
        "n_identifiers": len(identifiers),
        "n_direct_mapped": len(direct_mapped_items),
        "n_fallback_added": len(resolved_items - direct_mapped_items),
        "n_resolved_accessions": len(resolved_items),
        "n_unique_string_proteins": len(mapped),
        "n_unmapped": len(unmapped),
        "unmapped": unmapped,
        "collapsed_same_string_id": collapsed_originals,
        "required_score": required_score,
        "string_species": string_species,
        "add_nodes": add_nodes,
        "graph_n_nodes": G.number_of_nodes(),
        "graph_n_edges": G.number_of_edges(),
        "string_mapping": mapped.to_dict(orient="records"),
        "string_mapping_audit": mapping_audit_df.to_dict(orient="records"),
        "string_edges": edges_df.to_dict(orient="records"),
        "uniprot_records": uniprot_fallback_df.to_dict(orient="records") if not uniprot_fallback_df.empty else [],
        "added_string_ids": added_string_ids,
    }
    return G, meta


# =========================
# Evaluasi/modularity helpers
# =========================
def infer_random_seed(params: dict):
    for k in ("seed", "random_seed", "random_state", "rng_seed"):
        if k in params:
            return params[k]
    return None


def infer_weight_attr(G: nx.Graph):
    try:
        for _, _, d in G.edges(data=True):
            if not isinstance(d, dict):
                continue
            if "weight" in d:
                return "weight"
            if "score" in d:
                return "score"
    except Exception:
        pass
    return None


def crisp_partition_from_comms(G: nx.Graph, comms):
    covered = set()
    communities = []
    overlap = False
    for c in (comms or []):
        try:
            cset = {n for n in c if n in G}
        except TypeError:
            continue
        if not cset:
            continue
        if covered.intersection(cset):
            overlap = True
            break
        covered.update(cset)
        communities.append(cset)
    if overlap:
        return None
    for n in G.nodes():
        if n not in covered:
            communities.append({n})
    if not communities:
        communities = [{n} for n in G.nodes()]
    return communities


def crisp_partition_via_primary(G: nx.Graph, comms):
    comm_sets = []
    for c in (comms or []):
        try:
            cset = {n for n in c if n in G}
        except TypeError:
            continue
        if cset:
            comm_sets.append(cset)
    if not comm_sets:
        return [{n} for n in G.nodes()]
    sizes = [len(c) for c in comm_sets]
    node_to_cids = {n: [] for n in G.nodes()}
    for cid, cset in enumerate(comm_sets):
        for n in cset:
            node_to_cids[n].append(cid)
    assignment = {}
    for n, cids in node_to_cids.items():
        if not cids:
            assignment[n] = None
        else:
            assignment[n] = max(cids, key=lambda cid: (sizes[cid], -cid))
    cid_to_nodes = {}
    for n, cid in assignment.items():
        if cid is None:
            continue
        cid_to_nodes.setdefault(cid, set()).add(n)
    part = [cid_to_nodes[cid] for cid in sorted(cid_to_nodes.keys()) if cid_to_nodes[cid]]
    for n, cid in assignment.items():
        if cid is None:
            part.append({n})
    return part if part else [{n} for n in G.nodes()]


def compute_modularity_crisp(G: nx.Graph, comms, resolution: float = 1.0, allow_overlap_projection: bool = True):
    part = crisp_partition_from_comms(G, comms)
    if part is None:
        if not allow_overlap_projection:
            return None
        part = crisp_partition_via_primary(G, comms)
    H = G if not getattr(G, "is_directed", lambda: False)() else G.to_undirected()
    w = infer_weight_attr(H)
    try:
        return nx.algorithms.community.quality.modularity(H, part, weight=w, resolution=resolution)
    except Exception:
        return None


def get_components_info(G: nx.Graph):
    if G.number_of_nodes() == 0:
        return {"n_components": 0, "is_connected": True, "gc_nodes": set(), "gc_size": 0}
    UG = G.to_undirected(as_view=True)
    n_components = nx.number_connected_components(UG)
    is_connected = (n_components == 1)
    gc_nodes = max(nx.connected_components(UG), key=len) if not is_connected else set(G.nodes())
    gc_size = len(gc_nodes)
    return {"n_components": n_components, "is_connected": is_connected, "gc_nodes": gc_nodes, "gc_size": gc_size}


def decide_graph_scope(G_full: nx.Graph, use_giant: str):
    info = get_components_info(G_full)
    if use_giant not in {"auto", "true", "false"}:
        raise ValueError("use_giant harus salah satu dari: auto, true, false")
    if use_giant == "true":
        use_gc = True
    elif use_giant == "false":
        use_gc = False
    else:
        use_gc = not info["is_connected"]
    if use_gc and (not info["is_connected"]):
        G = G_full.subgraph(info["gc_nodes"]).copy()
        scope_label = "GC"
    else:
        G = G_full
        scope_label = "FULL"
    return G, info, scope_label, use_gc


def params_from_spec(spec) -> dict:
    params = {}
    for p in getattr(spec, "params", []) or []:
        params[p.key] = getattr(p, "default", None)
    if (getattr(spec, "key", "") or "").lower() == "cobs" and "sim_metric" not in params:
        params["sim_metric"] = "cosine"
    return params


def specs_overlap_only(include_cobs: bool = False, include_keys=None, exclude_keys=None):
    include_keys = set(include_keys or [])
    exclude_keys = {k.lower() for k in (exclude_keys or [])}
    specs = []
    for s in list_specs():
        key = (getattr(s, "key", "") or "").lower()
        name = (getattr(s, "name", "") or "").lower()
        desc = (getattr(s, "description", "") or "").lower()

        if key in {"dgm_crisp", "greedy_modularity"}:
            continue
        is_dgm = ("dgm" in key) or ("dgm" in name)
        is_overlap = ("overlap" in key) or ("overlap" in name) or ("overlap" in desc) or ("overlapping" in desc)
        if is_dgm and is_overlap:
            continue
        if key == "cobs" and not include_cobs:
            continue
        if include_keys and key not in include_keys:
            continue
        if key in exclude_keys:
            continue
        specs.append(s)
    return specs


def communities_to_rows(algo_key: str, algo_name: str, communities):
    rows = []
    for cid, comm in enumerate(communities or [], start=1):
        try:
            members = sorted(list(comm))
        except Exception:
            members = list(comm)
        for n in members:
            rows.append({"algo_key": algo_key, "algo_name": algo_name, "community_id": cid, "node": n})
    return rows


def memberships_to_rows(algo_key: str, algo_name: str, memberships):
    rows = []
    for node, mids in (memberships or {}).items():
        try:
            mids_list = list(mids)
        except Exception:
            mids_list = [mids]
        rows.append({
            "algo_key": algo_key,
            "algo_name": algo_name,
            "node": node,
            "n_memberships": len(mids_list),
            "memberships": json.dumps(mids_list, ensure_ascii=False),
        })
    return rows


# =========================
# Timeout execution helpers
# =========================
def _normalize_raw_communities(raw_comms):
    norm = []
    for c in (raw_comms or []):
        try:
            members = list(c)
        except Exception:
            continue
        norm.append(members)
    return norm




def _normalize_memberships(memberships):
    """
    Normalisasi memberships agar aman untuk berbagai bentuk output:
    - {node: [cid1, cid2]}
    - {node: set(...)}
    - {node: tuple(...)}
    - {node: int}
    - {node: numpy scalar}
    - list/tuple pasangan (node, membership)
    """
    if memberships is None:
        return {}

    # Bentuk dict adalah yang paling umum
    if isinstance(memberships, dict):
        out = {}
        for k, v in memberships.items():
            if isinstance(v, (list, tuple, set)):
                out[str(k)] = list(v)
            elif v is None:
                out[str(k)] = []
            else:
                # kasus v berupa int / numpy scalar / string tunggal
                out[str(k)] = [v]
        return out

    # fallback jika bentuknya iterable pasangan
    out = {}
    try:
        for item in memberships:
            if isinstance(item, (list, tuple)) and len(item) == 2:
                k, v = item
                if isinstance(v, (list, tuple, set)):
                    out[str(k)] = list(v)
                elif v is None:
                    out[str(k)] = []
                else:
                    out[str(k)] = [v]
    except Exception:
        return {}
    return out

def _sanitize_evals(evals: dict | None) -> dict:
    evals = evals or {}
    out = {}
    for k, v in evals.items():
        try:
            if pd.isna(v):
                out[k] = float("nan")
            elif hasattr(v, "item"):
                out[k] = v.item()
            else:
                out[k] = v
        except Exception:
            out[k] = safe_jsonable(v)
    return out


def _execute_algorithm_pipeline(G: nx.Graph, spec, params: dict, min_size: int):
    clustering = spec.run(G, params)
    raw_comms = getattr(clustering, "communities", None) or []

    covered = set()
    for c in raw_comms:
        try:
            covered.update([n for n in c if n in G])
        except TypeError:
            continue
    uncovered_nodes = [n for n in G.nodes() if n not in covered]

    comms, memberships, uncovered_after_min_size = postprocess_and_memberships(
        G, raw_comms, min_size=min_size
    )

    # Penting: evaluasi memakai object clustering asli, bukan object tiruan.
    clustering.communities = comms
    try:
        clustering.overlap = True
    except Exception:
        pass
    evals = _sanitize_evals(compute_evaluations(clustering) or {})

    return {
        "raw_communities": _normalize_raw_communities(raw_comms),
        "communities": _normalize_raw_communities(comms),
        "memberships": _normalize_memberships(memberships),
        "uncovered_nodes": list(uncovered_nodes),
        "uncovered_after_min_size": list(uncovered_after_min_size),
        "evals": evals,
    }


def _run_algo_worker(queue, algo_key: str, G: nx.Graph, params: dict, min_size: int):
    try:
        # import ulang di child process agar aman pada Windows/spawn
        from core.algos.registry import list_specs as _list_specs

        spec = None
        for s in _list_specs():
            if (getattr(s, "key", "") or "").lower() == (algo_key or "").lower():
                spec = s
                break
        if spec is None:
            raise RuntimeError(f"Spec algoritma tidak ditemukan: {algo_key}")

        result = _execute_algorithm_pipeline(G, spec, params, min_size=min_size)
        queue.put({
            "ok": True,
            **result,
        })
    except Exception as e:
        try:
            queue.put({
                "ok": False,
                "error_type": type(e).__name__,
                "error_message": str(e),
                "traceback": traceback.format_exc(),
            })
        except Exception:
            pass
    finally:
        try:
            queue.close()
        except Exception:
            pass
        try:
            queue.join_thread()
        except Exception:
            pass


def run_algorithm_with_optional_timeout(
    G: nx.Graph,
    spec,
    params: dict,
    timeout_sec: int | float | None = None,
    min_size: int = 1,
):
    timeout_sec = float(timeout_sec or 0)
    if timeout_sec <= 0:
        try:
            result = _execute_algorithm_pipeline(G, spec, params, min_size=min_size)
            return {
                "status": "ok",
                **result,
            }
        except Exception as e:
            return {
                "status": "error",
                "error_type": type(e).__name__,
                "error_message": str(e),
                "traceback": traceback.format_exc(),
            }

    ctx = mp.get_context("spawn")
    queue = ctx.Queue()
    proc = ctx.Process(
        target=_run_algo_worker,
        args=(queue, getattr(spec, "key", ""), G, params, min_size),
    )
    proc.start()

    result = None
    try:
        # Penting: baca hasil dari queue dulu. Ini mencegah deadlock/timeout palsu
        # saat child process sedang flush payload besar ke Queue pada Windows/spawn.
        result = queue.get(timeout=timeout_sec)
    except Exception:
        result = None

    if result is None:
        if proc.is_alive():
            proc.terminate()
            proc.join(5)
            if proc.is_alive():
                try:
                    proc.kill()
                except Exception:
                    pass
                proc.join(5)
        try:
            queue.close()
        except Exception:
            pass
        try:
            queue.join_thread()
        except Exception:
            pass
        return {
            "status": "timeout",
            "error_type": "TimeoutError",
            "error_message": f"Melebihi batas waktu {timeout_sec:.0f} detik",
            "traceback": "",
        }

    # Hasil sudah diterima; sekarang rapikan process/queue.
    proc.join(5)
    if proc.is_alive():
        proc.terminate()
        proc.join(5)
        if proc.is_alive():
            try:
                proc.kill()
            except Exception:
                pass
            proc.join(5)

    try:
        queue.close()
    except Exception:
        pass
    try:
        queue.join_thread()
    except Exception:
        pass

    if not result.get("ok"):
        return {
            "status": "error",
            "error_type": result.get("error_type", "RuntimeError"),
            "error_message": result.get("error_message", "Unknown worker error"),
            "traceback": result.get("traceback", ""),
        }

    return {
        "status": "ok",
        "raw_communities": [set(c) for c in result.get("raw_communities", [])],
        "communities": [set(c) for c in result.get("communities", [])],
        "memberships": result.get("memberships", {}) or {},
        "uncovered_nodes": result.get("uncovered_nodes", []) or [],
        "uncovered_after_min_size": result.get("uncovered_after_min_size", []) or [],
        "evals": result.get("evals", {}) or {},
    }


# =========================
# Main
# =========================
def main():
    ap = argparse.ArgumentParser(description="Batch run overlap community detection on a PPI graph from a gene or protein identifier file")
    ap.add_argument("--input-file", default="breastcancer.txt", help="Path file identifier, misalnya breastcancer.txt atau lungadenocarcinoma.txt")
    ap.add_argument("--dataset-label", default="breastcancer", help="Label dataset untuk log dan output, misalnya breastcancer atau lungadenocarcinoma")
    ap.add_argument("--species", type=int, default=9606, help="NCBI taxon / species")
    ap.add_argument("--required-score", type=int, default=400, help="STRING required score (0-999)")
    ap.add_argument("--add-nodes", type=int, default=0, help="Tambahkan top partner STRING di luar input")
    ap.add_argument("--use-giant", choices=["auto", "true", "false"], default="auto", help="Gunakan giant component")
    ap.add_argument("--min-size", type=int, default=1, help="Min ukuran komunitas untuk postprocess")
    ap.add_argument("--include-cobs", action="store_true", help="Ikut jalankan COBS juga")
    ap.add_argument("--include-keys", default="", help="Batasi ke key algoritma tertentu, dipisah koma")
    ap.add_argument("--exclude-keys", default="lfm", help="Kecualikan key algoritma tertentu, dipisah koma. Default: lfm")
    ap.add_argument("--timeout-sec", type=int, default=300, help="Batas waktu per algoritma (detik). 0=tanpa timeout")
    ap.add_argument("--outdir", default="outputs_breastcancer_batch", help="Folder output")
    args = ap.parse_args()

    include_keys = parse_csv_arg(args.include_keys)
    exclude_keys = parse_csv_arg(args.exclude_keys)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    dataset_slug = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in str(args.dataset_label).strip()) or "dataset"
    run_dir = outdir / f"{dataset_slug}_run_{run_id}"
    run_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print(f"MEMBANGUN GRAF {args.dataset_label.upper()}")
    print("=" * 80)
    print(f"Dataset label   : {args.dataset_label}")
    print(f"Input file      : {args.input_file}")
    print(f"Species         : {args.species}")
    print(f"required_score  : {args.required_score}")
    print(f"add_nodes       : {args.add_nodes}")
    print(f"use_giant       : {args.use_giant}")
    print(f"Output dir      : {run_dir}")
    print(f"timeout_sec     : {args.timeout_sec}")

    t0 = time.perf_counter()
    G_full, meta = build_graph_from_identifier_file_local(
        file_path=args.input_file,
        string_species=int(args.species),
        required_score=int(args.required_score),
        add_nodes=int(args.add_nodes),
    )
    graph_build_sec = time.perf_counter() - t0

    G, comp_info, scope_label, use_gc = decide_graph_scope(G_full, args.use_giant)
    print(f"Graph FULL       : |V|={G_full.number_of_nodes()} |E|={G_full.number_of_edges()}")
    print(f"Connected comps  : {comp_info['n_components']}")
    print(f"Giant component  : {comp_info['gc_size']} node")
    print(f"Scope dipakai    : {scope_label}")
    print(f"Graph aktif      : |V|={G.number_of_nodes()} |E|={G.number_of_edges()}")
    print(f"Build time       : {graph_build_sec:.2f} detik")

    nx.write_edgelist(G, run_dir / "graph_active.edgelist", data=True)
    with open(run_dir / "graph_meta.json", "w", encoding="utf-8") as f:
        json.dump(safe_jsonable(meta), f, ensure_ascii=False, indent=2)

    if isinstance(meta, dict):
        for key, filename in [
            ("uniprot_records", "meta_uniprot_records.csv"),
            ("string_mapping", "meta_string_mapping.csv"),
            ("string_mapping_audit", "meta_string_mapping_audit.csv"),
            ("string_edges", "meta_string_edges.csv"),
        ]:
            if key in meta:
                try:
                    pd.DataFrame(meta[key]).to_csv(run_dir / filename, index=False)
                except Exception:
                    pass

    specs = specs_overlap_only(
        include_cobs=args.include_cobs,
        include_keys=include_keys,
        exclude_keys=exclude_keys,
    )
    if not specs:
        raise RuntimeError("Tidak ada algoritma yang dipilih setelah filter include/exclude.")

    print("\n" + "=" * 80)
    print("ALGORITMA YANG AKAN DIJALANKAN")
    print("=" * 80)
    for i, s in enumerate(specs, start=1):
        print(f"{i:>2}. {s.key} -> {s.name}")

    summary_rows = []
    community_rows = []
    membership_rows = []
    error_rows = []

    total = len(specs)
    for idx, spec in enumerate(specs, start=1):
        algo_key = getattr(spec, "key", f"algo_{idx}")
        algo_name = getattr(spec, "name", algo_key)
        params = params_from_spec(spec)
        seed = infer_random_seed(params)

        print("\n" + "-" * 80)
        print(f"[{idx}/{total}] Menjalankan {algo_name} ({algo_key})")
        print(f"[{idx}/{total}] Params: {params}")

        t1 = time.perf_counter()
        try:
            run_result = run_algorithm_with_optional_timeout(
                G,
                spec,
                params,
                timeout_sec=args.timeout_sec,
                min_size=args.min_size,
            )
            if run_result.get("status") == "timeout":
                raise TimeoutError(run_result.get("error_message", "Timeout"))
            if run_result.get("status") != "ok":
                raise RuntimeError(run_result.get("error_message", "Unknown worker error"))

            raw_comms = run_result.get("raw_communities", []) or []
            comms = run_result.get("communities", []) or []
            memberships = run_result.get("memberships", {}) or {}
            uncovered_nodes = run_result.get("uncovered_nodes", []) or []
            uncovered_after_min_size = run_result.get("uncovered_after_min_size", []) or []
            evals = run_result.get("evals", {}) or {}

            gamma = float(params.get("resolution", 1.0))
            mod_crisp = compute_modularity_crisp(G, raw_comms, resolution=gamma, allow_overlap_projection=True)
            evals["modularity"] = float(mod_crisp) if mod_crisp is not None else float("nan")

            runtime_sec = time.perf_counter() - t1
            n_communities = len(comms)
            max_community_size = max((len(c) for c in comms), default=0)

            summary = {
                "run_id": run_id,
                "timestamp": now_jakarta_iso(),
                "status": "ok",
                "algo_key": algo_key,
                "algo_name": algo_name,
                "runtime_sec": runtime_sec,
                "random_seed": seed,
                "n_nodes": G.number_of_nodes(),
                "n_edges": G.number_of_edges(),
                "n_communities_raw": len(raw_comms),
                "n_communities": n_communities,
                "max_community_size_post": max_community_size,
                "n_uncovered_raw": len(uncovered_nodes),
                "n_uncovered_post": len(uncovered_after_min_size),
                "params_json": json.dumps(safe_jsonable(params), ensure_ascii=False),
            }
            summary.update({k: evals.get(k) for k in [
                "shen_modularity",
                "modularity_overlap",
                "modularity",
                "internal_edge_density",
                "conductance",
                "normalized_cut",
                "avg_transitivity",
                "n_communities",
                "max_community_size",
                "coverage_ratio",
                "uncovered_node_ratio",
                "overlap_node_ratio",
                "avg_memberships_per_covered_node",
                "overlap_intensity",
            ]})
            summary_rows.append(summary)

            community_rows.extend(communities_to_rows(algo_key, algo_name, comms))
            membership_rows.extend(memberships_to_rows(algo_key, algo_name, memberships))

            algo_slug = f"{idx:02d}_{algo_key}"
            with open(run_dir / f"communities_{algo_slug}.json", "w", encoding="utf-8") as f:
                json.dump(
                    {
                        "algo_key": algo_key,
                        "algo_name": algo_name,
                        "params": safe_jsonable(params),
                        "communities": [sorted(list(c)) for c in comms],
                    },
                    f,
                    ensure_ascii=False,
                    indent=2,
                )
            pd.DataFrame(communities_to_rows(algo_key, algo_name, comms)).to_csv(
                run_dir / f"communities_{algo_slug}.csv", index=False
            )
            pd.DataFrame(memberships_to_rows(algo_key, algo_name, memberships)).to_csv(
                run_dir / f"memberships_{algo_slug}.csv", index=False
            )

            print(
                f"[{idx}/{total}] Selesai {algo_name} | runtime={runtime_sec:.2f}s | "
                f"n_comms={n_communities} | modularity_overlap={evals.get('modularity_overlap')} | "
                f"coverage_ratio={evals.get('coverage_ratio')}"
            )

        except Exception as e:
            runtime_sec = time.perf_counter() - t1
            tb = traceback.format_exc()
            status_value = "timeout" if isinstance(e, TimeoutError) else "error"
            error_rows.append({
                "run_id": run_id,
                "timestamp": now_jakarta_iso(),
                "algo_key": algo_key,
                "algo_name": algo_name,
                "runtime_sec": runtime_sec,
                "params_json": json.dumps(safe_jsonable(params), ensure_ascii=False),
                "error_type": type(e).__name__,
                "error_message": str(e),
                "traceback": tb,
            })
            summary_rows.append({
                "run_id": run_id,
                "timestamp": now_jakarta_iso(),
                "status": status_value,
                "algo_key": algo_key,
                "algo_name": algo_name,
                "runtime_sec": runtime_sec,
                "random_seed": seed,
                "n_nodes": G.number_of_nodes(),
                "n_edges": G.number_of_edges(),
                "params_json": json.dumps(safe_jsonable(params), ensure_ascii=False),
                "error_type": type(e).__name__,
                "error_message": str(e),
            })
            print(f"[{idx}/{total}] GAGAL {algo_name} | {type(e).__name__}: {e}")
            continue

    summary_df = pd.DataFrame(summary_rows)
    errors_df = pd.DataFrame(error_rows)
    communities_df = pd.DataFrame(community_rows)
    memberships_df = pd.DataFrame(membership_rows)

    summary_csv = run_dir / "summary_metrics.csv"
    errors_csv = run_dir / "errors.csv"
    communities_csv = run_dir / "all_communities_long.csv"
    memberships_csv = run_dir / "all_memberships_long.csv"

    summary_df.to_csv(summary_csv, index=False)
    errors_df.to_csv(errors_csv, index=False)
    communities_df.to_csv(communities_csv, index=False)
    memberships_df.to_csv(memberships_csv, index=False)

    try:
        xlsx_path = run_dir / "batch_results.xlsx"
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            summary_df.to_excel(writer, sheet_name="summary", index=False)
            errors_df.to_excel(writer, sheet_name="errors", index=False)
            communities_df.to_excel(writer, sheet_name="communities", index=False)
            memberships_df.to_excel(writer, sheet_name="memberships", index=False)
    except Exception:
        pass

    metadata = {
        "run_id": run_id,
        "timestamp": now_jakarta_iso(),
        "dataset_label": args.dataset_label,
        "input_file": args.input_file,
        "species": args.species,
        "required_score": args.required_score,
        "add_nodes": args.add_nodes,
        "min_size": args.min_size,
        "use_giant_arg": args.use_giant,
        "use_gc_effective": use_gc,
        "scope_label": scope_label,
        "graph_full": {"n_nodes": G_full.number_of_nodes(), "n_edges": G_full.number_of_edges()},
        "graph_active": {"n_nodes": G.number_of_nodes(), "n_edges": G.number_of_edges()},
        "components": {
            "n_components": comp_info["n_components"],
            "is_connected": comp_info["is_connected"],
            "gc_size": comp_info["gc_size"],
        },
        "algorithms_requested": [getattr(s, "key", None) for s in specs],
        "n_algorithms": len(specs),
        "n_success": int((summary_df.get("status") == "ok").sum()) if not summary_df.empty else 0,
        "n_error": int((summary_df.get("status") == "error").sum()) if not summary_df.empty else 0,
        "files": {
            "summary_csv": str(summary_csv),
            "errors_csv": str(errors_csv),
            "communities_csv": str(communities_csv),
            "memberships_csv": str(memberships_csv),
        },
    }
    with open(run_dir / "run_metadata.json", "w", encoding="utf-8") as f:
        json.dump(metadata, f, ensure_ascii=False, indent=2)

    print("\n" + "=" * 80)
    print("SELESAI")
    print("=" * 80)
    print(f"Summary       : {summary_csv}")
    print(f"Errors        : {errors_csv}")
    print(f"Communities   : {communities_csv}")
    print(f"Memberships   : {memberships_csv}")
    print(f"Metadata      : {run_dir / 'run_metadata.json'}")

    if not summary_df.empty:
        cols_show = [
            c for c in [
                "status", "algo_key", "algo_name", "runtime_sec", "modularity_overlap",
                "modularity", "coverage_ratio", "overlap_node_ratio", "n_communities",
                "max_community_size", "error_type", "error_message",
            ] if c in summary_df.columns
        ]
        print("\nRingkasan hasil:")
        print(summary_df[cols_show].to_string(index=False))


if __name__ == "__main__":
    main()
