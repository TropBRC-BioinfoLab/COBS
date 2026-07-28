import math
import streamlit as st
import pandas as pd
import networkx as nx
import html
from core.data import (
    load_default_graph,
    build_graph_from_disease,
    build_graph_from_gene_list,
    build_graph_from_lungadenocarcinoma_file,
    build_graph_from_breastcancer_file,
    build_graph_from_breastcancerdisgenet_file,
)
from core.algos.registry import list_specs, get_spec
from core.metrics import postprocess_and_memberships, compute_evaluations
from core.ui_centrality import render_centrality, render_ccc_overlap
from core.ui_viz import render_visualization
from core.ui_communities import render_communities_paginated
from core.ui_enrichment import render_enrichment
from core.ui_hallmark import render_hallmark_jackknife
import json
from datetime import datetime
from zoneinfo import ZoneInfo

def _now_jakarta():
    return datetime.now(ZoneInfo("Asia/Jakarta"))

def _infer_random_seed(params: dict):
    # ambil seed kalau memang ada (tidak semua algoritma punya)
    for k in ("seed", "random_seed", "random_state", "rng_seed"):
        if k in params:
            return params[k]
    return None



def _infer_weight_attr(G: nx.Graph):
    # Deteksi atribut bobot edge yang umum dipakai
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


def _crisp_partition_from_comms(G: nx.Graph, comms):
    # Return list[set] jika crisp (tanpa overlap), tambah singleton utk node uncovered.
    # Jika overlap terdeteksi -> return None.
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


def _crisp_partition_via_primary(G: nx.Graph, comms):
    """
    Jika hasil komunitas overlap, kita proyeksikan menjadi partisi crisp:
    - Setiap node dipetakan ke 1 komunitas "utama" (pilih komunitas terbesar yang memuat node; tie -> id terkecil).
    - Node yang tidak ter-cover dijadikan singleton.
    Tujuan: bisa menghitung modularity crisp secara seragam untuk semua algoritma (meski output overlap).
    """
    # rapikan komunitas & buang node yang tidak ada di G
    comm_sets = []
    for c in (comms or []):
        try:
            cset = {n for n in c if n in G}
        except TypeError:
            continue
        if cset:
            comm_sets.append(cset)

    # fallback: kalau kosong, jadikan semua singleton
    if not comm_sets:
        return [{n} for n in G.nodes()]

    sizes = [len(c) for c in comm_sets]

    # membership map node -> list community_id (0-based)
    node_to_cids = {n: [] for n in G.nodes()}
    for cid, cset in enumerate(comm_sets):
        for n in cset:
            node_to_cids[n].append(cid)

    # assign ke komunitas utama
    assignment = {}
    for n, cids in node_to_cids.items():
        if not cids:
            assignment[n] = None
        else:
            # max by size; jika tie, pilih cid terkecil
            best = max(cids, key=lambda cid: (sizes[cid], -cid))
            assignment[n] = best

    # build partition sets
    cid_to_nodes = {}
    for n, cid in assignment.items():
        if cid is None:
            continue
        cid_to_nodes.setdefault(cid, set()).add(n)

    part = [cid_to_nodes[cid] for cid in sorted(cid_to_nodes.keys()) if cid_to_nodes[cid]]

    # uncovered nodes as singleton
    for n, cid in assignment.items():
        if cid is None:
            part.append({n})

    if not part:
        part = [{n} for n in G.nodes()]
    return part



def _compute_modularity_crisp(G: nx.Graph, comms, resolution: float = 1.0, allow_overlap_projection: bool = True):
    """
    Modularity Newman–Girvan (crisp) untuk graf (NetworkX).
    - Jika output komunitas sudah crisp (tanpa overlap): hitung modularity langsung.
    - Jika overlap dan allow_overlap_projection=True: proyeksikan terlebih dahulu menjadi partisi crisp via komunitas utama.
      (Ini bukan "modularity overlap", tapi berguna untuk tabel evaluasi yang seragam antar-algoritma.)
    """
    part = _crisp_partition_from_comms(G, comms)
    if part is None:
        if not allow_overlap_projection:
            return None
        part = _crisp_partition_via_primary(G, comms)

    H = G if not getattr(G, "is_directed", lambda: False)() else G.to_undirected()
    w = _infer_weight_attr(H)
    try:
        return nx.algorithms.community.quality.modularity(H, part, weight=w, resolution=resolution)
    except Exception:
        return None


def _reset_for_new_graph():
    # hapus cache centrality untuk dataset sebelumnya
    for k in list(st.session_state.keys()):
        if k.startswith("centrality_df::"):
            del st.session_state[k]
    # hapus hasil run sebelumnya
    st.session_state["result"] = None
    # reset pagination komunitas (jika ada)
    if "page" in st.session_state:
        st.session_state.page = 1

def _set_graph(G_new, label: str, meta=None):
    st.session_state["G"] = G_new
    st.session_state["GRAPH_LABEL"] = label
    st.session_state["GRAPH_META"] = meta or {}
    _reset_for_new_graph()


st.set_page_config(page_title="Overlapping Community Detection", layout="wide")
st.title("Overlapping Community Detection")


# =========================
# Load graph (session-based)
# =========================
if "GRAPH_LABEL" not in st.session_state:
    st.session_state["GRAPH_LABEL"] = "Zachary Karate (default)"
if "G" not in st.session_state:
    st.session_state["G"] = load_default_graph(st.session_state["GRAPH_LABEL"])
if "GRAPH_META" not in st.session_state:
    st.session_state["GRAPH_META"] = {}

GRAPH_LABEL_BASE = st.session_state["GRAPH_LABEL"]
G_FULL = st.session_state["G"]

# --- cek konektivitas (pakai view undirected) ---
if G_FULL.number_of_nodes() > 0:
    try:
        UG = G_FULL.to_undirected(as_view=True)
        _n_components = nx.number_connected_components(UG)
        _is_connected = (_n_components == 1)
        _gc_nodes = max(nx.connected_components(UG), key=len) if not _is_connected else None
        _gc_size = len(_gc_nodes) if _gc_nodes else G_FULL.number_of_nodes()
    except Exception:
        _n_components = None
        _is_connected = True
        _gc_nodes = None
        _gc_size = G_FULL.number_of_nodes()
else:
    _n_components = 0
    _is_connected = True
    _gc_nodes = None
    _gc_size = 0

# Scope graf diset dari sidebar (setelah tombol "Muat graf").
# Default: otomatis ON kalau graf tidak terhubung.
if "USE_GIANT_FOR_ALL" not in st.session_state:
    st.session_state["USE_GIANT_FOR_ALL"] = (not _is_connected)

use_giant_for_all = bool(st.session_state.get("USE_GIANT_FOR_ALL", (not _is_connected)))

if use_giant_for_all and (not _is_connected) and _gc_nodes:
    G = G_FULL.subgraph(_gc_nodes).copy()
    GRAPH_LABEL = f"{GRAPH_LABEL_BASE} [GC]"
else:
    G = G_FULL
    GRAPH_LABEL = f"{GRAPH_LABEL_BASE} [FULL]"

st.caption(f"Dataset aktif: {GRAPH_LABEL}")
st.write(f"Graph info: **|V|={G.number_of_nodes()}**, **|E|={G.number_of_edges()}**")

if use_giant_for_all and (not _is_connected) and _gc_nodes:
    st.info(
        f"Memakai giant component: {G.number_of_nodes()} node dari {G_FULL.number_of_nodes()} node "
        f"(komponen lain diabaikan untuk analisis & visualisasi)."
    )

with st.expander("Preview sumber data (opsional)", expanded=False):
    meta = st.session_state.get("GRAPH_META", {}) or {}
    if not meta:
        st.info("Belum ada metadata (dataset default / belum memuat dari API).")
    else:
        if "uniprot_records" in meta:
            st.write("UniProt (sample):")
            st.dataframe(pd.DataFrame(meta["uniprot_records"]).head(20), width='stretch')
        if "string_mapping" in meta:
            st.write("STRING mapping (sample):")
            st.dataframe(pd.DataFrame(meta["string_mapping"]).head(20), width='stretch')
        if "string_edges" in meta:
            st.write("STRING edges (sample):")
            st.dataframe(pd.DataFrame(meta["string_edges"]).head(20), width='stretch')


# =========================
# Centrality (Graph-level)
# =========================
render_centrality(G, GRAPH_LABEL)


# =========================
# Sidebar (auto dari spec)
# =========================

with st.sidebar:
    st.header("Kontrol")

    # -------------------------
    # Dataset / Graph Source
    # -------------------------
    st.subheader("Sumber graf")

    source = st.radio(
        "Pilih sumber dataset",
        [
            "Zachary Karate (default)",
            "Disease → UniProt → STRING",
            "Gene list → STRING",
            "Lung adenocarcinoma (lungadenocarcinoma.txt)",
            "Breast cancer (breastcancer.txt)",
            "Breast cancer DisGeNET (breastcancerdisgenet.txt)",
        ],
        index=0,
        label_visibility="collapsed",
    )

    # parameter umum untuk STRING
    st.caption("Parameter umum (STRING)")
    species = st.number_input("NCBI Taxon (species)", min_value=1, value=9606, step=1)
    required_score = st.slider("STRING required_score (0–999)", 0, 999, 700, 1)
    add_nodes = st.slider("Tambah neighbor (STRING add_nodes)", 0, 100, 0, 1)

    # sumber 2: disease → UniProt → STRING
    disease = None
    reviewed_only = True
    uniprot_limit = 200

    genes_text = ""
    lung_file_path = "lungadenocarcinoma.txt"
    breast_file_path = "breastcancer.txt"
    breast_disgenet_file_path = "breastcancerdisgenet.txt"

    if source == "Disease → UniProt → STRING":
        disease = st.text_input("Nama penyakit", value="lung cancer")
        reviewed_only = st.checkbox("UniProt reviewed saja", value=True)
        uniprot_limit = st.slider("Maks protein dari UniProt", 20, 500, 200, 10)

    # sumber 3: gene list → STRING
    if source == "Gene list → STRING":
        genes_text = st.text_area(
            "Daftar gen (pisahkan dengan koma / baris baru)",
            value="TP53, BRCA1, EGFR",
            height=120,
        )


    # sumber 4: lung adenocarcinoma file → STRING
    if source == "Lung adenocarcinoma (lungadenocarcinoma.txt)":
        lung_file_path = st.text_input(
            "Path file dataset lung adenocarcinoma",
            value=lung_file_path,
            help="File teks berisi UniProt accessions/identifiers; 1 per baris atau dipisah koma. Default: lungadenocarcinoma.txt (di folder proyek).",
        )

    # sumber 5: breast cancer file → STRING
    if source == "Breast cancer (breastcancer.txt)":
        breast_file_path = st.text_input(
            "Path file dataset breast cancer",
            value=breast_file_path,
            help="File teks berisi gene symbols/protein identifiers hasil filter Open Targets; 1 per baris atau dipisah koma. Default: breastcancer.txt (di folder proyek).",
        )

    # sumber 6: breast cancer DisGeNET file → STRING
    if source == "Breast cancer DisGeNET (breastcancerdisgenet.txt)":
        breast_disgenet_file_path = st.text_input(
            "Path file dataset breast cancer DisGeNET",
            value=breast_disgenet_file_path,
            help="File teks berisi gene symbols/protein identifiers hasil ekspor DisGeNET; 1 per baris atau dipisah koma. Default: breastcancerdisgenet.txt (di folder proyek).",
        )

    load_btn = st.button("Muat graf", type="secondary")

    if load_btn:
        try:
            with st.spinner("Mengambil data dari API dan membangun graf..."):
                if source == "Zachary Karate (default)":
                    G_new = load_default_graph("Zachary Karate (default)")
                    label = "Zachary Karate (default)"
                    meta = {"source": "zachary"}
                elif source == "Disease → UniProt → STRING":
                    G_new, meta = build_graph_from_disease(
                        disease=disease,
                        organism_id=int(species),
                        reviewed_only=bool(reviewed_only),
                        uniprot_size=int(uniprot_limit),
                        string_species=int(species),
                        required_score=int(required_score),
                        add_nodes=int(add_nodes),
                    )
                    label = f"Disease:{disease} | taxon:{species} | score≥{required_score} | add_nodes:{add_nodes}"
                elif source == "Gene list → STRING":
                    G_new, meta = build_graph_from_gene_list(
                        genes_text=genes_text,
                        string_species=int(species),
                        required_score=int(required_score),
                        add_nodes=int(add_nodes),
                    )
                    label = f"GeneList | taxon:{species} | score≥{required_score} | add_nodes:{add_nodes}"
                elif source == "Lung adenocarcinoma (lungadenocarcinoma.txt)":
                    G_new, meta = build_graph_from_lungadenocarcinoma_file(
                        file_path=lung_file_path,
                        string_species=int(species),
                        required_score=int(required_score),
                        add_nodes=int(add_nodes),
                    )
                    label = f"LungAdeno | file:{lung_file_path} | taxon:{species} | score≥{required_score} | add_nodes:{add_nodes}"
                elif source == "Breast cancer (breastcancer.txt)":
                    G_new, meta = build_graph_from_breastcancer_file(
                        file_path=breast_file_path,
                        string_species=int(species),
                        required_score=int(required_score),
                        add_nodes=int(add_nodes),
                    )
                    label = f"BreastCancer | file:{breast_file_path} | taxon:{species} | score≥{required_score} | add_nodes:{add_nodes}"
                elif source == "Breast cancer DisGeNET (breastcancerdisgenet.txt)":
                    G_new, meta = build_graph_from_breastcancerdisgenet_file(
                        file_path=breast_disgenet_file_path,
                        string_species=int(species),
                        required_score=int(required_score),
                        add_nodes=int(add_nodes),
                    )
                    label = f"BreastCancerDisGeNET | file:{breast_disgenet_file_path} | taxon:{species} | score≥{required_score} | add_nodes:{add_nodes}"
                else:
                    raise ValueError(f"Sumber dataset tidak dikenali: {source}")

            _set_graph(G_new, label, meta)
            st.success("Graf berhasil dimuat.")
            st.rerun()
        except Exception as e:
            st.error(f"Gagal memuat graf: {type(e).__name__}: {e}")

    
    # -------------------------
    # Scope graf (FULL vs GC)
    # -------------------------
    st.subheader("Scope graf")
    if _n_components is not None and _n_components > 1:
        st.caption(f"Graf terpecah: {_n_components} komponen | Giant component: {_gc_size} node")

    # init default (otomatis ON kalau graf tidak terhubung)
    if "USE_GIANT_FOR_ALL" not in st.session_state:
        st.session_state["USE_GIANT_FOR_ALL"] = (not _is_connected)
    if "_PREV_USE_GIANT_FOR_ALL" not in st.session_state:
        st.session_state["_PREV_USE_GIANT_FOR_ALL"] = st.session_state["USE_GIANT_FOR_ALL"]

    st.checkbox(
        "Gunakan giant component untuk analisis & visualisasi",
        key="USE_GIANT_FOR_ALL",
        help="Jika graf tidak terhubung, mode FULL akan menampilkan komponen terpisah. Mode GC hanya memakai komponen terbesar.",
    )

    # kalau scope berubah, reset hasil komputasi & rerun
    if st.session_state["_PREV_USE_GIANT_FOR_ALL"] != st.session_state["USE_GIANT_FOR_ALL"]:
        st.session_state["_PREV_USE_GIANT_FOR_ALL"] = st.session_state["USE_GIANT_FOR_ALL"]
        _reset_for_new_graph()
        st.rerun()

    st.divider()

    # -------------------------
    # Algorithm controls (existing)
    # -------------------------
    specs_all = list_specs()

    # Sembunyikan DGM (sementara) dari dropdown "Algoritma (overlap)"
    # Catatan: DGM_CRISP adalah versi crisp (views.py implementation) yang kita pakai untuk uji modularity.
    specs = []
    for s in specs_all:
        key = (getattr(s, "key", "") or "").lower()
        name = (getattr(s, "name", "") or "").lower()
        desc = (getattr(s, "description", "") or "").lower()

        # hide DGM crisp (key tetap DGM_CRISP)
        if key == "dgm_crisp":
            continue

        # hanya hide DGM yang overlap (jika suatu saat ada)
        is_dgm = ("dgm" in key) or ("dgm" in name)
        is_overlap = ("overlap" in key) or ("overlap" in name) or ("overlap" in desc) or ("overlapping" in desc)
        if is_dgm and is_overlap:
            continue

        specs.append(s)

    name_to_key = {s.name: s.key for s in specs}
    pick_name = st.selectbox("Algoritma (overlap)", list(name_to_key.keys()), index=0)
    algo_key = name_to_key[pick_name]
    spec = get_spec(algo_key)

    st.caption(spec.description)


with st.sidebar:
    st.subheader("Parameter")
    params = {}

    # COBS: tambah opsi similarity metric (selain Salton/cosine)
    _is_cobs = ((getattr(spec, "key", "") or "").lower() == "cobs")
    if _is_cobs:
        _metric_labels = {
            "cosine": "Cosine / Salton (Ochiai)",
            "jaccard": "Jaccard",
            "dice": "Sørensen–Dice",
            "overlap": "Overlap Coefficient (Szymkiewicz–Simpson)",
            "aa": "Adamic–Adar (weighted cosine)",
            "ra": "Resource Allocation (weighted cosine)",
        }
        _metric_keys = list(_metric_labels.keys())
        params["sim_metric"] = st.selectbox(
            "Similarity metric",
            _metric_keys,
            index=0,
            format_func=lambda k: _metric_labels.get(k, str(k)),
            key="cobs_sim_metric",
            help=(
                "Dipakai untuk: (i) similarity N(u) vs N(Cj) pada Step4, "
                "dan (ii) similarity Ci vs Cj pada Step5. "
                "AA/RA diimplementasikan sebagai *weighted cosine* berbasis derajat node."
            ),
        )

    for p in spec.params:
        # jika SPEC sudah punya sim_metric, jangan tampilkan dua kali (karena sudah pakai selectbox di atas)
        if _is_cobs and getattr(p, "key", None) == "sim_metric":
            continue
        params[p.key] = p.render(st)

    st.subheader("Postprocess")
    min_size = st.slider("Min ukuran komunitas (untuk tampilan)", 1, 20, 1, 1, key="post_min_comm_size")
    show_missing = st.checkbox("Tampilkan node yang tidak ter-cover", value=True, key="post_show_missing")

    _run_key = f"run_{(getattr(spec, 'key', '') or 'algo')}"
    run_btn = st.button(f"Jalankan {spec.name}", type="primary", key=_run_key)


# =========================
# Run & persist result
# =========================
if "result" not in st.session_state:
    st.session_state["result"] = None

if run_btn:
    ts = _now_jakarta()
    run_id = ts.strftime("%Y%m%d_%H%M%S")                 # aman untuk nama file
    timestamp = ts.isoformat(timespec="seconds")         # timestamp lengkap + offset timezone
    random_seed = _infer_random_seed(params)
    try:
        with st.spinner("Menjalankan algoritma..."):
            clustering = spec.run(G, params)
    except Exception as e:
        st.error(f"Gagal menjalankan {spec.name}: {type(e).__name__}: {e}")
        st.stop()

    raw_comms = clustering.communities or []
    covered = set()
    for c in raw_comms:
        for n in c:
            if n in G:   # guard kalau ada node aneh
                covered.add(n)

    uncovered_nodes = [n for n in G.nodes() if n not in covered]
    n_uncovered = len(uncovered_nodes)


    comms, memberships, uncovered_after_min_size = postprocess_and_memberships(
        G, raw_comms, min_size=min_size
    )

    clustering.communities = comms
    evals = compute_evaluations(clustering)

    # modularity crisp (Newman–Girvan): selalu ditampilkan.
    # Jika output overlap, modularity dihitung dari proyeksi crisp (primary community per node).
    gamma = float(params.get("resolution", 1.0))
    mod_crisp = _compute_modularity_crisp(G, raw_comms, resolution=gamma, allow_overlap_projection=True)
    evals["modularity"] = float(mod_crisp) if mod_crisp is not None else float("nan")

    # metadata eksperimen (reproducibility)
    meta = {
        "run_id": run_id,
        "timestamp": timestamp,
        "dataset": GRAPH_LABEL,
        "algo_key": spec.key,
        "algo_name": spec.name,
        "params": params,
        "random_seed": random_seed,
        "n_nodes": G.number_of_nodes(),
        "n_edges": G.number_of_edges(),
    }

    st.session_state["result"] = {
        # meta utama
        "meta": meta,

        "algo_name": spec.name,
        "algo_key": spec.key,
        "params": params,
        "min_size": min_size,
        #"ensure_cover": ensure_cover,
        "clustering": clustering,
        "communities": comms,
        "memberships": memberships,
        "evals": evals,

        "show_missing": show_missing,
        # raw vs display
        "communities_raw": raw_comms,
        "uncovered_nodes": uncovered_nodes,
        "n_uncovered": n_uncovered,

        "uncovered_nodes_after_min_size": uncovered_after_min_size,
        "n_uncovered_after_min_size": len(uncovered_after_min_size)
    }

res = st.session_state["result"]
if res is None:
    st.info("Klik tombol Jalankan untuk mulai.")
    st.stop()


selected_algo = spec.name
last_algo = res["algo_name"] if res else None

if res is None:
    st.markdown(f"**Algoritma dipilih:** {selected_algo}")
else:
    if selected_algo != last_algo:
        st.markdown(f"**Algoritma dipilih:** {selected_algo}")
        st.warning(f"Hasil di kanan masih dari **{last_algo}**. Klik **Jalankan** untuk menjalankan {selected_algo}.")
    else:
        st.markdown(f"**Algoritma terakhir dijalankan:** {last_algo}")


if res is not None and selected_algo == res["algo_name"]:
    p = res["params"]
    st.caption("Parameter: " + ", ".join([f"{k}={v}" for k, v in p.items()]))
else:
    # parameter yang sedang dipilih (belum dijalankan)
    st.caption("Parameter dipilih: " + ", ".join([f"{k}={v}" for k, v in params.items()]))

meta = res.get("meta", {})
if meta:
    # ringkasan kecil (tanpa json panjang)
    seed_txt = f", seed={meta['random_seed']}" if meta.get("random_seed") is not None else ""
    st.caption(
        f"Run: {meta.get('run_id')} | {meta.get('timestamp')} | "
        f"algo_key={meta.get('algo_key')} | |V|={meta.get('n_nodes')} | |E|={meta.get('n_edges')}{seed_txt}"
    )

    with st.expander("Metadata eksperimen", expanded=False):
        meta_json = json.dumps(meta, ensure_ascii=False, indent=2)
        st.download_button(
            "Download metadata (JSON)",
            data=meta_json,
            file_name=f"metadata_{meta.get('algo_key','ALGO')}_{meta.get('run_id','RUN')}.json",
            mime="application/json",
        )


if res.get("show_missing", True):
    uncovered = res.get("uncovered_nodes", [])
    n_uncovered = res.get("n_uncovered", len(uncovered))

    if n_uncovered == 0:
        st.success("Semua node ter-cover oleh komunitas (raw output algoritma).")
    else:
        st.warning(f"Ada **{n_uncovered}** node yang tidak ter-cover oleh komunitas (raw output algoritma).")

        with st.expander(f"Lihat daftar node tidak ter-cover (raw) ({n_uncovered})", expanded=False):
            # pagination khusus uncovered nodes (tidak bentrok dengan pagination komunitas)
            if "uc_page_size" not in st.session_state:
                st.session_state.uc_page_size = 50
            if "uc_page" not in st.session_state:
                st.session_state.uc_page = 1

            page_size = st.session_state.uc_page_size
            pages = max(1, math.ceil(n_uncovered / page_size))
            st.session_state.uc_page = max(1, min(st.session_state.uc_page, pages))
            page = st.session_state.uc_page

            start = (page - 1) * page_size
            end = min(start + page_size, n_uncovered)

            st.dataframe(
                pd.DataFrame({"node": uncovered[start:end]}),
                hide_index=True,
                width='stretch',
            )

            left, right = st.columns([3, 2], vertical_alignment="center")
            with left:
                st.caption(f"{start+1} to {end} of {n_uncovered}")
            with right:
                c1, c2, c3 = st.columns([1, 1, 2], vertical_alignment="center")
                if c1.button("◀", help="Prev", width='stretch', key="uc_prev"):
                    st.session_state.uc_page = max(1, st.session_state.uc_page - 1)
                    st.rerun()
                if c2.button("▶", help="Next", width='stretch', key="uc_next"):
                    st.session_state.uc_page = min(pages, st.session_state.uc_page + 1)
                    st.rerun()

                new_page = c3.number_input(
                    f"Page of {pages}",
                    min_value=1,
                    max_value=pages,
                    value=st.session_state.uc_page,
                    step=1,
                    label_visibility="collapsed",
                    key="uc_jump",
                )
                if new_page != st.session_state.uc_page:
                    st.session_state.uc_page = int(new_page)
                    st.rerun()

    # Transparansi dampak min_size untuk tampilan
    if res.get("min_size", 1) > 1:
        post = set(res.get("uncovered_nodes_after_min_size", []))
        raw = set(uncovered)
        extra_due_to_filter = sorted(list(post - raw))
        if len(extra_due_to_filter) > 0:
            st.caption(
                f"Catatan: karena **min_size={res['min_size']}** hanya untuk tampilan, "
                f"ada **{len(extra_due_to_filter)}** node yang terlihat 'tidak ter-cover' setelah filter (bukan dari raw)."
            )




# =========================
# Community-Consideration Centrality (Overlap-aware)
# =========================
render_ccc_overlap(G, res, GRAPH_LABEL)

# =========================
# Hallmark-of-Cancer (Jackknife curve)
# =========================
render_hallmark_jackknife(G, res, GRAPH_LABEL)

# =========================
# Visualisasi (FULL WIDTH)
# =========================
render_visualization(G, res)


# =========================
# Evaluasi 
# =========================
METRIC_ORDER = [
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
]

METRIC_INFO = {
    "modularity_overlap": ("↑", "Modularity versi overlap; makin besar umumnya struktur komunitas makin kuat."),
    "modularity": ("↑", "Modularity Newman–Girvan (crisp). Jika hasil overlap, dihitung dari proyeksi crisp (tiap node dipilih 1 komunitas utama; node uncovered jadi singleton) agar tabel evaluasi seragam."),
    "internal_edge_density": ("↑", "Kepadatan edge internal; makin besar komunitas makin kompak."),
    "conductance": ("↓", "Kebocoran edge keluar komunitas; makin kecil boundary makin rapat."),
    "normalized_cut": ("↓", "Cut relatif terhadap volume; makin kecil pemisahan komunitas makin baik."),
    "avg_transitivity": ("↑", "Kohesi triadik (segitiga); makin besar makin kohesif (jika graf kaya segitiga)."),
    "n_communities": ("", "Jumlah komunitas hasil (umumnya setelah filter min_size)."),
    "max_community_size": ("", "Ukuran komunitas terbesar; berguna untuk cek ‘komunitas raksasa’."),
    "coverage_ratio": ("↑", "Proporsi node yang masuk minimal 1 komunitas."),
    "uncovered_node_ratio": ("↓", "Proporsi node yang tidak masuk komunitas mana pun (membership=0)."),
    "overlap_node_ratio": ("", "Proporsi node yang masuk >1 komunitas (node overlap)."),
    "avg_memberships_per_covered_node": ("↑", "Rata-rata membership hanya untuk node yang ter-cover (>=1)."),
    "overlap_intensity": ("↑", "Rata-rata membership khusus node overlap saja (>1)."),
}


with st.expander("Evaluasi", expanded=True):
    evals = res.get("evals", {}) or {}

    # susun urutan: METRIC_ORDER dulu, sisanya di bawah
    ordered = [m for m in METRIC_ORDER if m in evals]
    rest = [m for m in evals.keys() if m not in ordered]
    metrics = ordered + rest

    rows = []
    for m in metrics:
        arrow, tip = METRIC_INFO.get(m, ("", ""))
        m_html = f'<span title="{html.escape(tip)}">{html.escape(m)}</span>'
        val = evals.get(m, None)

        # format angka
        if isinstance(val, (int, float)) and val == val:  # val==val untuk cek NaN
            score_str = f"{val:.4f}"
        else:
            score_str = "" if val is None else str(val)

        rows.append({"metric": m_html, "dir": arrow, "score": score_str})

    df = pd.DataFrame(rows)

    # render HTML table agar tooltip aktif
    table_html = df.to_html(index=False, escape=False)
    st.markdown(
        f'<div style="overflow-x:auto">{table_html}</div>',
        unsafe_allow_html=True
    )

    st.caption("Tip: arah ↑/↓ adalah kecenderungan umum untuk komunitas yang kompak; konteks graf/tujuan analisis tetap menentukan.")

    # ===== Export / Copy evaluasi (untuk Excel) =====
    try:
        export_rows = []
        for m in metrics:
            arrow, _tip = METRIC_INFO.get(m, ("", ""))
            export_rows.append({"metric": m, "dir": arrow, "score": evals.get(m, None)})
        df_export = pd.DataFrame(export_rows)
        # CSV download
        csv_bytes = df_export.to_csv(index=False).encode("utf-8")
        st.download_button(
            "Download evaluasi (CSV)",
            data=csv_bytes,
            file_name="evaluasi_metrics.csv",
            mime="text/csv",
        )
    except Exception:
        pass


    # ---- DGM: riwayat modularity (hanya jika tersedia) ----
    try:
        clustering_obj = res.get("clustering", None)
        mp = getattr(clustering_obj, "method_parameters", {}) or {}
    except Exception:
        mp = {}

    records = mp.get("records") if isinstance(mp, dict) else None
    q_hist = mp.get("q_history") if isinstance(mp, dict) else None
    best_hist = mp.get("best_history") if isinstance(mp, dict) else None

    has_records = isinstance(records, (list, tuple)) and len(records) > 0
    has_simple = (
        isinstance(q_hist, (list, tuple)) and isinstance(best_hist, (list, tuple))
        and len(q_hist) > 0 and len(best_hist) > 0
    )

    if has_records or has_simple:
        st.markdown("#### DGM: riwayat modularity")

        if has_records:
            hist_df = pd.DataFrame(records)
            # Basic stats (jika tersedia)
            stats = mp.get("stats", {}) if isinstance(mp, dict) else {}
            n_explore = int(stats.get("n_explore", int(hist_df["did_explore"].sum() if "did_explore" in hist_df else 0)))
            n_node = int(stats.get("n_random_node", int((hist_df["explore_type"]=="random_node").sum() if "explore_type" in hist_df else 0)))
            n_cond = int(stats.get("n_high_conductance", int((hist_df["explore_type"]=="high_conductance").sum() if "explore_type" in hist_df else 0)))
            n_improve = int(stats.get("n_improve", int(
                hist_df.get("accepted_exploit", pd.Series(dtype=int)).sum()
                + hist_df.get("accepted_explore", pd.Series(dtype=int)).sum()
            )))

            st.caption(
                f"Iterasi: {len(hist_df)} | exploration terjadi: {n_explore} kali "
                f"(random_node={n_node}, high_conductance={n_cond}) | improve accepted: {n_improve} kali"
            )

            # Build plot dataframe (kolom fleksibel)
            plot_cols = {}
            if "Q_used" in hist_df:
                plot_cols["modularity"] = hist_df["Q_used"]
            elif "modularity" in hist_df:
                plot_cols["modularity"] = hist_df["modularity"]

            if "Q_best" in hist_df:
                plot_cols["maximum_modularity"] = hist_df["Q_best"]
            elif "maximum_modularity" in hist_df:
                plot_cols["maximum_modularity"] = hist_df["maximum_modularity"]

            # Optional: modularity setelah disassembly (sebelum greedy)
            if "Q_trial_raw" in hist_df:
                plot_cols["after_disassembly_raw"] = hist_df["Q_trial_raw"]

            plot_df = pd.DataFrame(plot_cols)

            # downsample if too long
            L = len(plot_df)
            max_points = 5000
            if L > max_points:
                step = max(1, L // max_points)
                plot_df = plot_df.iloc[::step].reset_index(drop=True)
                L = len(plot_df)

            st.caption(f"Menampilkan {L} titik (1 titik per iterasi; jika terlalu panjang otomatis di-sampling).")
            st.line_chart(plot_df)

            with st.expander("Detail iterasi DGM (preview 50 baris)"):
                st.dataframe(hist_df.head(50), width='stretch')

            csv_bytes = hist_df.to_csv(index=False).encode("utf-8")
            st.download_button(
                "Download modularity history (CSV)",
                data=csv_bytes,
                file_name=f"dgm_history_{res.get('meta',{}).get('run_id','run')}.csv",
                mime="text/csv",
            )

        else:
            # fallback to simple histories (q_history & best_history)
            L = min(len(q_hist), len(best_hist))
            qh = list(q_hist)[:L]
            bh = list(best_hist)[:L]

            max_points = 5000
            if L > max_points:
                step = max(1, L // max_points)
                qh = qh[::step]
                bh = bh[::step]
                L = len(qh)

            st.caption(f"Menampilkan {L} titik (1 titik per iterasi; jika terlalu panjang otomatis di-sampling).")
            hist_df = pd.DataFrame({"modularity": qh, "maximum_modularity": bh})
            st.line_chart(hist_df)

            csv_bytes = hist_df.to_csv(index_label="iter").encode("utf-8")
            st.download_button(
                "Download modularity history (CSV)",
                data=csv_bytes,
                file_name=f"dgm_modularity_history_{res.get('meta',{}).get('run_id','run')}.csv",
                mime="text/csv",
            )

    # Jika tidak ada history, jangan apa-apa (hindari error len(None))
# =========================
# Communities with pagination
# =========================
render_communities_paginated(res)


# =========================
# Enrichment Analysis
# =========================
render_enrichment(res, GRAPH_LABEL)