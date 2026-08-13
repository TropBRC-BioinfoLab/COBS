import math
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

import networkx as nx
import pandas as pd
import streamlit as st


# =========================
# Hallmark parser
# =========================

def _canon_gene(x: str) -> str:
    """Normalisasi sederhana agar match ke simbol gen."""
    if x is None:
        return ""
    return str(x).strip().upper()


def load_hallmark_from_text(text: str) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Set[str]]:
    """Parse file hallmark: 1 baris = 1 hallmark, dipisah TAB.

    Format yang didukung:
      HALLMARK_NAME<TAB>GENE1<TAB>GENE2<...>

    Return:
      - hallmark_to_genes: {hallmark: set(gene)}
      - gene_to_hallmarks: {gene: set(hallmark)}
      - union_genes: set(gene) gabungan semua hallmark
    """
    hallmark_to_genes: Dict[str, Set[str]] = {}

    for raw in (text or "").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue

        parts = [p.strip() for p in line.split("\t") if p.strip()]
        if len(parts) < 2:
            parts = [p.strip() for p in line.split() if p.strip()]
        if len(parts) < 2:
            continue

        hallmark = parts[0]
        genes = {_canon_gene(g) for g in parts[1:] if _canon_gene(g)}
        if genes:
            hallmark_to_genes.setdefault(hallmark, set()).update(genes)

    gene_to_hallmarks: Dict[str, Set[str]] = {}
    for hm, genes in hallmark_to_genes.items():
        for g in genes:
            gene_to_hallmarks.setdefault(g, set()).add(hm)

    union = set().union(*hallmark_to_genes.values()) if hallmark_to_genes else set()
    return hallmark_to_genes, gene_to_hallmarks, union


def load_hallmark(path: str) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Set[str]]:
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        return load_hallmark_from_text(f.read())


# =========================
# Jackknife (cumulative hit) curve
# =========================

def jackknife_curve_genes(
    ranked_genes: Sequence[str],
    positive_set: Set[str],
    k_max: int,
) -> List[int]:
    """Kurva cumulative hits: y(k) = jumlah GEN hallmark yang ter-capture hingga top-k."""
    hits: List[int] = []
    c = 0
    for g in ranked_genes[: max(0, int(k_max))]:
        if _canon_gene(g) in positive_set:
            c += 1
        hits.append(c)
    return hits


def curve_auc(hits: Sequence[int]) -> float:
    """AUC trapezoid untuk kurva hit vs k (k=1..K)."""
    if not hits:
        return 0.0
    area = 0.0
    prev_y = 0.0
    for _, y in enumerate(hits, start=1):
        area += (prev_y + float(y)) / 2.0
        prev_y = float(y)
    return area


# =========================
# Centrality helpers
# =========================

def _infer_weight_attr(G: nx.Graph) -> Optional[str]:
    try:
        for _, _, d in G.edges(data=True):
            if isinstance(d, dict):
                if "weight" in d:
                    return "weight"
                if "score" in d:
                    return "score"
    except Exception:
        pass
    return None


def _safe_undirected(G: nx.Graph) -> nx.Graph:
    try:
        return G if not G.is_directed() else G.to_undirected()
    except Exception:
        return G


def _rank_from_scores(scores: Dict[str, float], descending: bool = True) -> List[str]:
    items = list(scores.items())
    items.sort(
        key=lambda kv: (float("-inf") if kv[1] is None else float(kv[1])),
        reverse=descending,
    )
    return [str(k) for k, _ in items]


@dataclass
class BuiltinCentralityConfig:
    degree: bool = True
    closeness: bool = True
    betweenness: bool = True
    eigenvector: bool = True
    katz: bool = True
    harmonic: bool = True
    load: bool = True
    information: bool = False
    communicability_betweenness: bool = False
    subgraph: bool = False


def compute_builtin_rankings(
    G: nx.Graph,
    cfg: BuiltinCentralityConfig,
    weight_attr: Optional[str] = None,
) -> Dict[str, List[str]]:
    """Return {method_name: ranked_genes}."""
    UG = _safe_undirected(G)
    w = weight_attr

    rankings: Dict[str, List[str]] = {}

    if cfg.degree:
        try:
            scores = nx.degree_centrality(UG)
        except Exception:
            scores = {n: float(d) for n, d in UG.degree()}
        rankings["degree"] = _rank_from_scores(scores)

    if cfg.closeness:
        try:
            scores = nx.closeness_centrality(UG)
            rankings["closeness"] = _rank_from_scores(scores)
        except Exception:
            pass

    if cfg.betweenness:
        try:
            scores = nx.betweenness_centrality(UG, normalized=True, weight=None)
            rankings["betweenness"] = _rank_from_scores(scores)
        except Exception:
            pass

    if cfg.harmonic:
        try:
            scores = nx.harmonic_centrality(UG)
            rankings["harmonic"] = _rank_from_scores(scores)
        except Exception:
            pass

    if cfg.load:
        try:
            scores = nx.load_centrality(UG, normalized=True, weight=None)
            rankings["load"] = _rank_from_scores(scores)
        except Exception:
            pass

    if cfg.eigenvector:
        try:
            scores = nx.eigenvector_centrality(UG, max_iter=2000, tol=1e-6, weight=w)
            rankings["eigenvector"] = _rank_from_scores(scores)
        except Exception:
            try:
                scores = nx.eigenvector_centrality_numpy(UG, weight=w)
                rankings["eigenvector"] = _rank_from_scores(scores)
            except Exception:
                pass

    if cfg.katz:
        try:
            alpha = 0.005
            try:
                ev = max(abs(x) for x in nx.adjacency_spectrum(UG))
                if ev and float(ev) > 0:
                    alpha = min(alpha, 0.9 / float(ev))
            except Exception:
                pass
            scores = nx.katz_centrality(UG, alpha=alpha, beta=1.0, max_iter=5000, tol=1e-8, weight=w)
            rankings["katz"] = _rank_from_scores(scores)
        except Exception:
            pass

    if cfg.information:
        try:
            scores = nx.information_centrality(UG, weight=None)
            rankings["information"] = _rank_from_scores(scores)
        except Exception:
            pass

    if cfg.communicability_betweenness:
        try:
            scores = nx.communicability_betweenness_centrality(UG)
            rankings["communicability_betweenness"] = _rank_from_scores(scores)
        except Exception:
            pass

    if cfg.subgraph:
        try:
            scores = nx.subgraph_centrality(UG)
            rankings["subgraph"] = _rank_from_scores(scores)
        except Exception:
            pass

    return rankings


# =========================
# CCC overlap (recompute from communities)
# =========================

def _memberships_from_communities(G: nx.Graph, communities: Sequence[Iterable[str]]) -> Dict[str, List[int]]:
    mem: Dict[str, List[int]] = {str(n): [] for n in G.nodes()}
    for idx, c in enumerate(communities or []):
        try:
            nodes = [str(n) for n in c if n in G]
        except TypeError:
            continue
        for n in nodes:
            mem.setdefault(n, []).append(idx)
    return mem


def _sum_shortest_paths_lengths(H: nx.Graph) -> Dict[str, float]:
    sums: Dict[str, float] = {str(n): 0.0 for n in H.nodes()}
    for i in H.nodes():
        lengths = nx.single_source_shortest_path_length(H, i)
        s = 0
        for _, d in lengths.items():
            s += d
        sums[str(i)] = float(s)
    return sums


def compute_ccc_overlap_rankings(
    G: nx.Graph,
    communities: Sequence[Iterable[str]],
    alphas: Sequence[float] = (0.1, 0.2),
    agg: str = "mean",  # mean|max
    compute_degree: bool = True,
    compute_closeness: bool = True,
    compute_betweenness: bool = True,
) -> Dict[str, List[str]]:
    """Hitung ulang CCC overlap untuk beberapa alpha dan return ranking.

    Catatan: betweenness selalu **exact** (tidak ada opsi approx) sesuai permintaan.
    """
    UG = _safe_undirected(G)
    nG = UG.number_of_nodes()
    if nG == 0:
        return {}

    _ = _memberships_from_communities(UG, communities)  # tetap dihitung (berguna untuk debugging/cek)

    dG = {str(n): float(UG.degree(n)) for n in UG.nodes()}
    sumDG = _sum_shortest_paths_lengths(UG) if compute_closeness else {}

    if compute_betweenness:
        try:
            bcG = nx.betweenness_centrality(UG, normalized=False)
            bcG = {str(n): float(v) for n, v in bcG.items()}
        except Exception:
            bcG = {str(n): 0.0 for n in UG.nodes()}
    else:
        bcG = {}

    PG = (nG - 1) * (nG - 2) / 2.0 if nG >= 3 else 0.0

    metrics: List[str] = []
    if compute_degree:
        metrics.append("ccdc")
    if compute_closeness:
        metrics.append("cccc")
    if compute_betweenness:
        metrics.append("ccbc")

    scores: Dict[str, Dict[float, Dict[str, List[float]]]] = {
        m: {float(a): {str(n): [] for n in UG.nodes()} for a in alphas} for m in metrics
    }

    for c in communities or []:
        try:
            c_nodes = [n for n in c if n in UG]
        except TypeError:
            continue
        if not c_nodes:
            continue

        C = UG.subgraph(c_nodes).copy()
        nC = C.number_of_nodes()
        if nC == 0:
            continue

        dC = {str(n): float(C.degree(n)) for n in C.nodes()} if compute_degree else {}
        sumDC = _sum_shortest_paths_lengths(C) if compute_closeness else {}

        if compute_betweenness:
            try:
                bcC = nx.betweenness_centrality(C, normalized=False)
                bcC = {str(n): float(v) for n, v in bcC.items()}
            except Exception:
                bcC = {str(n): 0.0 for n in C.nodes()}
        else:
            bcC = {}

        PC = (nC - 1) * (nC - 2) / 2.0 if nC >= 3 else 0.0

        for n in C.nodes():
            sn = str(n)
            for a in alphas:
                a = float(a)

                denom_deg = (nC + a * (nG - nC) - 1.0)
                denom_deg = denom_deg if denom_deg != 0 else float("nan")

                denom_cc = (sumDC.get(sn, 0.0) + a * (sumDG.get(sn, 0.0) - sumDC.get(sn, 0.0)))
                denom_cc = denom_cc if denom_cc != 0 else float("nan")

                denom_bc = (PC + a * (PG - PC))
                denom_bc = denom_bc if denom_bc != 0 else float("nan")

                if compute_degree:
                    val = (dC.get(sn, 0.0) + a * (dG.get(sn, 0.0) - dC.get(sn, 0.0))) / denom_deg
                    if math.isfinite(val):
                        scores["ccdc"][a][sn].append(float(val))

                if compute_closeness:
                    val = (nC + a * (nG - nC) - 1.0) / denom_cc
                    if math.isfinite(val):
                        scores["cccc"][a][sn].append(float(val))

                if compute_betweenness:
                    val = (bcC.get(sn, 0.0) + a * (bcG.get(sn, 0.0) - bcC.get(sn, 0.0))) / denom_bc
                    if math.isfinite(val):
                        scores["ccbc"][a][sn].append(float(val))

    def _agg(vals: List[float]) -> float:
        if not vals:
            return 0.0
        if agg == "max":
            return float(max(vals))
        return float(sum(vals) / len(vals))

    rankings: Dict[str, List[str]] = {}
    for metric in metrics:
        for a in alphas:
            a = float(a)
            sc = {n: _agg(scores[metric][a].get(str(n), [])) for n in UG.nodes()}
            rankings[f"{metric}_alpha={a:g}"] = _rank_from_scores(sc)

    return rankings


# =========================
# UI
# =========================

def render_hallmark_jackknife(
    G: nx.Graph,
    res: Optional[dict],
    graph_label: str,
    hallmark_path: str = "Hallmark of Cancer.txt",
):
    """Render evaluasi Hallmark-of-Cancer berbasis Jackknife curve.

    Perubahan sesuai permintaan:
    - File hallmark dianggap permanen (tanpa input path / upload).
    - Kurva selalu menghitung *Gen hallmark* (tanpa opsi kategori).
    - Betweenness CCC selalu exact (tanpa approx k).
    - Tidak ada pembanding dari CSV.
    """

    with st.expander("Hallmark-of-Cancer evaluation (Jackknife curve)", expanded=True):
        st.caption(
            "Membandingkan urutan protein (ranking centrality/CCC) dengan daftar Hallmark-of-Cancer. "
            "Kurva menunjukkan berapa banyak **GEN hallmark** yang ter-capture saat kita mengambil top-k node."
        )

        # --- Load hallmark (fixed) ---
        try:
            hallmark_to_genes, gene_to_hallmarks, union = load_hallmark(hallmark_path)
        except Exception as e:
            st.error(
                f"Gagal membaca hallmark file permanen '{hallmark_path}': {type(e).__name__}: {e}"
            )
            st.info("Pastikan file tersebut ada di folder proyek dan namanya persis sama.")
            return

        if not hallmark_to_genes:
            st.warning("Hallmark file kosong / format tidak terbaca.")
            return

        st.write(
            f"Hallmark file: **{hallmark_path}** | jumlah hallmark: **{len(hallmark_to_genes)}** | gen unik: **{len(union)}**"
        )

        # target set
        hallmarks = sorted(hallmark_to_genes.keys())
        target_mode = st.radio(
            "Target evaluasi",
            ["Union (semua hallmark)", "Satu hallmark tertentu"],
            horizontal=True,
            key=f"hm_target_mode::{graph_label}",
        )

        if target_mode == "Satu hallmark tertentu":
            pick_hm = st.selectbox(
                "Pilih hallmark",
                hallmarks,
                index=0,
                key=f"hm_pick::{graph_label}",
            )
            positive = set(hallmark_to_genes.get(pick_hm, set()))
            target_name = pick_hm
        else:
            positive = set(union)
            target_name = "UNION"

        # --- Methods ---
        st.subheader("Metode yang dibandingkan")
        w = _infer_weight_attr(_safe_undirected(G))
        st.caption(f"Edge weight terdeteksi: **{w or 'None'}**")

        left, right = st.columns([1, 1])

        with left:
            st.markdown("**Centrality (NetworkX)**")
            cfg = BuiltinCentralityConfig(
                degree=st.checkbox("degree", value=True, key=f"hm_m_deg::{graph_label}"),
                closeness=st.checkbox("closeness", value=True, key=f"hm_m_clo::{graph_label}"),
                betweenness=st.checkbox("betweenness", value=True, key=f"hm_m_bet::{graph_label}"),
                eigenvector=st.checkbox("eigenvector", value=True, key=f"hm_m_eig::{graph_label}"),
                katz=st.checkbox("katz", value=True, key=f"hm_m_katz::{graph_label}"),
                harmonic=st.checkbox("harmonic", value=False, key=f"hm_m_harm::{graph_label}"),
                load=st.checkbox("load", value=False, key=f"hm_m_load::{graph_label}"),
                information=st.checkbox("information", value=False, key=f"hm_m_info::{graph_label}"),
                communicability_betweenness=st.checkbox(
                    "communicability betweenness (berat)", value=False, key=f"hm_m_commb::{graph_label}"
                ),
                subgraph=st.checkbox("subgraph (berat)", value=False, key=f"hm_m_sub::{graph_label}"),
            )

        with right:
            st.markdown("**CCC overlap (dari komunitas)**")
            can_ccc = bool(res and res.get("communities"))
            if not can_ccc:
                st.info("Jalankan deteksi komunitas dulu agar CCC overlap tersedia.")

            use_ccc = st.checkbox(
                "Sertakan CCC overlap",
                value=can_ccc,
                disabled=not can_ccc,
                key=f"hm_use_ccc::{graph_label}",
            )

            ccc_deg = st.checkbox("CCDC", value=True, disabled=not use_ccc, key=f"hm_ccc_deg::{graph_label}")
            ccc_clo = st.checkbox("CCCC", value=True, disabled=not use_ccc, key=f"hm_ccc_clo::{graph_label}")
            ccc_bet = st.checkbox("CCBC", value=True, disabled=not use_ccc, key=f"hm_ccc_bet::{graph_label}")

            alphas = st.multiselect(
                "Alpha (multi)",
                options=[0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0],
                default=[0.1, 0.2],
                disabled=not use_ccc,
                key=f"hm_ccc_alpha::{graph_label}",
                help="Di paper CCC, alpha=0.1 dan 0.2 sering optimal.",
            )

            agg = st.selectbox(
                "Agregasi overlap",
                ["mean", "max"],
                index=0,
                disabled=not use_ccc,
                key=f"hm_ccc_agg::{graph_label}",
                help="mean = rata-rata skor antar membership; max = ambil yang terbesar.",
            )

        # --- Run evaluation ---
        st.subheader("Jalankan evaluasi")
        k_max = st.slider("Top-k maksimum", 5, 200, 50, 5, key=f"hm_kmax::{graph_label}")

        run = st.button("Hitung jackknife curve", type="primary", key=f"hm_run::{graph_label}")
        if not run:
            return

        with st.spinner("Menghitung ranking & jackknife curve..."):
            rankings: Dict[str, List[str]] = {}
            rankings.update(compute_builtin_rankings(G, cfg=cfg, weight_attr=w))

            if use_ccc and res and res.get("communities") and alphas:
                try:
                    ccc_rank = compute_ccc_overlap_rankings(
                        G,
                        communities=res.get("communities") or [],
                        alphas=alphas,
                        agg=str(agg),
                        compute_degree=bool(ccc_deg),
                        compute_closeness=bool(ccc_clo),
                        compute_betweenness=bool(ccc_bet),
                    )
                    rankings.update({f"CCC::{k}": v for k, v in ccc_rank.items()})
                except Exception as e:
                    st.warning(f"CCC overlap gagal dihitung: {type(e).__name__}: {e}")

            if not rankings:
                st.warning("Tidak ada metode yang dipilih / berhasil dihitung.")
                return

            curve_rows = []
            summary_rows = []
            for method, ranked in rankings.items():
                hits = jackknife_curve_genes(ranked, positive_set=positive, k_max=int(k_max))
                auc = curve_auc(hits)
                hit_k = hits[-1] if hits else 0
                hit_10 = hits[min(9, len(hits) - 1)] if hits else 0
                hit_20 = hits[min(19, len(hits) - 1)] if hits else 0

                summary_rows.append(
                    {
                        "method": method,
                        "target": target_name,
                        "k_max": int(k_max),
                        "hits@10": int(hit_10),
                        "hits@20": int(hit_20),
                        f"hits@{int(k_max)}": int(hit_k),
                        "auc": float(auc),
                    }
                )

                for k, h in enumerate(hits, start=1):
                    curve_rows.append({"k": k, "method": method, "hits": int(h)})

            df_curve = pd.DataFrame(curve_rows)
            df_sum = pd.DataFrame(summary_rows).sort_values(["hits@10", "auc"], ascending=[False, False])

        # --- Plot ---
        st.subheader("Kurva jackknife")
        st.caption("Sumbu-x: k (top-k node). Sumbu-y: jumlah **gen hallmark** yang ditemukan sampai k.")

        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(10, 6))
        ax = fig.add_subplot(111)
        for method in df_curve["method"].unique():
            d = df_curve[df_curve["method"] == method]
            ax.plot(d["k"].values, d["hits"].values, label=method)
        ax.set_xlabel("k (top-k)")
        ax.set_ylabel("# hallmark genes revealed")
        ax.set_title(f"Jackknife curve vs hallmark ({target_name})")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=8)
        st.pyplot(fig, use_container_width=True)

        # --- Summary table ---
        st.subheader("Ringkasan")
        st.dataframe(df_sum, use_container_width=True, hide_index=True)

        c1, c2 = st.columns(2)
        with c1:
            st.download_button(
                "Download summary (CSV)",
                data=df_sum.to_csv(index=False).encode("utf-8"),
                file_name="jackknife_summary.csv",
                mime="text/csv",
            )
        with c2:
            st.download_button(
                "Download curve (CSV)",
                data=df_curve.to_csv(index=False).encode("utf-8"),
                file_name="jackknife_curve.csv",
                mime="text/csv",
            )

        # --- Inspect top-k ---
        st.subheader("Detail top-k")
        pick_method = st.selectbox(
            "Pilih metode",
            list(df_sum["method"].values),
            index=0,
            key=f"hm_pickm::{graph_label}",
        )
        k_show = st.slider(
            "Tampilkan top-k",
            5,
            int(k_max),
            min(20, int(k_max)),
            5,
            key=f"hm_kshow::{graph_label}",
        )

        ranked = rankings.get(pick_method, [])
        top = ranked[: int(k_show)]
        rows = []
        for idx, g in enumerate(top, start=1):
            cg = _canon_gene(g)
            hms = sorted(list(gene_to_hallmarks.get(cg, set())))
            rows.append(
                {
                    "rank": idx,
                    "node": g,
                    "is_hallmark": bool(cg in positive),
                    "hallmarks": ", ".join(hms) if hms else "",
                }
            )
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
