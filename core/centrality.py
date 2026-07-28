# core/centrality.py
from __future__ import annotations

import html
from typing import Dict, Any, Optional, Iterable, List, Sequence, Tuple

import networkx as nx
import pandas as pd


CENTRALITY_ORDER = ["degree", "closeness", "betweenness"]

# (arrow, tooltip singkat)
CENTRALITY_INFO = {
    "degree": (
        "↑",
        "Degree centrality = k(v)/(n-1). Mengukur konektivitas lokal (hub).",
    ),
    "closeness": (
        "↑",
        "Closeness = (n-1)/Σ d(v,u). Mengukur kedekatan rata-rata ke semua node (akses global).",
    ),
    "betweenness": (
        "↑",
        "Betweenness = Σ_{s≠v≠t} σ_st(v)/σ_st (normalized). Mengukur peran broker/bottleneck pada shortest paths.",
    ),
}
 

def giant_component(G: nx.Graph) -> nx.Graph:
    """Ambil komponen terhubung terbesar (undirected)."""
    if G.number_of_nodes() == 0:
        return G.copy()
    UG = G.to_undirected(as_view=True)
    if nx.is_connected(UG):
        return G.copy()
    comp = max(nx.connected_components(UG), key=len)
    return G.subgraph(comp).copy()


def compute_basic_centralities(
    G: nx.Graph,
    use_giant_component: bool = False,
) -> Dict[str, Dict[Any, float]]:
    """
    Untuk graf tidak berarah, tidak berbobot.
    - Jika graf tidak terhubung dan use_giant_component=False -> raise.
    - Jika use_giant_component=True -> hitung pada giant component.
    """
    if G.number_of_nodes() == 0:
        return {"degree": {}, "closeness": {}, "betweenness": {}}

    H = G
    UG = G.to_undirected(as_view=True)
    if not nx.is_connected(UG):
        if not use_giant_component:
            raise ValueError("Graph tidak terhubung. Set use_giant_component=True untuk memakai giant component.")
        H = giant_component(G)

    deg = nx.degree_centrality(H)                          # normalized k/(n-1)
    clo = nx.closeness_centrality(H)                       # connected graph closeness
    bet = nx.betweenness_centrality(H, normalized=True)    # normalized

    return {"degree": deg, "closeness": clo, "betweenness": bet}


def centrality_dataframe(
    G: nx.Graph,
    use_giant_component: bool = False,
) -> pd.DataFrame:
    """
    Return df dengan index=node dan kolom: degree, closeness, betweenness.
    """
    c = compute_basic_centralities(G, use_giant_component=use_giant_component)

    # pastikan index lengkap sesuai node yang dihitung (H)
    df = pd.DataFrame({k: pd.Series(v, dtype="float64") for k, v in c.items()})
    df.index.name = "node"
    df = df.reindex(columns=CENTRALITY_ORDER)
    return df


def centrality_table_html(
    df: pd.DataFrame,
    sort_by: str = "betweenness",
    topk: int = 10,
    digits: int = 4,
    show_ranks: bool = True,
) -> str:
    """
    Buat HTML table dengan tooltip pada nama metrik (header).
    """
    if sort_by not in df.columns:
        sort_by = df.columns[0]

    view = df.copy()
    view = view.sort_values(sort_by, ascending=False).head(max(1, int(topk)))

    if show_ranks:
        for col in CENTRALITY_ORDER:
            if col in view.columns:
                view[f"rank_{col}"] = view[col].rank(ascending=False, method="min").astype("int64")

    # format angka
    for col in CENTRALITY_ORDER:
        if col in view.columns:
            view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{float(x):.{digits}f}")

    # pindahkan index node jadi kolom
    view = view.reset_index()

    # header html + tooltip
    headers = []
    for col in view.columns:
        if col in CENTRALITY_INFO:
            arrow, tip = CENTRALITY_INFO[col]
            headers.append(f"<th><span title='{html.escape(tip)}'>{html.escape(col)}</span> {html.escape(arrow)}</th>")
        else:
            headers.append(f"<th>{html.escape(str(col))}</th>")

    # body
    rows = []
    for _, r in view.iterrows():
        tds = []
        for col in view.columns:
            val = r[col]
            tds.append(f"<td>{html.escape(str(val))}</td>")
        rows.append("<tr>" + "".join(tds) + "</tr>")

    css = """
    <style>
      .cent-wrap{width:100%; overflow-x:auto;}
      table.cent{width:100%; border-collapse:collapse; table-layout:auto;}
      table.cent th, table.cent td{
        border:1px solid rgba(49,51,63,0.15);
        padding:8px 10px;
        font-size:14px;
        vertical-align:top;
      }
      table.cent th{
        background:rgba(49,51,63,0.05);
        font-weight:600;
        text-align:left;
        white-space:nowrap;
      }
      table.cent td{
        white-space:nowrap;
      }
    </style>
    """

    return f"""
    {css}
    <div class="cent-wrap">
      <table class="cent">
        <thead><tr>{''.join(headers)}</tr></thead>
        <tbody>{''.join(rows)}</tbody>
      </table>
    </div>
    """


# ============================
# Community-Consideration Centrality (Overlap-aware)
# ============================
# Implementasi ini menggeneralisasi konsep "community consideration centrality" untuk
# komunitas yang overlap (node bisa menjadi anggota >1 komunitas).
#
# Ide utama (membership-weighted):
# - Untuk setiap komunitas k (himpunan node C_k), kita hitung CCC relatif terhadap C_k
#   (inside = C_k, outside = V \ C_k).
# - Untuk node yang berada di banyak komunitas, skor akhirnya adalah ekspektasi / rata-rata
#   berbobot di semua komunitas yang ia ikuti (default bobot sama: 1/|M(i)|).
#
# Catatan:
# - Untuk graph yang tidak terhubung, fungsi ini mengikuti perilaku centrality_dataframe:
#   jika use_giant_component=False -> raise, jika True -> memakai giant component.

CCC_ORDER = ["ccdc", "cccc", "ccbc"]

CCC_INFO = {
    "ccdc": (
        "↑",
        "Community-Consideration Degree Centrality (overlap). alpha mengatur trade-off komunitas vs global.",
    ),
    "cccc": (
        "↑",
        "Community-Consideration Closeness Centrality (overlap). Berbasis shortest-path distance.",
    ),
    "ccbc": (
        "↑",
        "Community-Consideration Betweenness Centrality (overlap). Berbasis betweenness raw + normalisasi pasangan node.",
    ),
}


def _node_to_communities(communities: Sequence[Iterable[Any]]) -> Dict[Any, List[int]]:
    """Map node -> list community_id (0-based)."""
    m: Dict[Any, List[int]] = {}
    for k, com in enumerate(communities):
        for n in com:
            m.setdefault(n, []).append(k)
    return m


def _ensure_connected_or_gc(G: nx.Graph, use_giant_component: bool) -> nx.Graph:
    if G.number_of_nodes() == 0:
        return G.copy()
    UG = G.to_undirected(as_view=True)
    if nx.is_connected(UG):
        return G.copy()
    if not use_giant_component:
        raise ValueError("Graph tidak terhubung. Set use_giant_component=True untuk memakai giant component.")
    return giant_component(G)


def _normalize_betweenness_raw_to_01(raw: Dict[Any, float], n: int) -> Dict[Any, float]:
    """Konversi betweenness raw (normalized=False) menjadi [0,1] seperti normalized=True (undirected)."""
    if n <= 2:
        return {k: 0.0 for k in raw.keys()}
    scale = 2.0 / ((n - 1) * (n - 2))
    return {k: float(v) * scale for k, v in raw.items()}


def ccc_dataframe(
    G: nx.Graph,
    communities: Sequence[Iterable[Any]],
    *,
    alpha: float = 0.2,
    use_giant_component: bool = False,
    compute_ccdc: bool = True,
    compute_cccc: bool = True,
    compute_ccbc: bool = False,
    betweenness_k: Optional[int] = None,
    seed: Optional[int] = None,
) -> pd.DataFrame:
    """
    Return df dengan index=node dan kolom (subset dari): ccdc, cccc, ccbc (+ n_memberships).

    Parameter:
    - communities: list komunitas (iterable node). Boleh overlap.
    - alpha in [0,1]: 0 -> fokus komunitas, 1 -> kembali ke metrik global.
    - betweenness_k: jika diisi (int > 0), memakai estimasi betweenness (sampling) NetworkX.
    """
    if G.number_of_nodes() == 0:
        df = pd.DataFrame(index=pd.Index([], name="node"))
        for col in CCC_ORDER:
            df[col] = pd.Series(dtype="float64")
        df["n_memberships"] = pd.Series(dtype="int64")
        return df

    alpha = float(alpha)
    alpha = max(0.0, min(1.0, alpha))

    H = _ensure_connected_or_gc(G, use_giant_component=use_giant_component)
    n = H.number_of_nodes()

    # Filter komunitas: hanya node yang ada di H
    comm_nodes: List[List[Any]] = []
    comm_sets: List[set] = []
    comm_sizes: List[int] = []
    for com in communities or []:
        nodes = [x for x in com if x in H]
        if len(nodes) == 0:
            continue
        s = set(nodes)
        comm_nodes.append(nodes)
        comm_sets.append(s)
        comm_sizes.append(len(nodes))

    node_to_comms = _node_to_communities(comm_nodes)
    n_memberships = {node: len(ids) for node, ids in node_to_comms.items()}

    # Global centralities (fallback untuk node yang tidak punya komunitas)
    deg_global = nx.degree_centrality(H)
    clo_global = nx.closeness_centrality(H)

    bet_raw_global: Dict[Any, float] = {}
    bet_norm_global: Dict[Any, float] = {}
    if compute_ccbc:
        # raw betweenness (Brandes); k dan seed memungkinkan aproksimasi
        kw = {}
        if betweenness_k is not None and int(betweenness_k) > 0:
            kw["k"] = int(betweenness_k)
            kw["seed"] = seed
        bet_raw_global = nx.betweenness_centrality(H, normalized=False, **kw)  # type: ignore
        bet_norm_global = _normalize_betweenness_raw_to_01(bet_raw_global, n)

    # inisialisasi agregasi
    acc_ccdc: Dict[Any, float] = {node: 0.0 for node in H.nodes()}
    acc_cccc: Dict[Any, float] = {node: 0.0 for node in H.nodes()}
    acc_ccbc: Dict[Any, float] = {node: 0.0 for node in H.nodes()}
    wsum: Dict[Any, float] = {node: 0.0 for node in H.nodes()}

    # ---------- CCDC (Degree) ----------
    if compute_ccdc:
        for k, S in enumerate(comm_sets):
            n_k = comm_sizes[k]
            denom = (n_k + alpha * (n - n_k) - 1.0)
            if denom <= 0:
                continue

            # untuk setiap node dalam komunitas tsb
            for i in S:
                m_i = n_memberships.get(i, 0)
                if m_i <= 0:
                    continue
                w = 1.0 / float(m_i)

                # internal degree: tetangga i yang ada di S
                d_in = 0
                for nbr in H.neighbors(i):
                    if nbr in S:
                        d_in += 1
                d_out = int(H.degree(i)) - d_in

                val = (d_in + alpha * d_out) / denom
                acc_ccdc[i] += w * float(val)
                wsum[i] += w

    # ---------- CCCC (Closeness) ----------
    if compute_cccc:
        # precompute list komunitas per node untuk iterasi cepat
        # (node_to_comms sudah ada)
        for i in H.nodes():
            ids = node_to_comms.get(i, [])
            if len(ids) == 0:
                continue

            # shortest path lengths dari i ke semua node
            dist = nx.single_source_shortest_path_length(H, i)
            # karena H terhubung (atau GC), dist mencakup semua node
            S_total = float(sum(dist.values()))

            m_i = len(ids)
            w = 1.0 / float(m_i)

            for k in ids:
                # guard: k mungkin berubah jika komunitas kosong dibuang
                if k < 0 or k >= len(comm_nodes):
                    continue
                nodes_k = comm_nodes[k]
                n_k = comm_sizes[k]

                numer = (n_k + alpha * (n - n_k) - 1.0)

                S_in = 0.0
                # sum distances ke node dalam komunitas k
                for v in nodes_k:
                    S_in += float(dist.get(v, 0))
                S_out = S_total - S_in

                denom = S_in + alpha * S_out
                if denom <= 0:
                    val = float("nan")
                else:
                    val = numer / denom

                acc_cccc[i] += w * float(val)
                # wsum sudah dihitung di CCDC; jika CCDC tidak dihitung, kita butuh wsum
                if not compute_ccdc:
                    wsum[i] += w

    # ---------- CCBC (Betweenness) ----------
    if compute_ccbc:
        # pasangan node (undirected): P = (n-1)(n-2)/2
        P_G = ((n - 1) * (n - 2)) / 2.0 if n >= 3 else 0.0

        # betweenness komunitas dihitung pada subgraph induced C_k
        for k, S in enumerate(comm_sets):
            n_k = comm_sizes[k]
            # P_k untuk komunitas
            P_k = ((n_k - 1) * (n_k - 2)) / 2.0 if n_k >= 3 else 0.0
            denom = P_k + alpha * (P_G - P_k)
            if denom <= 0:
                continue

            # Consistent-domain OCCBC: both global and within-community
            # contributions use shortest paths in H. Restricting sources and
            # targets to S selects the endpoint-pair subset without changing
            # the graph on which shortest paths are computed.
            if len(S) < 3:
                bet_raw_k = {node: 0.0 for node in S}
            else:
                subset = nx.betweenness_centrality_subset(
                    H,
                    sources=list(S),
                    targets=list(S),
                    normalized=False,
                )
                bet_raw_k = {node: float(subset.get(node, 0.0)) for node in S}

            for i in S:
                m_i = n_memberships.get(i, 0)
                if m_i <= 0:
                    continue
                w = 1.0 / float(m_i)

                B_in = float(bet_raw_k.get(i, 0.0))
                B_Gi = float(bet_raw_global.get(i, 0.0))
                # Numerical guard: B_in is theoretically a subset of B_Gi.
                B_in = min(max(B_in, 0.0), max(B_Gi, 0.0))
                val = (B_in + alpha * (B_Gi - B_in)) / denom
                acc_ccbc[i] += w * float(val)
                wsum[i] += w

    # finalize: jika wsum(i)>0 pakai agregat, else fallback ke global
    out = pd.DataFrame(index=pd.Index(list(H.nodes()), name="node"))
    out["n_memberships"] = pd.Series({node: int(n_memberships.get(node, 0)) for node in H.nodes()}, dtype="int64")

    if compute_ccdc:
        out["ccdc"] = pd.Series(
            {node: (acc_ccdc[node] if wsum[node] > 0 else float(deg_global.get(node, 0.0))) for node in H.nodes()},
            dtype="float64",
        )
    else:
        out["ccdc"] = pd.Series(dtype="float64")

    if compute_cccc:
        # jika node tidak punya komunitas, fallback ke closeness global
        out["cccc"] = pd.Series(
            {node: (acc_cccc[node] if wsum[node] > 0 else float(clo_global.get(node, 0.0))) for node in H.nodes()},
            dtype="float64",
        )
    else:
        out["cccc"] = pd.Series(dtype="float64")

    if compute_ccbc:
        out["ccbc"] = pd.Series(
            {node: (acc_ccbc[node] if wsum[node] > 0 else float(bet_norm_global.get(node, 0.0))) for node in H.nodes()},
            dtype="float64",
        )
    else:
        out["ccbc"] = pd.Series(dtype="float64")

    # urutkan kolom: metrics dulu, lalu n_memberships
    cols = [c for c in CCC_ORDER if c in out.columns] + ["n_memberships"]
    out = out.reindex(columns=cols)
    return out


def ccc_table_html(
    df: pd.DataFrame,
    *,
    sort_by: str = "ccbc",
    topk: int = 10,
    digits: int = 4,
    show_ranks: bool = True,
) -> str:
    """HTML table untuk CCC (overlap-aware), dengan tooltip di header."""
    # tentukan metrik yang tersedia
    metric_cols = [c for c in CCC_ORDER if c in df.columns and df[c].notna().any()]
    if len(metric_cols) == 0:
        return "<div><em>Tidak ada metrik CCC yang dihitung.</em></div>"

    if sort_by not in metric_cols:
        sort_by = metric_cols[0]

    view = df.copy()
    view = view.sort_values(sort_by, ascending=False).head(max(1, int(topk)))

    if show_ranks:
        for col in metric_cols:
            view[f"rank_{col}"] = view[col].rank(ascending=False, method="min").astype("int64")

    # format angka
    for col in metric_cols:
        view[col] = view[col].map(lambda x: "" if pd.isna(x) else f"{float(x):.{digits}f}")

    # index -> kolom
    view = view.reset_index()

    headers = []
    for col in view.columns:
        if col in CCC_INFO:
            arrow, tip = CCC_INFO[col]
            headers.append(f"<th><span title='{html.escape(tip)}'>{html.escape(col)}</span> {html.escape(arrow)}</th>")
        else:
            headers.append(f"<th>{html.escape(str(col))}</th>")

    rows = []
    for _, r in view.iterrows():
        tds = []
        for col in view.columns:
            tds.append(f"<td>{html.escape(str(r[col]))}</td>")
        rows.append("<tr>" + "".join(tds) + "</tr>")

    css = """
    <style>
      .cent-wrap{width:100%; overflow-x:auto;}
      table.cent{width:100%; border-collapse:collapse; table-layout:auto;}
      table.cent th, table.cent td{
        border:1px solid rgba(49,51,63,0.15);
        padding:8px 10px;
        font-size:14px;
        vertical-align:top;
      }
      table.cent th{
        background:rgba(49,51,63,0.05);
        font-weight:600;
        text-align:left;
        white-space:nowrap;
      }
      table.cent td{white-space:nowrap;}
    </style>
    """
    return f"""
    {css}
    <div class="cent-wrap">
      <table class="cent">
        <thead><tr>{''.join(headers)}</tr></thead>
        <tbody>{''.join(rows)}</tbody>
      </table>
    </div>
    """
