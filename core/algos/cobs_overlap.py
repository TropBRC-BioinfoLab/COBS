# core/algos/cobs_overlap.py
# Crisp-to-Overlap via Boundary Similarity (COBS)
# Step 1: Greedy Modularity (crisp) -> pilih partisi modularitas bagus
# Step 2: Boundary node kandidat overlap
# Step 3: Profil komunitas berbasis union tetangga
# Step 4: Assign overlap via similarity node->community (Salton-like)
# Step 5: Merge komunitas bila overlap ratio + Scc melewati threshold
#
# Parameter UI: θ_sim, θ_overlap, dan pilihan similarity metric

import math
from collections import defaultdict, Counter
import networkx as nx

from .base import AlgoSpec, ParamSpec

# --- cdlib NodeClustering (fallback) ---
try:
    from cdlib import NodeClustering  # type: ignore
except Exception:
    try:
        from cdlib.classes.node_clustering import NodeClustering  # type: ignore
    except Exception:
        NodeClustering = None


def _infer_weight_attr(G: nx.Graph):
    """Ikuti pola app.py: prefer 'weight', fallback 'score'."""
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


def _normalize_metric(name: str | None) -> str:
    """
    Normalisasi nama metrik similarity (case-insensitive).
    Output canonical: cosine | jaccard | dice | overlap | aa | ra
    """
    if not name:
        return "cosine"
    s = str(name).strip().lower()

    aliases = {
        "salton": "cosine",
        "cos": "cosine",
        "cosine": "cosine",
        "ochiai": "cosine",

        "jaccard": "jaccard",
        "tanimoto": "jaccard",

        "dice": "dice",
        "sorensen": "dice",
        "sørensen": "dice",
        "sorensen-dice": "dice",
        "sorensen–dice": "dice",
        "sorensen_dice": "dice",

        "overlap": "overlap",
        "simpson": "overlap",
        "szymkiewicz-simpson": "overlap",
        "szymkiewicz–simpson": "overlap",
        "overlap_coefficient": "overlap",

        "aa": "aa",
        "adamic-adar": "aa",
        "adamic–adar": "aa",
        "adamic_adar": "aa",

        "ra": "ra",
        "resource-allocation": "ra",
        "resource_allocation": "ra",
    }
    return aliases.get(s, s)


def _intersection_count(A: set, B: set) -> int:
    """|A ∩ B| (iterasi set yang lebih kecil)."""
    if len(A) > len(B):
        A, B = B, A
    return sum(1 for x in A if x in B)


def _weight_value(deg: int, kind: str) -> float:
    """
    Bobot fitur untuk metrik AA/RA.
    - AA: 1 / log(1 + deg)
    - RA: 1 / deg
    """
    if deg <= 0:
        return 0.0
    if kind == "aa":
        # log1p(deg) aman untuk deg=1 (log(2))
        den = math.log1p(deg)
        return 0.0 if den <= 0 else (1.0 / den)
    if kind == "ra":
        return 1.0 / float(deg)
    return 0.0


def _norm_w2(S: set, deg_map: dict, kind: str) -> float:
    """∑ w(x)^2 untuk x ∈ S (dipakai sebagai ||v||^2 pada weighted cosine)."""
    s = 0.0
    for x in S:
        w = _weight_value(int(deg_map.get(x, 0)), kind)
        s += w * w
    return s


def _weighted_dot_w2(A: set, B: set, deg_map: dict, kind: str) -> float:
    """∑ w(x)^2 untuk x ∈ (A ∩ B)."""
    if len(A) > len(B):
        A, B = B, A
    dot = 0.0
    for x in A:
        if x in B:
            w = _weight_value(int(deg_map.get(x, 0)), kind)
            dot += w * w
    return dot

def _intersection_count_and_dot_w2(A: set, B: set, deg_map: dict, kind: str) -> tuple[int, float]:
    """Return (|A∩B|, ∑ w(x)^2 untuk x ∈ A∩B) dengan sekali iterasi."""
    if len(A) > len(B):
        A, B = B, A
    inter = 0
    dot = 0.0
    for x in A:
        if x in B:
            inter += 1
            w = _weight_value(int(deg_map.get(x, 0)), kind)
            dot += w * w
    return inter, dot


def _sim_from_stats(
    metric: str,
    nA: int,
    nB: int,
    inter: int,
    *,
    dot_w2: float | None = None,
    normA_w2: float | None = None,
    normB_w2: float | None = None,
) -> float:
    """
    Hitung similarity berdasarkan statistik sederhana.
    Semua output dibatasi di [0, 1] (untuk AA/RA kita pakai weighted cosine).
    """
    if inter <= 0 or nA <= 0 or nB <= 0:
        return 0.0

    if metric == "cosine":
        denom = math.sqrt(nA * nB)
        return 0.0 if denom == 0 else (inter / denom)

    if metric == "jaccard":
        denom = nA + nB - inter
        return 0.0 if denom == 0 else (inter / denom)

    if metric == "dice":
        denom = nA + nB
        return 0.0 if denom == 0 else ((2.0 * inter) / denom)

    if metric == "overlap":
        denom = min(nA, nB)
        return 0.0 if denom == 0 else (inter / denom)

    if metric in ("aa", "ra"):
        if dot_w2 is None or normA_w2 is None or normB_w2 is None:
            return 0.0
        denom = math.sqrt(normA_w2 * normB_w2)
        return 0.0 if denom == 0 else (dot_w2 / denom)

    # fallback: anggap cosine
    denom = math.sqrt(nA * nB)
    return 0.0 if denom == 0 else (inter / denom)


# Backward-compat (nama lama)
def _salton_sets(A: set, B: set) -> float:
    return _sim_from_stats("cosine", len(A), len(B), _intersection_count(A, B))


def _scc(Ci: set, Cj: set) -> float:
    # default lama: cosine pada node-set komunitas
    return _sim_from_stats("cosine", len(Ci), len(Cj), _intersection_count(Ci, Cj))


def _overlap_ratio(Ci: set, Cj: set) -> float:
    if not Ci:
        return 0.0
    return len(Ci & Cj) / len(Ci)


def _compute_modularity_crisp(G: nx.Graph, part: list[set], weight_attr: str | None):
    """Modularity hanya valid untuk partisi crisp; dipakai sebagai info."""
    try:
        H = G if not G.is_directed() else G.to_undirected()
        return nx.algorithms.community.quality.modularity(H, part, weight=weight_attr)
    except Exception:
        return None


def _run(G_in: nx.Graph, params: dict):
    # --- knobs ---
    theta_sim = float(params.get("theta_sim", 0.40))
    theta_overlap = float(params.get("theta_overlap", 0.30))
    sim_metric = _normalize_metric(params.get("sim_metric", "cosine"))

    # gunakan undirected untuk neighborhood
    G = G_in.to_undirected() if getattr(G_in, "is_directed", lambda: False)() else G_in
    G = G.copy()

    # precompute degree map (dipakai untuk AA/RA)
    deg_map = dict(G.degree()) if sim_metric in ("aa", "ra") else None

    # =========================
    # Step 1: Crisp (Greedy Modularity)
    # =========================
    w_attr = _infer_weight_attr(G)
    crisp_comms = list(nx.algorithms.community.greedy_modularity_communities(G, weight=w_attr))
    crisp_comms = [set(C) for C in crisp_comms if len(C) > 0]

    # map node -> crisp community id
    cid = {}
    for i, C in enumerate(crisp_comms):
        for u in C:
            cid[u] = i

    crisp_mod = _compute_modularity_crisp(G, crisp_comms, w_attr)

    # =========================
    # Step 2: Boundary nodes + candidate communities
    # =========================
    B = set()
    K = defaultdict(set)  # node -> set of candidate community ids (tetangga lintas komunitas)

    for u, v in G.edges():
        cu = cid.get(u, None)
        cv = cid.get(v, None)
        if cu is None or cv is None or cu == cv:
            continue
        B.add(u); B.add(v)
        K[u].add(cv)
        K[v].add(cu)

    # =========================
    # Step 3: Community profile N(Ci) = union_{x in Ci} N(x)
    # =========================
    neigh = {u: set(G.neighbors(u)) for u in G.nodes()}

    Ncomm = []
    Ncomm_size = []
    for C in crisp_comms:
        s = set()
        for x in C:
            s |= neigh.get(x, set())
        Ncomm.append(s)
        Ncomm_size.append(len(s))

    # precompute norma (||v||^2) untuk metrik AA/RA (weighted cosine)
    Ncomm_norm_w2 = None
    if sim_metric in ("aa", "ra") and deg_map is not None:
        Ncomm_norm_w2 = [_norm_w2(s, deg_map, sim_metric) for s in Ncomm]

    # cache norma N(u) untuk boundary nodes (biar tidak hitung berulang)
    Nu_norm_cache = {}

    # =========================
# Step 4: Overlap injection by node->community similarity
# sim(u, Cj) berbasis set similarity antara N(u) dan N(Cj).
# Opsi metrik:
#   - cosine  : |A∩B| / sqrt(|A||B|)   (Salton/Ochiai)
#   - jaccard : |A∩B| / |A∪B|
#   - dice    : 2|A∩B| / (|A|+|B|)
#   - overlap : |A∩B| / min(|A|,|B|)
#   - aa      : weighted cosine, w(x)=1/log(1+deg(x))
#   - ra      : weighted cosine, w(x)=1/deg(x)
# only for boundary nodes dan candidate communities K(u)
# =========================
    comms = [set(C) for C in crisp_comms]  # start from crisp

    for u in B:
        Nu = neigh.get(u, set())
        if not Nu:
            continue
        nNu = len(Nu)

        # cache norm untuk AA/RA (weighted cosine)
        Nu_norm = None
        if sim_metric in ("aa", "ra") and deg_map is not None:
            Nu_norm = Nu_norm_cache.get(u)
            if Nu_norm is None:
                Nu_norm = _norm_w2(Nu, deg_map, sim_metric)
                Nu_norm_cache[u] = Nu_norm

        for j in K[u]:
            NCj = Ncomm[j]
            nNC = Ncomm_size[j]
            if nNC == 0:
                continue

            if sim_metric in ("aa", "ra") and deg_map is not None and Ncomm_norm_w2 is not None:
                inter, dot_w2 = _intersection_count_and_dot_w2(Nu, NCj, deg_map, sim_metric)
                sim = _sim_from_stats(
                    sim_metric, nNu, nNC, inter,
                    dot_w2=dot_w2,
                    normA_w2=Nu_norm,
                    normB_w2=Ncomm_norm_w2[j],
                )
            else:
                inter = _intersection_count(Nu, NCj)
                sim = _sim_from_stats(sim_metric, nNu, nNC, inter)

            if sim >= theta_sim:
                comms[j].add(u)

    # =========================
    # Step 5: Merge communities if overlap ratio + Scc >= thresholds
        # merge rule (paper-like):
        #   |Ci∩Cj|/|Ci| >= theta_overlap  AND  Scc(Ci,Cj) >= theta_sim
        # =========================
    comms.sort(key=len, reverse=True)
    removed = [False] * len(comms)

    for i in range(len(comms)):
        if removed[i]:
            continue
        Ci = comms[i]
        for j in range(i):  # compare to larger communities
            if removed[j]:
                continue
            Cj = comms[j]
            inter_set = Ci & Cj
            if not inter_set:
                continue

            inter = len(inter_set)
            ov = inter / len(Ci) if len(Ci) > 0 else 0.0

            if sim_metric in ("aa", "ra") and deg_map is not None:
                # weighted cosine (AA/RA) pada node-set komunitas
                dot_w2 = 0.0
                for x in inter_set:
                    w = _weight_value(int(deg_map.get(x, 0)), sim_metric)
                    dot_w2 += w * w
                sim_cc = _sim_from_stats(
                    sim_metric, len(Ci), len(Cj), inter,
                    dot_w2=dot_w2,
                    normA_w2=_norm_w2(Ci, deg_map, sim_metric),
                    normB_w2=_norm_w2(Cj, deg_map, sim_metric),
                )
            else:
                sim_cc = _sim_from_stats(sim_metric, len(Ci), len(Cj), inter)

            if (ov >= theta_overlap) and (sim_cc >= theta_sim):
                Cj |= Ci
                removed[i] = True
                break

    comms_final = [comms[k] for k in range(len(comms)) if not removed[k] and len(comms[k]) > 0]

    # --- stats for UI/debug (no extra knobs) ---
    mem = defaultdict(int)
    for C in comms_final:
        for u in C:
            mem[u] += 1
    overlap_nodes = [u for u, c in mem.items() if c >= 2]
    dist = dict(sorted(Counter(mem.values()).items()))

    method_params = {
        "theta_sim": theta_sim,
        "theta_overlap": theta_overlap,
        "sim_metric": sim_metric,
        "base_crisp": "greedy_modularity",
        "crisp_k": len(crisp_comms),
        "crisp_modularity": crisp_mod,
        "boundary_nodes": len(B),
        "final_k": len(comms_final),
        "overlap_nodes": len(overlap_nodes),
        "membership_count_dist": dist,
        "weight_attr_for_crisp": w_attr,
    }

    comms_out = [sorted(list(C)) for C in comms_final]

    if NodeClustering is None:
        class _NC:
            def __init__(self, communities, graph, method_name, method_parameters, overlap):
                self.communities = communities
                self.graph = graph
                self.method_name = method_name
                self.method_parameters = method_parameters
                self.overlap = overlap
        return _NC(comms_out, G_in, "COBS", method_params, True)

    return NodeClustering(
        comms_out,
        G_in,
        method_name="COBS",
        method_parameters=method_params,
        overlap=True
    )


SPEC = AlgoSpec(
    key="COBS",
    name="Crisp→Overlap (Boundary Similarity)",
    description=(
        "Step1 greedy modularity (crisp) → Step2 boundary-node candidates → "
        "Step3 community neighbor profile → Step4 overlap by node→community similarity (Salton-like) → "
        "Step5 merge if overlap ratio & Scc pass thresholds. "
        "Parameter: sim_metric, θ_sim, dan θ_overlap."
    ),
    params=[
        ParamSpec(
            "sim_metric", "Similarity metric", "str", "cosine",
            help=(
                "Pilihan: cosine (Salton/Ochiai), jaccard, dice (Sørensen–Dice), "
                "overlap (Szymkiewicz–Simpson), aa (Adamic–Adar), ra (Resource Allocation). "
                "Dipakai untuk Step4 & Step5."
            )
        ),
        ParamSpec(
            "theta_sim", "θ_sim (similarity threshold)", "float", 0.40,
            min_value=0.00, max_value=1.00, step=0.05,
            help="Dipakai untuk (i) overlap assignment node→community dan (ii) community similarity Scc saat merge."
        ),
        ParamSpec(
            "theta_overlap", "θ_overlap (|Ci∩Cj|/|Ci|)", "float", 0.30,
            min_value=0.00, max_value=1.00, step=0.05,
            help="Ambang overlap ratio untuk merge komunitas pada Step5."
        ),
    ],
    run=_run
)