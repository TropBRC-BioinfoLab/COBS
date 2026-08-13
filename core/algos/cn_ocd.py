# core/algos/cn_ocd.py
import math
from collections import defaultdict
import networkx as nx

from .base import AlgoSpec, ParamSpec

# --- cdlib NodeClustering (fallback bila import berbeda) ---
try:
    from cdlib import NodeClustering  # type: ignore
except Exception:
    try:
        from cdlib.classes.node_clustering import NodeClustering  # type: ignore
    except Exception:
        NodeClustering = None


# =========================
# Helpers: Salton + Scc
# =========================

def _neighborsets(G: nx.Graph):
    return {u: set(G.neighbors(u)) for u in G.nodes()}

def _salton_from_sets(Nu: set, Nv: set) -> float:
    denom = math.sqrt(len(Nu) * len(Nv))
    if denom == 0:
        return 0.0
    # intersection size
    if len(Nu) > len(Nv):
        Nu, Nv = Nv, Nu
    inter = sum(1 for x in Nu if x in Nv)
    return inter / denom

def _weight_edges_salton(G: nx.Graph, weight_attr: str):
    nbr = _neighborsets(G)
    for u, v in G.edges():
        G[u][v][weight_attr] = float(_salton_from_sets(nbr[u], nbr[v]))

def _scc(Ci: set, Cj: set) -> float:
    denom = math.sqrt(len(Ci) * len(Cj))
    if denom == 0:
        return 0.0
    return len(Ci & Cj) / denom

def _overlap_ratio(Ci: set, Cj: set) -> float:
    if not Ci:
        return 0.0
    return len(Ci & Cj) / len(Ci)


# =========================
# Algorithm 2: community formation by merging edges
# (naskah: merge komunitas via edge weights high similarity)
# =========================

def _alg2_merge_by_edge_weight(G: nx.Graph, weight_attr: str, theta_sim: float):
    # init communities: C = {{v} | v in V}
    comm_of = {v: v for v in G.nodes()}         # node -> community id
    comm = {v: {v} for v in G.nodes()}          # community id -> set(nodes)

    # sort edges descending by weight
    edges = [(u, v, float(G[u][v].get(weight_attr, 0.0))) for u, v in G.edges()]
    edges.sort(key=lambda x: x[2], reverse=True)

    # merge only "high similarity edges" (pakai theta_sim, tanpa knob lain)
    for u, v, w in edges:
        if w < theta_sim:
            break

        cu = comm_of[u]
        cv = comm_of[v]
        if cu == cv:
            continue

        # merge smaller into larger
        if len(comm[cu]) < len(comm[cv]):
            cu, cv = cv, cu  # swap so cu is larger

        comm[cu] |= comm[cv]
        for x in comm[cv]:
            comm_of[x] = cu
        del comm[cv]

    communities = list(comm.values())
    return communities


# =========================
# Algorithm 3: refine communities + detect overlapping/hidden
# (naskah: merge jika overlap ratio & similarity memenuhi threshold)
# hidden community refinement untuk |Ci|=2 -> expand by common neighbors
# =========================

def _hidden_refinement_cn(G: nx.Graph, C: set):
    # versi minimal sesuai “common neighbor / hidden”
    if len(C) != 2:
        return C
    a, b = list(C)
    Na = set(G.neighbors(a))
    Nb = set(G.neighbors(b))
    return set(C) | (Na & Nb)

def _alg3_refine_merge(G: nx.Graph, communities: list, theta_overlap: float, theta_sim: float):
    comms = [set(c) for c in communities if len(c) > 0]
    comms.sort(key=len, reverse=True)

    removed = [False] * len(comms)

    for i in range(len(comms)):
        if removed[i]:
            continue
        Ci = comms[i]

        # compare with larger communities (j < i)
        for j in range(i):
            if removed[j]:
                continue
            Cj = comms[j]
            if not (Ci & Cj):
                continue

            r = _overlap_ratio(Ci, Cj)
            s = _scc(Ci, Cj)
            if (r >= theta_overlap) and (s >= theta_sim):
                # merge Ci into Cj
                Cj |= Ci
                removed[i] = True
                Ci = set()
                break

        comms[i] = Ci

    refined = [comms[k] for k in range(len(comms)) if (not removed[k]) and len(comms[k]) > 0]

    # hidden refinement (minimal, no extra knob)
    refined2 = []
    for C in refined:
        refined2.append(_hidden_refinement_cn(G, C) if len(C) == 2 else C)

    return refined2


# =========================
# Optimization: free nodes using CN(u)=max |Cv ∩ n_u|
# + final merge by Scc >= theta_sim  (naskah: threshold δ)
# =========================

def _assign_free_nodes_by_cn(G: nx.Graph, communities: list):
    comms = [set(c) for c in communities if len(c) > 0]

    # free nodes = komunitas singleton
    free_nodes = set()
    kept = []
    for C in comms:
        if len(C) == 1:
            free_nodes |= set(C)
        else:
            kept.append(C)
    comms = kept

    for u in list(free_nodes):
        scores = [len(set(G.neighbors(u)) & C) for C in comms]
        if not scores:
            comms.append({u})
            continue

        best = max(scores)
        if best == 0:
            comms.append({u})
            continue

        # CN(u) = max |Cv ∩ n_u|  -> assign to all communities that achieve max
        for idx, sc in enumerate(scores):
            if sc == best:
                comms[idx].add(u)

    return comms

def _final_merge_by_scc(communities: list, theta_sim: float):
    comms = [set(c) for c in communities if len(c) > 0]

    changed = True
    while changed:
        changed = False
        comms.sort(key=len, reverse=True)
        removed = [False] * len(comms)

        for i in range(len(comms)):
            if removed[i]:
                continue
            for j in range(i + 1, len(comms)):
                if removed[j]:
                    continue
                if _scc(comms[i], comms[j]) >= theta_sim:  # δ = theta_sim
                    comms[i] |= comms[j]
                    removed[j] = True
                    changed = True

        comms = [comms[k] for k in range(len(comms)) if not removed[k]]

    return comms


# =========================
# Streamlit Spec
# =========================

def _run(G_in: nx.Graph, params: dict):
    # undirected for neighborhood similarity
    try:
        is_directed = G_in.is_directed()
    except Exception:
        is_directed = False

    G = G_in.to_undirected() if is_directed else G_in.copy()

    weight_attr = "cn_salton"

    # ONLY 2 thresholds (sesuai permintaan)
    theta_overlap = float(params.get("theta_overlap", 0.30))
    theta_sim     = float(params.get("theta_sim", 0.40))

    # Alg 1
    _weight_edges_salton(G, weight_attr=weight_attr)

    # Alg 2
    comm2 = _alg2_merge_by_edge_weight(G, weight_attr=weight_attr, theta_sim=theta_sim)

    # Alg 3 (+ hidden)
    comm3 = _alg3_refine_merge(G, comm2, theta_overlap=theta_overlap, theta_sim=theta_sim)

    # Optimization: free nodes by CN(u), then final merge by Scc (δ=theta_sim)
    comm4 = _assign_free_nodes_by_cn(G, comm3)
    comm_final = _final_merge_by_scc(comm4, theta_sim=theta_sim)

    comms_out = [sorted(list(c)) for c in comm_final]

    method_params = {
        "theta_overlap": theta_overlap,
        "theta_sim": theta_sim,
        "delta_final": theta_sim,   # sesuai naskah (δ), tanpa knob baru
        "weight_attr": weight_attr,
        "n_communities_alg2": len(comm2),
        "n_communities_alg3": len(comm3),
        "n_communities_final": len(comm_final),
    }

    if NodeClustering is None:
        class _NC:
            def __init__(self, communities, graph, method_name, method_parameters, overlap):
                self.communities = communities
                self.graph = graph
                self.method_name = method_name
                self.method_parameters = method_parameters
                self.overlap = overlap
        return _NC(comms_out, G_in, "CN_OCD", method_params, True)

    return NodeClustering(
        comms_out,
        G_in,
        method_name="CN_OCD",
        method_parameters=method_params,
        overlap=True
    )


SPEC = AlgoSpec(
    key="CN_OCD",
    name="CN-OCD (paper)",
    description=(
        "Overlapping community detection berbasis neighborhood similarity (Salton) "
        "dengan refinement overlap + similarity sesuai naskah. "
        "Parameter hanya dua threshold: θ_overlap dan θ_sim."
    ),
    params=[
        ParamSpec("theta_overlap", "θ_overlap (|Ci∩Cj|/|Ci|)", "float", 0.30,
                  min_value=0.00, max_value=1.00, step=0.05),
        ParamSpec("theta_sim", "θ_sim (Similarity: Salton/Scc)", "float", 0.40,
                  min_value=0.00, max_value=1.00, step=0.05),
    ],
    run=_run
)