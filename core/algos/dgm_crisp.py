# core/algos/dgm_crisp.py
from __future__ import annotations

"""
DGM (Disassembly Greedy Modularity) - CRISP (views.py style) + history
=====================================================================

Ini menggantikan DGM lama (key tetap DGM_CRISP), dengan perilaku sama seperti views.py:

- Start: greedy_modularity_communities(G)
- Iterasi: iters = max(1, int(p * |V|))
- Per iterasi:
  1) Exploit: jalankan greedy dari partisi curr (via meta-graph contraction)
     -> diterima HANYA jika modularity > best-so-far (best-only acceptance)
  2) Explore (opsional, prob = exploration_prob):
       50% random node disassembly
       50% high-conductance community disassembly
     lalu greedy lagi dari hasil disassembly
     -> diterima HANYA jika modularity > best-so-far

Tambahan agar APP tidak error:
- q_history: modularity kandidat per iterasi
- best_history: running maximum (best modularity sepanjang iterasi)

Output: crisp (overlap=False)
"""

import random
from collections import defaultdict
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
from networkx.algorithms.community.quality import modularity as nx_modularity

from .base import AlgoSpec, ParamSpec


# -----------------------------
# Helpers
# -----------------------------
def _infer_weight_attr(G: nx.Graph) -> Optional[str]:
    """Infer attribute name for edge weights; return None if unweighted."""
    # quick sample
    for _, _, data in G.edges(data=True):
        if isinstance(data, dict):
            if "weight" in data:
                return "weight"
            if "score" in data:
                return "score"
        break
    # full scan fallback
    for _, _, data in G.edges(data=True):
        if isinstance(data, dict):
            if "weight" in data:
                return "weight"
            if "score" in data:
                return "score"
    return None


def _compute_modularity(
    G: nx.Graph,
    comms: Sequence[Set[Any]],
    weight: Optional[str] = None,
    resolution: float = 1.0,
) -> float:
    return float(nx_modularity(G, comms, weight=weight, resolution=resolution))


def _build_metagraph_from_partition(
    G: nx.Graph,
    comm_sets: Sequence[Set[Any]],
    weight: Optional[str] = None,
) -> nx.Graph:
    """
    Contract a partition into a meta-graph:
    - each community -> one meta node
    - meta edge weights = sum of original edge weights crossing communities
    - store within-community weight as node attr loop_weight (later added as self-loop)
    """
    node2c: Dict[Any, int] = {}
    for cid, S in enumerate(comm_sets):
        for u in S:
            node2c[u] = cid

    MG = nx.DiGraph() if G.is_directed() else nx.Graph()
    k = len(comm_sets)
    MG.add_nodes_from(range(k))

    def get_w(data: Dict[str, Any]) -> float:
        return float(data.get(weight, 1.0)) if weight else 1.0

    acc: Dict[int, Dict[int, float]] = defaultdict(lambda: defaultdict(float))

    if G.is_directed():
        for u, v, data in G.edges(data=True):
            cu, cv = node2c[u], node2c[v]
            acc[cu][cv] += get_w(data)
    else:
        for u, v, data in G.edges(data=True):
            cu, cv = node2c[u], node2c[v]
            w = get_w(data)
            if cu == cv:
                acc[cu][cu] += w
            else:
                acc[cu][cv] += w
                acc[cv][cu] += w

    # store internal weight as node attribute
    for i in range(k):
        MG.nodes[i]["loop_weight"] = acc[i][i] if (i in acc and i in acc[i]) else 0.0

    # add inter-community edges
    for i in range(k):
        for j, w in acc[i].items():
            if i == j:
                continue
            if MG.has_edge(i, j):
                MG[i][j]["weight"] += w
            else:
                MG.add_edge(i, j, weight=w)

    return MG


def _greedy_on_metagraph(MG: nx.Graph, resolution: float = 1.0) -> List[Set[int]]:
    """
    Run NetworkX greedy_modularity_communities on a meta-graph.
    Add self-loops from node attr loop_weight as edge (i,i) with weight.
    """
    MG2 = MG.copy()
    for i in MG2.nodes:
        wloop = float(MG2.nodes[i].get("loop_weight", 0.0))
        if wloop > 0:
            MG2.add_edge(i, i, weight=wloop)

    meta = greedy_modularity_communities(MG2, weight="weight", resolution=resolution)
    return [set(c) for c in meta]


def _expand_partition(comm_sets: Sequence[Set[Any]], meta_comms: Sequence[Set[int]]) -> List[Set[Any]]:
    out: List[Set[Any]] = []
    for meta in meta_comms:
        merged: Set[Any] = set()
        for cid in meta:
            merged |= set(comm_sets[cid])
        out.append(merged)
    return out


def _run_greedy_from_partition(
    G: nx.Graph,
    init_partition_sets: Sequence[Set[Any]],
    weight: Optional[str],
    resolution: float,
) -> List[Set[Any]]:
    """Greedy merge starting from an arbitrary partition via contraction."""
    MG = _build_metagraph_from_partition(G, init_partition_sets, weight=weight)
    meta = _greedy_on_metagraph(MG, resolution=resolution)
    return [set(c) for c in _expand_partition(init_partition_sets, meta)]


def _disassemble_random_node(G: nx.Graph, comm_sets: List[Set[Any]]) -> List[Set[Any]]:
    """
    Random node disassembly:
    - pick a community |C|>1
    - pick random node v* in it
    - v* becomes singleton
    - remaining nodes of C are split by connected components (undirected)
    """
    idx = [i for i, S in enumerate(comm_sets) if len(S) > 1]
    if not idx:
        return comm_sets

    i = random.choice(idx)
    v = random.choice(tuple(comm_sets[i]))

    S = set(comm_sets[i])
    S.remove(v)

    new_sets = comm_sets[:i] + comm_sets[i + 1 :]

    if S:
        H = G.subgraph(S)
        if G.is_directed():
            new_sets.append(S)
        else:
            new_sets.extend(list(nx.connected_components(H)))

    new_sets.append({v})
    return new_sets


def _internal_external_edges(G: nx.Graph, S: Set[Any]) -> Tuple[int, int]:
    inside = outside = 0
    Sset = set(S)
    for u in S:
        for v in G[u]:
            if v in Sset:
                inside += 1
            else:
                outside += 1
    if not G.is_directed():
        inside //= 2
    return inside, outside


def _disassemble_high_conductance(G: nx.Graph, comm_sets: List[Set[Any]]) -> List[Set[Any]]:
    """
    High-conductance community disassembly:
    - conductance(S) = m_out / (2*m_in + m_out)
    - pick community with highest conductance, then break it into singletons
    """
    def conductance(S: Set[Any]) -> float:
        inside, outside = _internal_external_edges(G, S)
        denom = 2 * inside + outside
        return (outside / denom) if denom > 0 else 1.0

    scores = [(conductance(S), i) for i, S in enumerate(comm_sets) if len(S) > 1]
    if not scores:
        return comm_sets

    _, i = max(scores)
    S = set(comm_sets[i])

    new_sets = comm_sets[:i] + comm_sets[i + 1 :]
    new_sets.extend([{v} for v in S])
    return new_sets


# -----------------------------
# DGM main (views.py style) + history
# -----------------------------
def dgm_views_style(
    G: nx.Graph,
    weight: Optional[str],
    resolution: float,
    p: float,
    exploration_prob: float,
    seed: Optional[int],
) -> Tuple[List[Set[Any]], float, List[float], List[float]]:
    if seed is not None:
        random.seed(int(seed))

    n = G.number_of_nodes()
    iters = max(1, int(float(p) * n))

    best = [set(c) for c in greedy_modularity_communities(G, weight=weight, resolution=resolution)]
    bestQ = _compute_modularity(G, best, weight, resolution)
    curr = [set(c) for c in best]

    q_history: List[float] = []
    best_history: List[float] = []

    for _ in range(iters):
        # Exploit
        comms = _run_greedy_from_partition(G, curr, weight, resolution)
        Q = _compute_modularity(G, comms, weight, resolution)
        q_iter = float(Q)

        if Q > bestQ:
            best = [set(c) for c in comms]
            bestQ = float(Q)
            curr = [set(c) for c in comms]

        # Explore (optional)
        if random.random() < float(exploration_prob):
            if random.random() < 0.5:
                trial = _disassemble_random_node(G, curr)
            else:
                trial = _disassemble_high_conductance(G, curr)

            comms2 = _run_greedy_from_partition(G, trial, weight, resolution)
            Q2 = _compute_modularity(G, comms2, weight, resolution)
            q_iter = float(Q2)

            if Q2 > bestQ:
                best = [set(c) for c in comms2]
                bestQ = float(Q2)
                curr = [set(c) for c in comms2]

        q_history.append(q_iter)
        best_history.append(float(bestQ))

    best_sorted = sorted([frozenset(c) for c in best], key=len, reverse=True)
    return [set(c) for c in best_sorted], float(bestQ), q_history, best_history


def run(G: nx.Graph, params: Dict[str, Any]):
    weight_attr = _infer_weight_attr(G)

    comms, bestQ, q_hist, best_hist = dgm_views_style(
        G,
        weight=weight_attr,
        resolution=float(params.get("resolution", 1.0)),
        p=float(params.get("p", 5.0)),
        exploration_prob=float(params.get("exploration_prob", 0.3)),
        seed=params.get("seed", 13579),
    )

    # return NodeClustering (cdlib)
    try:
        from cdlib.classes.node_clustering import NodeClustering
        return NodeClustering(
            communities=[list(c) for c in comms],
            graph=G,
            method_name="DGM_CRISP",
            method_parameters={
                "best_modularity": float(bestQ),
                "q_history": q_hist,
                "best_history": best_hist,
                "p": float(params.get("p", 5.0)),
                "exploration_prob": float(params.get("exploration_prob", 0.3)),
                "seed": params.get("seed", 13579),
                "resolution": float(params.get("resolution", 1.0)),
                "weight_attr_used": weight_attr,
                "iters": int(max(1, int(float(params.get("p", 5.0)) * G.number_of_nodes()))),
            },
            overlap=False,
        )
    except Exception:
        # fallback for environments without cdlib
        from types import SimpleNamespace
        return SimpleNamespace(
            communities=[list(c) for c in comms],
            method_name="DGM_CRISP",
            method_parameters={
                "best_modularity": float(bestQ),
                "q_history": q_hist,
                "best_history": best_hist,
            },
        )


# -----------------------------
# Streamlit AlgoSpec (key tetap DGM_CRISP)
# -----------------------------
SPEC = AlgoSpec(
    key="DGM_CRISP",
    name="DGM (crisp) - views.py implementation",
    description="Port 1:1 DGM dari views.py (best-only acceptance; random node / high conductance) + history.",
    params=[
        ParamSpec(
            key="p",
            label="p (multiplier): iters = p × |V|",
            ptype="float",
            default=5.0,
            min_value=0.1,
            max_value=30.0,
            step=0.1,
            help="Jumlah iterasi = p × jumlah node.",
        ),
        ParamSpec(
            key="exploration_prob",
            label="exploration_prob",
            ptype="float",
            default=0.3,
            min_value=0.0,
            max_value=1.0,
            step=0.01,
            help="Peluang menjalankan disassembly per iterasi.",
        ),
        ParamSpec(
            key="seed",
            label="seed (optional)",
            ptype="int",
            default=13579,
            min_value=0,
            max_value=10**9,
            step=1,
            help="Seed random untuk reproducibility.",
        ),
        ParamSpec(
            key="resolution",
            label="resolution (gamma)",
            ptype="float",
            default=1.0,
            min_value=0.1,
            max_value=5.0,
            step=0.05,
            help="Parameter resolution modularity (gamma).",
        ),
    ],
    run=run,
)