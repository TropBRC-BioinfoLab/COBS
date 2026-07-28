# core/algos/cover_refine.py
from __future__ import annotations

from typing import Dict, List, Set, Any, Tuple
import networkx as nx
from networkx.algorithms.community import louvain_communities, greedy_modularity_communities

from .base import AlgoSpec, ParamSpec

# cdlib NodeClustering (untuk konsisten dengan output algoritma lain)
try:
    from cdlib.classes import NodeClustering
except Exception:
    from cdlib.classes.node_clustering import NodeClustering


# =========================================================
# Helpers
# =========================================================
def _init_partition(G: nx.Graph, method: str, seed: int | None = None) -> List[Set[Any]]:
    """
    Menghasilkan komunitas disjoint awal (cover = semua node).
    """
    if method == "louvain":
        comms = louvain_communities(G, weight="weight", resolution=1, seed=seed)
        return [set(c) for c in comms]

    if method == "greedy_modularity":
        comms = greedy_modularity_communities(G, weight="weight", resolution=1)
        return [set(c) for c in comms]

    raise ValueError(f"Unknown init_partition: {method}")


def _build_memberships(comms: List[Set[Any]]) -> Dict[Any, Set[int]]:
    node2comms: Dict[Any, Set[int]] = {}
    for cid, nodes in enumerate(comms):
        for n in nodes:
            node2comms.setdefault(n, set()).add(cid)
    return node2comms


def _candidate_communities_neighbor_union(
    G: nx.Graph,
    node: Any,
    node2comms: Dict[Any, Set[int]],
) -> Set[int]:
    cand: Set[int] = set()
    for nb in G.neighbors(node):
        cand |= node2comms.get(nb, set())
    return cand


def _links_to_comm(G: nx.Graph, i: Any, comm_nodes: Set[Any]) -> int:
    """Jumlah tetangga i yang berada di komunitas comm_nodes."""
    cnt = 0
    for nb in G.neighbors(i):
        if nb in comm_nodes:
            cnt += 1
    return cnt


def _affinity_to_comm(G: nx.Graph, i: Any, comm_nodes: Set[Any]) -> float:
    """
    Proporsi bobot edge dari i yang menuju node dalam komunitas comm_nodes.
    Jika graf unweighted, bobot dianggap 1.0.
    """
    w_in = 0.0
    w_tot = 0.0
    for nb in G.neighbors(i):
        w = float(G[i][nb].get("weight", 1.0))
        w_tot += w
        if nb in comm_nodes:
            w_in += w
    return 0.0 if w_tot <= 0 else (w_in / w_tot)


def _objective_qov_penalized(
    G: nx.Graph,
    comms: List[Set[Any]],
    node2comms: Dict[Any, Set[int]],
    lambda_overlap: float = 0.0,
) -> float:
    """
    Overlapping modularity (varian belonging coefficient sederhana) + penalti overlap.

    Ringkas:
    - O_i = jumlah komunitas node i
    - Untuk edge (u,v): co(u,v) = |C(u) ∩ C(v)| / (O_u * O_v)
    - Qov = (2 * sum_{edges} w_uv * co(u,v)) / m2  -  (sum_c Kc^2) / (m2^2)
      dengan:
        m2 = sum_i k_i  (== 2m)
        Kc = sum_{i in c} k_i / O_i
    - objective = Qov - lambda_overlap * sum_i (O_i - 1)
    """
    if G.number_of_edges() == 0:
        return 0.0

    deg = dict(G.degree(weight="weight"))
    m2 = float(sum(deg.values()))
    if m2 <= 0:
        return 0.0

    # O_i
    O = {n: max(1, len(node2comms.get(n, set()))) for n in G.nodes()}

    # Term 1: edge contribution
    edge_sum = 0.0
    for u, v, data in G.edges(data=True):
        w = float(data.get("weight", 1.0))
        cu = node2comms.get(u, set())
        cv = node2comms.get(v, set())
        if not cu or not cv:
            continue
        inter = cu & cv
        if not inter:
            continue
        co = len(inter) / (O[u] * O[v])
        edge_sum += w * co

    term1 = (2.0 * edge_sum) / m2

    # Term 2: expected contribution via Kc
    sum_kc2 = 0.0
    for nodes in comms:
        if not nodes:
            continue
        kc = 0.0
        for i in nodes:
            kc += float(deg.get(i, 0.0)) / O[i]
        sum_kc2 += kc * kc

    term2 = sum_kc2 / (m2 * m2)
    qov = term1 - term2

    # Penalti overlap
    overlap_total = 0.0
    for n in G.nodes():
        overlap_total += max(0, O[n] - 1)

    return float(qov - lambda_overlap * overlap_total)


def _overlap_node_ratio_from_node2comms(G: nx.Graph, node2comms: Dict[Any, Set[int]]) -> float:
    n = G.number_of_nodes()
    if n <= 0:
        return 0.0
    ov = 0
    for u in G.nodes():
        if len(node2comms.get(u, set())) > 1:
            ov += 1
    return ov / n


def _clone_state(
    comms: List[Set[Any]],
    node2comms: Dict[Any, Set[int]],
) -> Tuple[List[Set[Any]], Dict[Any, Set[int]]]:
    comms2 = [set(c) for c in comms]
    node2comms2 = {k: set(v) for k, v in node2comms.items()}
    return comms2, node2comms2


# =========================================================
# Core refine loop (dipakai untuk run final + kalibrasi cepat)
# =========================================================
def _refine_loop(
    G: nx.Graph,
    comms: List[Set[Any]],
    node2comms: Dict[Any, Set[int]],
    params: dict,
    tau_add_value: float,
    max_iter_override: int | None = None,
) -> Tuple[List[Set[Any]], Dict[Any, Set[int]]]:
    max_iter = int(params["max_iter"]) if max_iter_override is None else int(max_iter_override)
    tau_remove = float(params["tau_remove"])
    lambda_overlap = float(params["lambda_overlap"])
    kmax = int(params["kmax"])
    neighbors_only = bool(params["neighbors_only"])
    remove_margin = float(params.get("remove_margin", 0.005))

    # Hard constraints agar tidak "nyasar"
    min_links_to_comm = int(params.get("min_links_to_comm", 1))
    min_affinity = float(params.get("min_affinity", 0.0))

    obj = _objective_qov_penalized(G, comms, node2comms, lambda_overlap=lambda_overlap)
    nodes_order = list(G.nodes())

    for _ in range(max_iter):
        changed = False

        # =========================
        # EXPANSION (ADD) - multi-add per node
        # =========================
        for i in nodes_order:
            while True:
                cur_comms = node2comms.get(i, set())
                if kmax > 0 and len(cur_comms) >= kmax:
                    break

                cand = (
                    _candidate_communities_neighbor_union(G, i, node2comms)
                    if neighbors_only
                    else set(range(len(comms)))
                )
                cand -= cur_comms
                if not cand:
                    break

                best_c = None
                best_delta = float("-inf")

                for c in cand:
                    # HARD FILTER: harus ada koneksi nyata ke komunitas c
                    links = _links_to_comm(G, i, comms[c])
                    if links < min_links_to_comm:
                        continue
                    if min_affinity > 0.0:
                        aff = _affinity_to_comm(G, i, comms[c])
                        if aff < min_affinity:
                            continue

                    # temp add
                    comms[c].add(i)
                    node2comms.setdefault(i, set()).add(c)

                    new_obj = _objective_qov_penalized(G, comms, node2comms, lambda_overlap=lambda_overlap)
                    delta = new_obj - obj

                    # rollback
                    node2comms[i].remove(c)
                    comms[c].remove(i)

                    if delta > best_delta:
                        best_delta = delta
                        best_c = c

                EPS = 1e-12
                if best_c is not None and best_delta >= (tau_add_value - EPS):
                    comms[best_c].add(i)
                    node2comms.setdefault(i, set()).add(best_c)
                    obj = _objective_qov_penalized(G, comms, node2comms, lambda_overlap=lambda_overlap)
                    changed = True
                else:
                    break

        # =========================
        # REDUCTION (REMOVE) - conservative (maks 1 remove per node per iter)
        # =========================
        for i in nodes_order:
            cur = node2comms.get(i, set())
            if len(cur) <= 1:
                continue

            best_c = None
            best_delta = float("-inf")

            for c in list(cur):
                comms[c].remove(i)
                node2comms[i].remove(c)

                new_obj = _objective_qov_penalized(G, comms, node2comms, lambda_overlap=lambda_overlap)
                delta = new_obj - obj

                node2comms[i].add(c)
                comms[c].add(i)

                if delta > best_delta:
                    best_delta = delta
                    best_c = c

            if best_c is not None and best_delta > (tau_remove + remove_margin):
                comms[best_c].remove(i)
                node2comms[i].remove(best_c)

                # safety
                if len(node2comms[i]) == 0:
                    node2comms[i].add(best_c)
                    comms[best_c].add(i)
                    continue

                obj = _objective_qov_penalized(G, comms, node2comms, lambda_overlap=lambda_overlap)
                changed = True

        if not changed:
            break

    return comms, node2comms


# =========================================================
# Auto-calibration tau_add (lock behavior)
# =========================================================
def _calibrate_tau_add(
    G: nx.Graph,
    comms0: List[Set[Any]],
    node2comms0: Dict[Any, Set[int]],
    params: dict,
) -> float:
    """
    Cari tau_add agar overlap_node_ratio mendekati target.
    Metode: binary search (trial terbatas) pada rentang [tau_add_min, tau_add_max].

    Catatan: Karena kita pakai hard constraint (min_links/min_affinity),
    monotonicity overlap_ratio terhadap tau_add biasanya lebih stabil.
    """
    target = float(params.get("target_overlap_node_ratio", 0.10))
    tol = float(params.get("calib_tol", 0.02))
    trials = int(params.get("calib_trials", 8))
    tau_lo = float(params.get("tau_add_min", -0.05))
    tau_hi = float(params.get("tau_add_max", 0.0))
    calib_max_iter = int(params.get("calib_max_iter", 4))

    best_tau = None
    best_gap = float("inf")

    lo, hi = tau_lo, tau_hi

    for _ in range(max(1, trials)):
        tau_mid = (lo + hi) / 2.0

        comms, node2comms = _clone_state(comms0, node2comms0)
        comms, node2comms = _refine_loop(
            G,
            comms,
            node2comms,
            params,
            tau_add_value=tau_mid,
            max_iter_override=calib_max_iter,
        )

        ratio = _overlap_node_ratio_from_node2comms(G, node2comms)
        gap = abs(ratio - target)

        if gap < best_gap:
            best_gap = gap
            best_tau = tau_mid

        if gap <= tol:
            break

        # arah pencarian:
        # ratio terlalu kecil => butuh lebih permisif => tau_add lebih negatif => geser ke kiri => hi=mid
        if ratio < target:
            hi = tau_mid
        else:
            lo = tau_mid

    return float(best_tau if best_tau is not None else float(params["tau_add"]))


# =========================================================
# Main run()
# =========================================================
def run(G: nx.Graph, params: dict):
    """
    CoverRefine (Modul 1, versi stabil):
    - Start dari partisi crisp (Louvain/Greedy) => cover 100% by construction.
    - Refinement jadi overlap via:
        * Expansion (ADD): terima jika ΔObjective >= tau_add
          + HARD CONSTRAINT: minimal koneksi ke komunitas target (min_links/min_affinity)
        * Reduction (REMOVE): terima jika ΔObjective > tau_remove + remove_margin
    - Optional "lock": auto_calibrate_tau_add untuk mencapai target overlap ratio.
    """
    if G.is_directed():
        G = G.to_undirected()

    init_partition = str(params["init_partition"])
    seed = int(params["seed"])
    tau_add = float(params["tau_add"])

    # 1) Start dari partisi crisp (cover 100% by construction)
    comms = _init_partition(G, init_partition, seed=seed)

    # ensure semua node masuk minimal 1 komunitas
    covered = set().union(*comms) if comms else set()
    for n in G.nodes():
        if n not in covered:
            comms.append({n})

    node2comms = _build_memberships(comms)

    # 2) Auto-calibrate (optional lock)
    params_used = dict(params)
    if bool(params.get("auto_calibrate_tau_add", False)):
        tau_add = _calibrate_tau_add(G, comms, node2comms, params_used)
        params_used["tau_add_calibrated"] = tau_add

    params_used["tau_add_used"] = tau_add

    # 3) Refinement final
    comms, node2comms = _refine_loop(G, comms, node2comms, params_used, tau_add_value=tau_add, max_iter_override=None)

    # stats (optional)
    params_used["overlap_node_ratio"] = _overlap_node_ratio_from_node2comms(G, node2comms)

    # 4) Final clean: buang komunitas kosong
    comms_final = [sorted(list(c)) for c in comms if len(c) > 0]

    try:
        return NodeClustering(
            communities=comms_final,
            graph=G,
            method_name="COVER_REFINE",
            method_parameters=params_used,
            overlap=True,
        )
    except TypeError:
        return NodeClustering(comms_final, G, "COVER_REFINE")


# =========================================================
# Algo spec (UI)
# =========================================================
SPEC = AlgoSpec(
    key="COVER_REFINE",
    name="CoverRefine (ΔQov + marginal + lock + hard-link)",
    description=(
        "Algoritma baru (Modul 1): mulai dari partisi crisp (Louvain/Greedy), "
        "lalu refinement menjadi overlap dengan constraint cover (setiap node minimal masuk 1 komunitas). "
        "Expansion memakai kriteria ΔObjective (Qov - penalti overlap) DAN hard constraint koneksi "
        "(min_links_to_comm / min_affinity) agar node tidak 'nyasar' ke komunitas yang jauh. "
        "Reduction memakai marginal removal. Opsional: auto-calibrate tau_add untuk mengunci tingkat overlap."
    ),
    params=[
        ParamSpec(
            "init_partition", "Init partition", "select", "louvain",
            options=["louvain", "greedy_modularity"],
            help="Komunitas awal (disjoint) sebelum dibuat overlap."
        ),
        ParamSpec(
            "seed", "seed", "int", 42,
            min_value=0, max_value=9999, step=1,
            help="Seed untuk init partition (khusus Louvain)."
        ),
        ParamSpec(
            "max_iter", "max_iter", "int", 10,
            min_value=1, max_value=50, step=1,
            help="Maksimum iterasi refinement final."
        ),

        # ===== manual tau_add (tetap ada) =====
        ParamSpec(
            "tau_add", "tau_add (min Δ untuk ADD)", "float", 0.0,
            min_value=-0.05, max_value=0.05, step=0.001,
            help="ADD(i,c) diterima jika ΔObjective >= tau_add. Negatif untuk memancing overlap."
        ),

        ParamSpec(
            "tau_remove", "tau_remove (min Δ untuk REMOVE)", "float", 0.0,
            min_value=0.0, max_value=1.0, step=0.01,
            help="REMOVE(i,c) diterima jika ΔObjective > tau_remove + remove_margin."
        ),
        ParamSpec(
            "remove_margin", "remove_margin (buffer REMOVE)", "float", 0.005,
            min_value=0.0, max_value=0.05, step=0.001,
            help="Buffer tambahan agar REMOVE lebih konservatif."
        ),

        ParamSpec(
            "lambda_overlap", "lambda_overlap (penalti overlap)", "float", 0.01,
            min_value=0.0, max_value=0.20, step=0.005,
            help="Penalti agar overlap tidak berlebihan (semakin besar → overlap makin dibatasi)."
        ),
        ParamSpec(
            "kmax", "kmax (maks komunitas per node)", "int", 3,
            min_value=1, max_value=10, step=1,
            help="Batas jumlah komunitas yang boleh dimiliki node."
        ),
        ParamSpec(
            "neighbors_only", "Kandidat komunitas dari neighbors saja", "bool", True,
            help="Jika True, ADD hanya dipertimbangkan ke komunitas yang dimiliki tetangga node."
        ),

        # ===== hard constraint: mencegah membership 'nyasar' =====
        ParamSpec(
            "min_links_to_comm", "min_links_to_comm", "int", 1,
            min_value=0, max_value=20, step=1,
            help="Hard constraint: minimal jumlah tetangga node yang berada di komunitas target."
        ),
        ParamSpec(
            "min_affinity", "min_affinity (0=off)", "float", 0.0,
            min_value=0.0, max_value=1.0, step=0.01,
            help="Hard constraint: minimal proporsi bobot edge node yang menuju komunitas target. 0=off."
        ),

        # ===== lock mode =====
        ParamSpec(
            "auto_calibrate_tau_add", "Lock: auto-calibrate tau_add", "bool", False,
            help="Jika True, tau_add dicari otomatis agar overlap_node_ratio mendekati target."
        ),
        ParamSpec(
            "target_overlap_node_ratio", "Target overlap node ratio", "float", 0.10,
            min_value=0.0, max_value=0.50, step=0.01,
            help="Target proporsi node overlap (membership>1)."
        ),
        ParamSpec(
            "calib_tol", "Toleransi target overlap", "float", 0.02,
            min_value=0.0, max_value=0.20, step=0.01,
            help="Toleransi absolut |ratio-target|."
        ),
        ParamSpec(
            "calib_trials", "Calib trials", "int", 8,
            min_value=1, max_value=20, step=1,
            help="Jumlah percobaan pencarian tau_add (binary search)."
        ),
        ParamSpec(
            "tau_add_min", "tau_add min (lock)", "float", -0.05,
            min_value=-0.20, max_value=0.0, step=0.01,
            help="Batas bawah pencarian tau_add (lebih negatif = lebih permisif)."
        ),
        ParamSpec(
            "tau_add_max", "tau_add max (lock)", "float", 0.0,
            min_value=-0.20, max_value=0.20, step=0.01,
            help="Batas atas pencarian tau_add."
        ),
        ParamSpec(
            "calib_max_iter", "Calib max_iter", "int", 4,
            min_value=1, max_value=15, step=1,
            help="Iterasi kecil saat kalibrasi (biar cepat)."
        ),
    ],
    run=run,
)
