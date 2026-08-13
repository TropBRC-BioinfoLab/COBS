import numpy as np


def _infer_weight_attr(G):
    """Deteksi atribut bobot edge yang tersedia (mis. 'weight' atau 'score')."""
    try:
        for _, _, d in G.edges(data=True):
            if not isinstance(d, dict):
                continue
            if 'weight' in d:
                return 'weight'
            if 'score' in d:
                return 'score'
    except Exception:
        pass
    return 'weight'


def _objective_qov_penalized(G, comms, lambda_overlap=0.0, weight="weight"):
    """
    Overlapping modularity sederhana + penalti overlap.
    (Selaras dengan objective di modul CoverRefine)
    """
    if G.number_of_edges() == 0:
        return 0.0

    deg = dict(G.degree(weight=weight))
    m2 = float(sum(deg.values()))
    if m2 <= 0:
        return 0.0

    # membership map
    node_to_cids = {n: [] for n in G.nodes()}
    for cid, c in enumerate(comms):
        for n in c:
            if n in node_to_cids:
                node_to_cids[n].append(cid)

    # O_i
    O = {n: max(1, len(node_to_cids.get(n, []))) for n in G.nodes()}

    # Term 1: edge contribution
    edge_sum = 0.0
    for u, v, data in G.edges(data=True):
        w = float(data.get(weight, 1.0))
        cu = set(node_to_cids.get(u, []))
        cv = set(node_to_cids.get(v, []))
        inter = cu & cv
        if not inter:
            continue
        co = len(inter) / (O[u] * O[v])
        edge_sum += w * co
    term1 = (2.0 * edge_sum) / m2

    # Term 2: expected contribution via Kc
    sum_kc2 = 0.0
    for c in comms:
        if not c:
            continue
        kc = 0.0
        for i in c:
            kc += float(deg.get(i, 0.0)) / O[i]
        sum_kc2 += kc * kc
    term2 = sum_kc2 / (m2 * m2)

    qov = term1 - term2

    # penalti overlap
    overlap_total = 0.0
    for n in G.nodes():
        overlap_total += max(0, O[n] - 1)

    return float(qov - lambda_overlap * overlap_total)


def postprocess_and_memberships(G, communities, min_size=1):
    # 1) raw communities (apa adanya, hanya rapikan duplikat persis)
    raw = []
    seen = set()
    for c in communities:
        c = list(set(c))
        key = frozenset(c)
        if key not in seen:
            seen.add(key)
            raw.append(c)

    # missing nodes dihitung dari RAW (bukan setelah min_size)
    covered_raw = set()
    for c in raw:
        covered_raw.update(c)
    missing_nodes = [n for n in G.nodes() if n not in covered_raw]

    # 2) communities untuk DISPLAY (boleh difilter min_size)
    disp = [c for c in raw if len(c) >= min_size]
    disp = sorted(disp, key=lambda x: len(x), reverse=True)

    # memberships berdasarkan DISPLAY (untuk highlight overlap di tabel)
    mem = {n: 0 for n in G.nodes()}
    for c in disp:
        for n in c:
            mem[n] += 1

    return disp, mem, sorted(missing_nodes)


def compute_evaluations(clustering):
    def _try(fn):
        try:
            r = fn()
            # cdlib biasanya mengembalikan object dengan atribut .score
            return float(getattr(r, "score", r))
        except Exception:
            return float("nan")

    out = {
        # ===== kualitas struktur (cdlib) =====
        "modularity_overlap": _try(lambda: clustering.modularity_overlap()),
        "internal_edge_density": _try(lambda: clustering.internal_edge_density()),
        "conductance": _try(lambda: clustering.conductance()),
        "normalized_cut": _try(lambda: clustering.normalized_cut()),
        "avg_transitivity": _try(lambda: clustering.avg_transitivity()),
    }

    # akses graph + communities
    try:
        G = clustering.graph
    except Exception:
        G = None

    try:
        comms = clustering.communities
    except Exception:
        comms = None

    if G is None or comms is None:
        return out

    # ===== Shen modularity (Q_ov; Shen et al.) =====
    # Menggunakan formula dengan faktor 1/(O_i O_j); selaras dengan paper Shen.
    try:
        w_attr = _infer_weight_attr(G)
        out["shen_modularity"] = float(_objective_qov_penalized(G, comms, lambda_overlap=0.0, weight=w_attr))
    except Exception:
        out["shen_modularity"] = float("nan")

    # ===== ringkasan hasil =====
    out["n_communities"] = int(len(comms))
    try:
        sizes = [len(c) for c in comms if c]
        out["max_community_size"] = int(max(sizes)) if sizes else 0
    except Exception:
        out["max_community_size"] = float("nan")

    # ===== coverage & overlap stats =====
    n = G.number_of_nodes()
    if n <= 0:
        out["coverage_ratio"] = float("nan")
        out["uncovered_node_ratio"] = float("nan")
        out["overlap_node_ratio"] = float("nan")
        out["avg_memberships_per_covered_node"] = float("nan")
        out["overlap_intensity"] = float("nan")
        return out

    # membership map: node -> list community_id
    try:
        node_to_comms, _ = build_membership_maps(G, comms)
        memberships = np.array([len(v) for v in node_to_comms.values()], dtype=float)

        covered_mask = memberships > 0
        overlap_mask = memberships > 1

        out["coverage_ratio"] = float(np.mean(covered_mask))
        out["uncovered_node_ratio"] = float(np.mean(~covered_mask))
        out["overlap_node_ratio"] = float(np.mean(overlap_mask))

        out["avg_memberships_per_covered_node"] = (
            float(np.mean(memberships[covered_mask])) if np.any(covered_mask) else float("nan")
        )
        out["overlap_intensity"] = (
            float(np.mean(memberships[overlap_mask])) if np.any(overlap_mask) else float("nan")
        )
    except Exception:
        out["coverage_ratio"] = float("nan")
        out["uncovered_node_ratio"] = float("nan")
        out["overlap_node_ratio"] = float("nan")
        out["avg_memberships_per_covered_node"] = float("nan")
        out["overlap_intensity"] = float("nan")

    return out



def build_membership_maps(G, comms):
    """
    comms: list of communities (list/set of nodes), sudah untuk DISPLAY atau raw—sesuai kebutuhan mas Heru.
    Return:
      node_to_comms: dict node -> list community_id (1-based)
      primary_comm: dict node -> community_id utama (pilih komunitas terbesar yang memuat node)
    """
    node_to_comms = {n: [] for n in G.nodes()}

    for cid, c in enumerate(comms, start=1):
        for n in c:
            if n in node_to_comms:
                node_to_comms[n].append(cid)

    # primary = komunitas terbesar (agar stabil untuk graf besar)
    sizes = {cid: len(comms[cid - 1]) for cid in range(1, len(comms) + 1)}
    primary_comm = {}
    for n, cids in node_to_comms.items():
        if not cids:
            primary_comm[n] = None
        else:
            primary_comm[n] = max(cids, key=lambda x: sizes.get(x, 0))

    return node_to_comms, primary_comm


