import matplotlib.pyplot as plt
import networkx as nx
from itertools import combinations
import numpy as np
from matplotlib.offsetbox import AnnotationBbox, DrawingArea
from matplotlib.patches import Wedge, Circle




def draw_graph_with_overlap(G, memberships: dict):
    """
    memberships: dict node -> number of memberships
    Overlap node: memberships > 1
    """
    fig = plt.figure(figsize=(10, 7))
    pos = nx.spring_layout(G, seed=42)

    overlap_nodes = [n for n, k in memberships.items() if k > 1]
    normal_nodes = [n for n, k in memberships.items() if k <= 1]

    nx.draw_networkx_edges(G, pos, alpha=0.4)

    # normal nodes (default)
    nx.draw_networkx_nodes(G, pos, nodelist=normal_nodes, node_size=300, alpha=0.9)

    # overlap nodes (lebih besar + bentuk beda)
    nx.draw_networkx_nodes(G, pos, nodelist=overlap_nodes, node_size=550, node_shape="s", alpha=0.95)

    nx.draw_networkx_labels(G, pos, font_size=9)
    plt.axis("off")
    return fig


def _topk_by_degree(G, k):
    if k <= 0:
        return []
    deg = dict(G.degree())
    return [n for n, _ in sorted(deg.items(), key=lambda x: x[1], reverse=True)[:k]]

def draw_overview_community_map(G, primary_comm, node_to_comms, seed=42, node_size=520, font_size=10):
    pos = nx.spring_layout(G, seed=seed)

    overlap_nodes = [n for n, cids in node_to_comms.items() if len(cids) > 1]
    non_overlap_nodes = [n for n in G.nodes() if n not in overlap_nodes]

    fig = plt.figure(figsize=(11, 7))
    ax = plt.gca()

    nx.draw_networkx_edges(G, pos, alpha=0.25, ax=ax)

    # node non-overlap (warna komunitas utama)
    base_colors = []
    for n in non_overlap_nodes:
        base_colors.append(_cid_color(primary_comm.get(n)))

    nx.draw_networkx_nodes(
        G, pos,
        nodelist=non_overlap_nodes,
        node_size=node_size,
        node_color=base_colors,
        alpha=0.95,
        ax=ax
    )

    # node overlap sebagai donut/pie
    for n in overlap_nodes:
        x, y = pos[n]
        cids = node_to_comms.get(n, [])
        _add_pie_node(ax, x, y, cids, node_size=node_size, donut=0.55)

    # label SEMUA node
    labels = {n: str(n) for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=font_size, ax=ax)

    ax.set_title("Overview: community map (warna=komunitas utama, pie=overlap membership)")
    ax.axis("off")
    return fig


def draw_subgraph_for_community(G, comms, selected_cid, node_to_comms, seed=42,
                               node_size=520, font_size=10, add_1hop=True):
    if selected_cid < 1 or selected_cid > len(comms):
        selected_cid = 1

    C = set(comms[selected_cid - 1])
    nodes = set(C)

    if add_1hop:
        for n in list(C):
            nodes.update(G.neighbors(n))

    SG = G.subgraph(nodes).copy()
    pos = nx.spring_layout(SG, seed=seed)

    overlap_nodes = [n for n in SG.nodes() if len(node_to_comms.get(n, [])) > 1]
    in_comm_nodes = [n for n in SG.nodes() if n in C]
    out_nodes = [n for n in SG.nodes() if n not in C]

    fig = plt.figure(figsize=(11, 7))
    ax = plt.gca()

    nx.draw_networkx_edges(SG, pos, alpha=0.25, ax=ax)

    # node luar komunitas (abu-abu)
    if out_nodes:
        nx.draw_networkx_nodes(
            SG, pos,
            nodelist=out_nodes,
            node_size=node_size,
            node_color="lightgray",
            alpha=0.55,
            ax=ax
        )

    # node dalam komunitas (warna default)
    nx.draw_networkx_nodes(
        SG, pos,
        nodelist=in_comm_nodes,
        node_size=node_size,
        alpha=0.95,
        ax=ax
    )

    # overlap outline (biar tetap terlihat overlap meskipun warnanya sama)
    if overlap_nodes:
        nx.draw_networkx_nodes(
            SG, pos,
            nodelist=overlap_nodes,
            node_size=node_size * 1.10,
            node_color="none",
            edgecolors="black",
            linewidths=2.0,
            ax=ax
        )

    # label SEMUA node di subgraph
    labels = {n: str(n) for n in SG.nodes()}
    nx.draw_networkx_labels(SG, pos, labels=labels, font_size=font_size, ax=ax)

    ax.set_title(f"Komunitas {selected_cid}: subgraph (abu=luar, outline=overlap)")
    ax.axis("off")
    return fig


def draw_community_interaction_graph(G, comms, seed=42, node_size=900, font_size=11):
    comm_sets = [set(c) for c in comms]
    nC = len(comm_sets)

    CG = nx.Graph()
    for i in range(nC):
        CG.add_node(i + 1)

    for i, j in combinations(range(nC), 2):
        w = len(comm_sets[i].intersection(comm_sets[j]))
        if w > 0:
            CG.add_edge(i + 1, j + 1, weight=w)

    pos = nx.spring_layout(CG, seed=seed)

    weights = [CG.edges[e]["weight"] for e in CG.edges()]
    widths = [1.0 + 0.8 * w for w in weights] if weights else []

    fig = plt.figure(figsize=(9, 7))
    ax = plt.gca()

    nx.draw_networkx_edges(CG, pos, width=widths, alpha=0.35, ax=ax)
    nx.draw_networkx_nodes(CG, pos, node_size=node_size, alpha=0.95, ax=ax)

    # label SEMUA supernode (komunitas)
    labels = {n: f"C{n}" for n in CG.nodes()}
    nx.draw_networkx_labels(CG, pos, labels=labels, font_size=font_size, ax=ax)

    ax.set_title("Interaksi antar komunitas (edge weight = jumlah overlap nodes)")
    ax.axis("off")
    return fig


def _cid_color(cid):
    if cid is None:
        return (0.75, 0.75, 0.75, 1.0)
    return plt.cm.tab20((cid - 1) % 20)

def _add_pie_node(ax, x, y, cids, node_size, donut=0.55, edgecolor="black", lw=1.6):
    fig = ax.figure
    dpi = fig.dpi

    dia_pt = max(10.0, np.sqrt(float(node_size)))        # points
    dia_px = int(dia_pt * dpi / 72.0)                    # px

    da = DrawingArea(dia_px, dia_px, 0, 0)
    r = dia_px / 2.0
    cx = cy = r

    k = max(1, len(cids))
    for i, cid in enumerate(cids):
        theta1 = 360.0 * i / k
        theta2 = 360.0 * (i + 1) / k
        da.add_artist(Wedge((cx, cy), r, theta1, theta2, facecolor=_cid_color(cid), edgecolor="none"))

    if donut and donut > 0:
        da.add_artist(Circle((cx, cy), r * donut, facecolor="white", edgecolor="none"))

    da.add_artist(Circle((cx, cy), r, facecolor="none", edgecolor=edgecolor, linewidth=lw))

    ab = AnnotationBbox(da, (x, y), frameon=False, box_alignment=(0.5, 0.5))
    ax.add_artist(ab)
