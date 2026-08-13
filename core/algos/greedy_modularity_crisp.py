"""Greedy Modularity (CNM) - crisp community detection.

This is the standard greedy modularity maximization algorithm (Clauset-Newman-Moore)
from NetworkX, exposed as an AlgoSpec so it can be run & compared against DGM.

Output is CRISP (non-overlapping) communities.
"""

from __future__ import annotations

import networkx as nx

from .base import AlgoSpec, ParamSpec


def _infer_weight_attr(G: nx.Graph):
    """Heuristik: pakai 'weight' jika ada, else 'score', else None."""
    try:
        for _, _, data in G.edges(data=True):
            if isinstance(data, dict):
                if "weight" in data:
                    return "weight"
                if "score" in data:
                    return "score"
            break
    except Exception:
        pass
    return None


def run(G, params):
    # NetworkX greedy modularity bekerja pada graf undirected.
    UG = G.to_undirected(as_view=False) if getattr(G, "is_directed", lambda: False)() else G

    weight_attr = _infer_weight_attr(UG)
    resolution = float(params.get("resolution", 1.0))

    # Import dari modul yang sama dengan referensi mas Heru
    from networkx.algorithms.community.modularity_max import greedy_modularity_communities

    comms = greedy_modularity_communities(
        UG,
        weight=weight_attr,
        resolution=resolution,
    )

    comms_list = [list(c) for c in comms]

    # Return CDLIB-like clustering object (agar kompatibel dengan pipeline app)
    try:
        from cdlib.classes import NodeClustering

        clustering = NodeClustering(
            communities=comms_list,
            graph=UG,
            method_name="GREEDY_MODULARITY",
            method_parameters={**params, "weight_attr": weight_attr},
            overlap=False,
        )
        return clustering
    except Exception:
        class _Simple:
            def __init__(self, communities, graph):
                self.communities = communities
                self.graph = graph
                self.overlap = False
                self.method_name = "GREEDY_MODULARITY"
                self.method_parameters = dict(params)

        return _Simple(comms_list, UG)


SPEC = AlgoSpec(
    key="GREEDY_MODULARITY",
    name="Greedy Modularity (crisp)",
    description=(
        "Greedy modularity maximization (CNM) dari NetworkX. Output crisp (non-overlap). "
        "Cocok untuk baseline pembanding DGM."
    ),
    params=[
        ParamSpec(
            key="resolution",
            label="resolution (gamma)",
            ptype="float",
            default=1.0,
            min_value=0.1,
            max_value=5.0,
            step=0.1,
            help="Parameter resolusi modularity (gamma). >1 cenderung komunitas lebih kecil.",
        ),
    ],
    run=run,
)
