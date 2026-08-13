from cdlib import algorithms
from .base import AlgoSpec, ParamSpec


def _ensure_weight_attr(G, attr="weight", default=1.0):
    # jika user menyalakan mode weighted tetapi graph belum punya atribut weight
    for u, v, d in G.edges(data=True):
        if attr not in d:
            d[attr] = default


def run(G, params):
    use_weights = bool(params["use_weights"])
    if use_weights:
        _ensure_weight_attr(G, attr="weight", default=1.0)
        return algorithms.graph_entropy(G, weights="weight")

    # unweighted
    return algorithms.graph_entropy(G)


SPEC = AlgoSpec(
    key="GRAPH_ENT",
    name="GraphEntropy",
    description="Graph Entropy: cluster lokal dengan minimal graph entropy.",
    params=[
        ParamSpec("use_weights", "use_weights (pakai edge weight?)", "bool", False),
    ],
    run=run,
)
