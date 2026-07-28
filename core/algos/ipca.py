from cdlib import algorithms
from .base import AlgoSpec, ParamSpec


def _ensure_weight_attr(G, attr="weight", default=1.0):
    for u, v, d in G.edges(data=True):
        if attr not in d:
            d[attr] = default


def run(G, params):
    use_weights = bool(params["use_weights"])
    t_in = float(params["t_in"])

    if use_weights:
        _ensure_weight_attr(G, attr="weight", default=1.0)
        return algorithms.ipca(G, weights="weight", t_in=t_in)

    return algorithms.ipca(G, t_in=t_in)


SPEC = AlgoSpec(
    key="IPCA",
    name="IPCA",
    description="IPCA: varian DPClus untuk mendeteksi protein complexes/komunitas overlap.",
    params=[
        ParamSpec("use_weights", "use_weights (pakai edge weight?)", "bool", False),
        ParamSpec("t_in", "t_in", "float", 0.5, min_value=0.0, max_value=1.0, step=0.05),
    ],
    run=run,
)
