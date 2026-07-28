from cdlib import algorithms
from .base import AlgoSpec, ParamSpec


def _ensure_weight_attr(G, attr="weight", default=1.0):
    for u, v, d in G.edges(data=True):
        if attr not in d:
            d[attr] = default


def run(G, params):
    use_weights = bool(params["use_weights"])
    if use_weights:
        _ensure_weight_attr(G, attr="weight", default=1.0)
        return algorithms.dpclus(
            G,
            weights="weight",
            d_threshold=float(params["d_threshold"]),
            cp_threshold=float(params["cp_threshold"]),
        )

    return algorithms.dpclus(
        G,
        d_threshold=float(params["d_threshold"]),
        cp_threshold=float(params["cp_threshold"]),
    )


SPEC = AlgoSpec(
    key="DPCLUS",
    name="DPClus",
    description="DPClus: density-peak style untuk komunitas overlap (sering dipakai konteks protein complex).",
    params=[
        ParamSpec("use_weights", "use_weights (pakai edge weight?)", "bool", False),
        ParamSpec("d_threshold", "d_threshold", "float", 0.9, min_value=0.0, max_value=1.0, step=0.05),
        ParamSpec("cp_threshold", "cp_threshold", "float", 0.5, min_value=0.0, max_value=1.0, step=0.05),
    ],
    run=run,
)
