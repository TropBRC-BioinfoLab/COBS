from cdlib import algorithms
from .base import AlgoSpec, ParamSpec

def run(G, params):
    return algorithms.demon(
        G,
        epsilon=float(params["epsilon"]),
        min_com_size=int(params["min_com_size"]),
    )

SPEC = AlgoSpec(
    key="DEMON",
    name="DEMON",
    description="Overlapping community detection berbasis ego-network label propagation.",
    params=[
        ParamSpec("epsilon", "epsilon", "float", 0.25, min_value=0.0, max_value=1.0, step=0.05),
        ParamSpec("min_com_size", "min_com_size", "int", 3, min_value=2, max_value=10, step=1),
    ],
    run=run,
)
