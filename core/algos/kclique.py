from cdlib import algorithms
from .base import AlgoSpec, ParamSpec

def run(G, params):
    return algorithms.kclique(G, k=int(params["k"]))

SPEC = AlgoSpec(
    key="KCLIQUE",
    name="kClique",
    description="Overlapping community berbasis clique percolation.",
    params=[
        ParamSpec("k", "k", "int", 3, min_value=3, max_value=8, step=1),
    ],
    run=run,
)
