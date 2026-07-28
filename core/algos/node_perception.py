from cdlib import algorithms
from .base import AlgoSpec, ParamSpec

def run(G, params):
    return algorithms.node_perception(
        G,
        threshold=float(params["threshold"]),
        overlap_threshold=float(params["overlap_threshold"]),
    )

SPEC = AlgoSpec(
    key="NODEP",
    name="NodePerception",
    description="Overlapping community berbasis persepsi node; cenderung menghasilkan banyak komunitas kecil.",
    params=[
        ParamSpec("threshold", "threshold", "float", 0.25, min_value=0.0, max_value=1.0, step=0.05),
        ParamSpec("overlap_threshold", "overlap_threshold", "float", 0.25, min_value=0.0, max_value=1.0, step=0.05),
    ],
    run=run,
)
