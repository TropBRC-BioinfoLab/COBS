from cdlib import algorithms
from .base import AlgoSpec, ParamSpec

def run(G, params):
    return algorithms.slpa(G, t=int(params["t"]), r=float(params["r"]))

SPEC = AlgoSpec(
    key="SLPA",
    name="SLPA",
    description="Speaker-Listener Label Propagation untuk overlapping communities.",
    params=[
        ParamSpec("t", "t (iterations)", "int", 50, min_value=10, max_value=200, step=10),
        ParamSpec("r", "r (threshold)", "float", 0.40, min_value=0.05, max_value=0.90, step=0.05),
    ],
    run=run,
)
