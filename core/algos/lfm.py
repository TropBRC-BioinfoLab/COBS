from cdlib import algorithms
from .base import AlgoSpec, ParamSpec

def run(G, params):
    return algorithms.lfm(G, alpha=float(params["alpha"]))

SPEC = AlgoSpec(
    key="LFM",
    name="LFM",
    description="Local Fitness Measure: komunitas dibangun dari seed dan fitness lokal.",
    params=[
        ParamSpec("alpha", "alpha", "float", 1.0, min_value=0.1, max_value=5.0, step=0.1),
    ],
    run=run,
)
