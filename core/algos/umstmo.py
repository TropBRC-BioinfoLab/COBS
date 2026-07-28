from cdlib import algorithms
from .base import AlgoSpec

def run(G, params):
    return algorithms.umstmo(G)

SPEC = AlgoSpec(
    key="UMSTMO",
    name="UMSTMO",
    description="UMSTMO: overlapping communities berbasis union of maximum spanning trees.",
    params=[],
    run=run,
)
