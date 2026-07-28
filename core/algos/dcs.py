from cdlib import algorithms
from .base import AlgoSpec


def run(G, params):
    return algorithms.dcs(G)


SPEC = AlgoSpec(
    key="DCS",
    name="DCS",
    description="DCS (Divide and Conquer Strategy) untuk komunitas overlap.",
    params=[],
    run=run,
)
