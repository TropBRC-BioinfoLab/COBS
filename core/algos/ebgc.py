from cdlib import algorithms
from .base import AlgoSpec


def run(G, params):
    return algorithms.ebgc(G)


SPEC = AlgoSpec(
    key="EBGC",
    name="EBGC",
    description="EBGC (Entropy-based clustering): tumbuhkan cluster dari seed untuk minimisasi entropy.",
    params=[],
    run=run,
)
