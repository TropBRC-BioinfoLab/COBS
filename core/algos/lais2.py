from cdlib import algorithms
from .base import AlgoSpec


def run(G, params):
    # LAIS2 tidak memerlukan parameter tambahan (sesuai API CDlib)
    return algorithms.lais2(G)


SPEC = AlgoSpec(
    key="LAIS2",
    name="LAIS2",
    description="LAIS2: overlapping community detection berbasis local-structure.",
    params=[],
    run=run,
)
