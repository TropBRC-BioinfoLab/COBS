from cdlib import algorithms
from .base import AlgoSpec, ParamSpec


def run(G, params):
    return algorithms.core_expansion(
        G,
        tolerance=float(params["tolerance"]),
    )


SPEC = AlgoSpec(
    key="CORE_EXP",
    name="CoreExpansion",
    description="Core Expansion: deteksi core tiap komunitas lalu ekspansi (overlap).",
    params=[
        ParamSpec("tolerance", "tolerance", "float", 0.5, min_value=0.0, max_value=1.0, step=0.05),
    ],
    run=run,
)
