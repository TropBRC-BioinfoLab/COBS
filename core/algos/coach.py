from cdlib import algorithms
from .base import AlgoSpec, ParamSpec

def run(G, params):
    return algorithms.coach(
        G,
        density_threshold=float(params["density_threshold"]),
        affinity_threshold=float(params["affinity_threshold"]),
        closeness_threshold=float(params["closeness_threshold"]),
    )

SPEC = AlgoSpec(
    key="COACH",
    name="CoAch",
    description="Core-Attachment (CoAch). Awalnya populer untuk deteksi protein complex: core padat + attachment.",
    params=[
        ParamSpec(
            "density_threshold", "density_threshold", "float", 0.70,
            min_value=0.0, max_value=1.0, step=0.05,
            help="Minimum core density."
        ),
        ParamSpec(
            "affinity_threshold", "affinity_threshold", "float", 0.225,
            min_value=0.0, max_value=1.0, step=0.025,
            help="Maximum core affinity."
        ),
        ParamSpec(
            "closeness_threshold", "closeness_threshold", "float", 0.50,
            min_value=0.0, max_value=1.0, step=0.05,
            help="Minimum neighbor closeness."
        ),
    ],
    run=run,
)
