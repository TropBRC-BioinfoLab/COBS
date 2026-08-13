from cdlib import algorithms
from .base import AlgoSpec, ParamSpec


def run(G, params):
    return algorithms.lpanni(G, threshold=params["threshold"])


SPEC = AlgoSpec(
    key="LPANNI",
    name="LPANNI",
    description="LPANNI: label propagation dengan neighbor-node influence (overlap).",
    params=[
        ParamSpec(
            key="threshold",
            label="threshold",
            ptype="float",
            default=0.1,
            min_value=0.0,
            max_value=1.0,
            step=0.01,
            help="Ambang post-processing membership pada LPANNI.",
        ),
    ],
    run=run,
)
