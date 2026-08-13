from cdlib import algorithms
from .base import AlgoSpec, ParamSpec


def run(G, params):
    return algorithms.angel(
        G,
        threshold=params["threshold"],
        min_community_size=params["min_community_size"],
    )


SPEC = AlgoSpec(
    key="ANGEL",
    name="ANGEL",
    description="ANGEL (node-centric) untuk overlapping community detection.",
    params=[
        ParamSpec(
            key="threshold",
            label="threshold",
            ptype="float",
            default=0.25,
            min_value=0.0,
            max_value=1.0,
            step=0.01,
            help="Ambang kedekatan/kemiripan untuk pembentukan komunitas.",
        ),
        ParamSpec(
            key="min_community_size",
            label="min_community_size",
            ptype="int",
            default=3,
            min_value=2,
            max_value=50,
            step=1,
            help="Ukuran minimal komunitas.",
        ),
    ],
    run=run,
)
