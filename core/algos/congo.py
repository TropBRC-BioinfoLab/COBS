from cdlib import algorithms
from .base import AlgoSpec, ParamSpec


def run(G, params):
    return algorithms.congo(
        G,
        number_communities=params["number_communities"],
        height=params["height"],
    )


SPEC = AlgoSpec(
    key="CONGO",
    name="CONGO",
    description="CONGO: varian CONGA yang lebih efisien (graph-contraction).",
    params=[
        ParamSpec(
            key="number_communities",
            label="number_communities",
            ptype="int",
            default=3,
            min_value=2,
            max_value=50,
            step=1,
            help="Target jumlah komunitas.",
        ),
        ParamSpec(
            key="height",
            label="height",
            ptype="int",
            default=2,
            min_value=1,
            max_value=10,
            step=1,
            help="Parameter 'height'.",
        ),
    ],
    run=run,
)
