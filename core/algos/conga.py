from cdlib import algorithms
from .base import AlgoSpec, ParamSpec


def run(G, params):
    return algorithms.conga(
        G,
        number_communities=int(params.get("number_communities", 3)),
    )


SPEC = AlgoSpec(
    key="CONGA",
    name="CONGA",
    description="CONGA: deteksi overlap dengan split node/edge betweenness (tanpa height pada versi CDlib Anda).",
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
    ],
    run=run,
)
