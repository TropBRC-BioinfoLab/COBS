from cdlib import algorithms
from .base import AlgoSpec, ParamSpec

def run(G, params):
    return algorithms.walkscan(
        G,
        nb_steps=int(params["nb_steps"]),
        eps=float(params["eps"]),
        min_samples=int(params["min_samples"]),
        # init_vector dibiarkan None (default CDlib)
    )

SPEC = AlgoSpec(
    key="WALKSCAN",
    name="WalkSCAN",
    description="WalkSCAN: random-walk + PageRank scoring, lalu clustering (DBSCAN) pada embedding.",
    params=[
        ParamSpec(
            "nb_steps", "nb_steps", "int", 2,
            min_value=1, max_value=10, step=1,
            help="Panjang random walk."
        ),
        ParamSpec(
            "eps", "eps", "float", 0.10,
            min_value=0.01, max_value=1.0, step=0.01,
            help="DBSCAN eps."
        ),
        ParamSpec(
            "min_samples", "min_samples", "int", 3,
            min_value=1, max_value=30, step=1,
            help="DBSCAN min_samples."
        ),
    ],
    run=run,
)
