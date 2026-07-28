from cdlib import algorithms
from .base import AlgoSpec, ParamSpec

def run(G, params):
    return algorithms.ego_networks(G, level=int(params["level"]))

SPEC = AlgoSpec(
    key="EGONET",
    name="EgoNetworks",
    description="Komunitas dari ego-network; biasanya overlap sangat besar.",
    params=[
        ParamSpec("level", "level", "int", 1, min_value=1, max_value=3, step=1),
    ],
    run=run,
)
