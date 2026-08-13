from cdlib import algorithms
from .base import AlgoSpec, ParamSpec

def run(G, params):
    seed_node = int(params["seed_node"])
    if seed_node not in G:
        seed_node = list(G.nodes())[0]

    return algorithms.overlapping_seed_set_expansion(
        G,
        seeds=[seed_node],
        ninf=bool(params["ninf"]),
        expansion=str(params["expansion"]),
        stopping=str(params["stopping"]),
        nworkers=int(params["nworkers"]),
        nruns=int(params["nruns"]),
        alpha=float(params["alpha"]),
        delta=float(params["delta"]),
        # maxexpand tidak saya expose (default INF)
    )

SPEC = AlgoSpec(
    key="OSSE",
    name="OSSE (Seed Expansion)",
    description="OSSE: seed set expansion (Personalized PageRank) untuk overlapping communities.",
    params=[
        ParamSpec(
            "seed_node", "seed_node", "int", 0,
            min_value=0, max_value=33, step=1,
            help="Satu seed node (untuk Karate: 0..33)."
        ),
        ParamSpec("ninf", "ninf", "bool", False, help="Neighbourhood inflation."),
        ParamSpec(
            "expansion", "expansion", "select", "ppr",
            options=["ppr", "vppr"],
            help="Metode ekspansi seed."
        ),
        ParamSpec(
            "stopping", "stopping", "select", "cond",
            options=["cond"],
            help="Stopping criteria."
        ),
        ParamSpec(
            "nworkers", "nworkers", "int", 1,
            min_value=1, max_value=8, step=1,
            help="Jumlah worker."
        ),
        ParamSpec(
            "nruns", "nruns", "int", 13,
            min_value=1, max_value=50, step=1,
            help="Jumlah run."
        ),
        ParamSpec(
            "alpha", "alpha", "float", 0.99,
            min_value=0.50, max_value=0.999, step=0.001,
            help="Alpha untuk Personalized PageRank."
        ),
        ParamSpec(
            "delta", "delta", "float", 0.20,
            min_value=0.0, max_value=1.0, step=0.05,
            help="Minimum distance untuk near-duplicate communities."
        ),
    ],
    run=run,
)
