from cdlib import algorithms
from .base import AlgoSpec, ParamSpec

def run(G, params):
    seed_node = int(params["seed_node"])

    # jika seed_node tidak ada (misal nanti graph bukan 0..N), fallback ke node pertama
    if seed_node not in G:
        seed_node = list(G.nodes())[0]

    return algorithms.multicom(G, seed_node=seed_node)

SPEC = AlgoSpec(
    key="MULTICOM",
    name="MULTICOM",
    description="MULTICOM: multiple local communities (bisa overlap) dengan ekspansi dari seed.",
    params=[
        ParamSpec(
            "seed_node", "seed_node", "int", 0,
            min_value=0, max_value=33, step=1,
            help="Seed node (untuk Karate: 0..33)."
        ),
    ],
    run=run,
)
