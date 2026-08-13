from cdlib import algorithms
from .base import AlgoSpec, ParamSpec

def run(G, params):
    weight_name = str(params["weightName"])

    # wCommunity butuh edge weights; kalau belum ada, set default 1.0
    for u, v, data in G.edges(data=True):
        if weight_name not in data:
            data[weight_name] = 1.0

    return algorithms.wCommunity(
        G,
        min_bel_degree=float(params["min_bel_degree"]),
        threshold_bel_degree=float(params["threshold_bel_degree"]),
        weightName=weight_name,
    )

SPEC = AlgoSpec(
    key="WCOMMUNITY",
    name="WCommunity",
    description="wCommunity: local algorithm untuk overlapping communities pada weighted graph.",
    params=[
        ParamSpec(
            "min_bel_degree", "min_bel_degree", "float", 0.70,
            min_value=0.0, max_value=1.0, step=0.05,
            help="Tolerance (belonging degree) untuk memasukkan node ke community."
        ),
        ParamSpec(
            "threshold_bel_degree", "threshold_bel_degree", "float", 0.70,
            min_value=0.0, max_value=1.0, step=0.05,
            help="Tolerance (belonging degree) untuk memasukkan node ke 'NLU' community."
        ),
        ParamSpec(
            "weightName", "weightName", "select", "weight",
            options=["weight"],
            help="Nama atribut bobot edge. Default: weight."
        ),
    ],
    run=run,
)
