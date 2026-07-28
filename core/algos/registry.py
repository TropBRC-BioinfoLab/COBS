"""Algoritma registry untuk Streamlit app.

File ini adalah versi update dari registry.py mas Heru, ditambah:
- DGM_CRISP (disassembly greedy modularity, crisp)
- GREEDY_MODULARITY (baseline greedy modularity, crisp)

Silakan merge manual jika registry.py lokal mas Heru sudah punya tambahan lain.
"""

from .demon import SPEC as DEMON
from .lfm import SPEC as LFM
from .kclique import SPEC as KCLIQUE
from .slpa import SPEC as SLPA
from .ego_networks import SPEC as EGONET
from .node_perception import SPEC as NODEP

from .angel import SPEC as ANGEL
from .conga import SPEC as CONGA
from .congo import SPEC as CONGO
from .lais2 import SPEC as LAIS2
from .lpanni import SPEC as LPANNI

from .core_expansion import SPEC as CORE_EXP
from .dcs import SPEC as DCS
from .dpclus import SPEC as DPCLUS
from .ipca import SPEC as IPCA
from .ebgc import SPEC as EBGC
from .graph_entropy import SPEC as GRAPH_ENT

from .coach import SPEC as COACH
from .umstmo import SPEC as UMSTMO
from .walkscan import SPEC as WALKSCAN
from .wcommunity import SPEC as WCOMMUNITY
from .multicom import SPEC as MULTICOM
from .osse import SPEC as OSSE

from .nocd_gnn import SPEC as NOCD_GNN
from .cn_ocd import SPEC as CN_OCD
from .cobs_overlap import SPEC as COBS

# --- crisp additions ---
from .dgm_crisp import SPEC as DGM_CRISP
from .greedy_modularity_crisp import SPEC as GREEDY_MODULARITY


SPECS = [
    DEMON,
    LFM,
    KCLIQUE,
    SLPA,
    EGONET,
    NODEP,
    ANGEL,
    CONGA,
    CONGO,
    LAIS2,
    LPANNI,
    CORE_EXP,
    DCS,
    DPCLUS,
    IPCA,
    EBGC,
    GRAPH_ENT,
    COACH,
    UMSTMO,
    WALKSCAN,
    WCOMMUNITY,
    MULTICOM,
    OSSE,
    CN_OCD,
    COBS,
    NOCD_GNN,

    # crisp
    GREEDY_MODULARITY,
    DGM_CRISP,
]


def list_specs():
    return SPECS


def get_spec(key: str):
    for s in SPECS:
        if s.key == key:
            return s
    raise ValueError(f"Unknown algo key: {key}")
