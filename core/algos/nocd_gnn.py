from __future__ import annotations

from .base import AlgoSpec, ParamSpec
import torch
import torch.nn as nn
import torch.nn.functional as F



def run(G, params):
    """
    NOCD-style overlapping community detection (GNN encoder + Bernoulli-Poisson decoder).
    Menghasilkan komunitas overlap (list of list).
    """
    # Import heavy deps lazily supaya startup Streamlit tetap ringan
    import random
    import numpy as np
    import scipy.sparse as sp
    import networkx as nx


    # cdlib NodeClustering (fallback stub kalau cdlib tidak ada)
    try:
        from cdlib.classes.node_clustering import NodeClustering  # type: ignore
    except Exception:  # pragma: no cover
        NodeClustering = None  # type: ignore

    # ----------------
    # Params
    # ----------------
    k = int(params.get("k", 3))  # DEFAULT K=3 (karate)
    hidden_dim = int(params.get("hidden_dim", 64))
    n_layers = int(params.get("n_layers", 2))
    dropout = float(params.get("dropout", 0.5))
    lr = float(params.get("lr", 0.01))
    epochs = int(params.get("epochs", 200))
    weight_decay = float(params.get("weight_decay", 0.0))
    features = params.get("features", "degree")  # degree | onehot | random
    rand_dim = int(params.get("rand_dim", 16))
    device_pref = params.get("device", "auto")  # auto | cpu | cuda
    seed = int(params.get("seed", 42))
    overlap_ratio = float(params.get("overlap_ratio", 0.8))
    neg_ratio = float(params.get("neg_ratio", 1.0))

    # Reproducible
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)

    # Device
    if device_pref == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    elif device_pref == "cuda" and not torch.cuda.is_available():
        device = "cpu"
    else:
        device = device_pref

    # Untuk implementasi ini: dibuat undirected dulu
    G0 = G.to_undirected() if getattr(G, "is_directed", lambda: False)() else G
    nodes = list(G0.nodes())
    n = len(nodes)
    if n == 0:
        return _make_clustering(NodeClustering, [], G0, "NOCD (GNN)")

    # Adjacency binary
    A = nx.to_scipy_sparse_array(G0, nodelist=nodes, format="csr", dtype=np.float32)
    if A.nnz > 0:
        A.data = np.ones_like(A.data, dtype=np.float32)

    # Normalize adjacency with self-loop: D^{-1/2}(A+I)D^{-1/2}
    I = sp.identity(n, dtype=np.float32, format="csr")
    A_hat = A + I
    deg = np.asarray(A_hat.sum(axis=1)).reshape(-1).astype(np.float32)
    deg_inv_sqrt = np.zeros_like(deg)
    mask = deg > 0
    deg_inv_sqrt[mask] = np.power(deg[mask], -0.5)
    D_inv_sqrt = sp.diags(deg_inv_sqrt)
    A_norm = (D_inv_sqrt @ A_hat @ D_inv_sqrt).tocsr()
    A_norm_t = _sp_to_torch_sparse(A_norm).to(device)

    # ----------------
    # Features X
    # ----------------
    if features == "onehot":
        X_sp = sp.identity(n, dtype=np.float32, format="csr")
        X = _sp_to_torch_sparse(X_sp).to(device)
        in_dim = n
        x_is_sparse = True
    elif features == "random":
        rng = np.random.default_rng(seed)
        X_np = rng.normal(size=(n, rand_dim)).astype(np.float32)
        X = torch.from_numpy(X_np).to(device)
        in_dim = rand_dim
        x_is_sparse = False
    else:  # "degree" default
        deg_feat = np.asarray(A.sum(axis=1)).reshape(-1, 1).astype(np.float32)
        X = torch.from_numpy(deg_feat).to(device)
        in_dim = 1
        x_is_sparse = False

    # ----------------
    # Model
    # ----------------
    model = _GCNEncoder(
        in_dim, hidden_dim, k,
        n_layers=n_layers,
        dropout=dropout,
        x_is_sparse=x_is_sparse
    ).to(device)

    opt = torch.optim.Adam(model.parameters(), lr=lr)

    # Edge list (ambil unik i<j)
    A_csr = A.tocsr()
    rows, cols = A_csr.nonzero()
    edges = np.vstack([rows, cols]).T
    if edges.size > 0:
        edges = edges[edges[:, 0] < edges[:, 1]]
    m = int(edges.shape[0])

    total_pairs = n * (n - 1) // 2
    num_non_edges = total_pairs - m

    # Untuk karate (kecil) aman pakai exact non-edge; untuk besar pakai negative sampling
    exact_non_edges = (n <= 500)

    for _ in range(max(1, epochs)):
        model.train()
        Z = model(X, A_norm_t)  # (n,k), non-negative

        loss = _berpo_loss(
            Z, edges, A_csr,
            exact_non_edges=exact_non_edges,
            neg_ratio=neg_ratio,
            total_pairs=total_pairs,
            num_non_edges=num_non_edges,
        )
        if weight_decay > 0:
            l2 = sum((p ** 2).sum() for p in model.parameters())
            loss = loss + weight_decay * l2

        opt.zero_grad()
        loss.backward()
        opt.step()

    # ----------------
    # Hard memberships -> communities
    # ----------------
    model.eval()
    with torch.no_grad():
        Z = model(X, A_norm_t).detach().cpu().numpy()

    # Membership overlap: selalu include argmax; tambah komunitas lain kalau >= overlap_ratio * max
    maxv = Z.max(axis=1, keepdims=True) + 1e-12
    membership = Z >= (overlap_ratio * maxv)
    primary = Z.argmax(axis=1)
    membership[np.arange(n), primary] = True

    comms = []
    for c in range(k):
        members = [nodes[i] for i in range(n) if membership[i, c]]
        if members:
            comms.append(members)

    return _make_clustering(NodeClustering, comms, G0, "NOCD (GNN)")


# ----------------------------
# Helpers
# ----------------------------
class _GCNEncoder(nn.Module):
    """GCN encoder sederhana (sparse-aware) -> Softplus output (non-negative memberships)."""
    def __init__(self, in_dim: int, hidden_dim: int, out_dim: int,
                 n_layers: int = 2, dropout: float = 0.5, x_is_sparse: bool = False):
        super().__init__()
        self.n_layers = max(1, int(n_layers))
        self.dropout = float(dropout)
        self.x_is_sparse = x_is_sparse

        if self.n_layers == 1:
            self.W_out = nn.Linear(in_dim, out_dim, bias=False)
            self.W_hid = nn.ModuleList()
        else:
            self.W_hid = nn.ModuleList()
            self.W_hid.append(nn.Linear(in_dim, hidden_dim, bias=False))
            for _ in range(self.n_layers - 2):
                self.W_hid.append(nn.Linear(hidden_dim, hidden_dim, bias=False))
            self.W_out = nn.Linear(hidden_dim, out_dim, bias=False)

    def forward(self, X, A_norm):
        H = X
        if self.n_layers == 1:
            H = self._apply_linear(H, self.W_out)
            H = torch.sparse.mm(A_norm, H)
            return F.softplus(H)

        for W in self.W_hid:
            H = self._apply_linear(H, W)
            H = torch.sparse.mm(A_norm, H)
            H = F.relu(H)
            H = F.dropout(H, p=self.dropout, training=self.training)

        H = self._apply_linear(H, self.W_out)
        H = torch.sparse.mm(A_norm, H)
        return F.softplus(H)

    def _apply_linear(self, H, layer: nn.Linear):
        if getattr(H, "is_sparse", False):
            W = layer.weight.t()  # (d_in x d_out)
            return torch.sparse.mm(H, W)
        return layer(H)


def _sp_to_torch_sparse(mat):
    import numpy as np
    import torch
    mat = mat.tocoo()
    idx = np.vstack((mat.row, mat.col)).astype(np.int64)
    indices = torch.from_numpy(idx)
    values = torch.from_numpy(mat.data.astype(np.float32))
    shape = torch.Size(mat.shape)
    return torch.sparse_coo_tensor(indices, values, shape).coalesce()


def _berpo_loss(Z, edges, A_csr, exact_non_edges: bool,
                neg_ratio: float, total_pairs: int, num_non_edges: int):
    """
    Bernoulli-Poisson model:
      P(A_ij=1) = 1 - exp(- Z_i^T Z_j), Z>=0
      log P(A_ij=0) = - Z_i^T Z_j
    """
    import numpy as np
    import torch

    eps = 1e-12
    n = Z.shape[0]
    device = Z.device

    if edges.size == 0:
        # Kalau tak ada edge: cukup minimisasi rate rata-rata
        if n < 2:
            return torch.tensor(0.0, device=device)
        i = torch.randint(0, n, (min(1000, n * 10),), device=device)
        j = torch.randint(0, n, (min(1000, n * 10),), device=device)
        mask = i != j
        i, j = i[mask], j[mask]
        r = (Z[i] * Z[j]).sum(dim=1)
        return r.mean()

    e = torch.from_numpy(edges.astype(np.int64)).to(device)
    u, v = e[:, 0], e[:, 1]
    r_e = (Z[u] * Z[v]).sum(dim=1)

    # -log(1-exp(-r)) stable: log(-expm1(-r))
    edge_nll = -torch.log(-torch.expm1(-r_e) + eps).mean()

    if exact_non_edges and n <= 500:
        # exact O(n^2) untuk kecil (karate aman)
        R = Z @ Z.t()
        tri = torch.triu_indices(n, n, offset=1, device=device)
        r_all = R[tri[0], tri[1]]

        A_dense = torch.from_numpy(A_csr.toarray().astype(np.float32)).to(device)
        edge_mask = (A_dense[tri[0], tri[1]] > 0)
        r_non = r_all[~edge_mask]
        non_nll = r_non.mean() if r_non.numel() > 0 else torch.tensor(0.0, device=device)
        return edge_nll + non_nll

    # Negative sampling untuk non-edges (lebih cocok graph besar)
    n_edges = edges.shape[0]
    n_neg = max(1, int(neg_ratio * n_edges))

    neg_u, neg_v = [], []
    tries = 0
    max_tries = n_neg * 50
    while len(neg_u) < n_neg and tries < max_tries:
        tries += 1
        a = np.random.randint(0, n)
        b = np.random.randint(0, n)
        if a == b:
            continue
        if a > b:
            a, b = b, a
        if A_csr[a, b] != 0:
            continue
        neg_u.append(a)
        neg_v.append(b)

    if len(neg_u) == 0:
        return edge_nll

    neg = torch.tensor(np.vstack([neg_u, neg_v]).T, dtype=torch.long, device=device)
    nu, nv = neg[:, 0], neg[:, 1]
    r_neg = (Z[nu] * Z[nv]).sum(dim=1)

    non_mean = r_neg.mean()
    return edge_nll + non_mean


def _make_clustering(NodeClustering, comms, G, method_name: str):
    if NodeClustering is None:  # pragma: no cover
        class _Stub:
            def __init__(self, communities, graph, method_name):
                self.communities = communities
                self.graph = graph
                self.method_name = method_name
                self.overlap = True
        return _Stub(comms, G, method_name)

    return NodeClustering(comms, G, method_name=method_name, overlap=True)


SPEC = AlgoSpec(
    key="NOCD_GNN",
    name="NOCD (GNN)",
    description="Deteksi komunitas overlap berbasis GNN (GCN) + Bernoulli-Poisson decoder (NOCD-style).",
    params=[
    ParamSpec("k", "Jumlah komunitas (K)", "int", 3, min_value=2, max_value=20, step=1),
    ParamSpec("features", "Fitur node", "select", "degree", options=["degree", "onehot", "random"]),

    ParamSpec("rand_dim", "Dimensi fitur random", "int", 16, min_value=4, max_value=256, step=4),
    ParamSpec("hidden_dim", "Hidden dim", "int", 64, min_value=8, max_value=256, step=8),
    ParamSpec("n_layers", "Jumlah layer GCN", "int", 2, min_value=1, max_value=4, step=1),

    ParamSpec("dropout", "Dropout", "float", 0.5, min_value=0.0, max_value=0.9, step=0.05),
    ParamSpec("lr", "Learning rate", "float", 0.01, min_value=0.0001, max_value=0.1, step=0.001),
    ParamSpec("epochs", "Epoch", "int", 200, min_value=10, max_value=2000, step=10),
    ParamSpec("weight_decay", "L2 regularization", "float", 0.0, min_value=0.0, max_value=0.01, step=0.0005),

    ParamSpec("overlap_ratio", "Rasio overlap (>= ratio * max)", "float", 0.8, min_value=0.0, max_value=1.0, step=0.05),
    ParamSpec("neg_ratio", "Neg sampling ratio", "float", 1.0, min_value=0.1, max_value=10.0, step=0.1),

    ParamSpec("device", "Device", "select", "auto", options=["auto", "cpu", "cuda"]),
    ParamSpec("seed", "Random seed", "int", 42, min_value=0, max_value=9999, step=1),
],

    run=run,
)
