import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.data import Data
from torch_geometric.nn import GCNConv, GAE
from torch_geometric.utils import add_self_loops, negative_sampling


# ── Hilfsfunktion ─────────────────────────────────────────────────────────────

def standardize_np(X: np.ndarray, eps: float = 1e-8) -> np.ndarray:
    """Z-score standardization per column."""
    mu = X.mean(axis=0, keepdims=True)
    sd = X.std(axis=0, keepdims=True)
    return (X - mu) / (sd + eps)


# ── Modell ────────────────────────────────────────────────────────────────────

class Encoder(nn.Module):
    """2-layer GCN encoder: node features + graph → embeddings Z."""
    def __init__(self, in_channels: int, hidden: int = 64, out_channels: int = 8, dropout: float = 0.1):
        super().__init__()
        self.conv1 = GCNConv(in_channels, hidden)
        self.conv2 = GCNConv(hidden, out_channels)
        self.bn1 = nn.BatchNorm1d(hidden)
        self.drop = nn.Dropout(p=dropout)

    def forward(self, x, edge_index, edge_weight=None):
        x = F.relu(self.bn1(self.conv1(x, edge_index, edge_weight=edge_weight)))
        x = self.drop(x)
        return self.conv2(x, edge_index, edge_weight=edge_weight)


class FeatureDecoder(nn.Module):
    """Reconstructs input features from embeddings Z (prevents Z from ignoring feature signal)."""
    def __init__(self, in_channels: int, out_channels: int, hidden: int = 64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_channels, hidden),
            nn.ReLU(),
            nn.Linear(hidden, out_channels),
        )

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.net(z)


# ── Training ──────────────────────────────────────────────────────────────────

def train_gae_embeddings(
    X: np.ndarray,
    edge_index: np.ndarray,
    edge_weight: np.ndarray = None,
    out_dim: int = 8,
    hidden: int = 64,
    epochs: int = 1000,
    lr: float = 5e-4,
    add_loops: bool = True,
    seed: int = 0,
    alpha: float = 0.5,  # 0 = only feature loss | 1 = only graph loss
    dropout: float = 0.1,
    kmer_weight: float = 1.0,
    n_tiara_features: int = 0,
    n_length_features: int = 0,
) -> np.ndarray:
    """Trains GAE, returns node embeddings Z of shape (n_nodes, out_dim).

    Feature layout in X: [coverage(1), length(0|1), kmer(n_kmer), tiara(0|4)]
      - coverage: standardized, no downweighting
      - length: standardized, no downweighting (continuous)
      - kmer: standardized, then multiplied by kmer_weight
      - tiara: binary (0/1), not standardized, no downweighting
    """
    torch.manual_seed(seed)
    np.random.seed(seed)

    # ── Standardisierung + k-mer Downweighting ───────────────────────────────
    Xs = standardize_np(X).astype(np.float32)
    if n_tiara_features > 0:
        # Tiara is binary — restore raw values after standardization
        Xs[:, -n_tiara_features:] = X[:, -n_tiara_features:].astype(np.float32)
    n_kmer = Xs.shape[1] - 1 - n_length_features - n_tiara_features
    kmer_start = 1 + n_length_features
    if kmer_weight != 1.0 and n_kmer > 0:
        Xs[:, kmer_start:kmer_start + n_kmer] *= kmer_weight
        print(f"[GNN] k-mer downweighting: x{kmer_weight} on {n_kmer} columns")

    x = torch.tensor(Xs, dtype=torch.float32)
    ei = torch.tensor(edge_index, dtype=torch.long)
    ew = torch.tensor(edge_weight, dtype=torch.float32) if edge_weight is not None else None

    if add_loops:
        ei, ew = add_self_loops(ei, edge_attr=ew, num_nodes=x.size(0), fill_value=1.0)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    data = Data(x=x, edge_index=ei, edge_weight=ew).to(device)
    print(f"[GNN] Device: {device} | Features: {x.size(1)} | Hidden: {hidden} | Out: {out_dim}")

    gae_model = GAE(Encoder(x.size(1), hidden=hidden, out_channels=out_dim, dropout=dropout)).to(device)
    feat_decoder = FeatureDecoder(out_dim, x.size(1), hidden=hidden).to(device)

    opt = torch.optim.Adam(
        list(gae_model.parameters()) + list(feat_decoder.parameters()),
        lr=lr, weight_decay=1e-4,
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)

    print(f"[GNN] Training {epochs} epochs | alpha={alpha}")
    for epoch in range(1, epochs + 1):
        gae_model.train()
        feat_decoder.train()
        opt.zero_grad()

        z = gae_model.encode(data.x, data.edge_index, data.edge_weight)

        loss_link = torch.tensor(0.0, device=device)
        if alpha > 0.0 and data.edge_index.size(1) > 0:
            neg_ei = negative_sampling(data.edge_index, num_nodes=data.num_nodes,
                                       num_neg_samples=data.edge_index.size(1))
            loss_link = gae_model.recon_loss(z, data.edge_index, neg_ei)

        loss_feat = torch.tensor(0.0, device=device)
        if alpha < 1.0:
            loss_feat = F.mse_loss(feat_decoder(z), data.x)

        loss = alpha * loss_link + (1.0 - alpha) * loss_feat
        loss.backward()
        torch.nn.utils.clip_grad_norm_(
            list(gae_model.parameters()) + list(feat_decoder.parameters()), max_norm=1.0
        )
        opt.step()
        scheduler.step()

        if epoch % 100 == 0 or epoch == 1:
            print(f"[GNN] Epoch {epoch:4d}/{epochs} | "
                  f"loss={loss.item():.4f}  (link={loss_link.item():.4f}, feat={loss_feat.item():.4f})")

    # ── Embeddings extrahieren ────────────────────────────────────────────────
    gae_model.eval()
    with torch.no_grad():
        Z = gae_model.encode(data.x, data.edge_index, data.edge_weight).cpu().numpy()

    if np.isnan(Z).any():
        print(f"[GNN] WARNUNG: {int(np.isnan(Z).sum())} NaN in Z — replaced with 0")
        Z = np.nan_to_num(Z, nan=0.0)

    return Z
