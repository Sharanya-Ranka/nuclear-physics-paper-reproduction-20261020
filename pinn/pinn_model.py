import io
import math
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, Dataset
import pytorch_lightning as pl


# -----------------------------------------------------------------------------
# 1. Dataset & Collocation Generator
# -----------------------------------------------------------------------------
class SingleNuclideDataset(Dataset):
    """Holds empirical data points and samples unsupervised collocation points."""

    def __init__(
        self,
        df: pd.DataFrame,
        e_min: float = 0.0,
        e_max: float = 30.0,
        n_collocation: int = 500,
    ):
        # Extract energy (MeV) and cross-section (barns -> mb)
        energy_col = "('Energy', 'MeV')"
        data_col = "('Data', 'barns')"

        self.energies = torch.tensor(df[energy_col].values, dtype=torch.float32).view(
            -1, 1
        )
        self.cross_sections = torch.tensor(
            df[data_col].values * 1000.0, dtype=torch.float32
        ).view(
            -1, 1
        )  # barns to mb

        # Extract fixed nuclide constants
        self.A = int(df["TargetA"].iloc[0])
        self.Z = int(df["TargetZ"].iloc[0])
        self.N = self.A - self.Z
        self.Sn = float(df["Sn"].iloc[0]) / 1000.0  # keV -> MeV

        # SLO parameters (mb and MeV)
        self.E1 = float(df["E1"].iloc[0])
        self.CS1 = float(df["CS1"].iloc[0])
        self.W1 = float(df["W1"].iloc[0])

        # Collocation bounds for physics residual evaluation
        self.e_min = max(0.0, e_min)
        self.e_max = e_max
        self.n_collocation = n_collocation

    def __len__(self):
        return len(self.energies)

    def __getitem__(self, idx):
        return self.energies[idx], self.cross_sections[idx]

    def sample_collocation_points(self) -> torch.Tensor:
        """Uniformly samples points along [0, 30] MeV with gradient tracking."""
        pts = torch.linspace(
            self.e_min, self.e_max, self.n_collocation, dtype=torch.float32
        ).view(-1, 1)
        return pts


# -----------------------------------------------------------------------------
# 2. PINN Lightning Module
# -----------------------------------------------------------------------------
class NuclearPINN(pl.LightningModule):

    def __init__(
        self,
        A: int,
        Z: int,
        Sn: float,
        E1: float,
        CS1: float,
        W1: float,
        e_min: float = 0.0,
        e_max: float = 30.0,
        hidden_dim: int = 64,
        num_layers: int = 4,
        lr: float = 1e-3,
        lambda_slo: float = 0.05,
        lambda_smooth: float = 1e-4,
        lambda_trk: float = 1e-4,
    ):
        super().__init__()
        self.save_hyperparameters()

        self.A = A
        self.Z = Z
        self.N = A - Z
        self.Sn = Sn
        self.E1 = E1
        self.CS1 = CS1
        self.W1 = W1
        self.e_min = e_min
        self.e_max = e_max
        self.lr = lr

        # Scaling: Max expected height is bounded by roughly ~1.5x of peak Lorentzian CS1
        self.sigma_scale = max(self.CS1 * 1.5, 100.0)

        # Physics loss weights
        self.lambda_slo = lambda_slo
        self.lambda_smooth = lambda_smooth
        self.lambda_trk = lambda_trk

        # TRK classical dipole sum rule limit in MeV*mb: 60 * (N * Z / A)
        self.trk_limit = 60.0 * (self.N * self.Z) / float(self.A)

        # MLP architecture with residual connections
        layers = []
        layers.append(nn.Linear(1, hidden_dim))
        layers.append(nn.Tanh())
        for _ in range(num_layers - 1):
            layers.append(nn.Linear(hidden_dim, hidden_dim))
            layers.append(nn.Tanh())
        layers.append(nn.Linear(hidden_dim, 1))

        self.net = nn.Sequential(*layers)
        self.gate_tau = 0.1  # Transition sharpness at threshold Sn

    def standard_lorentzian(self, E: torch.Tensor) -> torch.Tensor:
        """Computes reference SLO cross section in mb."""
        E2 = E**2
        E1_2 = self.E1**2
        W1_2 = self.W1**2
        denom = (E2 - E1_2) ** 2 + E2 * W1_2
        num = self.CS1 * E2 * W1_2
        return torch.where(denom > 0, num / denom, torch.zeros_like(E))

    def forward(self, E: torch.Tensor) -> torch.Tensor:
        """
        Forward pass with normalized coordinate input and hard physical gating:
        1. E is normalized to [-1, 1] for neural network stability.
        2. Network output is constrained positive via Softplus.
        3. Hard gate Sigmoid((E - Sn) / tau) drives predictions to 0 below threshold.
        """
        # Normalize input domain [e_min, e_max] -> [-1, 1]
        E_norm = 2.0 * (E - self.e_min) / (self.e_max - self.e_min) - 1.0

        raw_out = self.net(E_norm)
        pos_out = nn.functional.softplus(raw_out) * self.sigma_scale

        # Hard threshold gating: zero below Sn
        gate = torch.sigmoid((E - self.Sn) / self.gate_tau)
        return pos_out * gate

    def training_step(self, batch, batch_idx):
        e_exp, sigma_exp = batch

        # 1. Experimental Data Loss (Normalized Huber Loss)
        sigma_pred_exp = self.forward(e_exp)
        loss_data = nn.functional.smooth_l1_loss(
            sigma_pred_exp / self.sigma_scale, sigma_exp / self.sigma_scale
        )

        # 2. Unsupervised Collocation Grid
        e_coll = torch.linspace(self.e_min, self.e_max, 300, device=self.device).view(
            -1, 1
        )
        e_coll.requires_grad_(True)

        sigma_coll = self.forward(e_coll)

        # 3. Standard Lorentzian Prior Loss (weak anchor in the GDR region)
        with torch.no_grad():
            sigma_slo = self.standard_lorentzian(e_coll)
        # Apply prior primarily above Sn
        above_sn_mask = (e_coll >= self.Sn).float()
        loss_slo = torch.mean(
            above_sn_mask * ((sigma_coll - sigma_slo) / self.sigma_scale) ** 2
        )

        # 4. Smoothness / Curvature Penalty (via autograd)
        d_sigma_dE = torch.autograd.grad(
            outputs=sigma_coll,
            inputs=e_coll,
            grad_outputs=torch.ones_like(sigma_coll),
            create_graph=True,
            retain_graph=True,
        )[0]

        d2_sigma_dE2 = torch.autograd.grad(
            outputs=d_sigma_dE,
            inputs=e_coll,
            grad_outputs=torch.ones_like(d_sigma_dE),
            create_graph=True,
            retain_graph=True,
        )[0]

        loss_smooth = torch.mean((d2_sigma_dE2 / self.sigma_scale) ** 2)

        # 5. TRK Sum Rule Penalty (Upper Bound Constraint)
        # Numerical integration using trapezoidal rule
        delta_e = (self.e_max - self.e_min) / (len(e_coll) - 1)
        integrated_area = torch.sum(sigma_coll) * delta_e
        loss_trk = nn.functional.relu(integrated_area - self.trk_limit) / (
            self.trk_limit + 1e-6
        )

        # Total Loss
        total_loss = (
            loss_data
            + self.lambda_slo * loss_slo
            + self.lambda_smooth * loss_smooth
            + self.lambda_trk * loss_trk
        )

        self.log_dict(
            {
                "loss/total": total_loss,
                "loss/data": loss_data,
                "loss/slo": loss_slo,
                "loss/smooth": loss_smooth,
                "loss/trk": loss_trk,
            },
            prog_bar=True,
        )

        return total_loss

    def configure_optimizers(self):
        optimizer = torch.optim.AdamW(self.parameters(), lr=self.lr, weight_decay=1e-5)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
            optimizer, T_max=500, eta_min=1e-5
        )
        return {
            "optimizer": optimizer,
            "lr_scheduler": {"scheduler": scheduler, "monitor": "loss/total"},
        }

    # -------------------------------------------------------------------------
    # 3. Clean Predict Function
    # -------------------------------------------------------------------------
    @torch.no_grad()
    def predict_cross_section(
        self, energies: np.ndarray, return_units: str = "mb"
    ) -> np.ndarray:
        """
        Inference function.
        :param energies: 1D numpy array of photon energies in MeV.
        :param return_units: 'mb' (millibarns) or 'barns'.
        :return: Predicted cross sections as numpy array.
        """
        self.eval()
        e_tensor = torch.tensor(energies, dtype=torch.float32, device=self.device).view(
            -1, 1
        )
        pred_mb = self.forward(e_tensor).cpu().numpy().flatten()

        # Hard zero out anything strictly <= Sn
        pred_mb[energies <= self.Sn] = 0.0

        if return_units.lower() == "barns":
            return pred_mb / 1000.0
        return pred_mb


# -----------------------------------------------------------------------------
# 4. Training Pipeline Demonstration
# -----------------------------------------------------------------------------
def run_pinn_training():
    # Load raw data snippet
    csv_data = """TargetA,TargetZ,"('Energy', 'MeV')","('Data', 'barns')",Element,Sn,E1,CS1,W1
237,93,6.73,0.0327,Np,6577.81306,11.01,253.8,2.91
237,93,7.23,0.0164,Np,6577.81306,11.01,253.8,2.91
237,93,7.38,0.0278,Np,6577.81306,11.01,253.8,2.91
"""
    df = pd.read_csv(io.StringIO(csv_data))

    # Initialize Dataset and Model
    dataset = SingleNuclideDataset(df, e_min=0.0, e_max=30.0)
    train_loader = DataLoader(
        dataset, batch_size=len(dataset), shuffle=False
    )  # Small dataset: full-batch

    model = NuclearPINN(
        A=dataset.A,
        Z=dataset.Z,
        Sn=dataset.Sn,
        E1=dataset.E1,
        CS1=dataset.CS1,
        W1=dataset.W1,
        e_min=0.0,
        e_max=30.0,
        hidden_dim=32,
        num_layers=3,
        lr=2e-3,
    )

    # Train with PyTorch Lightning
    trainer = pl.Trainer(
        max_epochs=300,
        accelerator="auto",
        devices=1,
        enable_checkpointing=False,
        logger=False,
    )
    trainer.fit(model, train_loader)

    # Evaluate on full uniform [0, 30] MeV grid at 0.1 MeV steps
    test_grid = np.arange(0.0, 30.1, 0.1)
    predictions_mb = model.predict_cross_section(test_grid, return_units="mb")
    predictions_barns = model.predict_cross_section(test_grid, return_units="barns")

    print("\n--- Predictions for 237Np(gamma, n) ---")
    print(f"Reaction Threshold (Sn): {dataset.Sn:.3f} MeV")
    print(
        f"E = 5.0 MeV (Below Sn) -> Cross Section: {predictions_mb[test_grid == 5.0][0]:.4f} mb"
    )
    print(
        f"E = 6.6 MeV (At Threshold) -> Cross Section: {predictions_mb[np.isclose(test_grid, 6.6)][0]:.4f} mb"
    )
    print(
        f"E = 11.0 MeV (Resonance Peak) -> Cross Section: {predictions_mb[np.isclose(test_grid, 11.0)][0]:.2f} mb ({predictions_barns[np.isclose(test_grid, 11.0)][0]:.4f} barns)"
    )


if __name__ == "__main__":
    run_pinn_training()
