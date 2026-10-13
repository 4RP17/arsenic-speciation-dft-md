"""
arsenic_langevin_md.py

Purpose
-------
Minimal, runnable molecular-dynamics-style simulation for arsenic species in water
under a uniform magnetic field. This is NOT a full atomistic water MD simulation.
It is an implicit-water Langevin/Lorentz model for the center-of-mass motion of
environmentally relevant arsenic species.

Scientific use in thesis
------------------------
Use this as a preliminary computational support model:
- It tests whether charged arsenate species in water would experience a measurable
  Lorentz-force perturbation at room temperature and Tesla-scale fields.
- It calculates diffusion/MSD and a dimensionless Lorentz perturbation number.
- It supports the conclusion that direct magnetic deflection of dissolved ions is
  likely very weak, so future sensing requires preconcentration/immobilization
  and blank-corrected magnetometric readout.

Requirements
------------
pip install numpy pandas matplotlib

Run
---
python arsenic_langevin_md.py

Outputs
-------
md_summary.csv
md_lorentz_number.csv
md_msd_curves.csv
md_msd_curves.png
md_lorentz_number.png
md_magnetic_perturbation_index.png
"""

from dataclasses import dataclass
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# -----------------------------
# Physical constants (SI units)
# -----------------------------
KB = 1.380649e-23        # Boltzmann constant, J/K
E_CHARGE = 1.602176634e-19  # elementary charge, C
NA = 6.02214076e23       # Avogadro constant, mol^-1
AMU = 1.66053906660e-27  # kg

OUTPUT_DIR = Path("md_outputs")
OUTPUT_DIR.mkdir(exist_ok=True)

@dataclass
class Species:
    name: str
    formula: str
    charge_e: int          # net molecular charge in units of e
    molar_mass_g_mol: float
    diffusion_m2_s: float  # approximate aqueous diffusion coefficient

    @property
    def mass_kg(self) -> float:
        return (self.molar_mass_g_mol / 1000.0) / NA

    @property
    def charge_C(self) -> float:
        return self.charge_e * E_CHARGE

def simulate_langevin_lorentz(
    species: Species,
    B_T: float,
    temp_K: float = 298.15,
    n_particles: int = 200,
    n_steps: int = 50000,
    dt_s: float = 1.0e-15,
    flow_velocity_m_s=(1.0e-3, 0.0, 0.0),
    seed: int = 7,
):
    """
    Underdamped Langevin dynamics with Lorentz force:
        m dv/dt = q(v x B) - zeta(v - u_flow) + random_force
        dr/dt = v

    zeta is estimated from Einstein relation:
        D = kBT / zeta  =>  zeta = kBT / D

    Magnetic field is applied along +z:
        B = (0, 0, B_T)

    Notes:
    - This is an implicit-solvent center-of-mass model, not atomistic water MD.
    - It is useful for estimating whether direct Lorentz perturbation is large or small.
    """

    rng = np.random.default_rng(seed)

    m = species.mass_kg
    q = species.charge_C
    D = species.diffusion_m2_s
    zeta = KB * temp_K / D  # kg/s

    B = np.array([0.0, 0.0, B_T], dtype=float)
    u = np.array(flow_velocity_m_s, dtype=float)

    # Thermal initial velocities from Maxwell-Boltzmann distribution.
    v_sigma = np.sqrt(KB * temp_K / m)
    vel = rng.normal(0.0, v_sigma, size=(n_particles, 3))

    pos = np.zeros((n_particles, 3), dtype=float)
    pos0 = pos.copy()

    # Random force standard deviation for underdamped Langevin.
    # <F_random(t)F_random(t')> = 2 zeta kBT delta(t-t')
    noise_std = np.sqrt(2.0 * zeta * KB * temp_K / dt_s)

    sample_every = max(1, n_steps // 500)
    times = []
    msd = []
    mean_y = []
    mean_speed = []

    for step in range(n_steps + 1):
        if step % sample_every == 0:
            disp = pos - pos0
            times.append(step * dt_s)
            msd.append(np.mean(np.sum(disp * disp, axis=1)))
            mean_y.append(float(np.mean(pos[:, 1])))
            mean_speed.append(float(np.mean(np.linalg.norm(vel, axis=1))))

        if step == n_steps:
            break

        lorentz = q * np.cross(vel, B)                         # q(v x B)
        drag = -zeta * (vel - u)                               # viscous relaxation toward flow
        random_force = rng.normal(0.0, noise_std, size=(n_particles, 3))

        force = lorentz + drag + random_force
        vel += (force / m) * dt_s
        pos += vel * dt_s

    df = pd.DataFrame({
        "time_s": np.array(times),
        "time_ps": np.array(times) * 1e12,
        "MSD_m2": np.array(msd),
        "mean_y_m": np.array(mean_y),
        "mean_speed_m_s": np.array(mean_speed),
        "species": species.name,
        "formula": species.formula,
        "B_T": B_T,
    })

    # Estimate diffusion from final 50% of MSD curve: MSD = 6Dt.
    fit_start = len(df) // 2
    slope, intercept = np.polyfit(df["time_s"].iloc[fit_start:], df["MSD_m2"].iloc[fit_start:], 1)
    D_est = slope / 6.0

    # Dimensionless Lorentz perturbation number: qB/zeta = qBD/kBT.
    # In viscous water, this is expected to be << 1 at laboratory fields.
    lorentz_number = abs(q) * abs(B_T) / zeta

    summary = {
        "species": species.name,
        "formula": species.formula,
        "charge_e": species.charge_e,
        "molar_mass_g_mol": species.molar_mass_g_mol,
        "input_D_m2_s": species.diffusion_m2_s,
        "estimated_D_from_MSD_m2_s": D_est,
        "B_T": B_T,
        "temperature_K": temp_K,
        "zeta_kg_s": zeta,
        "lorentz_number_abs_qB_over_zeta": lorentz_number,
        "final_MSD_m2": float(df["MSD_m2"].iloc[-1]),
        "final_mean_y_m": float(df["mean_y_m"].iloc[-1]),
        "mean_speed_final_m_s": float(df["mean_speed_m_s"].iloc[-1]),
        "n_particles": n_particles,
        "n_steps": n_steps,
        "dt_s": dt_s,
    }

    return df, summary

def main():
    # Approximate room-temperature aqueous diffusion coefficients.
    # These are order-of-magnitude values suitable for a preliminary implicit-water model.
    species_list = [
        Species("As(III), neutral arsenous acid", "H3AsO3", 0, 125.943, 5.0e-10),
        Species("As(V), dihydrogen arsenate", "H2AsO4-", -1, 140.927, 6.0e-10),
        Species("As(V), hydrogen arsenate", "HAsO4^2-", -2, 139.919, 6.0e-10),
    ]

    fields_T = [0.0, 1.5, 5.0, 10.0]

    all_curves = []
    summaries = []

    for sp in species_list:
        for B in fields_T:
            curve, summary = simulate_langevin_lorentz(
                species=sp,
                B_T=B,
                temp_K=298.15,
                n_particles=200,
                n_steps=50000,
                dt_s=1.0e-15,
                flow_velocity_m_s=(1.0e-3, 0.0, 0.0),
                seed=1000 + int(abs(sp.charge_e) * 100),  # same seed across B for cleaner comparison
            )
            all_curves.append(curve)
            summaries.append(summary)

    curves = pd.concat(all_curves, ignore_index=True)
    summary_df = pd.DataFrame(summaries)

    # Calculate magnetic perturbation index relative to B=0 for each species.
    rows = []
    for sp_name, group in summary_df.groupby("species"):
        D0 = group.loc[group["B_T"] == 0.0, "estimated_D_from_MSD_m2_s"].iloc[0]
        for _, row in group.iterrows():
            mpi = abs(row["estimated_D_from_MSD_m2_s"] - D0) / abs(D0) if D0 != 0 else np.nan
            rows.append({
                "species": sp_name,
                "formula": row["formula"],
                "B_T": row["B_T"],
                "magnetic_perturbation_index_D": mpi,
                "lorentz_number_abs_qB_over_zeta": row["lorentz_number_abs_qB_over_zeta"],
                "estimated_D_m2_s": row["estimated_D_from_MSD_m2_s"],
            })
    perturb_df = pd.DataFrame(rows)

    curves.to_csv(OUTPUT_DIR / "md_msd_curves.csv", index=False)
    summary_df.to_csv(OUTPUT_DIR / "md_summary.csv", index=False)
    perturb_df.to_csv(OUTPUT_DIR / "md_lorentz_number.csv", index=False)

    # Plot 1: MSD curves
    fig, ax = plt.subplots(figsize=(10, 6))
    for (sp, B), g in curves.groupby(["formula", "B_T"]):
        label = f"{sp}, {B:g} T"
        ax.plot(g["time_ps"], g["MSD_m2"] * 1e18, linewidth=1.3, label=label)
    ax.set_xlabel("Time (ps)")
    ax.set_ylabel("MSD (nm$^2$)")
    ax.set_title("Implicit-water Langevin MD: MSD of arsenic species under magnetic field")
    ax.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "md_msd_curves.png", dpi=600)
    plt.close(fig)

    # Plot 2: Lorentz perturbation number
    fig, ax = plt.subplots(figsize=(8, 5))
    for formula, g in perturb_df.groupby("formula"):
        ax.plot(g["B_T"], g["lorentz_number_abs_qB_over_zeta"], marker="o", linewidth=1.8, label=formula)
    ax.set_xlabel("Magnetic field B (T)")
    ax.set_ylabel(r"Lorentz number, $|q|B/\zeta$")
    ax.set_title("Dimensionless Lorentz perturbation in water")
    ax.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "md_lorentz_number.png", dpi=600)
    plt.close(fig)

    # Plot 3: Perturbation index from simulated diffusion
    fig, ax = plt.subplots(figsize=(8, 5))
    for formula, g in perturb_df.groupby("formula"):
        ax.plot(g["B_T"], g["magnetic_perturbation_index_D"], marker="o", linewidth=1.8, label=formula)
    ax.set_xlabel("Magnetic field B (T)")
    ax.set_ylabel(r"$|D(B)-D(0)|/D(0)$")
    ax.set_title("Simulated magnetic perturbation index from diffusion")
    ax.legend()
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "md_magnetic_perturbation_index.png", dpi=600)
    plt.close(fig)

    print("\nCompleted successfully.")
    print(f"Outputs saved in: {OUTPUT_DIR.resolve()}")
    print("\nKey interpretation:")
    print("1. H3AsO3 is neutral, so the direct Lorentz term is zero in this model.")
    print("2. H2AsO4- and HAsO4^2- experience Lorentz force, but the dimensionless")
    print("   Lorentz number is expected to be very small in water at Tesla-scale fields.")
    print("3. Therefore, direct magnetic deflection of dissolved arsenic ions is not")
    print("   a practical standalone sensing route. Preconcentration or immobilization")
    print("   is required for a realistic magnetometric approach.\n")
    print(summary_df[[
        "formula", "charge_e", "B_T", "estimated_D_from_MSD_m2_s",
        "lorentz_number_abs_qB_over_zeta"
    ]].to_string(index=False))

if __name__ == "__main__":
    main()
