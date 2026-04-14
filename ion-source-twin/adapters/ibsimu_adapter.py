"""
Adapter for the IBSIMU C++ beam extraction binary.
Writes config, runs simulation, parses all CSV/PNG outputs.
"""
import os, subprocess, numpy as np
from config import IBSIMU_BINARY, DATA_DIR, EV_TO_J, AMU_TO_KG, ECHARGE


def write_config(params: dict, electrodes: list, path: str):
    """Write beam_config.txt for beam_sim."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for k in ["grid_points", "particles", "current_density", "beam_energy",
                   "beam_mass", "beam_charge"]:
            if k in params:
                f.write(f"{k}={params[k]}\n")
        f.write(f"solver={params.get('solver','bicgstab')}\n")
        f.write(f"mg_levels={params.get('mg_levels',4)}\n")
        f.write(f"beam_mode={params.get('beam_mode','energy')}\n")
        if params.get("beam_mode") == "twiss":
            for k in ["beam_alpha", "beam_beta", "beam_emittance"]:
                f.write(f"{k}={params.get(k,0)}\n")
        f.write(f"bfield={params.get('bfield','none')}\n")
        if params.get("bfield") == "solenoid":
            for k in ["sol_B0", "sol_z1", "sol_z2"]:
                f.write(f"{k}={params.get(k,0)}\n")
        f.write(f"electrodes={len(electrodes)}\n")
        for e in electrodes:
            f.write(f"{e['dist']} {e['apt']} {e['volt']} {e['thick']} {e['chamfer']}\n")


def run_extraction(params: dict, electrodes: list, outdir: str = None) -> dict:
    """Run beam_sim and return parsed results."""
    outdir = outdir or os.path.join(DATA_DIR, "ibsimu")
    os.makedirs(outdir, exist_ok=True)
    cfg_path = os.path.join(outdir, "beam_config.txt")
    write_config(params, electrodes, cfg_path)
    cmd = [IBSIMU_BINARY, "--config", cfg_path, "--outdir", outdir]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if proc.returncode != 0:
        raise RuntimeError(f"beam_sim failed:\n{proc.stderr}")
    return parse_results(outdir)


def parse_results(outdir: str) -> dict:
    """Parse all CSV outputs from an IBSIMU run."""
    results = {}
    for name in ["emittance_profile", "phase_space", "field_along_axis",
                  "convergence", "electrode_currents", "energy_distribution",
                  "exit_particles"]:
        path = os.path.join(outdir, f"{name}.csv")
        if os.path.exists(path):
            results[name] = _read_csv(path)
    traj = os.path.join(outdir, "trajectory.png")
    if os.path.exists(traj):
        results["trajectory_png"] = traj
    return results


def extract_exit_particles(outdir: str) -> np.ndarray:
    """Load exit-plane particle data: (N, 5) array of [y_m, vx, vy, vz, ek_eV]."""
    path = os.path.join(outdir, "exit_particles.csv")
    if not os.path.exists(path):
        return np.empty((0, 5))
    data = _read_csv(path)
    if not data or not data.get("y_m"):
        return np.empty((0, 5))
    return np.column_stack([
        np.array(data["y_m"]),
        np.array(data["vx_ms"]),
        np.array(data["vy_ms"]),
        np.array(data["vz_ms"]),
        np.array(data["ek_eV"]),
    ])


def last_electrode_exit_mm(electrodes: list) -> float:
    """X position of the last electrode's exit face (mm)."""
    return max(float(e["dist"]) + float(e["thick"]) for e in electrodes)


def _read_csv(path):
    data = {}
    with open(path) as f:
        hdr = f.readline().strip().split(",")
        if not hdr or not hdr[0]:
            return {}
        for h in hdr:
            data[h] = []
        for line in f:
            parts = line.strip().split(",")
            if len(parts) != len(hdr):
                continue
            for h, v in zip(hdr, parts):
                try:
                    data[h].append(float(v))
                except ValueError:
                    data[h].append(0.0)
    return data
