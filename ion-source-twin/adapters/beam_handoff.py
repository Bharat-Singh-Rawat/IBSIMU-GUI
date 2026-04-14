"""
Coordinate transform and particle format conversion between IBSIMU and PY-BEMCS.

IBSIMU:   x=axial (m), r=radial (m), velocities in m/s
PY-BEMCS: x=axial (mm), y=radial (mm), velocities in m/s
"""
import numpy as np
from config import M_TO_MM


def ibsimu_to_pybemcs(particles: np.ndarray,
                      inject_x_mm: float = 0.1) -> dict:
    """
    Transform IBSIMU exit-plane particles to PY-BEMCS injection format.

    Input:  (N, 5) array: [y_m, vx_m/s, vy_m/s, vz_m/s, ek_eV]
    Output: dict with arrays: x_mm, y_mm, vx, vy, vz (ready for _add_ions)
    """
    if len(particles) == 0:
        return {"x_mm": np.array([]), "y_mm": np.array([]),
                "vx": np.array([]), "vy": np.array([]), "vz": np.array([]),
                "ek_eV": np.array([])}

    y_m = particles[:, 0]
    vx = particles[:, 1]
    vy = particles[:, 2]
    vz = particles[:, 3]
    ek = particles[:, 4]

    y_mm = np.abs(y_m) * M_TO_MM   # radial coordinate is always positive in PY-BEMCS
    x_mm = np.full_like(y_mm, inject_x_mm)

    return {
        "x_mm": x_mm.astype(np.float32),
        "y_mm": y_mm.astype(np.float32),
        "vx": vx.astype(np.float32),
        "vy": vy.astype(np.float32),
        "vz": vz.astype(np.float32),
        "ek_eV": ek.astype(np.float32),
    }


def pybemcs_aperture_to_ibsimu(aperture_mm: float) -> float:
    """Convert PY-BEMCS aperture (mm) to IBSIMU electrode config value (mm)."""
    return max(0.1, aperture_mm)


def estimate_beam_current(particles: np.ndarray, beam_charge: float) -> float:
    """
    Estimate total beam current from particle data.
    This is approximate; the actual current comes from emittance_profile.csv.
    """
    # For a rough estimate: I = N * q * v_avg / L (not very accurate)
    # Better to use the current_A from IBSIMU's emittance profile
    return 0.0


def validate_handoff(beam: dict) -> list:
    """Sanity-check the handed-off beam data. Returns list of warnings."""
    warnings = []
    if len(beam["y_mm"]) == 0:
        warnings.append("No particles in handoff data")
        return warnings
    if np.any(beam["y_mm"] < 0):
        warnings.append("Negative y positions found (should be >= 0)")
    if np.any(beam["ek_eV"] < 0):
        warnings.append("Negative kinetic energies found")
    vmax = np.max(np.abs(beam["vx"]))
    if vmax > 3e8:
        warnings.append(f"Superluminal vx detected: {vmax:.2e} m/s")
    return warnings
