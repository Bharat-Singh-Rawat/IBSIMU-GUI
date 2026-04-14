"""
Adapter for the PY-BEMCS DigitalTwinSimulator.
Handles import, beam injection override, CEX backflow extraction, erosion readout.
"""
import sys, os, numpy as np
from config import PYBEMCS_DIR, ECHARGE, AMU_TO_KG

# Add PY-BEMCS to path so we can import physics_engine
if PYBEMCS_DIR not in sys.path:
    sys.path.insert(0, PYBEMCS_DIR)


def create_simulator(mass_amu: float, charge_state: int, material: dict = None):
    """Create and configure a DigitalTwinSimulator instance."""
    from physics_engine import DigitalTwinSimulator
    sim = DigitalTwinSimulator()
    sim.m_ion = mass_amu * AMU_TO_KG
    sim.Z_ion = abs(charge_state)
    sim.q_ion = sim.Z_ion * ECHARGE
    if material:
        sim.set_material(props=material)
    return sim


def build_domain(sim, grids: list, plasma_params: dict):
    """Build PY-BEMCS simulation domain."""
    params = dict(plasma_params)
    params["grids"] = grids
    params.setdefault("Lx", 20.0)
    params.setdefault("Ly", 3.0)
    params.setdefault("beam_mass_amu", sim.m_ion / AMU_TO_KG)
    params.setdefault("beam_charge_state", sim.Z_ion)
    sim.build_domain(params)
    return params


def make_suppressed_params(base_params: dict) -> dict:
    """Return params copy with near-zero plasma density to suppress default injection."""
    p = dict(base_params)
    p["n0_plasma"] = 1.0  # effectively zero -> ~0 particles injected per step
    return p


def inject_beam(sim, beam: dict, fraction: float = 1.0):
    """
    Inject externally provided particles into the simulator.
    beam: dict with x_mm, y_mm, vx, vy, vz arrays.
    fraction: inject only this fraction of particles (for rate matching).
    """
    n = len(beam["x_mm"])
    if n == 0:
        return
    if fraction < 1.0:
        mask = np.random.rand(n) < fraction
        x = beam["x_mm"][mask]
        y = beam["y_mm"][mask]
        vx = beam["vx"][mask]
        vy = beam["vy"][mask]
        vz = beam["vz"][mask]
    else:
        x, y = beam["x_mm"], beam["y_mm"]
        vx, vy, vz = beam["vx"], beam["vy"], beam["vz"]

    is_cex = np.zeros(len(x), dtype=bool)
    sim._add_ions(x.astype(np.float32), y.astype(np.float32),
                  vx.astype(np.float32), vy.astype(np.float32),
                  vz.astype(np.float32), is_cex)


def run_steps(sim, params: dict, n_steps: int, beam: dict = None,
              inject_every: int = 1, callback=None) -> dict:
    """
    Run PY-BEMCS for n_steps iterations.
    If beam is provided, injects particles periodically instead of default injection.
    callback(step, n_steps, div, t_grids) called each step for progress.
    """
    if beam is not None:
        params = make_suppressed_params(params)
        # Compute injection fraction per step to match beam current
        inject_frac = 1.0 / max(1, inject_every)

    history = {"div": [], "min_pot": [], "t_grids": [], "remeshed": []}

    for i in range(n_steps):
        if beam is not None and (i % inject_every == 0):
            inject_beam(sim, beam, fraction=inject_frac)

        remeshed, min_pot, div, t_grids = sim.step(params)
        history["div"].append(div)
        history["min_pot"].append(min_pot)
        history["t_grids"].append(list(t_grids) if t_grids else [])
        history["remeshed"].append(remeshed)

        if callback and i % 5 == 0:
            callback(i, n_steps, div, t_grids)

    return history


def extract_cex_backflow(sim) -> dict:
    """Extract backward-moving CEX ions for backflow analysis."""
    n = sim.num_p
    if n == 0:
        return {"count": 0, "current_A": 0, "energies_eV": np.array([]),
                "y_mm": np.array([])}

    cex = sim.p_isCEX[:n]
    vx = sim.p_vx[:n]
    backward = cex & (vx < 0)
    count = int(np.sum(backward))

    if count == 0:
        return {"count": 0, "current_A": 0, "energies_eV": np.array([]),
                "y_mm": np.array([])}

    v_sq = vx[backward]**2 + sim.p_vy[:n][backward]**2 + sim.p_vz[:n][backward]**2
    energies = 0.5 * sim.m_ion * v_sq / ECHARGE

    # Rough current estimate
    current = sim.q_ion * sim.macro_weight * count / (sim.dt * 1000)

    return {
        "count": count,
        "current_A": current,
        "energies_eV": energies,
        "y_mm": sim.p_y[:n][backward].copy(),
    }


def extract_eroded_aperture(sim, grid_index: int, grids: list) -> float:
    """
    Find effective aperture after erosion by scanning isBound at grid center.
    Returns aperture radius in mm.
    """
    # Compute grid x-center position
    current_x = 1.0  # PY-BEMCS grids start at x=1mm
    for i, g in enumerate(grids):
        g_start = current_x
        g_end = g_start + g["t"]
        if i == grid_index:
            x_center = (g_start + g_end) / 2.0
            break
        current_x = g_end + g["gap"]
    else:
        return grids[grid_index]["r"]

    # Find the grid column index closest to x_center
    ix = int(round(x_center / sim.dx))
    ix = min(ix, sim.isBound.shape[1] - 1)

    # Scan from y=0 upward to find first solid cell
    for jy in range(sim.isBound.shape[0]):
        if sim.isBound[jy, ix]:
            return jy * sim.dy  # aperture radius in mm

    # No solid found at all - grid fully eroded
    return sim.Ly


def get_damage_map(sim) -> np.ndarray:
    """Return a copy of the cumulative damage map."""
    return sim.damage_map.copy() if hasattr(sim, "damage_map") else np.zeros((1, 1))


def get_temperature_map(sim) -> np.ndarray:
    """Return a copy of the temperature map."""
    return sim.T_map.copy() if hasattr(sim, "T_map") else np.full((1, 1), 300.0)
