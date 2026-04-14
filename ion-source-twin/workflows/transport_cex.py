"""
Workflow: PY-BEMCS transport simulation with externally injected beam.
"""
from adapters import pybemcs_adapter
from config import DEFAULT_TRANSPORT, DEFAULT_MATERIAL


def run(beam: dict, grids: list, plasma_params: dict = None,
        material: dict = None, mass_amu: float = 1.0, charge_state: int = 1,
        n_steps: int = 500, callback=None) -> dict:
    """
    Run PY-BEMCS transport with an externally provided beam.

    beam: dict from beam_handoff.ibsimu_to_pybemcs() with current_A field
    grids: list of dicts with keys: V, t, gap, r, cham
    plasma_params: PY-BEMCS plasma parameters (defaults from config)
    material: electrode material properties
    n_steps: number of PIC iterations
    callback: progress callback(step, total, div, t_grids)

    Returns dict with: history, backflow, damage_map, T_map, sim reference
    """
    plasma = dict(DEFAULT_TRANSPORT)
    if plasma_params:
        plasma.update(plasma_params)
    mat = material or DEFAULT_MATERIAL

    sim = pybemcs_adapter.create_simulator(mass_amu, charge_state, mat)
    params = pybemcs_adapter.build_domain(sim, grids, plasma)

    # Compute injection rate: how often to inject a batch of particles
    # PY-BEMCS dt = 1e-9s, beam_current in A
    # Particles per step = I * dt / (q * macro_weight)
    if beam.get("current_A", 0) > 0 and beam["n_particles"] > 0:
        from config import ECHARGE
        parts_per_step = (beam["current_A"] * sim.dt) / (sim.q_ion * sim.macro_weight)
        # inject_every = how many steps between injections of the full batch
        if parts_per_step > 0:
            inject_every = max(1, int(beam["n_particles"] / parts_per_step))
        else:
            inject_every = 10
    else:
        inject_every = 10

    history = pybemcs_adapter.run_steps(
        sim, params, n_steps, beam=beam,
        inject_every=inject_every, callback=callback)

    backflow = pybemcs_adapter.extract_cex_backflow(sim)
    damage = pybemcs_adapter.get_damage_map(sim)
    temp = pybemcs_adapter.get_temperature_map(sim)

    return {
        "history": history,
        "backflow": backflow,
        "damage_map": damage,
        "T_map": temp,
        "sim": sim,
        "grids": grids,
        "params": params,
    }
