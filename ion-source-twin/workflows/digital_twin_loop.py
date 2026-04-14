"""
Workflow: Long-duration digital twin.
Iterative loop: extraction -> transport/erosion -> geometry update -> repeat.
Predicts grid lifetime and performance degradation over time.
"""
import numpy as np
from workflows import extraction, transport_cex
from adapters import pybemcs_adapter, beam_handoff


def run(initial_params: dict, initial_electrodes: list,
        bemcs_grids: list, plasma_params: dict = None,
        material: dict = None, mass_amu: float = 1.0, charge_state: int = 1,
        steps_per_cycle: int = 500, max_cycles: int = 10,
        hours_per_cycle: float = 100.0,
        failure_div_mrad: float = 50.0, failure_temp_K: float = 2000.0,
        callback=None) -> dict:
    """
    Run the iterative digital twin loop.

    Each cycle:
    1. Run IBSIMU extraction with current geometry
    2. Hand off beam to PY-BEMCS
    3. Run PY-BEMCS for steps_per_cycle iterations
    4. Extract eroded aperture from PY-BEMCS
    5. Update electrode geometry
    6. Check failure criteria
    7. Repeat

    callback: (cycle, max_cycles, cycle_result) called after each cycle

    Returns dict with time-series arrays:
        hours, aperture_mm (per grid), divergence_mrad, grid_temp_K,
        erosion_damage, beam_current_A, backflow_current_A, failed_at_cycle
    """
    history = {
        "hours": [],
        "apertures": [],        # list of lists (one per grid)
        "divergence": [],
        "grid_temps": [],       # list of lists
        "total_damage": [],
        "beam_current": [],
        "backflow_current": [],
    }

    current_electrodes = [dict(e) for e in initial_electrodes]
    current_bemcs_grids = [dict(g) for g in bemcs_grids]
    failed_at = None
    sim_hours = 0.0

    for cycle in range(max_cycles):
        sim_hours += hours_per_cycle

        # 1. Run IBSIMU extraction
        try:
            ibsimu_results, beam = extraction.run_with_handoff(
                initial_params, current_electrodes)
        except Exception as e:
            print(f"Cycle {cycle+1}: Extraction failed: {e}")
            failed_at = cycle
            break

        # Get divergence and current
        ep = ibsimu_results.get("emittance_profile", {})
        div_list = ep.get("divergence_mrad", [0])
        div = div_list[-1] if div_list else 0
        curr = beam.get("current_A", 0)

        # 2-3. Run PY-BEMCS transport + erosion
        try:
            tr = transport_cex.run(
                beam, current_bemcs_grids, plasma_params, material,
                mass_amu, charge_state, n_steps=steps_per_cycle)
        except Exception as e:
            print(f"Cycle {cycle+1}: Transport failed: {e}")
            failed_at = cycle
            break

        sim = tr["sim"]
        backflow = tr["backflow"]

        # 4. Extract eroded apertures
        apertures = []
        for gi in range(len(current_bemcs_grids)):
            new_r = pybemcs_adapter.extract_eroded_aperture(
                sim, gi, current_bemcs_grids)
            apertures.append(new_r)

        # Grid temperatures
        t_grids = []
        if tr["history"]["t_grids"]:
            t_grids = tr["history"]["t_grids"][-1]

        # Total damage
        total_dmg = float(np.sum(tr["damage_map"]))

        # Record history
        history["hours"].append(sim_hours)
        history["apertures"].append(apertures)
        history["divergence"].append(div)
        history["grid_temps"].append(t_grids)
        history["total_damage"].append(total_dmg)
        history["beam_current"].append(curr)
        history["backflow_current"].append(backflow["current_A"])

        # 5. Update geometry for next cycle
        # Map PY-BEMCS eroded apertures back to IBSIMU electrodes
        # Only update electrodes that correspond to BEMCS grids
        for gi, new_r in enumerate(apertures):
            if gi + 1 < len(current_electrodes):  # grid 0 in BEMCS -> electrode 1 in IBSIMU
                old_r = float(current_electrodes[gi + 1]["apt"])
                updated_r = beam_handoff.pybemcs_aperture_to_ibsimu(new_r)
                if updated_r > old_r:
                    current_electrodes[gi + 1]["apt"] = str(updated_r)
                    current_bemcs_grids[gi]["r"] = updated_r

        cycle_result = {
            "cycle": cycle + 1,
            "hours": sim_hours,
            "divergence": div,
            "apertures": apertures,
            "grid_temps": t_grids,
            "damage": total_dmg,
            "current": curr,
            "backflow": backflow["current_A"],
        }

        if callback:
            callback(cycle, max_cycles, cycle_result)

        # 6. Check failure criteria
        if div > failure_div_mrad:
            print(f"Cycle {cycle+1}: Divergence {div:.1f} mrad exceeds limit {failure_div_mrad}")
            failed_at = cycle + 1
            break
        if t_grids and max(t_grids) > failure_temp_K:
            print(f"Cycle {cycle+1}: Grid temp {max(t_grids):.0f} K exceeds limit {failure_temp_K}")
            failed_at = cycle + 1
            break

    history["failed_at_cycle"] = failed_at
    history["total_hours"] = sim_hours
    return history
