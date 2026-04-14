"""
Workflow: Perveance scan coupled with erosion prediction.
At each operating point, run extraction + transport to get both
beam quality and grid degradation rate.
"""
import numpy as np
from workflows import extraction, transport_cex
from config import ECHARGE


def run_scan(base_params: dict, electrodes: list, bemcs_grids: list,
             scan_type: str, scan_electrode: int,
             scan_min: float, scan_max: float, scan_steps: int,
             plasma_params: dict = None, material: dict = None,
             erosion_steps: int = 200, callback=None) -> dict:
    """
    Run a perveance scan where each point also runs PY-BEMCS for erosion.

    scan_type: "voltage" or "current"
    scan_electrode: electrode index (1-based) for voltage scan
    callback: (step, total, scan_val, results_so_far) called after each point

    Returns dict with arrays:
        scan_values, perveance, divergence, erosion_rate, grid_ratio,
        backflow_current, labels
    """
    values = np.linspace(scan_min, scan_max, scan_steps)
    results = {
        "scan_values": [], "perveance": [], "divergence": [],
        "erosion_rate": [], "grid_ratio": [], "backflow_current": [],
        "labels": [],
    }

    for i, val in enumerate(values):
        # Modify parameters for this scan point
        if scan_type == "voltage":
            elecs = [dict(e) for e in electrodes]
            elecs[scan_electrode - 1]["volt"] = val
            label = f"{val:+.1f} kV"
            params = dict(base_params)
        else:
            params = dict(base_params)
            params["current_density"] = val
            elecs = electrodes
            label = f"{val:.0f} A/m\u00b2"

        try:
            # Run IBSIMU extraction + handoff
            ibsimu_results, beam = extraction.run_with_handoff(params, elecs)

            # Compute perveance
            if scan_type == "voltage":
                v_volts = abs(val) * 1e3
            else:
                v_kv = float(electrodes[1]["volt"]) if len(electrodes) > 1 else 8.0
                v_volts = abs(v_kv) * 1e3

            if v_volts > 0 and beam["current_A"] > 0:
                perv = (abs(beam["current_A"]) / (v_volts ** 1.5)) * 1e6
            else:
                perv = 0

            # Get divergence at exit
            ep = ibsimu_results.get("emittance_profile", {})
            div = ep.get("divergence_mrad", [0])[-1] if ep.get("divergence_mrad") else 0

            # Run PY-BEMCS transport for erosion estimate
            mass_amu = float(params.get("beam_mass", 1.0))
            charge = int(abs(float(params.get("beam_charge", 1))))

            tr = transport_cex.run(
                beam, bemcs_grids, plasma_params, material,
                mass_amu, charge, n_steps=erosion_steps)

            # Erosion rate from damage map
            dmg = tr["damage_map"]
            erosion = float(np.sum(dmg)) / max(1, erosion_steps)

            # Grid ratio from IBSIMU
            curr_list = ep.get("current_A", [])
            if len(curr_list) >= 2:
                gr = 1.0 - abs(curr_list[-1]) / abs(curr_list[0])
            else:
                gr = 0

            # Backflow
            bf_current = tr["backflow"]["current_A"]

        except Exception as e:
            print(f"  Scan point {i+1} failed: {e}")
            perv, div, erosion, gr, bf_current = 0, 0, 0, 0, 0

        results["scan_values"].append(val)
        results["perveance"].append(perv)
        results["divergence"].append(div)
        results["erosion_rate"].append(erosion)
        results["grid_ratio"].append(gr * 100)
        results["backflow_current"].append(bf_current)
        results["labels"].append(label)

        if callback:
            callback(i + 1, scan_steps, label, results)

    return results
