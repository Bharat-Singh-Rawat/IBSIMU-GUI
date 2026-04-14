"""
Workflow: IBSIMU beam extraction + particle export for handoff.
"""
import os
from adapters import ibsimu_adapter, beam_handoff
from config import DATA_DIR


def run(params: dict, electrodes: list, outdir: str = None) -> dict:
    """Run IBSIMU extraction and return parsed results."""
    outdir = outdir or os.path.join(DATA_DIR, "ibsimu")
    return ibsimu_adapter.run_extraction(params, electrodes, outdir)


def run_with_handoff(params: dict, electrodes: list,
                     inject_x_mm: float = 0.1) -> tuple:
    """
    Run extraction and prepare beam handoff data for PY-BEMCS.
    Returns (results_dict, beam_dict_for_pybemcs).
    """
    outdir = os.path.join(DATA_DIR, "ibsimu")
    results = ibsimu_adapter.run_extraction(params, electrodes, outdir)

    # Load exit-plane particles
    particles = ibsimu_adapter.extract_exit_particles(outdir)

    # Transform to PY-BEMCS coordinates
    beam = beam_handoff.ibsimu_to_pybemcs(particles, inject_x_mm)

    # Get beam current from emittance profile (at exit)
    beam_current = 0.0
    if "emittance_profile" in results:
        ep = results["emittance_profile"]
        if ep.get("current_A"):
            x_exit = ibsimu_adapter.last_electrode_exit_mm(electrodes)
            # Find closest x
            x_list = ep["x_mm"]
            idx = min(range(len(x_list)), key=lambda i: abs(x_list[i] - x_exit))
            beam_current = abs(ep["current_A"][idx])

    beam["current_A"] = beam_current
    beam["n_particles"] = len(particles)

    warnings = beam_handoff.validate_handoff(beam)
    if warnings:
        for w in warnings:
            print(f"  Handoff warning: {w}")

    return results, beam
