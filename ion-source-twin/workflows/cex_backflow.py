"""
Workflow: CEX backflow analysis.
After PY-BEMCS transport, analyze backward-moving charge-exchange ions.
"""
import numpy as np
from adapters import pybemcs_adapter


def analyze(sim) -> dict:
    """
    Extract and analyze CEX backflow from a completed PY-BEMCS simulation.

    Returns dict with:
        count: number of backward-moving CEX macroparticles
        current_A: estimated backflow current
        energies_eV: array of backflow ion energies
        y_mm: radial positions of backflow ions
        mean_energy_eV: average energy
        max_energy_eV: maximum energy
        energy_histogram: (counts, bin_edges) for energy distribution
    """
    bf = pybemcs_adapter.extract_cex_backflow(sim)

    result = dict(bf)
    if bf["count"] > 0:
        result["mean_energy_eV"] = float(np.mean(bf["energies_eV"]))
        result["max_energy_eV"] = float(np.max(bf["energies_eV"]))
        result["energy_histogram"] = np.histogram(bf["energies_eV"], bins=30)
    else:
        result["mean_energy_eV"] = 0.0
        result["max_energy_eV"] = 0.0
        result["energy_histogram"] = (np.array([]), np.array([]))

    return result


def report(results: dict) -> str:
    """Format a human-readable backflow report."""
    if results["count"] == 0:
        return "No CEX backflow detected."
    return (
        f"CEX Backflow Analysis:\n"
        f"  Backward CEX ions: {results['count']}\n"
        f"  Backflow current:  {results['current_A']:.4e} A\n"
        f"  Mean energy:       {results['mean_energy_eV']:.1f} eV\n"
        f"  Max energy:        {results['max_energy_eV']:.1f} eV\n"
    )
