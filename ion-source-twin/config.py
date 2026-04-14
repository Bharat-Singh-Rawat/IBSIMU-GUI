"""
Central configuration: paths, constants, defaults.
"""
import os

# ---- Paths ----
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(PROJECT_DIR, "data")

IBSIMU_BINARY = os.path.expanduser("~/ibsimu/gui/beam_sim")
PYBEMCS_DIR = os.path.expanduser("~/PYBEMCS/PY-BEMCS/Python")

# ---- Unit Conversions ----
M_TO_MM = 1e3
MM_TO_M = 1e-3
MRAD_TO_RAD = 1e-3
EV_TO_J = 1.602176634e-19
AMU_TO_KG = 1.66053906660e-27
ECHARGE = 1.602176634e-19

# ---- Default Simulation Parameters ----
DEFAULT_EXTRACTION = {
    "grid_points": 241,
    "particles": 5000,
    "current_density": 600.0,
    "beam_energy": 5.0,
    "beam_mass": 1.00728,
    "beam_charge": 1.0,
    "solver": "bicgstab",
    "beam_mode": "energy",
    "bfield": "none",
}

DEFAULT_TRANSPORT = {
    "Lx": 20.0,          # domain length (mm)
    "Ly": 3.0,           # domain height (mm)
    "n0_plasma": 1e17,   # plasma density (m^-3)
    "Te_up": 3.0,        # electron temp (eV)
    "Ti": 2.0,           # ion temp (eV)
    "Tn": 300,           # neutral temp (K)
    "n0": 1e20,          # neutral density (m^-3)
    "Accel": 1.0,        # erosion acceleration factor
    "Thresh": 10000.0,   # cell failure threshold
    "sim_mode": "Both",
}

DEFAULT_MATERIAL = {
    "name": "Molybdenum",
    "k": 138.0,
    "rho": 10280.0,
    "cp": 250.0,
    "emissivity": 0.80,
    "alpha": 4.8e-6,
    "E_mod": 329e9,
    "Y_coeff": 1.05e-4,
    "E_th": 30.0,
}

# Species presets: (mass_amu, charge_state)
SPECIES = {
    "Xe+": (131.293, 1), "Kr+": (83.798, 1), "Ar+": (39.948, 1),
    "N2+": (28.014, 1), "O2+": (31.998, 1), "H+": (1.00728, 1),
    "H2+": (2.014, 1), "D+": (2.014, 1), "He+": (4.003, 1),
    "He2+": (4.003, 2), "e-": (0.000548579909, -1),
}
