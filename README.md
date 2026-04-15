
# IBSIMU-GUI

<img width="800" height="800" alt="Screenshot from 2026-04-14 22-20-09" src="https://github.com/user-attachments/assets/8c7253a5-390f-4014-aae2-9bdcbce22bad" />

Interactive GUI for axisymmetric ion/electron beam extraction simulation using [IBSimu](http://ibsimu.sourceforge.net/) (Ion Beam Simulator).


## An example of a perveance scan of a D+ beam between 2-300kV range.

![Perveancescan_300_2kV_d+](https://github.com/user-attachments/assets/a99f80aa-bc2c-4cca-85ab-66468298a432)


## Features

### Geometry & Electrodes
- **Multi-electrode geometry** - up to 10 electrodes with individual gap, aperture radius, voltage, thickness, chamfer angle, and wall radius
- **Gap-based positioning** - specify gap between electrodes rather than absolute positions
- **Independent wall radius** - electrode wall height decoupled from aperture (0 = auto 3x aperture)
- Input validation with overlap detection

### Beam Species & Definition
- **Beam species selector** - presets for e-, H+, D+, He+, He2+, Ar+, Xe+, or custom mass/charge
- **Beam definition modes**:
  - **Energy/Temperature** - define beam by kinetic energy and thermal spread (default)
  - **Twiss parameters** - define beam by alpha, beta, emittance (Gaussian distribution) for beam matching
- Electron beam support with automatic plasma model handling

### Solvers
- **BiCGSTAB** - iterative solver (default, good for most cases)
- **Multigrid** - configurable levels, faster for large grids with compatible dimensions
- **Gauss-Seidel** - simple iterative solver for validation

### Magnetic Field
- **Solenoid field** - uniform axial B-field with configurable strength (Tesla), start/end position
- Properly handled in particle trajectory integration

### Diagnostics & Plots
- **Particle trajectories** - axisymmetric geometry plot with trajectory density
- **Phase-space plot** (y, y') with density-based coloring (KDE) and colorbar, RMS emittance ellipse overlay, adjustable via slider
- **Beam profile** - spatial histogram of particle positions at any axial location with RMS radius markers
- **Beam envelope, divergence & transmission** profile along the beam path
- **Field diagnostics** - electric potential, E-field, and B-field along the beam axis
- **Convergence history** - log-scale plot of potential and space charge convergence per iteration
- **Energy distribution** - histogram of particle kinetic energy at the exit
- **Per-electrode current** - current intercepted by each electrode and domain boundary
- **Grid current analysis** - transmission % profile, grid current ratio in beam info panel

### Perveance Scan
- Sweep **voltage** or **current density** with real-time plot updates
- Signed voltage values: negative for ions, positive for electrons
- **Configurable density override** for voltage scans
- **Configurable divergence measurement location** - specify distance (mm) past last electrode exit
- Divergence displayed in **degrees**
- **Dual y-axis**: divergence + grid current % vs perveance
- All plots update live during scan (trajectory, phase space, field diagnostics, convergence, energy)
- Optimal point (minimum divergence) marked with red star

### Matched Beam Finder
- **Golden section search** to automatically find the electrode voltage that minimizes beam divergence
- Configurable electrode, voltage range (min/max), and convergence tolerance
- Typically converges in ~10 simulation evaluations
- Updates all plots with the optimal configuration when complete

### Auto-Optimizer
- **Multi-parameter grid search** over any combination of voltage, gap, and aperture
- Enable/disable each parameter independently with individual min/max ranges
- Configurable number of steps per axis
- Tracks divergence, transmission, and current for every evaluated combination
- Automatically restores the best configuration found
- Results stored for inclusion in PDF reports

### Export & Reporting
- **Save All Plots (PNG)** - 8 high-resolution plots (trajectory, phase space, beam profile, envelope, perveance scan, field diagnostics, convergence, energy distribution)
- **Save Scan GIF** - animated GIF of beam evolution during perveance scan
- **PDF Report Generator** - multi-page report containing:
  - Simulation configuration summary (beam parameters, electrode table)
  - Key results at electrode exit (divergence, RMS radius, current, transmission)
  - All diagnostic plots (trajectory, phase space, beam profile, envelope, field, convergence, energy)
  - Perveance scan results (if available)
  - Optimizer results table with all evaluated parameter combinations (if available)

## Prerequisites

### System dependencies

```bash
# Ubuntu/Debian
sudo apt install build-essential pkg-config \
    libgsl-dev libcairo2-dev libgtk-3-dev libpng-dev \
    libfreetype-dev libfontconfig-dev zlib1g-dev \
    python3 python3-tk python3-matplotlib python3-pil python3-numpy python3-scipy
```

### IBSimu library

```bash
git clone git://ibsimu.git.sourceforge.net/gitroot/ibsimu/ibsimu
cd ibsimu
./reconf
./configure --prefix=$HOME/ibsimu-install
make -j$(nproc)
make install
```

## Building

```bash
git clone https://github.com/Bharat-Singh-Rawat/IBSIMU-GUI.git
cd IBSIMU-GUI
# Edit Makefile: set IBSIMU_SRC and IBSIMU_LIB to your ibsimu build paths
make
```

## Running

```bash
python3 beam_gui.py
```

### Quick start

1. Select beam species and definition mode (Energy or Twiss)
2. Choose solver (BiCGSTAB recommended)
3. Optionally enable solenoid magnetic field
4. Set electrode parameters
5. Click **Run Simulation**
6. Explore results with the slider and tabs:
   - **Beam Profile** - spatial distribution histogram at selected axial position
   - **Envelope** - beam size, divergence, transmission along axis
   - **Field Diag** - potential and E-field profile
   - **Convergence** - iteration convergence history
   - **Perveance Scan** - sweep voltage/current for optimization
   - **Energy Dist.** - particle energy histogram at exit
7. Use **Find Matched Beam** for automatic voltage optimization
8. Use **Run Optimizer** for multi-parameter sweeps
9. Click **Generate Report (PDF)** to export a complete summary

### Voltage conventions

| Particle | Extraction electrode | Scan range |
|----------|---------------------|------------|
| Positive ions (H+, Ar+) | Negative (e.g. -8 kV) | -2 to -15 kV |
| Electrons | Positive (e.g. +8 kV) | +2 to +15 kV |

## Packaging as standalone application

```bash
pip install pyinstaller matplotlib pillow numpy
chmod +x build_package.sh
./build_package.sh
# Output: dist/IBSimuGUI/run_ibsimu.sh
```

## File structure

| File | Description |
|------|-------------|
| `beam_gui.py` | Python/tkinter GUI with matplotlib plots |
| `beam_sim.cpp` | C++ simulation backend using IBSimu |
| `Makefile` | Builds `beam_sim` |
| `build_package.sh` | PyInstaller packaging script |

## How it works

1. GUI writes config file with all parameters (geometry, beam, solver, B-field)
2. `beam_sim` reads config, builds geometry, runs Vlasov iteration
3. Outputs: trajectory PNG, emittance/phase-space CSV, field diagnostics CSV, convergence history, electrode currents, energy distribution
4. GUI loads and displays results interactively

## License

GUI code provided as-is. IBSimu is licensed under GNU GPL v2+ by Taneli Kalvas.

## Acknowledgements

- [IBSimu](http://ibsimu.sourceforge.net/) by Taneli Kalvas (University of Jyvaskyla)
- Built with [matplotlib](https://matplotlib.org/), [tkinter](https://docs.python.org/3/library/tkinter.html), [PyInstaller](https://pyinstaller.org/)
