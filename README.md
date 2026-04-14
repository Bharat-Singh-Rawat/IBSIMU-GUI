
# IBSIMU-GUI


<img width="800" height="800" alt="Screenshot from 2026-04-14 22-20-09" src="https://github.com/user-attachments/assets/8c7253a5-390f-4014-aae2-9bdcbce22bad" />


Interactive GUI for axisymmetric ion beam extraction simulation using [IBSimu](http://ibsimu.sourceforge.net/) (Ion Beam Simulator).

## Features

- **Multi-electrode geometry** - define any number of electrodes with individual position, aperture radius, voltage, and chamfer angle
- **Axisymmetric (cylindrical) particle trajectory** simulation with space charge (Vlasov iteration)
- **Phase-space plot** (y, y') with RMS emittance ellipse overlay, adjustable via slider along the beam axis
- **Beam envelope and divergence** profile along the beam path
- **Perveance scan** - sweep voltage or current density and plot perveance vs. divergence in real time to find the optimal operating point
- **Standalone packaging** - build as a distributable application with `build_package.sh`

## Prerequisites

### System dependencies

```bash
# Ubuntu/Debian
sudo apt install build-essential pkg-config \
    libgsl-dev libcairo2-dev libgtk-3-dev libpng-dev \
    libfreetype-dev libfontconfig-dev zlib1g-dev \
    python3 python3-tk python3-matplotlib python3-pil python3-numpy
```

### IBSimu library

IBSIMU must be built and available on your system. Clone and build:

```bash
git clone git://ibsimu.git.sourceforge.net/gitroot/ibsimu/ibsimu
cd ibsimu
./reconf
./configure --prefix=$HOME/ibsimu-install
make -j$(nproc)
make install
```

## Building

1. Clone this repository:
```bash
git clone https://github.com/Bharat-Singh-Rawat/IBSIMU-GUI.git
cd IBSIMU-GUI
```

2. Edit the `Makefile` to point `IBSIMU_SRC` and `IBSIMU_LIB` to your IBSIMU build:
```makefile
IBSIMU_SRC = /path/to/ibsimu/src
IBSIMU_LIB = /path/to/ibsimu/src/.libs
```

3. Build the C++ simulation backend:
```bash
make
```

## Running

```bash
python3 beam_gui.py
```

### Quick start

1. Set electrode parameters (defaults: plasma electrode at 0 V, puller at -8 kV)
2. Click **Run Simulation**
3. Use the **slider** to explore the beam emittance at different axial positions
4. Switch to the **Perveance Scan** tab to sweep voltage/current and find the optimal perveance

### Perveance scan

1. Select scan type: **Voltage scan** (vary extraction voltage) or **Current scan** (vary beam current density)
2. Set the electrode number to scan, min/max range, and number of steps
3. Click **Run Scan** - the perveance vs. divergence curve updates in real time
4. The optimal point (minimum divergence) is marked with a red star

## Packaging as standalone application

To create a distributable package (bundles Python + C++ binary + all libraries):

```bash
pip install pyinstaller matplotlib pillow numpy
chmod +x build_package.sh
./build_package.sh
```

The output is in `dist/IBSimuGUI/`. Distribute the folder as a tar.gz:
```bash
cd dist && tar czf IBSimuGUI-linux-x86_64.tar.gz IBSimuGUI/
```

Users run it with:
```bash
./IBSimuGUI/run_ibsimu.sh
```

## File structure

| File | Description |
|------|-------------|
| `beam_gui.py` | Python/tkinter GUI with matplotlib plots |
| `beam_sim.cpp` | C++ simulation backend (multi-electrode, cylindrical geometry) |
| `Makefile` | Builds `beam_sim` from `beam_sim.cpp` |
| `build_package.sh` | Creates a distributable PyInstaller bundle |

## How it works

1. The GUI writes electrode parameters to a config file
2. `beam_sim` reads the config, builds the IBSIMU geometry, runs 5 Vlasov iterations with space-charge self-consistency
3. Outputs: trajectory plot (PNG), emittance profile (CSV), phase-space scatter data (CSV)
4. The GUI displays the results interactively

## License

The GUI code in this repository is provided as-is. IBSimu itself is licensed under the GNU GPL v2+ by Taneli Kalvas.

## Acknowledgements

- [IBSimu](http://ibsimu.sourceforge.net/) by Taneli Kalvas (University of Jyvaskyla)
- Built with [matplotlib](https://matplotlib.org/), [tkinter](https://docs.python.org/3/library/tkinter.html), and [PyInstaller](https://pyinstaller.org/)
