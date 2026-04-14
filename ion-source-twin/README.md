# Ion Source Digital Twin (Under test phase)

<img width="1920" height="1080" alt="Screenshot from 2026-04-15 00-43-18" src="https://github.com/user-attachments/assets/296c7579-3ac6-4447-8878-70dfae22f27e" />


Unified simulation tool combining [IBSIMU](http://ibsimu.sourceforge.net/) beam extraction with [PY-BEMCS](https://github.com/Bharat-Singh-Rawat/PY-BEMCS) transport, CEX collisions, erosion, and thermal modeling.

## What it does

**Stage 1 (IBSIMU):** Design electrode geometry, optimize extraction optics, compute beam emittance and divergence with Vlasov self-consistency.

**Stage 2 (PY-BEMCS):** Take the extracted beam, propagate through downstream grids with charge-exchange (CEX) collisions, compute sputtering erosion, track grid temperatures.

**Stage 3 (Digital Twin):** Iterate: extraction -> erosion -> geometry update -> re-extract. Predicts grid lifetime and performance degradation.

## Features

### 1. Beam Handoff
- IBSIMU exports full particle state (y, vx, vy, vz, energy) at the exit plane
- PY-BEMCS imports these particles instead of default Bohm injection
- Coordinate transform handled automatically (m -> mm, r -> y)

### 2. Perveance Scan + Erosion
- Sweep voltage or current density
- At each operating point, run both extraction AND transport
- Plot: perveance vs divergence vs erosion rate on the same chart

### 3. Combined 4-Tab GUI
- **Tab 1: Extraction** - IBSIMU electrode design with trajectory, phase space, envelope plots
- **Tab 2: Transport & CEX** - PY-BEMCS particle scatter, divergence history, CEX backflow report
- **Tab 3: Erosion & Thermal** - Damage map, temperature map, grid temperature history
- **Tab 4: Digital Twin** - Automated loop controller + perveance-erosion scan

### 4. CEX Backflow Analysis
- Identifies backward-moving charge-exchange ions
- Reports backflow current, energy spectrum, radial distribution

### 5. Long-Duration Digital Twin Loop
- Iterative: extract beam -> transport/erode -> update aperture -> repeat
- Time-series plots: aperture, divergence, current, damage vs operating hours
- Failure detection: divergence or temperature threshold exceeded

## Prerequisites

1. **IBSIMU** built at `~/ibsimu/gui/beam_sim` (see [IBSIMU-GUI](https://github.com/Bharat-Singh-Rawat/IBSIMU-GUI))
2. **PY-BEMCS** at `~/PYBEMCS/PY-BEMCS/Python/` (see [PY-BEMCS](https://github.com/Bharat-Singh-Rawat/PY-BEMCS))
3. Python 3.10+ with: `numpy matplotlib pillow tkinter`
4. Optional: `taichi` for GPU-accelerated PY-BEMCS transport

## Running

```bash
cd ion-source-twin
python3 main.py
```

### Workflow
1. **Tab 1:** Set electrodes, species, run extraction
2. Click **"Send Beam to Transport Tab >>"**
3. **Tab 2:** Set transport grids and plasma params, run transport
4. **Tab 3:** View erosion damage and temperature maps
5. **Tab 4:** Run multi-cycle digital twin or perveance+erosion scan

## Project Structure

```
ion-source-twin/
├── main.py                    # Entry point
├── config.py                  # Paths, constants, defaults
├── adapters/
│   ├── ibsimu_adapter.py      # Wraps IBSIMU C++ binary
│   ├── pybemcs_adapter.py     # Wraps PY-BEMCS DigitalTwinSimulator
│   └── beam_handoff.py        # Coordinate transform (m->mm)
├── workflows/
│   ├── extraction.py          # IBSIMU run + particle export
│   ├── transport_cex.py       # PY-BEMCS with external beam
│   ├── cex_backflow.py        # Backward CEX ion analysis
│   ├── perveance_erosion.py   # Coupled scan
│   └── digital_twin_loop.py   # Iterative erosion loop
├── gui/
│   ├── app.py                 # Main window with 4 tabs
│   ├── tab_extraction.py      # Tab 1: IBSIMU
│   ├── tab_transport.py       # Tab 2: PY-BEMCS
│   ├── tab_erosion.py         # Tab 3: Erosion/thermal
│   ├── tab_twin.py            # Tab 4: Digital twin
│   └── plot_helpers.py        # Shared matplotlib utilities
└── data/                      # Runtime outputs
```

## License

Provided as-is. IBSIMU: GNU GPL v2+ (Taneli Kalvas). PY-BEMCS: see its repository.
