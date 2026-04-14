#!/usr/bin/env python3
"""
Ion Source Digital Twin
Combines IBSIMU beam extraction with PY-BEMCS transport, CEX, erosion, and thermal modeling.
"""
import sys
import os
import tkinter as tk

# Ensure project root is in path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config import IBSIMU_BINARY, PYBEMCS_DIR, DATA_DIR

def check_dependencies():
    errors = []
    if not os.path.exists(IBSIMU_BINARY):
        errors.append(f"IBSIMU binary not found: {IBSIMU_BINARY}\n"
                      f"  Build it: cd ~/ibsimu/gui && make")
    if not os.path.exists(os.path.join(PYBEMCS_DIR, "physics_engine.py")):
        errors.append(f"PY-BEMCS not found: {PYBEMCS_DIR}\n"
                      f"  Clone: git clone https://github.com/Bharat-Singh-Rawat/PY-BEMCS.git")
    return errors

def main():
    errors = check_dependencies()
    if errors:
        print("Missing dependencies:")
        for e in errors:
            print(f"  {e}")
        sys.exit(1)

    os.makedirs(DATA_DIR, exist_ok=True)

    from gui.app import IonSourceTwinApp
    root = tk.Tk()
    app = IonSourceTwinApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
