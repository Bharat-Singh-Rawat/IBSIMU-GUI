#!/usr/bin/env python3
"""
IBSIMU Beam Extraction GUI - Multi-electrode with perveance scan
"""

import os
import sys
import subprocess
import threading
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from PIL import Image

def _get_base_dir():
    """Return base directory - works both in dev and PyInstaller bundle."""
    if getattr(sys, 'frozen', False):
        # Running inside PyInstaller bundle
        return sys._MEIPASS
    return os.path.dirname(os.path.abspath(__file__))

def _get_output_dir():
    """Writable output directory for simulation results."""
    if getattr(sys, 'frozen', False):
        # When bundled, use a folder next to the executable
        return os.path.join(os.path.dirname(sys.executable), "sim_output")
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "sim_output")

BASE_DIR = _get_base_dir()
SIM_BINARY = os.path.join(BASE_DIR, "simbin", "beam_sim") if getattr(sys, 'frozen', False) \
    else os.path.join(os.path.dirname(os.path.abspath(__file__)), "beam_sim")
OUTPUT_DIR = _get_output_dir()
CONFIG_FILE = os.path.join(OUTPUT_DIR, "beam_config.txt")


def compute_twiss(y, yp):
    n = len(y)
    if n < 3:
        return 0, 0, 0, 0
    y_c, yp_c = y - np.mean(y), yp - np.mean(yp)
    syy = np.mean(y_c**2)
    spp = np.mean(yp_c**2)
    syp = np.mean(y_c * yp_c)
    det = syy * spp - syp**2
    if det <= 0:
        return 0, 0, 0, 0
    eps = np.sqrt(det)
    return eps, -syp / eps, syy / eps, spp / eps


class BeamExtractionGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("IBSIMU - Multi-Electrode Beam Extraction")
        self.root.geometry("1500x950")
        self.root.minsize(1200, 700)

        self.emittance_data = None
        self.phase_space = None
        self.traj_image = None

        # Scan state
        self.scan_running = False
        self.scan_stop = False
        self.scan_data = {"perveance": [], "divergence": [],
                          "voltage": [], "current": []}

        self._build_ui()

    # ==================================================================
    # UI construction
    # ==================================================================
    def _build_ui(self):
        main_pw = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_pw.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        # ---- LEFT: scrollable control panel ----
        left_outer = ttk.Frame(main_pw, width=330)
        main_pw.add(left_outer, weight=0)
        left_canvas = tk.Canvas(left_outer, highlightthickness=0, width=320)
        left_sb = ttk.Scrollbar(left_outer, orient=tk.VERTICAL,
                                command=left_canvas.yview)
        left = ttk.Frame(left_canvas)
        left.bind("<Configure>",
                  lambda e: left_canvas.configure(
                      scrollregion=left_canvas.bbox("all")))
        left_canvas.create_window((0, 0), window=left, anchor="nw")
        left_canvas.configure(yscrollcommand=left_sb.set)
        left_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        left_sb.pack(side=tk.RIGHT, fill=tk.Y)

        # --- General Parameters ---
        ttk.Label(left, text="General Parameters",
                  font=("Helvetica", 12, "bold")).pack(pady=(8, 8))
        self.general = {}
        for lbl, key, dflt in [
            ("Grid points (x)", "grid", "241"),
            ("Particles", "particles", "5000"),
            ("Current density (A/m\u00b2)", "current", "600"),
            ("Beam energy (eV)", "energy", "5.0"),
        ]:
            f = ttk.Frame(left); f.pack(fill=tk.X, padx=12, pady=2)
            ttk.Label(f, text=lbl, width=22, anchor="w").pack(side=tk.LEFT)
            v = tk.StringVar(value=dflt)
            ttk.Entry(f, textvariable=v, width=10).pack(side=tk.RIGHT)
            self.general[key] = v

        # --- Beam Species ---
        ttk.Separator(left).pack(fill=tk.X, padx=12, pady=8)
        ttk.Label(left, text="Beam Species",
                  font=("Helvetica", 12, "bold")).pack(pady=(0, 4))

        # Presets: name -> (mass_amu, charge_e)
        self._species_presets = {
            "e\u207b (electron)":       (0.000548579909, -1.0),
            "H\u207a (proton)":         (1.00728, 1.0),
            "H\u2082\u207a":            (2.01410, 1.0),
            "H\u2083\u207a":            (3.02140, 1.0),
            "D\u207a (deuteron)":       (2.01410, 1.0),
            "He\u207a":                 (4.00260, 1.0),
            "He\u00b2\u207a":           (4.00260, 2.0),
            "N\u207a":                  (14.003, 1.0),
            "O\u207a":                  (15.999, 1.0),
            "Ar\u207a":                 (39.948, 1.0),
            "Ar\u00b2\u207a":           (39.948, 2.0),
            "Xe\u207a":                 (131.29, 1.0),
            "Custom":                    (None, None),
        }
        sf = ttk.Frame(left); sf.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(sf, text="Species:", width=10, anchor="w").pack(side=tk.LEFT)
        self.species_var = tk.StringVar(value="H\u207a (proton)")
        species_cb = ttk.Combobox(sf, textvariable=self.species_var, width=16,
                                  values=list(self._species_presets.keys()),
                                  state="readonly")
        species_cb.pack(side=tk.LEFT)
        species_cb.bind("<<ComboboxSelected>>", self._on_species_changed)

        mf = ttk.Frame(left); mf.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(mf, text="Mass (amu):", width=14, anchor="w").pack(side=tk.LEFT)
        self.beam_mass_var = tk.StringVar(value="1.00728")
        self.mass_entry = ttk.Entry(mf, textvariable=self.beam_mass_var, width=12)
        self.mass_entry.pack(side=tk.LEFT)

        cf = ttk.Frame(left); cf.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(cf, text="Charge (e):", width=14, anchor="w").pack(side=tk.LEFT)
        self.beam_charge_var = tk.StringVar(value="1.0")
        self.charge_entry = ttk.Entry(cf, textvariable=self.beam_charge_var, width=12)
        self.charge_entry.pack(side=tk.LEFT)

        # --- Electrodes ---
        ttk.Separator(left).pack(fill=tk.X, padx=12, pady=8)
        ttk.Label(left, text="Electrodes",
                  font=("Helvetica", 12, "bold")).pack(pady=(0, 4))
        hdr = ttk.Frame(left); hdr.pack(fill=tk.X, padx=12)
        for txt, w in [("#", 3), ("Dist(mm)", 9), ("Apt(mm)", 8),
                        ("V(kV)", 7), ("Chm(\u00b0)", 7)]:
            ttk.Label(hdr, text=txt, width=w, anchor="center",
                      font=("Helvetica", 9, "bold")).pack(side=tk.LEFT, padx=1)

        self.elec_frame = ttk.Frame(left)
        self.elec_frame.pack(fill=tk.X, padx=12)
        self.electrode_rows = []
        self._add_electrode_row("0.0", "0.5", "0.0", "63.4")
        self._add_electrode_row("10.0", "1.5", "-8.0", "45.0")

        bf = ttk.Frame(left); bf.pack(fill=tk.X, padx=12, pady=4)
        ttk.Button(bf, text="+ Add", command=self._on_add_electrode).pack(
            side=tk.LEFT, padx=2)
        ttk.Button(bf, text="- Remove", command=self._on_remove_electrode).pack(
            side=tk.LEFT, padx=2)

        # --- Run Simulation ---
        ttk.Separator(left).pack(fill=tk.X, padx=12, pady=8)
        self.run_btn = ttk.Button(left, text="Run Simulation",
                                  command=self._on_run)
        self.run_btn.pack(fill=tk.X, padx=12, ipady=5)
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(left, textvariable=self.status_var,
                  foreground="gray").pack(padx=12, pady=2)
        self.progress = ttk.Progressbar(left, mode="indeterminate")
        self.progress.pack(fill=tk.X, padx=12, pady=2)

        # --- Perveance Scan ---
        ttk.Separator(left).pack(fill=tk.X, padx=12, pady=8)
        ttk.Label(left, text="Perveance Scan",
                  font=("Helvetica", 12, "bold")).pack(pady=(0, 4))

        sf = ttk.Frame(left); sf.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(sf, text="Scan type:", width=12, anchor="w").pack(side=tk.LEFT)
        self.scan_type = tk.StringVar(value="Voltage scan")
        ttk.Combobox(sf, textvariable=self.scan_type, width=14,
                     values=["Voltage scan", "Current scan"],
                     state="readonly").pack(side=tk.LEFT)

        sf2 = ttk.Frame(left); sf2.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(sf2, text="Electrode #:", width=12, anchor="w").pack(side=tk.LEFT)
        self.scan_electrode = tk.StringVar(value="2")
        ttk.Entry(sf2, textvariable=self.scan_electrode, width=5).pack(side=tk.LEFT)
        ttk.Label(sf2, text="(for V-scan)", foreground="gray").pack(side=tk.LEFT, padx=4)

        self.scan_params = {}
        self.scan_labels = {}
        for lbl, key, dflt, unit in [
            ("Min", "min", "2.0", "kV"),
            ("Max", "max", "15.0", "kV"),
            ("Steps", "steps", "8", ""),
        ]:
            f = ttk.Frame(left); f.pack(fill=tk.X, padx=12, pady=2)
            ttk.Label(f, text=lbl, width=6, anchor="w").pack(side=tk.LEFT)
            v = tk.StringVar(value=dflt)
            ttk.Entry(f, textvariable=v, width=10).pack(side=tk.LEFT)
            self.scan_params[key] = v
            if unit:
                ulbl = ttk.Label(f, text=unit, width=8, foreground="gray")
                ulbl.pack(side=tk.LEFT, padx=4)
                self.scan_labels[key] = ulbl

        # Update units when scan type changes
        self.scan_type.trace_add("write", self._on_scan_type_changed)

        sbf = ttk.Frame(left); sbf.pack(fill=tk.X, padx=12, pady=4)
        self.scan_btn = ttk.Button(sbf, text="Run Scan",
                                   command=self._on_run_scan)
        self.scan_btn.pack(side=tk.LEFT, padx=2, fill=tk.X, expand=True)
        self.scan_stop_btn = ttk.Button(sbf, text="Stop",
                                        command=self._on_stop_scan,
                                        state="disabled")
        self.scan_stop_btn.pack(side=tk.LEFT, padx=2)

        self.scan_status = tk.StringVar(value="")
        ttk.Label(left, textvariable=self.scan_status,
                  foreground="gray").pack(padx=12)

        # --- Save / Export ---
        ttk.Separator(left).pack(fill=tk.X, padx=12, pady=8)
        ttk.Label(left, text="Save / Export",
                  font=("Helvetica", 12, "bold")).pack(pady=(0, 4))
        save_f1 = ttk.Frame(left); save_f1.pack(fill=tk.X, padx=12, pady=2)
        ttk.Button(save_f1, text="Save All Plots (PNG)",
                   command=self._on_save_plots).pack(fill=tk.X)
        save_f2 = ttk.Frame(left); save_f2.pack(fill=tk.X, padx=12, pady=2)
        ttk.Button(save_f2, text="Save Scan GIF",
                   command=self._on_save_gif).pack(fill=tk.X)
        self.save_status = tk.StringVar(value="")
        ttk.Label(left, textvariable=self.save_status,
                  foreground="gray", font=("Helvetica", 8)).pack(padx=12)

        # --- Beam Info ---
        ttk.Separator(left).pack(fill=tk.X, padx=12, pady=8)
        ttk.Label(left, text="Beam Info at Slider",
                  font=("Helvetica", 11, "bold")).pack(padx=12, anchor="w")
        self.info_var = tk.StringVar(value="Run simulation first")
        ttk.Label(left, textvariable=self.info_var, justify=tk.LEFT,
                  wraplength=280, font=("Courier", 9)).pack(
            padx=12, pady=4, anchor="w")

        # ---- RIGHT PANEL ----
        right = ttk.Frame(main_pw)
        main_pw.add(right, weight=1)
        right_pw = ttk.PanedWindow(right, orient=tk.VERTICAL)
        right_pw.pack(fill=tk.BOTH, expand=True)

        # Top: trajectory
        top_f = ttk.LabelFrame(right_pw, text="Particle Trajectories (Axisymmetric)")
        right_pw.add(top_f, weight=1)
        self.traj_fig = Figure(figsize=(10, 3.2), dpi=100, facecolor="#f5f5f5")
        self.traj_ax = self.traj_fig.add_subplot(111)
        self.traj_ax.text(0.5, 0.5, "No simulation data", ha="center",
                          va="center", transform=self.traj_ax.transAxes,
                          color="gray")
        self.traj_ax.axis("off")
        self.traj_canvas = FigureCanvasTkAgg(self.traj_fig, master=top_f)
        self.traj_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Bottom row
        bot_f = ttk.Frame(right_pw)
        right_pw.add(bot_f, weight=1)

        # Slider
        sl_f = ttk.Frame(bot_f)
        sl_f.pack(fill=tk.X, padx=10, pady=(4, 0))
        ttk.Label(sl_f, text="Axial position:").pack(side=tk.LEFT)
        self.pos_var = tk.StringVar(value="--")
        ttk.Label(sl_f, textvariable=self.pos_var,
                  font=("Courier", 10)).pack(side=tk.LEFT, padx=8)
        self.slider = ttk.Scale(sl_f, from_=0, to=100, orient=tk.HORIZONTAL,
                                command=self._on_slider)
        self.slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=8)

        plots_f = ttk.Frame(bot_f)
        plots_f.pack(fill=tk.BOTH, expand=True)

        # Phase-space (left)
        self.emit_fig = Figure(figsize=(5, 3.5), dpi=100)
        self.emit_ax = self.emit_fig.add_subplot(111)
        self.emit_ax.set_title("Phase Space (y, y')")
        self._placeholder(self.emit_ax)
        self.emit_canvas = FigureCanvasTkAgg(self.emit_fig, master=plots_f)
        self.emit_canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH,
                                              expand=True)

        # Right: notebook with Envelope and Scan tabs
        self.right_nb = ttk.Notebook(plots_f)
        self.right_nb.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        env_f = ttk.Frame(self.right_nb)
        self.right_nb.add(env_f, text="Envelope & Divergence")
        self.div_fig = Figure(figsize=(5, 3.5), dpi=100)
        self.div_ax = self.div_fig.add_subplot(111)
        self._placeholder(self.div_ax, "Beam Envelope & Divergence")
        self.div_canvas = FigureCanvasTkAgg(self.div_fig, master=env_f)
        self.div_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        scan_f = ttk.Frame(self.right_nb)
        self.right_nb.add(scan_f, text="Perveance Scan")
        self.scan_fig = Figure(figsize=(5, 3.5), dpi=100)
        self.scan_ax = self.scan_fig.add_subplot(111)
        self._placeholder(self.scan_ax, "Perveance Scan")
        self.scan_canvas = FigureCanvasTkAgg(self.scan_fig, master=scan_f)
        self.scan_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    @staticmethod
    def _placeholder(ax, title=""):
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, color="gray")
        if title:
            ax.set_title(title)

    # ==================================================================
    # Species selection
    # ==================================================================
    def _on_species_changed(self, _event=None):
        name = self.species_var.get()
        preset = self._species_presets.get(name)
        if preset and preset[0] is not None:
            self.beam_mass_var.set(str(preset[0]))
            self.beam_charge_var.set(str(preset[1]))
            self.mass_entry.config(state="disabled")
            self.charge_entry.config(state="disabled")
        else:
            self.mass_entry.config(state="normal")
            self.charge_entry.config(state="normal")

    def _on_scan_type_changed(self, *_args):
        is_v = self.scan_type.get() == "Voltage scan"
        unit = "kV" if is_v else "A/m\u00b2"
        self.scan_labels["min"].config(text=unit)
        self.scan_labels["max"].config(text=unit)

    # ==================================================================
    # Electrode helpers
    # ==================================================================
    def _add_electrode_row(self, dist="0.0", apt="1.0", volt="0.0",
                           chamfer="45.0"):
        idx = len(self.electrode_rows) + 1
        row = ttk.Frame(self.elec_frame)
        row.pack(fill=tk.X, pady=1)
        ttk.Label(row, text=str(idx), width=3, anchor="center").pack(
            side=tk.LEFT, padx=1)
        vars_ = {}
        for key, dflt, w in [("dist", dist, 9), ("apt", apt, 8),
                              ("volt", volt, 7), ("chamfer", chamfer, 7)]:
            v = tk.StringVar(value=dflt)
            ttk.Entry(row, textvariable=v, width=w).pack(side=tk.LEFT, padx=1)
            vars_[key] = v
        self.electrode_rows.append({"frame": row, "vars": vars_})

    def _on_add_electrode(self):
        if len(self.electrode_rows) >= 10:
            return
        self._add_electrode_row()

    def _on_remove_electrode(self):
        if len(self.electrode_rows) <= 1:
            return
        self.electrode_rows.pop()["frame"].destroy()

    # ==================================================================
    # Config + validation
    # ==================================================================
    def _write_config(self, voltage_override=None, current_override=None):
        """Write config file. Optional overrides for scan."""
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        cur = current_override or self.general["current"].get()
        with open(CONFIG_FILE, "w") as f:
            f.write(f"grid_points={self.general['grid'].get()}\n")
            f.write(f"particles={self.general['particles'].get()}\n")
            f.write(f"current_density={cur}\n")
            f.write(f"beam_energy={self.general['energy'].get()}\n")
            f.write(f"beam_mass={self.beam_mass_var.get()}\n")
            f.write(f"beam_charge={self.beam_charge_var.get()}\n")
            n = len(self.electrode_rows)
            f.write(f"electrodes={n}\n")
            for i, r in enumerate(self.electrode_rows):
                v = r["vars"]
                volt = v["volt"].get()
                if voltage_override and (i + 1) == voltage_override[0]:
                    volt = str(voltage_override[1])
                f.write(f"{v['dist'].get()} {v['apt'].get()} "
                        f"{volt} {v['chamfer'].get()}\n")

    def _validate_inputs(self):
        try:
            g = int(self.general["grid"].get())
            if g < 20: return "Grid points must be >= 20"
        except ValueError:
            return "Grid points must be an integer"
        try:
            p = int(self.general["particles"].get())
            if p < 100: return "Particles must be >= 100"
        except ValueError:
            return "Particles must be an integer"
        try:
            c = float(self.general["current"].get())
            if c <= 0: return "Current density must be > 0"
        except ValueError:
            return "Current density must be a number"

        positions = []
        for i, row in enumerate(self.electrode_rows):
            v = row["vars"]
            try:
                dist = float(v["dist"].get())
                apt = float(v["apt"].get())
                chamfer = float(v["chamfer"].get())
            except ValueError:
                return f"Electrode #{i+1}: fields must be numbers"
            if apt <= 0:
                return f"Electrode #{i+1}: aperture must be > 0"
            if chamfer < 5:
                return f"Electrode #{i+1}: chamfer must be >= 5 deg"
            thickness = max(2.0, 4.0 * apt)
            positions.append((dist, dist + thickness, i + 1))
        positions.sort()
        for j in range(len(positions) - 1):
            _, e1, i1 = positions[j]
            s2, _, i2 = positions[j + 1]
            if e1 > s2:
                return (f"Electrodes #{i1} and #{i2} overlap! "
                        f"#{i1} ends at {e1:.1f} mm, #{i2} starts at {s2:.1f} mm")
        return None

    # ==================================================================
    # Run single simulation
    # ==================================================================
    def _on_run(self):
        err = self._validate_inputs()
        if err:
            messagebox.showwarning("Invalid Input", err)
            return
        self.run_btn.config(state="disabled")
        self.scan_btn.config(state="disabled")
        self.status_var.set("Running simulation...")
        self.progress.start(15)
        threading.Thread(target=self._run_single, daemon=True).start()

    def _run_single(self):
        """Thread target for single simulation run."""
        ok = self._run_sim()
        if ok:
            self.root.after(0, self._sim_done)

    def _run_sim(self, voltage_override=None, current_override=None):
        """Run one simulation. Returns True on success."""
        try:
            self._write_config(voltage_override, current_override)
            cmd = [SIM_BINARY, "--config", CONFIG_FILE, "--outdir", OUTPUT_DIR]
            env = os.environ.copy()
            # When bundled, shared libs for beam_sim are in simbin/
            if getattr(sys, 'frozen', False):
                simbin_dir = os.path.dirname(SIM_BINARY)
                env["LD_LIBRARY_PATH"] = simbin_dir + ":" + env.get("LD_LIBRARY_PATH", "")
            proc = subprocess.run(cmd, capture_output=True, text=True,
                                  timeout=300, env=env)
            if proc.returncode != 0:
                self.root.after(0, lambda: self._sim_error(proc.stderr))
                return False
            self._load_results()
            return True
        except subprocess.TimeoutExpired:
            self.root.after(0, lambda: self._sim_error("Timed out (5 min)"))
            return False
        except Exception as e:
            self.root.after(0, lambda: self._sim_error(str(e)))
            return False

    def _sim_error(self, msg):
        self.progress.stop()
        self.run_btn.config(state="normal")
        self.scan_btn.config(state="normal")
        self.status_var.set("Error")
        messagebox.showerror("Simulation Error", msg)

    def _sim_done(self):
        self.progress.stop()
        self.run_btn.config(state="normal")
        self.scan_btn.config(state="normal")
        self.status_var.set("Done")
        self._update_trajectory_plot()
        self._update_divergence_plot()
        self._setup_slider()
        self._on_slider(None)

    # ==================================================================
    # Perveance scan
    # ==================================================================
    def _on_run_scan(self):
        err = self._validate_inputs()
        if err:
            messagebox.showwarning("Invalid Input", err)
            return
        try:
            smin = float(self.scan_params["min"].get())
            smax = float(self.scan_params["max"].get())
            steps = int(self.scan_params["steps"].get())
            if steps < 2:
                raise ValueError
        except ValueError:
            messagebox.showwarning("Invalid Scan", "Check min/max/steps values")
            return

        is_vscan = self.scan_type.get() == "Voltage scan"
        if is_vscan:
            try:
                elec_idx = int(self.scan_electrode.get())
                if elec_idx < 1 or elec_idx > len(self.electrode_rows):
                    raise ValueError
            except ValueError:
                messagebox.showwarning("Invalid",
                                       f"Electrode # must be 1-{len(self.electrode_rows)}")
                return

        self.scan_running = True
        self.scan_stop = False
        self.scan_is_vscan = is_vscan
        self.scan_frames = []  # PIL Images for GIF
        self.scan_data = {"perveance": [], "divergence": [],
                          "voltage": [], "current": [],
                          "scan_val": [], "scan_unit": "kV" if is_vscan else "A/m\u00b2"}
        self.run_btn.config(state="disabled")
        self.scan_btn.config(state="disabled")
        self.scan_stop_btn.config(state="normal")

        # Switch to scan tab
        self.right_nb.select(1)

        # Clear scan plot
        self.scan_ax.clear()
        self.scan_ax.set_xlabel("Perveance (\u00b5A / V\u00b3\u02f2\u00b2)")
        self.scan_ax.set_ylabel("RMS Divergence (mrad)")
        self.scan_ax.set_title("Perveance Scan")
        self.scan_ax.grid(True, alpha=0.3)
        self.scan_fig.tight_layout()
        self.scan_canvas.draw()

        threading.Thread(target=self._scan_thread,
                         args=(smin, smax, steps, is_vscan),
                         daemon=True).start()

    def _on_stop_scan(self):
        self.scan_stop = True

    def _scan_thread(self, smin, smax, steps, is_vscan):
        values = np.linspace(smin, smax, steps)
        elec_idx = int(self.scan_electrode.get()) if is_vscan else 0

        for i, val in enumerate(values):
            if self.scan_stop:
                break

            self.root.after(0, lambda i=i, s=steps, v=val:
                            self.scan_status.set(
                                f"Scan {i+1}/{s}  "
                                f"({'V' if is_vscan else 'J'}={v:.2f})..."))

            if is_vscan:
                # Use sign from the electrode: negative for ions, positive for electrons
                try:
                    orig_v = float(self.electrode_rows[elec_idx - 1]["vars"]["volt"].get())
                    sign = -1.0 if orig_v <= 0 else 1.0
                except (IndexError, ValueError):
                    sign = -1.0
                v_kv = sign * abs(val)
                ok = self._run_sim(voltage_override=(elec_idx, v_kv))
                v_volts = abs(val) * 1e3
                scan_label = f"{abs(val):.1f} kV"
            else:
                ok = self._run_sim(current_override=str(val))
                try:
                    v_kv = float(self.electrode_rows[1]["vars"]["volt"].get())
                except (IndexError, ValueError):
                    v_kv = -8.0
                v_volts = abs(v_kv) * 1e3
                scan_label = f"{val:.0f} A/m\u00b2"

            if not ok:
                break

            # Extract divergence at last sampling point
            if self.emittance_data and len(self.emittance_data["divergence_mrad"]) > 0:
                div = self.emittance_data["divergence_mrad"][-1]
                curr = self.emittance_data["current_A"][-1]
            else:
                continue

            # Perveance: P = I / V^(3/2), in microperv
            if v_volts > 0:
                perv = (curr / (v_volts ** 1.5)) * 1e6
            else:
                continue

            self.scan_data["perveance"].append(perv)
            self.scan_data["divergence"].append(div)
            self.scan_data["voltage"].append(v_volts)
            self.scan_data["current"].append(curr)
            self.scan_data["scan_val"].append(scan_label)

            # Update ALL plots in real-time after each scan point
            self.root.after(0, self._update_all_scan_plots)

        self.root.after(0, self._scan_finished)

    def _update_all_scan_plots(self):
        """Refresh trajectory, phase space, envelope, and scan plot."""
        self._update_trajectory_plot()
        self._update_divergence_plot()
        self._setup_slider()
        self._on_slider(None)
        self._redraw_scan_plot()
        self._capture_scan_frame()

    def _redraw_scan_plot(self):
        self.scan_ax.clear()
        p = self.scan_data["perveance"]
        d = self.scan_data["divergence"]
        labels = self.scan_data["scan_val"]

        self.scan_ax.plot(p, d, "o-", color="#1565C0", markersize=8, lw=2)
        for i in range(len(p)):
            self.scan_ax.annotate(labels[i], (p[i], d[i]),
                                  textcoords="offset points", xytext=(6, 6),
                                  fontsize=7, color="gray")

        self.scan_ax.set_xlabel("Perveance (\u00b5A / V\u00b3\u02f2\u00b2)")
        self.scan_ax.set_ylabel("RMS Divergence (mrad)")
        scan_type = "Voltage" if self.scan_is_vscan else "Current density"
        self.scan_ax.set_title(f"Perveance vs Divergence ({scan_type} scan)")
        self.scan_ax.grid(True, alpha=0.3)

        if len(d) > 1:
            imin = int(np.argmin(d))
            self.scan_ax.plot(p[imin], d[imin], "*", color="red",
                              markersize=15, zorder=5,
                              label=f"Min div: {d[imin]:.2f} mrad @ {labels[imin]}")
            self.scan_ax.legend(fontsize=9)

        self.scan_fig.tight_layout()
        self.scan_canvas.draw()

    def _scan_finished(self):
        self.scan_running = False
        self.run_btn.config(state="normal")
        self.scan_btn.config(state="normal")
        self.scan_stop_btn.config(state="disabled")
        if self.scan_stop:
            self.scan_status.set("Scan stopped")
        else:
            self.scan_status.set(f"Scan done ({len(self.scan_data['perveance'])} pts)")
        # Also update the single-sim plots with the last run
        self.root.after(0, self._sim_done)

    # ==================================================================
    # Frame capture for GIF
    # ==================================================================
    def _capture_scan_frame(self):
        """Capture a composite frame of trajectory + phase space + scan plot."""
        try:
            # Render each figure to a PIL image
            imgs = []
            for fig in [self.traj_fig, self.emit_fig, self.scan_fig]:
                fig.canvas.draw()
                buf = fig.canvas.buffer_rgba()
                w, h = fig.canvas.get_width_height()
                img = Image.frombuffer("RGBA", (w, h), buf).copy()
                imgs.append(img)

            # Stack: trajectory on top, phase space + scan side by side below
            traj_img = imgs[0]
            emit_img = imgs[1]
            scan_img = imgs[2]

            # Resize bottom row to same width as top
            top_w = traj_img.width
            bot_h = max(emit_img.height, scan_img.height)
            half_w = top_w // 2
            emit_r = emit_img.resize((half_w, bot_h), Image.LANCZOS)
            scan_r = scan_img.resize((top_w - half_w, bot_h), Image.LANCZOS)

            composite = Image.new("RGBA",
                                  (top_w, traj_img.height + bot_h),
                                  (255, 255, 255, 255))
            composite.paste(traj_img, (0, 0))
            composite.paste(emit_r, (0, traj_img.height))
            composite.paste(scan_r, (half_w, traj_img.height))

            self.scan_frames.append(composite.convert("RGB"))
        except Exception:
            pass  # skip frame on error

    # ==================================================================
    # Save plots and GIF
    # ==================================================================
    def _on_save_plots(self):
        from tkinter import filedialog
        folder = filedialog.askdirectory(title="Select folder to save plots")
        if not folder:
            return
        try:
            self.traj_fig.savefig(os.path.join(folder, "trajectory.png"),
                                  dpi=150, bbox_inches="tight")
            self.emit_fig.savefig(os.path.join(folder, "phase_space.png"),
                                  dpi=150, bbox_inches="tight")
            self.div_fig.savefig(os.path.join(folder, "envelope_divergence.png"),
                                 dpi=150, bbox_inches="tight")
            self.scan_fig.savefig(os.path.join(folder, "perveance_scan.png"),
                                   dpi=150, bbox_inches="tight")
            self.save_status.set(f"Saved 4 plots to {folder}")
        except Exception as e:
            messagebox.showerror("Save Error", str(e))

    def _on_save_gif(self):
        if not self.scan_frames:
            messagebox.showinfo("No frames",
                                "Run a perveance scan first to generate frames.")
            return
        from tkinter import filedialog
        path = filedialog.asksaveasfilename(
            title="Save scan animation",
            defaultextension=".gif",
            filetypes=[("GIF", "*.gif")])
        if not path:
            return
        try:
            self.save_status.set("Saving GIF...")
            self.root.update_idletasks()
            # duration per frame in ms; last frame lingers longer
            durations = [800] * len(self.scan_frames)
            durations[-1] = 2000
            self.scan_frames[0].save(
                path, save_all=True,
                append_images=self.scan_frames[1:],
                duration=durations, loop=0, optimize=True)
            self.save_status.set(f"GIF saved: {os.path.basename(path)}")
        except Exception as e:
            messagebox.showerror("GIF Error", str(e))

    # ==================================================================
    # Load results
    # ==================================================================
    def _load_results(self):
        self.emittance_data = self._read_csv(
            os.path.join(OUTPUT_DIR, "emittance_profile.csv"))

        raw = self._read_csv(os.path.join(OUTPUT_DIR, "phase_space.csv"))
        self.phase_space = {}
        if raw and len(raw.get("x_mm", [])) > 0:
            x_v = np.array(raw["x_mm"])
            y_v = np.array(raw["y_mm"])
            yp_v = np.array(raw["yp_mrad"])
            for xv in np.unique(x_v):
                m = x_v == xv
                yh, yph = y_v[m], yp_v[m]
                self.phase_space[xv] = (np.concatenate([yh, -yh]),
                                        np.concatenate([yph, -yph]))

        tfile = os.path.join(OUTPUT_DIR, "trajectory.png")
        self.traj_image = Image.open(tfile) if os.path.exists(tfile) else None

    def _read_csv(self, path):
        if not os.path.exists(path):
            return {}
        data = {}
        with open(path) as f:
            header = f.readline().strip().split(",")
            for h in header:
                data[h] = []
            for line in f:
                parts = line.strip().split(",")
                if len(parts) != len(header):
                    continue
                for h, v in zip(header, parts):
                    try:
                        data[h].append(float(v))
                    except ValueError:
                        data[h].append(0.0)
        return data

    # ==================================================================
    # Plot updates
    # ==================================================================
    def _update_trajectory_plot(self):
        self.traj_ax.clear()
        if self.traj_image is not None:
            self.traj_ax.imshow(np.array(self.traj_image), aspect="auto")
        self.traj_ax.axis("off")
        self.traj_fig.tight_layout(pad=0.5)
        self.traj_canvas.draw()

    def _update_divergence_plot(self):
        # Fully recreate axes to avoid stale twin-axis issues
        self.div_fig.clear()
        self.div_ax = self.div_fig.add_subplot(111)

        if self.emittance_data is None or not self.emittance_data.get("x_mm"):
            self._placeholder(self.div_ax, "Beam Envelope & Divergence")
            self.div_canvas.draw()
            return

        x = np.array(self.emittance_data["x_mm"])
        r_rms = np.array(self.emittance_data["r_rms_mm"])
        div = np.array(self.emittance_data["divergence_mrad"])

        c1, c2 = "#2196F3", "#FF5722"
        self.div_ax.set_xlabel("x (mm)")
        self.div_ax.set_ylabel("RMS beam size (mm)", color=c1)
        ln1 = self.div_ax.plot(x, r_rms, color=c1, lw=2, label="r_rms")
        self.div_ax.tick_params(axis="y", labelcolor=c1)

        self._div_ax2 = self.div_ax.twinx()
        self._div_ax2.set_ylabel("Divergence (mrad)", color=c2)
        ln2 = self._div_ax2.plot(x, div, color=c2, lw=2, ls="--",
                                  label="Divergence")
        self._div_ax2.tick_params(axis="y", labelcolor=c2)

        lns = ln1 + ln2
        self.div_ax.legend(lns, [l.get_label() for l in lns],
                           loc="upper right", fontsize=9)
        self.div_ax.set_title("Beam Envelope & Divergence vs Position")
        self.div_ax.grid(True, alpha=0.3)
        self._div_vline = None
        self.div_fig.tight_layout()
        self.div_canvas.draw()

    def _setup_slider(self):
        if not self.emittance_data or not self.emittance_data.get("x_mm"):
            return
        n = len(self.emittance_data["x_mm"])
        self.slider.config(from_=0, to=n - 1)
        self.slider.set(n // 2)

    # ==================================================================
    # Slider
    # ==================================================================
    def _on_slider(self, _val):
        if not self.emittance_data or not self.phase_space:
            return

        x_list = self.emittance_data["x_mm"]
        if not x_list:
            return
        idx = max(0, min(int(float(self.slider.get())), len(x_list) - 1))
        x_mm = x_list[idx]

        ps_keys = sorted(self.phase_space.keys())
        if not ps_keys:
            return
        closest_x = min(ps_keys, key=lambda k: abs(k - x_mm))
        y_data, yp_data = self.phase_space[closest_x]

        eps_mm, alpha, beta_mm, _ = compute_twiss(y_data, yp_data)

        r_rms = self.emittance_data["r_rms_mm"][idx]
        div_mrad = self.emittance_data["divergence_mrad"][idx]
        curr = self.emittance_data["current_A"][idx]

        self.pos_var.set(f"x = {x_mm:.2f} mm")
        species = self.species_var.get()
        self.info_var.set(
            f"Species  : {species}\n"
            f"x = {x_mm:.3f} mm\n"
            f"eps(y,y')= {eps_mm:.4f} mm-mrad\n"
            f"alpha    = {alpha:.3f}\n"
            f"beta     = {beta_mm:.4f} mm/mrad\n"
            f"r_rms    = {r_rms:.4f} mm\n"
            f"div(rms) = {div_mrad:.3f} mrad\n"
            f"current  = {curr:.4e} A"
        )

        # Phase-space scatter
        self.emit_ax.clear()
        self.emit_ax.scatter(y_data, yp_data, s=1.5, alpha=0.4, c="#1565C0",
                             edgecolors="none")
        self.emit_ax.set_xlabel("y (mm)")
        self.emit_ax.set_ylabel("y' (mrad)")
        self.emit_ax.set_title(
            f"Phase Space x={x_mm:.2f}mm  "
            f"\u03b5={eps_mm:.4f} mm-mrad", fontsize=10)
        self.emit_ax.axhline(0, color="gray", lw=0.5)
        self.emit_ax.axvline(0, color="gray", lw=0.5)
        self.emit_ax.grid(True, alpha=0.3)
        if eps_mm > 0 and beta_mm > 0:
            self._draw_rms_ellipse(self.emit_ax, eps_mm, alpha, beta_mm)
        self.emit_fig.tight_layout()
        self.emit_canvas.draw()

        # Vertical marker on divergence plot
        if self._div_vline is not None:
            try:
                self._div_vline.remove()
            except Exception:
                pass
            self._div_vline = None
        self._div_vline = self.div_ax.axvline(x_mm, color="red", lw=1.5,
                                               ls=":", alpha=0.7)
        self.div_canvas.draw()

    @staticmethod
    def _draw_rms_ellipse(ax, eps, alpha, beta):
        t = np.linspace(0, 2 * np.pi, 200)
        sb = np.sqrt(beta)
        y_e = np.sqrt(eps) * sb * np.cos(t)
        yp_e = np.sqrt(eps) * (-alpha / sb * np.cos(t)
                                + 1.0 / sb * np.sin(t))
        ax.plot(y_e, yp_e, "r-", lw=1.5, alpha=0.8, label="RMS ellipse")
        ax.legend(fontsize=8, loc="upper left")


def main():
    if not os.path.exists(SIM_BINARY):
        print(f"Error: {SIM_BINARY} not found. Build: cd gui && make")
        sys.exit(1)
    root = tk.Tk()
    BeamExtractionGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
