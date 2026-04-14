#!/usr/bin/env python3
"""
IBSIMU Beam Extraction GUI - Full-featured version
Solver selection, beam modes, B-field, field diagnostics, convergence, electrode currents.
"""

import os, sys, subprocess, threading
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from PIL import Image

def _base():
    return sys._MEIPASS if getattr(sys,'frozen',False) else os.path.dirname(os.path.abspath(__file__))
def _outdir():
    return os.path.join(os.path.dirname(sys.executable),"sim_output") if getattr(sys,'frozen',False) \
        else os.path.join(os.path.dirname(os.path.abspath(__file__)),"sim_output")

BASE_DIR = _base()
SIM_BINARY = os.path.join(BASE_DIR,"simbin","beam_sim") if getattr(sys,'frozen',False) \
    else os.path.join(os.path.dirname(os.path.abspath(__file__)),"beam_sim")
OUTPUT_DIR = _outdir()
CONFIG_FILE = os.path.join(OUTPUT_DIR,"beam_config.txt")

def compute_twiss(y, yp):
    if len(y) < 3: return 0,0,0,0
    yc, ypc = y-np.mean(y), yp-np.mean(yp)
    syy, spp, syp = np.mean(yc**2), np.mean(ypc**2), np.mean(yc*ypc)
    det = syy*spp - syp**2
    if det <= 0: return 0,0,0,0
    eps = np.sqrt(det)
    return eps, -syp/eps, syy/eps, spp/eps


class BeamGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("IBSIMU - Beam Extraction Simulator")
        self.root.geometry("1550x950")
        self.root.minsize(1200, 700)

        self.emittance_data = None
        self.phase_space = None
        self.traj_image = None
        self.field_data = None
        self.convergence_data = None
        self.electrode_currents = None
        self.energy_data = None

        self.scan_running = False
        self.scan_stop = False
        self.scan_is_vscan = True
        self.scan_frames = []
        self.scan_data = {"perveance":[],"divergence":[],"voltage":[],"current":[],
                          "grid_ratio":[],"scan_val":[]}

        self._species_presets = {
            "e\u207b (electron)":(0.000548579909,-1.0),
            "H\u207a (proton)":(1.00728,1.0), "H\u2082\u207a":(2.014,1.0),
            "D\u207a (deuteron)":(2.014,1.0), "He\u207a":(4.003,1.0),
            "He\u00b2\u207a":(4.003,2.0), "N\u207a":(14.003,1.0),
            "O\u207a":(15.999,1.0), "Ar\u207a":(39.948,1.0),
            "Ar\u00b2\u207a":(39.948,2.0), "Xe\u207a":(131.29,1.0),
            "Custom":(None,None),
        }
        self._build_ui()

    # ==================================================================
    # UI
    # ==================================================================
    def _build_ui(self):
        pw = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        pw.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

        # ---- LEFT: scrollable ----
        lo = ttk.Frame(pw, width=340); pw.add(lo, weight=0)
        lc = tk.Canvas(lo, highlightthickness=0, width=330)
        lsb = ttk.Scrollbar(lo, orient=tk.VERTICAL, command=lc.yview)
        left = ttk.Frame(lc)
        left.bind("<Configure>", lambda e: lc.configure(scrollregion=lc.bbox("all")))
        lc.create_window((0,0), window=left, anchor="nw")
        lc.configure(yscrollcommand=lsb.set)
        lc.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        lsb.pack(side=tk.RIGHT, fill=tk.Y)

        # --- General ---
        self._section(left, "General Parameters")
        self.general = {}
        for l,k,d in [("Grid points (x)","grid","241"),("Particles","particles","5000"),
                       ("Current density (A/m\u00b2)","current","600"),("Beam energy (eV)","energy","5.0")]:
            self._param_row(left, l, k, d, self.general)

        # --- Solver ---
        sf = ttk.Frame(left); sf.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(sf, text="Solver:", width=16, anchor="w").pack(side=tk.LEFT)
        self.solver_var = tk.StringVar(value="BiCGSTAB")
        ttk.Combobox(sf, textvariable=self.solver_var, width=12,
                     values=["BiCGSTAB","Multigrid","Gauss-Seidel"],
                     state="readonly").pack(side=tk.LEFT)

        mgf = ttk.Frame(left); mgf.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(mgf, text="MG levels:", width=16, anchor="w").pack(side=tk.LEFT)
        self.mg_levels_var = tk.StringVar(value="4")
        ttk.Entry(mgf, textvariable=self.mg_levels_var, width=5).pack(side=tk.LEFT)
        ttk.Label(mgf, text="(multigrid only)", foreground="gray").pack(side=tk.LEFT, padx=4)

        # --- Species ---
        self._section(left, "Beam Species")
        sf = ttk.Frame(left); sf.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(sf, text="Species:", width=10, anchor="w").pack(side=tk.LEFT)
        self.species_var = tk.StringVar(value="H\u207a (proton)")
        cb = ttk.Combobox(sf, textvariable=self.species_var, width=16,
                          values=list(self._species_presets.keys()), state="readonly")
        cb.pack(side=tk.LEFT); cb.bind("<<ComboboxSelected>>", self._on_species)

        self.beam_mass_var = tk.StringVar(value="1.00728")
        self.beam_charge_var = tk.StringVar(value="1.0")
        mf = ttk.Frame(left); mf.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(mf, text="Mass (amu):", width=14, anchor="w").pack(side=tk.LEFT)
        self.mass_entry = ttk.Entry(mf, textvariable=self.beam_mass_var, width=12)
        self.mass_entry.pack(side=tk.LEFT)
        cf = ttk.Frame(left); cf.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(cf, text="Charge (e):", width=14, anchor="w").pack(side=tk.LEFT)
        self.charge_entry = ttk.Entry(cf, textvariable=self.beam_charge_var, width=12)
        self.charge_entry.pack(side=tk.LEFT)

        # --- Beam mode ---
        self._section(left, "Beam Definition")
        bmf = ttk.Frame(left); bmf.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(bmf, text="Mode:", width=10, anchor="w").pack(side=tk.LEFT)
        self.beam_mode_var = tk.StringVar(value="Energy/Temp")
        bm_cb = ttk.Combobox(bmf, textvariable=self.beam_mode_var, width=14,
                             values=["Energy/Temp","Twiss (KV)"], state="readonly")
        bm_cb.pack(side=tk.LEFT)
        bm_cb.bind("<<ComboboxSelected>>", self._on_beam_mode)

        self.twiss_frame = ttk.Frame(left)
        self.twiss_frame.pack(fill=tk.X, padx=12)
        self.twiss_vars = {}
        for l,k,d in [("Alpha:","tw_alpha","0.0"),("Beta (m/rad):","tw_beta","0.01"),
                       ("Emittance (m-rad):","tw_emit","1e-6")]:
            f = ttk.Frame(self.twiss_frame); f.pack(fill=tk.X, pady=1)
            ttk.Label(f, text=l, width=18, anchor="w").pack(side=tk.LEFT)
            v = tk.StringVar(value=d)
            ttk.Entry(f, textvariable=v, width=12).pack(side=tk.LEFT)
            self.twiss_vars[k] = v
        self.twiss_frame.pack_forget()  # hidden by default

        # --- Magnetic field ---
        self._section(left, "Magnetic Field")
        bff = ttk.Frame(left); bff.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(bff, text="Type:", width=10, anchor="w").pack(side=tk.LEFT)
        self.bfield_var = tk.StringVar(value="None")
        bf_cb = ttk.Combobox(bff, textvariable=self.bfield_var, width=12,
                             values=["None","Solenoid"], state="readonly")
        bf_cb.pack(side=tk.LEFT)
        bf_cb.bind("<<ComboboxSelected>>", self._on_bfield_mode)

        self.sol_frame = ttk.Frame(left)
        self.sol_vars = {}
        for l,k,d in [("B0 (T):","sol_B0","0.1"),("Z start (mm):","sol_z1","5.0"),
                       ("Z end (mm):","sol_z2","15.0")]:
            f = ttk.Frame(self.sol_frame); f.pack(fill=tk.X, pady=1)
            ttk.Label(f, text=l, width=14, anchor="w").pack(side=tk.LEFT)
            v = tk.StringVar(value=d)
            ttk.Entry(f, textvariable=v, width=10).pack(side=tk.LEFT)
            self.sol_vars[k] = v
        # hidden by default

        # --- Electrodes ---
        self._section(left, "Electrodes")
        hdr = ttk.Frame(left); hdr.pack(fill=tk.X, padx=12)
        for t,w in [("#",3),("Dist",8),("Apt",7),("V(kV)",7),("Chm\u00b0",6)]:
            ttk.Label(hdr, text=t, width=w, anchor="center",
                      font=("Helvetica",9,"bold")).pack(side=tk.LEFT, padx=1)
        self.elec_frame = ttk.Frame(left); self.elec_frame.pack(fill=tk.X, padx=12)
        self.electrode_rows = []
        self._add_elec("0.0","0.5","0.0","63.4")
        self._add_elec("10.0","1.5","-8.0","45.0")
        bf = ttk.Frame(left); bf.pack(fill=tk.X, padx=12, pady=4)
        ttk.Button(bf, text="+ Add", command=lambda: self._add_elec()).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="- Remove", command=self._rm_elec).pack(side=tk.LEFT, padx=2)

        # --- Run ---
        self._section(left, "")
        self.run_btn = ttk.Button(left, text="Run Simulation", command=self._on_run)
        self.run_btn.pack(fill=tk.X, padx=12, ipady=5)
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(left, textvariable=self.status_var, foreground="gray").pack(padx=12, pady=2)
        self.progress = ttk.Progressbar(left, mode="indeterminate")
        self.progress.pack(fill=tk.X, padx=12, pady=2)

        # --- Perveance Scan ---
        self._section(left, "Perveance Scan")
        sf = ttk.Frame(left); sf.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(sf, text="Type:", width=8, anchor="w").pack(side=tk.LEFT)
        self.scan_type = tk.StringVar(value="Voltage scan")
        ttk.Combobox(sf, textvariable=self.scan_type, width=14,
                     values=["Voltage scan","Current scan"], state="readonly").pack(side=tk.LEFT)
        sf2 = ttk.Frame(left); sf2.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(sf2, text="Elec #:", width=8, anchor="w").pack(side=tk.LEFT)
        self.scan_electrode = tk.StringVar(value="2")
        ttk.Entry(sf2, textvariable=self.scan_electrode, width=5).pack(side=tk.LEFT)

        self.scan_params = {}; self.scan_labels = {}
        for l,k,d,u in [("Min","min","-2.0","kV"),("Max","max","-15.0","kV"),("Steps","steps","8","")]:
            f = ttk.Frame(left); f.pack(fill=tk.X, padx=12, pady=2)
            ttk.Label(f, text=l, width=6, anchor="w").pack(side=tk.LEFT)
            v = tk.StringVar(value=d); ttk.Entry(f, textvariable=v, width=10).pack(side=tk.LEFT)
            self.scan_params[k] = v
            if u:
                ul = ttk.Label(f, text=u, width=8, foreground="gray"); ul.pack(side=tk.LEFT, padx=4)
                self.scan_labels[k] = ul
        self.scan_type.trace_add("write", self._on_scan_type)

        sbf = ttk.Frame(left); sbf.pack(fill=tk.X, padx=12, pady=4)
        self.scan_btn = ttk.Button(sbf, text="Run Scan", command=self._on_run_scan)
        self.scan_btn.pack(side=tk.LEFT, padx=2, fill=tk.X, expand=True)
        self.scan_stop_btn = ttk.Button(sbf, text="Stop", command=lambda: setattr(self,'scan_stop',True), state="disabled")
        self.scan_stop_btn.pack(side=tk.LEFT, padx=2)
        self.scan_status = tk.StringVar(value="")
        ttk.Label(left, textvariable=self.scan_status, foreground="gray").pack(padx=12)

        # --- Save ---
        self._section(left, "Save / Export")
        ttk.Button(left, text="Save All Plots (PNG)", command=self._save_plots).pack(fill=tk.X, padx=12, pady=2)
        ttk.Button(left, text="Save Scan GIF", command=self._save_gif).pack(fill=tk.X, padx=12, pady=2)
        self.save_status = tk.StringVar(value="")
        ttk.Label(left, textvariable=self.save_status, foreground="gray", font=("Helvetica",8)).pack(padx=12)

        # --- Beam Info ---
        self._section(left, "Beam Info at Slider")
        self.info_var = tk.StringVar(value="Run simulation first")
        ttk.Label(left, textvariable=self.info_var, justify=tk.LEFT, wraplength=300,
                  font=("Courier",9)).pack(padx=12, pady=4, anchor="w")

        # ---- RIGHT PANEL ----
        right = ttk.Frame(pw); pw.add(right, weight=1)
        rpw = ttk.PanedWindow(right, orient=tk.VERTICAL)
        rpw.pack(fill=tk.BOTH, expand=True)

        # Top: trajectory
        tf = ttk.LabelFrame(rpw, text="Particle Trajectories (Axisymmetric)")
        rpw.add(tf, weight=1)
        self.traj_fig = Figure(figsize=(10,3), dpi=100, facecolor="#f5f5f5")
        self.traj_ax = self.traj_fig.add_subplot(111); self.traj_ax.axis("off")
        self.traj_canvas = FigureCanvasTkAgg(self.traj_fig, master=tf)
        self.traj_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Bottom
        botf = ttk.Frame(rpw); rpw.add(botf, weight=1)
        # Slider
        slf = ttk.Frame(botf); slf.pack(fill=tk.X, padx=10, pady=(4,0))
        ttk.Label(slf, text="Axial position:").pack(side=tk.LEFT)
        self.pos_var = tk.StringVar(value="--")
        ttk.Label(slf, textvariable=self.pos_var, font=("Courier",10)).pack(side=tk.LEFT, padx=8)
        self.slider = ttk.Scale(slf, from_=0, to=100, orient=tk.HORIZONTAL, command=self._on_slider)
        self.slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=8)

        pf = ttk.Frame(botf); pf.pack(fill=tk.BOTH, expand=True)

        # Phase space (left)
        self.emit_fig = Figure(figsize=(5,3.3), dpi=100)
        self.emit_ax = self.emit_fig.add_subplot(111)
        self.emit_canvas = FigureCanvasTkAgg(self.emit_fig, master=pf)
        self.emit_canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Right: notebook with tabs
        self.right_nb = ttk.Notebook(pf)
        self.right_nb.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        # Tab: Envelope
        ef = ttk.Frame(self.right_nb); self.right_nb.add(ef, text="Envelope")
        self.div_fig = Figure(figsize=(5,3.3), dpi=100)
        self.div_ax = self.div_fig.add_subplot(111)
        self.div_canvas = FigureCanvasTkAgg(self.div_fig, master=ef)
        self.div_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Tab: Field Diagnostics
        fdf = ttk.Frame(self.right_nb); self.right_nb.add(fdf, text="Field Diag")
        self.field_fig = Figure(figsize=(5,3.3), dpi=100)
        self.field_ax = self.field_fig.add_subplot(111)
        self.field_canvas = FigureCanvasTkAgg(self.field_fig, master=fdf)
        self.field_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Tab: Convergence
        cvf = ttk.Frame(self.right_nb); self.right_nb.add(cvf, text="Convergence")
        self.conv_fig = Figure(figsize=(5,3.3), dpi=100)
        self.conv_ax = self.conv_fig.add_subplot(111)
        self.conv_canvas = FigureCanvasTkAgg(self.conv_fig, master=cvf)
        self.conv_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Tab: Perveance Scan
        scf = ttk.Frame(self.right_nb); self.right_nb.add(scf, text="Perveance Scan")
        self.scan_fig = Figure(figsize=(5,3.3), dpi=100)
        self.scan_ax = self.scan_fig.add_subplot(111)
        self.scan_canvas = FigureCanvasTkAgg(self.scan_fig, master=scf)
        self.scan_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Tab: Energy
        ekf = ttk.Frame(self.right_nb); self.right_nb.add(ekf, text="Energy Dist.")
        self.ek_fig = Figure(figsize=(5,3.3), dpi=100)
        self.ek_ax = self.ek_fig.add_subplot(111)
        self.ek_canvas = FigureCanvasTkAgg(self.ek_fig, master=ekf)
        self.ek_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self._div_vline = None

    # -- UI helpers --
    def _section(self, parent, title):
        ttk.Separator(parent).pack(fill=tk.X, padx=12, pady=6)
        if title:
            ttk.Label(parent, text=title, font=("Helvetica",11,"bold")).pack(pady=(0,4))

    def _param_row(self, parent, label, key, default, store):
        f = ttk.Frame(parent); f.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(f, text=label, width=22, anchor="w").pack(side=tk.LEFT)
        v = tk.StringVar(value=default)
        ttk.Entry(f, textvariable=v, width=10).pack(side=tk.RIGHT)
        store[key] = v

    def _add_elec(self, dist="0.0", apt="1.0", volt="0.0", chamfer="45.0"):
        idx = len(self.electrode_rows)+1
        row = ttk.Frame(self.elec_frame); row.pack(fill=tk.X, pady=1)
        ttk.Label(row, text=str(idx), width=3, anchor="center").pack(side=tk.LEFT, padx=1)
        vs = {}
        for k,d,w in [("dist",dist,8),("apt",apt,7),("volt",volt,7),("chamfer",chamfer,6)]:
            v = tk.StringVar(value=d); ttk.Entry(row, textvariable=v, width=w).pack(side=tk.LEFT, padx=1)
            vs[k] = v
        self.electrode_rows.append({"frame":row, "vars":vs})

    def _rm_elec(self):
        if len(self.electrode_rows) <= 1: return
        self.electrode_rows.pop()["frame"].destroy()

    # -- Callbacks --
    def _on_species(self, _e=None):
        p = self._species_presets.get(self.species_var.get())
        if p and p[0] is not None:
            self.beam_mass_var.set(str(p[0])); self.beam_charge_var.set(str(p[1]))
            self.mass_entry.config(state="disabled"); self.charge_entry.config(state="disabled")
        else:
            self.mass_entry.config(state="normal"); self.charge_entry.config(state="normal")

    def _on_beam_mode(self, _e=None):
        if self.beam_mode_var.get() == "Twiss (KV)":
            self.twiss_frame.pack(fill=tk.X, padx=12)
        else:
            self.twiss_frame.pack_forget()

    def _on_bfield_mode(self, _e=None):
        if self.bfield_var.get() == "Solenoid":
            self.sol_frame.pack(fill=tk.X, padx=12)
        else:
            self.sol_frame.pack_forget()

    def _on_scan_type(self, *_):
        is_v = self.scan_type.get() == "Voltage scan"
        u = "kV" if is_v else "A/m\u00b2"
        self.scan_labels["min"].config(text=u); self.scan_labels["max"].config(text=u)
        if is_v: self.scan_params["min"].set("-2.0"); self.scan_params["max"].set("-15.0")
        else: self.scan_params["min"].set("100"); self.scan_params["max"].set("1000")

    # ==================================================================
    # Config + validation
    # ==================================================================
    def _write_config(self, voltage_override=None, current_override=None):
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        cur = current_override or self.general["current"].get()
        with open(CONFIG_FILE,"w") as f:
            f.write(f"grid_points={self.general['grid'].get()}\n")
            f.write(f"particles={self.general['particles'].get()}\n")
            f.write(f"current_density={cur}\n")
            f.write(f"beam_energy={self.general['energy'].get()}\n")
            f.write(f"beam_mass={self.beam_mass_var.get()}\n")
            f.write(f"beam_charge={self.beam_charge_var.get()}\n")
            # Solver
            sm = {"BiCGSTAB":"bicgstab","Multigrid":"multigrid","Gauss-Seidel":"gaussseidel"}
            f.write(f"solver={sm.get(self.solver_var.get(),'bicgstab')}\n")
            f.write(f"mg_levels={self.mg_levels_var.get()}\n")
            # Beam mode
            if self.beam_mode_var.get() == "Twiss (KV)":
                f.write("beam_mode=twiss\n")
                f.write(f"beam_alpha={self.twiss_vars['tw_alpha'].get()}\n")
                f.write(f"beam_beta={self.twiss_vars['tw_beta'].get()}\n")
                f.write(f"beam_emittance={self.twiss_vars['tw_emit'].get()}\n")
            else:
                f.write("beam_mode=energy\n")
            # B-field
            if self.bfield_var.get() == "Solenoid":
                f.write("bfield=solenoid\n")
                f.write(f"sol_B0={self.sol_vars['sol_B0'].get()}\n")
                f.write(f"sol_z1={self.sol_vars['sol_z1'].get()}\n")
                f.write(f"sol_z2={self.sol_vars['sol_z2'].get()}\n")
            else:
                f.write("bfield=none\n")
            # Electrodes
            n = len(self.electrode_rows); f.write(f"electrodes={n}\n")
            for i,r in enumerate(self.electrode_rows):
                v = r["vars"]; volt = v["volt"].get()
                if voltage_override and (i+1)==voltage_override[0]: volt = str(voltage_override[1])
                f.write(f"{v['dist'].get()} {v['apt'].get()} {volt} {v['chamfer'].get()}\n")

    def _validate(self):
        try:
            if int(self.general["grid"].get()) < 20: return "Grid >= 20"
        except: return "Grid must be integer"
        try:
            if int(self.general["particles"].get()) < 100: return "Particles >= 100"
        except: return "Particles must be integer"
        positions = []
        for i,r in enumerate(self.electrode_rows):
            v = r["vars"]
            try:
                d,a,ch = float(v["dist"].get()), float(v["apt"].get()), float(v["chamfer"].get())
            except: return f"Electrode #{i+1}: invalid numbers"
            if a <= 0: return f"Electrode #{i+1}: aperture > 0"
            if ch < 5: return f"Electrode #{i+1}: chamfer >= 5"
            t = max(2.0, 4.0*a); positions.append((d, d+t, i+1))
        positions.sort()
        for j in range(len(positions)-1):
            if positions[j][1] > positions[j+1][0]:
                return f"Electrodes #{positions[j][2]} and #{positions[j+1][2]} overlap"
        return None

    # ==================================================================
    # Run simulation
    # ==================================================================
    def _on_run(self):
        err = self._validate()
        if err: messagebox.showwarning("Invalid", err); return
        self.run_btn.config(state="disabled"); self.scan_btn.config(state="disabled")
        self.status_var.set("Running..."); self.progress.start(15)
        threading.Thread(target=self._run_single, daemon=True).start()

    def _run_single(self):
        if self._run_sim(): self.root.after(0, self._sim_done)

    def _run_sim(self, voltage_override=None, current_override=None):
        try:
            self._write_config(voltage_override, current_override)
            cmd = [SIM_BINARY, "--config", CONFIG_FILE, "--outdir", OUTPUT_DIR]
            env = os.environ.copy()
            if getattr(sys,'frozen',False):
                env["LD_LIBRARY_PATH"] = os.path.dirname(SIM_BINARY)+":"+env.get("LD_LIBRARY_PATH","")
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=env)
            if proc.returncode != 0:
                self.root.after(0, lambda: self._sim_error(proc.stderr)); return False
            self._load_results(); return True
        except subprocess.TimeoutExpired:
            self.root.after(0, lambda: self._sim_error("Timeout")); return False
        except Exception as e:
            self.root.after(0, lambda: self._sim_error(str(e))); return False

    def _sim_error(self, msg):
        self.progress.stop(); self.run_btn.config(state="normal"); self.scan_btn.config(state="normal")
        self.status_var.set("Error"); messagebox.showerror("Error", msg)

    def _sim_done(self):
        self.progress.stop(); self.run_btn.config(state="normal"); self.scan_btn.config(state="normal")
        self.status_var.set("Done")
        self._update_trajectory_plot(); self._update_divergence_plot()
        self._update_field_plot(); self._update_convergence_plot(); self._update_energy_plot()
        self._setup_slider(); self._on_slider(None)

    # ==================================================================
    # Analysis helpers
    # ==================================================================
    def _last_elec_exit_mm(self):
        exits = []
        for r in self.electrode_rows:
            v = r["vars"]; d = float(v["dist"].get()); a = float(v["apt"].get())
            exits.append(d + max(2.0, 4.0*a))
        return max(exits)

    def _idx_at_x(self, xt):
        xl = self.emittance_data["x_mm"]
        return min(range(len(xl)), key=lambda i: abs(xl[i]-xt))

    def _grid_ratio(self):
        c = self.emittance_data["current_A"]
        if not c: return 0.0
        i0 = abs(c[0]); ie = abs(c[self._idx_at_x(self._last_elec_exit_mm())])
        return 1.0 - ie/i0 if i0 > 0 else 0.0

    def _transmission(self):
        c = np.array(self.emittance_data["current_A"])
        i0 = abs(c[0]) if abs(c[0]) > 0 else 1.0
        return np.abs(c) / i0 * 100.0

    # ==================================================================
    # Load results
    # ==================================================================
    def _load_results(self):
        self.emittance_data = self._csv(os.path.join(OUTPUT_DIR,"emittance_profile.csv"))
        raw = self._csv(os.path.join(OUTPUT_DIR,"phase_space.csv"))
        self.phase_space = {}
        if raw and raw.get("x_mm"):
            xv,yv,ypv = np.array(raw["x_mm"]),np.array(raw["y_mm"]),np.array(raw["yp_mrad"])
            for x in np.unique(xv):
                m = xv==x; yh,yph = yv[m],ypv[m]
                self.phase_space[x] = (np.concatenate([yh,-yh]), np.concatenate([yph,-yph]))
        tf = os.path.join(OUTPUT_DIR,"trajectory.png")
        self.traj_image = Image.open(tf) if os.path.exists(tf) else None
        self.field_data = self._csv(os.path.join(OUTPUT_DIR,"field_along_axis.csv"))
        self.convergence_data = self._load_convergence()
        self.electrode_currents = self._csv(os.path.join(OUTPUT_DIR,"electrode_currents.csv"))
        self.energy_data = self._csv(os.path.join(OUTPUT_DIR,"energy_distribution.csv"))

    def _csv(self, path):
        if not os.path.exists(path): return {}
        data = {}
        with open(path) as f:
            hdr = f.readline().strip().split(",")
            for h in hdr: data[h] = []
            for line in f:
                parts = line.strip().split(",")
                if len(parts) != len(hdr): continue
                for h,v in zip(hdr,parts):
                    try: data[h].append(float(v))
                    except: data[h].append(0.0)
        return data

    def _load_convergence(self):
        path = os.path.join(OUTPUT_DIR,"convergence.csv")
        if not os.path.exists(path): return None
        iters, epot_d, sch_d = [], [], []
        with open(path) as f:
            for line in f:
                if line.startswith("#") or not line.strip(): continue
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        iters.append(int(parts[0]))
                        epot_d.append(float(parts[1]))
                        sch_d.append(float(parts[2]))
                    except: pass
        if iters: return {"iter":iters, "epot":epot_d, "scharge":sch_d}
        return None

    # ==================================================================
    # Plot updates
    # ==================================================================
    def _update_trajectory_plot(self):
        self.traj_ax.clear()
        if self.traj_image: self.traj_ax.imshow(np.array(self.traj_image), aspect="auto")
        self.traj_ax.axis("off"); self.traj_fig.tight_layout(pad=0.5); self.traj_canvas.draw()

    def _update_divergence_plot(self):
        self.div_fig.clear(); self.div_ax = self.div_fig.add_subplot(111)
        if not self.emittance_data or not self.emittance_data.get("x_mm"):
            self.div_canvas.draw(); return
        x = np.array(self.emittance_data["x_mm"])
        r_rms = np.array(self.emittance_data["r_rms_mm"])
        div = np.array(self.emittance_data["divergence_mrad"])
        trans = self._transmission()
        c1,c2,c3 = "#2196F3","#FF5722","#4CAF50"
        ln1 = self.div_ax.plot(x, r_rms, color=c1, lw=2, label="r_rms")
        self.div_ax.set_xlabel("x (mm)"); self.div_ax.set_ylabel("RMS size (mm)", color=c1)
        self.div_ax.tick_params(axis="y", labelcolor=c1)
        ax2 = self.div_ax.twinx()
        ln2 = ax2.plot(x, div, color=c2, lw=2, ls="--", label="Divergence")
        ln3 = ax2.plot(x, trans, color=c3, lw=2, ls=":", label="Transmission %")
        ax2.set_ylabel("Div (mrad) / Trans (%)", color=c2)
        try:
            xe = self._last_elec_exit_mm()
            self.div_ax.axvline(xe, color="gray", lw=1, ls="-.", alpha=0.5)
        except: pass
        lns = ln1+ln2+ln3
        self.div_ax.legend(lns, [l.get_label() for l in lns], loc="upper right", fontsize=7)
        self.div_ax.set_title("Envelope, Divergence & Transmission"); self.div_ax.grid(True, alpha=0.3)
        self._div_vline = None; self.div_fig.tight_layout(); self.div_canvas.draw()

    def _update_field_plot(self):
        self.field_fig.clear(); self.field_ax = self.field_fig.add_subplot(111)
        if not self.field_data or not self.field_data.get("x_mm"):
            self.field_canvas.draw(); return
        x = np.array(self.field_data["x_mm"])
        pot = np.array(self.field_data["potential_V"])
        ex = np.array(self.field_data["Ex_Vm"])
        c1,c2 = "#1565C0","#E65100"
        ln1 = self.field_ax.plot(x, pot, color=c1, lw=2, label="Potential (V)")
        self.field_ax.set_xlabel("x (mm)"); self.field_ax.set_ylabel("Potential (V)", color=c1)
        self.field_ax.tick_params(axis="y", labelcolor=c1)
        ax2 = self.field_ax.twinx()
        ln2 = ax2.plot(x, ex, color=c2, lw=1.5, ls="--", label="E_x (V/m)")
        ax2.set_ylabel("E-field (V/m)", color=c2); ax2.tick_params(axis="y", labelcolor=c2)
        # B-field if present
        if self.field_data.get("Bx_T"):
            bx = np.array(self.field_data["Bx_T"])
            if np.any(np.abs(bx) > 1e-10):
                ax3 = self.field_ax.twinx()
                ax3.spines["right"].set_position(("axes", 1.15))
                ln3 = ax3.plot(x, bx, color="#4CAF50", lw=2, ls=":", label="B_x (T)")
                ax3.set_ylabel("B-field (T)", color="#4CAF50")
                ln2 = ln2 + ln3
        lns = ln1+ln2
        self.field_ax.legend(lns, [l.get_label() for l in lns], loc="upper right", fontsize=7)
        self.field_ax.set_title("Field Diagnostics Along Axis"); self.field_ax.grid(True, alpha=0.3)
        self.field_fig.tight_layout(); self.field_canvas.draw()

    def _update_convergence_plot(self):
        self.conv_fig.clear(); self.conv_ax = self.conv_fig.add_subplot(111)
        if not self.convergence_data:
            self.conv_canvas.draw(); return
        it = self.convergence_data["iter"]
        ep = self.convergence_data["epot"]
        sc = self.convergence_data["scharge"]
        self.conv_ax.semilogy(it, ep, "o-", color="#1565C0", label="Epot diff", markersize=6)
        self.conv_ax.semilogy(it, sc, "s-", color="#E65100", label="Scharge diff", markersize=6)
        self.conv_ax.set_xlabel("Iteration"); self.conv_ax.set_ylabel("Max difference")
        self.conv_ax.set_title("Convergence History"); self.conv_ax.legend(fontsize=9)
        self.conv_ax.grid(True, alpha=0.3, which="both")
        self.conv_fig.tight_layout(); self.conv_canvas.draw()

    def _update_energy_plot(self):
        self.ek_fig.clear(); self.ek_ax = self.ek_fig.add_subplot(111)
        if not self.energy_data or not self.energy_data.get("x_mm"):
            self.ek_canvas.draw(); return
        x = np.array(self.energy_data["x_mm"])
        ek = np.array(self.energy_data["ek_eV"])
        unique_x = np.unique(x)
        # Plot energy histogram at exit
        if len(unique_x) > 0:
            xe = unique_x[-1]  # last sampled position
            mask = x == xe
            ek_exit = ek[mask]
            if len(ek_exit) > 0:
                self.ek_ax.hist(ek_exit, bins=30, color="#1565C0", alpha=0.7, edgecolor="black")
                self.ek_ax.set_xlabel("Kinetic Energy (eV)")
                self.ek_ax.set_ylabel("Particle count")
                mean_e = np.mean(ek_exit); spread = np.std(ek_exit)
                self.ek_ax.set_title(f"Energy at x={xe:.1f}mm  "
                                     f"E={mean_e:.1f}\u00b1{spread:.2f} eV")
                self.ek_ax.axvline(mean_e, color="red", ls="--", lw=1.5)
        self.ek_ax.grid(True, alpha=0.3)
        self.ek_fig.tight_layout(); self.ek_canvas.draw()

    def _setup_slider(self):
        if not self.emittance_data or not self.emittance_data.get("x_mm"): return
        n = len(self.emittance_data["x_mm"])
        self.slider.config(from_=0, to=n-1); self.slider.set(n//2)

    # ==================================================================
    # Slider
    # ==================================================================
    def _on_slider(self, _v):
        if not self.emittance_data or not self.phase_space: return
        xl = self.emittance_data["x_mm"]
        if not xl: return
        idx = max(0, min(int(float(self.slider.get())), len(xl)-1))
        xm = xl[idx]
        psk = sorted(self.phase_space.keys())
        if not psk: return
        cx = min(psk, key=lambda k: abs(k-xm))
        yd, ypd = self.phase_space[cx]
        em, al, be, _ = compute_twiss(yd, ypd)
        r_rms = self.emittance_data["r_rms_mm"][idx]
        div = self.emittance_data["divergence_mrad"][idx]
        cur = self.emittance_data["current_A"][idx]
        try: gr = self._grid_ratio()*100
        except: gr = 0
        # Electrode currents
        elec_str = ""
        if self.electrode_currents and self.electrode_currents.get("electrode"):
            for i in range(len(self.electrode_currents["electrode"])):
                eid = self.electrode_currents["electrode"][i]
                ec = self.electrode_currents["current_A"][i]
                try: eid_s = f"#{int(eid)}"
                except: eid_s = str(eid)
                elec_str += f"  {eid_s}: {ec:.3e}A\n"

        self.pos_var.set(f"x = {xm:.2f} mm")
        self.info_var.set(
            f"Species : {self.species_var.get()}\n"
            f"x = {xm:.3f} mm\n"
            f"eps     = {em:.4f} mm-mrad\n"
            f"alpha   = {al:.3f}\n"
            f"beta    = {be:.4f} mm/mrad\n"
            f"r_rms   = {r_rms:.4f} mm\n"
            f"div     = {div:.3f} mrad\n"
            f"current = {cur:.4e} A\n"
            f"grid I% = {gr:.1f}%\n"
            f"transmit= {100-gr:.1f}%\n"
            + (f"Electrode currents:\n{elec_str}" if elec_str else ""))

        # Phase space
        self.emit_ax.clear()
        self.emit_ax.scatter(yd, ypd, s=1.5, alpha=0.4, c="#1565C0", edgecolors="none")
        self.emit_ax.set_xlabel("y (mm)"); self.emit_ax.set_ylabel("y' (mrad)")
        self.emit_ax.set_title(f"x={xm:.2f}mm  \u03b5={em:.4f} mm-mrad", fontsize=10)
        self.emit_ax.axhline(0,color="gray",lw=0.5); self.emit_ax.axvline(0,color="gray",lw=0.5)
        self.emit_ax.grid(True, alpha=0.3)
        if em > 0 and be > 0:
            t = np.linspace(0,2*np.pi,200); sb = np.sqrt(be)
            ye = np.sqrt(em)*sb*np.cos(t)
            ype = np.sqrt(em)*(-al/sb*np.cos(t)+1/sb*np.sin(t))
            self.emit_ax.plot(ye, ype, "r-", lw=1.5, alpha=0.8, label="RMS ellipse")
            self.emit_ax.legend(fontsize=8, loc="upper left")
        self.emit_fig.tight_layout(); self.emit_canvas.draw()

        # Div vline
        if self._div_vline:
            try: self._div_vline.remove()
            except: pass
            self._div_vline = None
        self._div_vline = self.div_ax.axvline(xm, color="red", lw=1.5, ls=":", alpha=0.7)
        self.div_canvas.draw()

    # ==================================================================
    # Perveance scan
    # ==================================================================
    def _on_run_scan(self):
        err = self._validate()
        if err: messagebox.showwarning("Invalid", err); return
        try:
            smin,smax = float(self.scan_params["min"].get()),float(self.scan_params["max"].get())
            steps = int(self.scan_params["steps"].get())
            if steps < 2: raise ValueError
        except: messagebox.showwarning("Invalid","Check min/max/steps"); return
        is_v = self.scan_type.get() == "Voltage scan"
        self.scan_running = True; self.scan_stop = False; self.scan_is_vscan = is_v
        self.scan_frames = []
        self.scan_data = {"perveance":[],"divergence":[],"voltage":[],"current":[],
                          "grid_ratio":[],"scan_val":[]}
        self.run_btn.config(state="disabled"); self.scan_btn.config(state="disabled")
        self.scan_stop_btn.config(state="normal")
        self.right_nb.select(3)  # Perveance Scan tab
        threading.Thread(target=self._scan_thread, args=(smin,smax,steps,is_v), daemon=True).start()

    def _scan_thread(self, smin, smax, steps, is_v):
        vals = np.linspace(smin, smax, steps)
        eidx = int(self.scan_electrode.get()) if is_v else 0
        for i,val in enumerate(vals):
            if self.scan_stop: break
            self.root.after(0, lambda i=i,s=steps,v=val:
                self.scan_status.set(f"Scan {i+1}/{s} ({'V' if is_v else 'J'}={v:.2f})..."))
            if is_v:
                ok = self._run_sim(voltage_override=(eidx, val))
                v_volts = abs(val)*1e3; label = f"{val:+.1f} kV"
            else:
                ok = self._run_sim(current_override=str(val))
                try: vk = float(self.electrode_rows[1]["vars"]["volt"].get())
                except: vk = -8.0
                v_volts = abs(vk)*1e3; label = f"{val:.0f} A/m\u00b2"
            if not ok: break
            if not self.emittance_data or not self.emittance_data.get("divergence_mrad"): continue
            xe = self._last_elec_exit_mm(); ie = self._idx_at_x(xe)
            div = self.emittance_data["divergence_mrad"][ie]
            cur = self.emittance_data["current_A"][0]
            gr = self._grid_ratio()
            if v_volts > 0:
                perv = (abs(cur)/(v_volts**1.5))*1e6
            else: continue
            self.scan_data["perveance"].append(perv); self.scan_data["divergence"].append(div)
            self.scan_data["voltage"].append(v_volts); self.scan_data["current"].append(cur)
            self.scan_data["grid_ratio"].append(gr*100); self.scan_data["scan_val"].append(label)
            self.root.after(0, self._update_all_scan)
        self.root.after(0, self._scan_done)

    def _update_all_scan(self):
        self._update_trajectory_plot(); self._update_divergence_plot()
        self._update_field_plot(); self._update_convergence_plot(); self._update_energy_plot()
        self._setup_slider(); self._on_slider(None)
        self._redraw_scan(); self._capture_frame()

    def _redraw_scan(self):
        self.scan_fig.clear(); ax1 = self.scan_fig.add_subplot(111); self.scan_ax = ax1
        p,d,gr,lb = self.scan_data["perveance"],self.scan_data["divergence"],\
                     self.scan_data["grid_ratio"],self.scan_data["scan_val"]
        c1,c2 = "#1565C0","#E65100"
        ln1 = ax1.plot(p,d,"o-",color=c1,markersize=7,lw=2,label="Divergence")
        ax1.set_xlabel("Perveance ($\\mu$A / V$^{3/2}$)")
        ax1.set_ylabel("RMS Divergence (mrad)",color=c1); ax1.tick_params(axis="y",labelcolor=c1)
        ax2 = ax1.twinx()
        ln2 = ax2.plot(p,gr,"s--",color=c2,markersize=5,lw=1.5,label="Grid I %")
        ax2.set_ylabel("Grid current (%)",color=c2); ax2.tick_params(axis="y",labelcolor=c2)
        ax2.set_ylim(bottom=0)
        for i in range(len(p)):
            ax1.annotate(lb[i],(p[i],d[i]),textcoords="offset points",xytext=(6,6),fontsize=7,color="gray")
        st = "Voltage" if self.scan_is_vscan else "Current"
        ax1.set_title(f"Perveance Scan ({st})"); ax1.grid(True,alpha=0.3)
        if len(d) > 1:
            im = int(np.argmin(d))
            ax1.plot(p[im],d[im],"*",color="red",markersize=15,zorder=5,
                     label=f"Min: {d[im]:.2f}mrad @ {lb[im]}")
        lns = ln1+ln2; ax1.legend(lns,[l.get_label() for l in lns],loc="upper left",fontsize=7)
        self.scan_fig.tight_layout(); self.scan_canvas.draw()

    def _scan_done(self):
        self.scan_running = False; self.run_btn.config(state="normal")
        self.scan_btn.config(state="normal"); self.scan_stop_btn.config(state="disabled")
        self.scan_status.set(f"Done ({len(self.scan_data['perveance'])} pts)" if not self.scan_stop else "Stopped")
        self.root.after(0, self._sim_done)

    # ==================================================================
    # GIF capture + save
    # ==================================================================
    def _capture_frame(self):
        try:
            imgs = []
            for fig in [self.traj_fig, self.emit_fig, self.scan_fig]:
                fig.canvas.draw(); buf = fig.canvas.buffer_rgba()
                w,h = fig.canvas.get_width_height()
                imgs.append(Image.frombuffer("RGBA",(w,h),buf).copy())
            tw = imgs[0].width; bh = max(imgs[1].height, imgs[2].height); hw = tw//2
            comp = Image.new("RGBA",(tw, imgs[0].height+bh),(255,255,255,255))
            comp.paste(imgs[0],(0,0))
            comp.paste(imgs[1].resize((hw,bh),Image.LANCZOS),(0,imgs[0].height))
            comp.paste(imgs[2].resize((tw-hw,bh),Image.LANCZOS),(hw,imgs[0].height))
            self.scan_frames.append(comp.convert("RGB"))
        except: pass

    def _save_plots(self):
        folder = filedialog.askdirectory(title="Save plots to")
        if not folder: return
        try:
            for fig,name in [(self.traj_fig,"trajectory"),(self.emit_fig,"phase_space"),
                              (self.div_fig,"envelope"),(self.scan_fig,"perveance_scan"),
                              (self.field_fig,"field_diagnostics"),(self.conv_fig,"convergence"),
                              (self.ek_fig,"energy_distribution")]:
                fig.savefig(os.path.join(folder,f"{name}.png"), dpi=150, bbox_inches="tight")
            self.save_status.set(f"Saved 7 plots to {folder}")
        except Exception as e: messagebox.showerror("Error",str(e))

    def _save_gif(self):
        if not self.scan_frames:
            messagebox.showinfo("No frames","Run a scan first."); return
        path = filedialog.asksaveasfilename(title="Save GIF",defaultextension=".gif",
                                            filetypes=[("GIF","*.gif")])
        if not path: return
        try:
            dur = [800]*len(self.scan_frames); dur[-1] = 2000
            self.scan_frames[0].save(path, save_all=True, append_images=self.scan_frames[1:],
                                     duration=dur, loop=0, optimize=True)
            self.save_status.set(f"GIF saved: {os.path.basename(path)}")
        except Exception as e: messagebox.showerror("Error",str(e))


def main():
    if not os.path.exists(SIM_BINARY):
        print(f"Error: {SIM_BINARY} not found. Build: cd gui && make"); sys.exit(1)
    root = tk.Tk(); BeamGUI(root); root.mainloop()

if __name__ == "__main__":
    main()
