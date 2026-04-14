"""
Tab 1: IBSIMU Beam Extraction Design
Electrode geometry, beam species, solver, B-field, run + results.
"""
import os, threading
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from PIL import Image

from config import SPECIES, DATA_DIR
from workflows import extraction


class ExtractionTab:
    def __init__(self, parent, app):
        self.parent = parent
        self.app = app
        self.results = None
        self._build_ui()

    def _build_ui(self):
        pw = ttk.PanedWindow(self.parent, orient=tk.HORIZONTAL)
        pw.pack(fill=tk.BOTH, expand=True)

        # --- LEFT: controls ---
        lo = ttk.Frame(pw, width=320); pw.add(lo, weight=0)
        lc = tk.Canvas(lo, highlightthickness=0, width=310)
        lsb = ttk.Scrollbar(lo, orient=tk.VERTICAL, command=lc.yview)
        left = ttk.Frame(lc)
        left.bind("<Configure>", lambda e: lc.configure(scrollregion=lc.bbox("all")))
        lc.create_window((0, 0), window=left, anchor="nw")
        lc.configure(yscrollcommand=lsb.set)
        lc.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        lsb.pack(side=tk.RIGHT, fill=tk.Y)

        # General
        self._section(left, "General")
        self.gen = {}
        for l, k, d in [("Grid points", "grid_points", "241"),
                         ("Particles", "particles", "5000"),
                         ("Current (A/m\u00b2)", "current_density", "600"),
                         ("Beam energy (eV)", "beam_energy", "5.0")]:
            self._row(left, l, k, d, self.gen)

        # Solver
        sf = ttk.Frame(left); sf.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(sf, text="Solver:", width=16, anchor="w").pack(side=tk.LEFT)
        self.solver_var = tk.StringVar(value="BiCGSTAB")
        ttk.Combobox(sf, textvariable=self.solver_var, width=12,
                     values=["BiCGSTAB", "Multigrid", "Gauss-Seidel"],
                     state="readonly").pack(side=tk.LEFT)

        # Species
        self._section(left, "Beam Species")
        sf = ttk.Frame(left); sf.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(sf, text="Species:", width=10, anchor="w").pack(side=tk.LEFT)
        self.species_var = tk.StringVar(value="H+")
        cb = ttk.Combobox(sf, textvariable=self.species_var, width=14,
                          values=list(SPECIES.keys()) + ["Custom"], state="readonly")
        cb.pack(side=tk.LEFT); cb.bind("<<ComboboxSelected>>", self._on_species)
        self.mass_var = tk.StringVar(value="1.00728")
        self.charge_var = tk.StringVar(value="1.0")
        for l, v in [("Mass (amu):", self.mass_var), ("Charge (e):", self.charge_var)]:
            f = ttk.Frame(left); f.pack(fill=tk.X, padx=12, pady=1)
            ttk.Label(f, text=l, width=14, anchor="w").pack(side=tk.LEFT)
            ttk.Entry(f, textvariable=v, width=12).pack(side=tk.LEFT)

        # Electrodes
        self._section(left, "Electrodes")
        hdr = ttk.Frame(left); hdr.pack(fill=tk.X, padx=12)
        for t, w in [("#", 3), ("Dist", 7), ("Apt", 6), ("V(kV)", 6), ("Thk", 5), ("Chm", 5)]:
            ttk.Label(hdr, text=t, width=w, anchor="center",
                      font=("Helvetica", 9, "bold")).pack(side=tk.LEFT, padx=1)
        self.elec_frame = ttk.Frame(left); self.elec_frame.pack(fill=tk.X, padx=12)
        self.electrode_rows = []
        self._add_elec("0.0", "0.5", "0.0", "2.0", "0")
        self._add_elec("10.0", "1.5", "-8.0", "2.0", "0")
        bf = ttk.Frame(left); bf.pack(fill=tk.X, padx=12, pady=4)
        ttk.Button(bf, text="+ Add", command=lambda: self._add_elec()).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="- Remove", command=self._rm_elec).pack(side=tk.LEFT, padx=2)

        # Run + handoff
        self._section(left, "")
        self.run_btn = ttk.Button(left, text="Run Extraction", command=self._on_run)
        self.run_btn.pack(fill=tk.X, padx=12, ipady=5)
        self.handoff_btn = ttk.Button(left, text="Send Beam to Transport Tab >>",
                                      command=self._on_handoff, state="disabled")
        self.handoff_btn.pack(fill=tk.X, padx=12, pady=4, ipady=3)
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(left, textvariable=self.status_var, foreground="gray").pack(padx=12)
        self.progress = ttk.Progressbar(left, mode="indeterminate")
        self.progress.pack(fill=tk.X, padx=12, pady=2)

        # Info
        self._section(left, "Beam Info")
        self.info_var = tk.StringVar(value="Run extraction first")
        ttk.Label(left, textvariable=self.info_var, justify=tk.LEFT,
                  wraplength=280, font=("Courier", 9)).pack(padx=12, anchor="w")

        # --- RIGHT: plots ---
        right = ttk.Frame(pw); pw.add(right, weight=1)
        rpw = ttk.PanedWindow(right, orient=tk.VERTICAL)
        rpw.pack(fill=tk.BOTH, expand=True)

        # Trajectory
        tf = ttk.LabelFrame(rpw, text="Particle Trajectories")
        rpw.add(tf, weight=1)
        self.traj_fig = Figure(figsize=(10, 3), dpi=100, facecolor="#f5f5f5")
        self.traj_ax = self.traj_fig.add_subplot(111); self.traj_ax.axis("off")
        self.traj_canvas = FigureCanvasTkAgg(self.traj_fig, master=tf)
        self.traj_canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Bottom: phase space + envelope
        botf = ttk.Frame(rpw); rpw.add(botf, weight=1)
        slf = ttk.Frame(botf); slf.pack(fill=tk.X, padx=10, pady=(4, 0))
        ttk.Label(slf, text="Axial pos:").pack(side=tk.LEFT)
        self.pos_var = tk.StringVar(value="--")
        ttk.Label(slf, textvariable=self.pos_var, font=("Courier", 10)).pack(side=tk.LEFT, padx=8)
        self.slider = ttk.Scale(slf, from_=0, to=100, orient=tk.HORIZONTAL, command=self._on_slider)
        self.slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=8)

        pf = ttk.Frame(botf); pf.pack(fill=tk.BOTH, expand=True)
        self.emit_fig = Figure(figsize=(5, 3.3), dpi=100)
        self.emit_ax = self.emit_fig.add_subplot(111)
        self.emit_canvas = FigureCanvasTkAgg(self.emit_fig, master=pf)
        self.emit_canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.div_fig = Figure(figsize=(5, 3.3), dpi=100)
        self.div_ax = self.div_fig.add_subplot(111)
        self.div_canvas = FigureCanvasTkAgg(self.div_fig, master=pf)
        self.div_canvas.get_tk_widget().pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self._div_vline = None

    # -- helpers --
    def _section(self, p, title):
        ttk.Separator(p).pack(fill=tk.X, padx=12, pady=6)
        if title:
            ttk.Label(p, text=title, font=("Helvetica", 11, "bold")).pack(pady=(0, 4))

    def _row(self, p, label, key, default, store):
        f = ttk.Frame(p); f.pack(fill=tk.X, padx=12, pady=2)
        ttk.Label(f, text=label, width=18, anchor="w").pack(side=tk.LEFT)
        v = tk.StringVar(value=default)
        ttk.Entry(f, textvariable=v, width=10).pack(side=tk.RIGHT)
        store[key] = v

    def _add_elec(self, dist="0.0", apt="1.0", volt="0.0", thick="2.0", chamfer="0"):
        idx = len(self.electrode_rows) + 1
        row = ttk.Frame(self.elec_frame); row.pack(fill=tk.X, pady=1)
        ttk.Label(row, text=str(idx), width=3, anchor="center").pack(side=tk.LEFT, padx=1)
        vs = {}
        for k, d, w in [("dist", dist, 7), ("apt", apt, 6), ("volt", volt, 6),
                         ("thick", thick, 5), ("chamfer", chamfer, 5)]:
            v = tk.StringVar(value=d)
            ttk.Entry(row, textvariable=v, width=w).pack(side=tk.LEFT, padx=1)
            vs[k] = v
        self.electrode_rows.append({"frame": row, "vars": vs})

    def _rm_elec(self):
        if len(self.electrode_rows) <= 1: return
        self.electrode_rows.pop()["frame"].destroy()

    def _on_species(self, _e=None):
        s = self.species_var.get()
        if s in SPECIES:
            m, q = SPECIES[s]
            self.mass_var.set(str(m)); self.charge_var.set(str(q))

    def _get_params(self):
        sm = {"BiCGSTAB": "bicgstab", "Multigrid": "multigrid", "Gauss-Seidel": "gaussseidel"}
        return {k: v.get() for k, v in self.gen.items()} | {
            "beam_mass": self.mass_var.get(),
            "beam_charge": self.charge_var.get(),
            "solver": sm.get(self.solver_var.get(), "bicgstab"),
            "beam_mode": "energy", "bfield": "none",
        }

    def _get_electrodes(self):
        return [{k: v.get() for k, v in r["vars"].items()} for r in self.electrode_rows]

    # -- run --
    def _on_run(self):
        self.run_btn.config(state="disabled"); self.handoff_btn.config(state="disabled")
        self.status_var.set("Running extraction..."); self.progress.start(15)
        threading.Thread(target=self._run_thread, daemon=True).start()

    def _run_thread(self):
        try:
            params = self._get_params()
            electrodes = self._get_electrodes()
            results, beam = extraction.run_with_handoff(params, electrodes)
            self.results = results
            self.app.shared["ibsimu_results"] = results
            self.app.shared["beam_handoff"] = beam
            self.app.shared["extraction_params"] = params
            self.app.shared["extraction_electrodes"] = electrodes
            self._load_phase_space()
            self.parent.after(0, self._done)
        except Exception as e:
            self.parent.after(0, lambda: self._error(str(e)))

    def _error(self, msg):
        self.progress.stop(); self.run_btn.config(state="normal")
        self.status_var.set("Error"); messagebox.showerror("Error", msg)

    def _done(self):
        self.progress.stop(); self.run_btn.config(state="normal")
        self.handoff_btn.config(state="normal")
        n = self.app.shared["beam_handoff"]["n_particles"]
        cur = self.app.shared["beam_handoff"]["current_A"]
        self.status_var.set(f"Done - {n} particles, I={cur:.4e} A")
        self._update_plots()

    def _on_handoff(self):
        if self.app.shared["beam_handoff"] is None:
            messagebox.showinfo("No beam", "Run extraction first"); return
        self.app.tab2.receive_beam(self.app.shared["beam_handoff"])
        self.app.switch_to_tab(1)

    # -- data --
    def _load_phase_space(self):
        self.phase_space = {}
        r = self.results
        if not r or "phase_space" not in r:
            return
        ps = r["phase_space"]
        if not ps.get("x_mm"):
            return
        xv = np.array(ps["x_mm"]); yv = np.array(ps["y_mm"]); ypv = np.array(ps["yp_mrad"])
        for x in np.unique(xv):
            m = xv == x; yh, yph = yv[m], ypv[m]
            self.phase_space[x] = (np.concatenate([yh, -yh]), np.concatenate([yph, -yph]))

    # -- plots --
    def _update_plots(self):
        # Trajectory
        self.traj_ax.clear()
        tpath = self.results.get("trajectory_png")
        if tpath and os.path.exists(tpath):
            self.traj_ax.imshow(np.array(Image.open(tpath)), aspect="auto")
        self.traj_ax.axis("off"); self.traj_fig.tight_layout(pad=0.5); self.traj_canvas.draw()

        # Envelope
        self.div_fig.clear(); self.div_ax = self.div_fig.add_subplot(111)
        ep = self.results.get("emittance_profile", {})
        if ep.get("x_mm"):
            x = np.array(ep["x_mm"])
            self.div_ax.plot(x, ep["r_rms_mm"], color="#2196F3", lw=2, label="r_rms")
            ax2 = self.div_ax.twinx()
            ax2.plot(x, ep["divergence_mrad"], color="#FF5722", lw=2, ls="--", label="Div")
            self.div_ax.set_xlabel("x (mm)"); self.div_ax.set_ylabel("r_rms (mm)")
            ax2.set_ylabel("Div (mrad)")
            self.div_ax.set_title("Envelope & Divergence"); self.div_ax.grid(True, alpha=0.3)
        self._div_vline = None; self.div_fig.tight_layout(); self.div_canvas.draw()

        # Slider
        if ep.get("x_mm"):
            n = len(ep["x_mm"])
            self.slider.config(from_=0, to=n - 1); self.slider.set(n // 2)
        self._on_slider(None)

    def _on_slider(self, _v):
        if not self.results or not hasattr(self, "phase_space"):
            return
        ep = self.results.get("emittance_profile", {})
        if not ep.get("x_mm"):
            return
        idx = max(0, min(int(float(self.slider.get())), len(ep["x_mm"]) - 1))
        xm = ep["x_mm"][idx]
        self.pos_var.set(f"x = {xm:.2f} mm")

        # Phase space
        self.emit_ax.clear()
        psk = sorted(self.phase_space.keys())
        if psk:
            cx = min(psk, key=lambda k: abs(k - xm))
            yd, ypd = self.phase_space[cx]
            self.emit_ax.scatter(yd, ypd, s=1.5, alpha=0.4, c="#1565C0", edgecolors="none")
            self.emit_ax.set_title(f"Phase Space x={xm:.2f}mm", fontsize=10)
        self.emit_ax.set_xlabel("y (mm)"); self.emit_ax.set_ylabel("y' (mrad)")
        self.emit_ax.grid(True, alpha=0.3)
        self.emit_fig.tight_layout(); self.emit_canvas.draw()

        # Info
        div = ep["divergence_mrad"][idx] if ep.get("divergence_mrad") else 0
        cur = ep["current_A"][idx] if ep.get("current_A") else 0
        self.info_var.set(
            f"x = {xm:.3f} mm\n"
            f"div = {div:.3f} mrad\n"
            f"current = {cur:.4e} A\n"
            f"species = {self.species_var.get()}")

        # Div vline
        if self._div_vline:
            try: self._div_vline.remove()
            except: pass
        self._div_vline = self.div_ax.axvline(xm, color="red", lw=1.5, ls=":", alpha=0.7)
        self.div_canvas.draw()
