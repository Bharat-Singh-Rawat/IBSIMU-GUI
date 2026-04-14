"""
Tab 4: Digital Twin - Automated Loop Controller
Runs iterative extraction -> erosion -> geometry update cycles.
Also supports perveance scan + erosion coupling.
"""
import threading
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from gui.plot_helpers import make_canvas, placeholder
from workflows import digital_twin_loop, perveance_erosion


class DigitalTwinTab:
    def __init__(self, parent, app):
        self.parent = parent
        self.app = app
        self.running = False
        self.stop_flag = False
        self._build_ui()

    def _build_ui(self):
        pw = ttk.PanedWindow(self.parent, orient=tk.HORIZONTAL)
        pw.pack(fill=tk.BOTH, expand=True)

        # --- LEFT: controls ---
        lo = ttk.Frame(pw, width=300); pw.add(lo, weight=0)
        lc = tk.Canvas(lo, highlightthickness=0, width=290)
        lsb = ttk.Scrollbar(lo, orient=tk.VERTICAL, command=lc.yview)
        left = ttk.Frame(lc)
        left.bind("<Configure>", lambda e: lc.configure(scrollregion=lc.bbox("all")))
        lc.create_window((0,0), window=left, anchor="nw")
        lc.configure(yscrollcommand=lsb.set)
        lc.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        lsb.pack(side=tk.RIGHT, fill=tk.Y)

        # -- Digital Twin Loop --
        ttk.Label(left, text="Digital Twin Loop", font=("Helvetica",12,"bold")).pack(pady=(8,4))
        ttk.Label(left, text="Uses geometry from Tab 1\nand grids from Tab 2",
                  foreground="gray", wraplength=260).pack(padx=12)

        self.twin_params = {}
        for l,k,d in [("Max cycles","max_cycles","10"),
                       ("Steps/cycle","steps_per_cycle","500"),
                       ("Hours/cycle","hours_per_cycle","100"),
                       ("Fail div (mrad)","failure_div","50"),
                       ("Fail temp (K)","failure_temp","2000")]:
            f = ttk.Frame(left); f.pack(fill=tk.X, padx=12, pady=2)
            ttk.Label(f, text=l, width=16, anchor="w").pack(side=tk.LEFT)
            v = tk.StringVar(value=d); ttk.Entry(f, textvariable=v, width=8).pack(side=tk.LEFT)
            self.twin_params[k] = v

        ttk.Separator(left).pack(fill=tk.X, padx=12, pady=8)
        self.twin_btn = ttk.Button(left, text="Run Digital Twin", command=self._on_run_twin)
        self.twin_btn.pack(fill=tk.X, padx=12, ipady=5)

        # -- Perveance + Erosion Scan --
        ttk.Separator(left).pack(fill=tk.X, padx=12, pady=8)
        ttk.Label(left, text="Perveance + Erosion Scan", font=("Helvetica",12,"bold")).pack(pady=(0,4))

        self.scan_params = {}
        for l,k,d in [("Scan type","scan_type","Voltage"),
                       ("Electrode #","scan_elec","2"),
                       ("Min","scan_min","-2"),("Max","scan_max","-15"),
                       ("Steps","scan_steps","6"),("Erosion iters","erosion_steps","200")]:
            f = ttk.Frame(left); f.pack(fill=tk.X, padx=12, pady=2)
            ttk.Label(f, text=l, width=14, anchor="w").pack(side=tk.LEFT)
            v = tk.StringVar(value=d); ttk.Entry(f, textvariable=v, width=10).pack(side=tk.LEFT)
            self.scan_params[k] = v

        self.scan_btn = ttk.Button(left, text="Run Perveance+Erosion Scan",
                                   command=self._on_run_scan)
        self.scan_btn.pack(fill=tk.X, padx=12, ipady=5, pady=4)

        # Stop + status
        self.stop_btn = ttk.Button(left, text="Stop",
                                   command=lambda: setattr(self, 'stop_flag', True),
                                   state="disabled")
        self.stop_btn.pack(fill=tk.X, padx=12)
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(left, textvariable=self.status_var, foreground="gray").pack(padx=12, pady=4)

        # --- RIGHT: plots ---
        right = ttk.Frame(pw); pw.add(right, weight=1)
        nb = ttk.Notebook(right)
        nb.pack(fill=tk.BOTH, expand=True)

        # Twin time-series tab
        tsf = ttk.Frame(nb); nb.add(tsf, text="Twin Time-Series")
        self.ts_fig, self.ts_ax, self.ts_canvas = make_canvas(tsf, figsize=(10, 7))
        placeholder(self.ts_ax, "Run digital twin to see time-series")

        # Perveance+Erosion tab
        pef = ttk.Frame(nb); nb.add(pef, text="Perveance+Erosion")
        self.pe_fig, self.pe_ax, self.pe_canvas = make_canvas(pef, figsize=(10, 7))
        placeholder(self.pe_ax, "Run scan to see results")

    # ==================================================================
    # Digital Twin Loop
    # ==================================================================
    def _on_run_twin(self):
        params = self.app.shared.get("extraction_params")
        elecs = self.app.shared.get("extraction_electrodes")
        if not params or not elecs:
            messagebox.showinfo("Setup", "Run extraction in Tab 1 first"); return

        self.running = True; self.stop_flag = False
        self.twin_btn.config(state="disabled"); self.scan_btn.config(state="disabled")
        self.stop_btn.config(state="normal")
        threading.Thread(target=self._twin_thread, daemon=True).start()

    def _twin_thread(self):
        try:
            params = self.app.shared["extraction_params"]
            elecs = self.app.shared["extraction_electrodes"]
            bemcs_grids = self.app.tab2._get_grids()
            plasma = self.app.tab2._get_plasma()
            mass = float(params.get("beam_mass", 1.0))
            charge = int(abs(float(params.get("beam_charge", 1))))
            tp = self.twin_params

            def cb(cycle, total, result):
                self.parent.after(0, lambda: self._update_twin_plot(cycle, total, result))

            history = digital_twin_loop.run(
                params, elecs, bemcs_grids, plasma,
                mass_amu=mass, charge_state=charge,
                steps_per_cycle=int(tp["steps_per_cycle"].get()),
                max_cycles=int(tp["max_cycles"].get()),
                hours_per_cycle=float(tp["hours_per_cycle"].get()),
                failure_div_mrad=float(tp["failure_div"].get()),
                failure_temp_K=float(tp["failure_temp"].get()),
                callback=cb)

            self.parent.after(0, lambda: self._twin_done(history))
        except Exception as e:
            self.parent.after(0, lambda: self._error(str(e)))

    def _update_twin_plot(self, cycle, total, result):
        self.status_var.set(f"Cycle {cycle+1}/{total}  div={result['divergence']:.2f}mrad")

    def _twin_done(self, history):
        self.running = False
        self.twin_btn.config(state="normal"); self.scan_btn.config(state="normal")
        self.stop_btn.config(state="disabled")

        failed = history.get("failed_at_cycle")
        self.status_var.set(
            f"Done - {len(history['hours'])} cycles, {history['total_hours']:.0f}h"
            + (f" (failed at cycle {failed})" if failed else ""))

        # Plot time-series
        self.ts_fig.clear()
        h = history
        if not h["hours"]:
            return

        ax1 = self.ts_fig.add_subplot(221)
        ax1.plot(h["hours"], h["divergence"], "o-", color="#1565C0", lw=2)
        ax1.set_xlabel("Hours"); ax1.set_ylabel("Divergence (mrad)")
        ax1.set_title("Divergence vs Time"); ax1.grid(True, alpha=0.3)

        ax2 = self.ts_fig.add_subplot(222)
        if h["apertures"] and len(h["apertures"][0]) > 0:
            ng = len(h["apertures"][0])
            for gi in range(ng):
                apt = [a[gi] for a in h["apertures"]]
                ax2.plot(h["hours"], apt, "o-", lw=2, label=f"Grid {gi+1}")
            ax2.legend(fontsize=8)
        ax2.set_xlabel("Hours"); ax2.set_ylabel("Aperture (mm)")
        ax2.set_title("Aperture vs Time"); ax2.grid(True, alpha=0.3)

        ax3 = self.ts_fig.add_subplot(223)
        ax3.plot(h["hours"], h["beam_current"], "o-", color="#4CAF50", lw=2, label="Beam")
        ax3.plot(h["hours"], h["backflow_current"], "s--", color="#E65100", lw=1.5, label="Backflow")
        ax3.set_xlabel("Hours"); ax3.set_ylabel("Current (A)")
        ax3.set_title("Currents vs Time"); ax3.legend(fontsize=8); ax3.grid(True, alpha=0.3)

        ax4 = self.ts_fig.add_subplot(224)
        ax4.plot(h["hours"], h["total_damage"], "o-", color="#9C27B0", lw=2)
        ax4.set_xlabel("Hours"); ax4.set_ylabel("Total Damage")
        ax4.set_title("Cumulative Erosion"); ax4.grid(True, alpha=0.3)

        self.ts_fig.tight_layout(); self.ts_canvas.draw()

    # ==================================================================
    # Perveance + Erosion Scan
    # ==================================================================
    def _on_run_scan(self):
        params = self.app.shared.get("extraction_params")
        elecs = self.app.shared.get("extraction_electrodes")
        if not params or not elecs:
            messagebox.showinfo("Setup", "Run extraction in Tab 1 first"); return

        self.running = True; self.stop_flag = False
        self.twin_btn.config(state="disabled"); self.scan_btn.config(state="disabled")
        self.stop_btn.config(state="normal")
        threading.Thread(target=self._scan_thread, daemon=True).start()

    def _scan_thread(self):
        try:
            params = self.app.shared["extraction_params"]
            elecs = self.app.shared["extraction_electrodes"]
            bemcs_grids = self.app.tab2._get_grids()
            plasma = self.app.tab2._get_plasma()
            sp = self.scan_params

            def cb(step, total, label, results):
                self.parent.after(0, lambda: self._update_scan_plot(results))
                self.parent.after(0, lambda: self.status_var.set(
                    f"Scan {step}/{total} {label}"))

            results = perveance_erosion.run_scan(
                params, elecs, bemcs_grids,
                scan_type=sp["scan_type"].get().lower(),
                scan_electrode=int(sp["scan_elec"].get()),
                scan_min=float(sp["scan_min"].get()),
                scan_max=float(sp["scan_max"].get()),
                scan_steps=int(sp["scan_steps"].get()),
                plasma_params=plasma,
                erosion_steps=int(sp["erosion_steps"].get()),
                callback=cb)

            self.parent.after(0, lambda: self._scan_done(results))
        except Exception as e:
            self.parent.after(0, lambda: self._error(str(e)))

    def _update_scan_plot(self, results):
        self.pe_fig.clear()
        p, d = results["perveance"], results["divergence"]
        er, lb = results["erosion_rate"], results["labels"]
        if not p: return

        ax1 = self.pe_fig.add_subplot(111)
        c1, c2 = "#1565C0", "#E65100"
        ln1 = ax1.plot(p, d, "o-", color=c1, markersize=7, lw=2, label="Divergence")
        ax1.set_xlabel("Perveance ($\\mu$A/V$^{3/2}$)")
        ax1.set_ylabel("Divergence (mrad)", color=c1)
        ax1.tick_params(axis="y", labelcolor=c1)

        ax2 = ax1.twinx()
        ln2 = ax2.plot(p, er, "s--", color=c2, markersize=5, lw=1.5, label="Erosion rate")
        ax2.set_ylabel("Erosion rate (damage/step)", color=c2)
        ax2.tick_params(axis="y", labelcolor=c2)

        for i in range(len(p)):
            ax1.annotate(lb[i], (p[i], d[i]), textcoords="offset points",
                         xytext=(6, 6), fontsize=7, color="gray")

        if len(d) > 1:
            im = int(np.argmin(d))
            ax1.plot(p[im], d[im], "*", color="red", markersize=15, zorder=5,
                     label=f"Min: {d[im]:.2f}mrad")

        lns = ln1 + ln2
        ax1.legend(lns, [l.get_label() for l in lns], loc="upper left", fontsize=8)
        ax1.set_title("Perveance vs Divergence + Erosion Rate")
        ax1.grid(True, alpha=0.3)
        self.pe_fig.tight_layout(); self.pe_canvas.draw()

    def _scan_done(self, results):
        self.running = False
        self.twin_btn.config(state="normal"); self.scan_btn.config(state="normal")
        self.stop_btn.config(state="disabled")
        self.status_var.set(f"Scan done ({len(results['perveance'])} points)")
        self._update_scan_plot(results)

    def _error(self, msg):
        self.running = False
        self.twin_btn.config(state="normal"); self.scan_btn.config(state="normal")
        self.stop_btn.config(state="disabled")
        self.status_var.set("Error"); messagebox.showerror("Error", msg)
