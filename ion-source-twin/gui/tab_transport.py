"""
Tab 2: PY-BEMCS Transport & CEX Simulation
Receives beam from Tab 1, runs PY-BEMCS, shows particle trajectories + CEX backflow.
"""
import threading
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from gui.plot_helpers import make_canvas, placeholder
from workflows import transport_cex, cex_backflow
from config import DEFAULT_TRANSPORT, DEFAULT_MATERIAL


class TransportTab:
    def __init__(self, parent, app):
        self.parent = parent
        self.app = app
        self.beam = None
        self.transport_result = None
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

        # Beam status
        ttk.Label(left, text="Beam from Extraction", font=("Helvetica",11,"bold")).pack(pady=(8,4))
        self.beam_status = tk.StringVar(value="No beam received yet")
        ttk.Label(left, textvariable=self.beam_status, foreground="gray",
                  wraplength=260).pack(padx=12)

        # Transport grid definition
        ttk.Separator(left).pack(fill=tk.X, padx=12, pady=8)
        ttk.Label(left, text="Transport Grids", font=("Helvetica",11,"bold")).pack(pady=(0,4))
        hdr = ttk.Frame(left); hdr.pack(fill=tk.X, padx=12)
        for t,w in [("#",3),("V(V)",7),("t(mm)",6),("gap",5),("r(mm)",6),("chm",5)]:
            ttk.Label(hdr, text=t, width=w, anchor="center",
                      font=("Helvetica",8,"bold")).pack(side=tk.LEFT, padx=1)
        self.grid_frame = ttk.Frame(left); self.grid_frame.pack(fill=tk.X, padx=12)
        self.grid_rows = []
        self._add_grid("1650","1.0","1.0","1.0","0")
        self._add_grid("-350","1.0","1.0","0.6","0")

        bf = ttk.Frame(left); bf.pack(fill=tk.X, padx=12, pady=4)
        ttk.Button(bf, text="+ Grid", command=lambda: self._add_grid()).pack(side=tk.LEFT, padx=2)
        ttk.Button(bf, text="- Grid", command=self._rm_grid).pack(side=tk.LEFT, padx=2)

        # Plasma params
        ttk.Separator(left).pack(fill=tk.X, padx=12, pady=8)
        ttk.Label(left, text="Plasma Parameters", font=("Helvetica",11,"bold")).pack(pady=(0,4))
        self.plasma = {}
        for l,k,d in [("Plasma n0 (m\u207b\u00b3)","n0_plasma","1e17"),
                       ("Te (eV)","Te_up","3.0"), ("Ti (eV)","Ti","2.0"),
                       ("Neutral n0 (m\u207b\u00b3)","n0","1e20"),
                       ("Steps","n_steps","500")]:
            f = ttk.Frame(left); f.pack(fill=tk.X, padx=12, pady=1)
            ttk.Label(f, text=l, width=18, anchor="w").pack(side=tk.LEFT)
            v = tk.StringVar(value=d); ttk.Entry(f, textvariable=v, width=10).pack(side=tk.LEFT)
            self.plasma[k] = v

        # Run + Stop
        ttk.Separator(left).pack(fill=tk.X, padx=12, pady=8)
        rbf = ttk.Frame(left); rbf.pack(fill=tk.X, padx=12)
        self.run_btn = ttk.Button(rbf, text="Run Transport", command=self._on_run)
        self.run_btn.pack(side=tk.LEFT, fill=tk.X, expand=True, ipady=5, padx=(0,2))
        self.stop_btn = ttk.Button(rbf, text="Stop", command=self._on_stop, state="disabled")
        self.stop_btn.pack(side=tk.LEFT, ipady=5)
        self.stop_flag = False
        self.status_var = tk.StringVar(value="Ready")
        ttk.Label(left, textvariable=self.status_var, foreground="gray").pack(padx=12, pady=2)
        self.progress = ttk.Progressbar(left, mode="indeterminate")
        self.progress.pack(fill=tk.X, padx=12, pady=2)

        # Backflow report
        ttk.Separator(left).pack(fill=tk.X, padx=12, pady=8)
        ttk.Label(left, text="CEX Backflow", font=("Helvetica",11,"bold")).pack(pady=(0,4))
        self.bf_var = tk.StringVar(value="Run transport first")
        ttk.Label(left, textvariable=self.bf_var, justify=tk.LEFT,
                  wraplength=260, font=("Courier",9)).pack(padx=12, anchor="w")

        # --- RIGHT: plots ---
        right = ttk.Frame(pw); pw.add(right, weight=1)
        rpw = ttk.PanedWindow(right, orient=tk.VERTICAL)
        rpw.pack(fill=tk.BOTH, expand=True)

        # Particle scatter
        tf = ttk.LabelFrame(rpw, text="Particle Positions")
        rpw.add(tf, weight=1)
        self.part_fig, self.part_ax, self.part_canvas = make_canvas(tf, figsize=(10,3.5))
        placeholder(self.part_ax, "Run transport to see particles")

        # Bottom: divergence history + backflow energy
        botf = ttk.Frame(rpw); rpw.add(botf, weight=1)
        self.div_fig, self.div_ax, self.div_canvas = make_canvas(botf, figsize=(5,3))
        placeholder(self.div_ax, "Divergence history")

    # -- grid helpers --
    def _add_grid(self, v="0", t="1.0", gap="1.0", r="1.0", cham="0"):
        idx = len(self.grid_rows)+1
        row = ttk.Frame(self.grid_frame); row.pack(fill=tk.X, pady=1)
        ttk.Label(row, text=str(idx), width=3, anchor="center").pack(side=tk.LEFT, padx=1)
        vs = {}
        for k,d,w in [("V",v,7),("t",t,6),("gap",gap,5),("r",r,6),("cham",cham,5)]:
            sv = tk.StringVar(value=d); ttk.Entry(row, textvariable=sv, width=w).pack(side=tk.LEFT, padx=1)
            vs[k] = sv
        self.grid_rows.append({"frame":row, "vars":vs})

    def _rm_grid(self):
        if len(self.grid_rows) <= 1: return
        self.grid_rows.pop()["frame"].destroy()

    def _get_grids(self):
        return [{k: float(v.get()) for k,v in r["vars"].items()} for r in self.grid_rows]

    def _get_plasma(self):
        p = {}
        for k,v in self.plasma.items():
            try: p[k] = float(v.get())
            except: p[k] = DEFAULT_TRANSPORT.get(k, 0)
        return p

    def receive_beam(self, beam):
        """Called by Tab 1 when user clicks 'Send Beam'."""
        self.beam = beam
        n = beam.get("n_particles", 0)
        cur = beam.get("current_A", 0)
        self.beam_status.set(f"Received {n} particles, I={cur:.4e} A")

    # -- run --
    def _on_stop(self):
        self.stop_flag = True

    def _on_run(self):
        if self.beam is None:
            messagebox.showinfo("No beam", "Go to Tab 1, run extraction, then click 'Send Beam'")
            return
        self.stop_flag = False
        self.run_btn.config(state="disabled"); self.stop_btn.config(state="normal")
        self.status_var.set("Running transport..."); self.progress.start(15)
        threading.Thread(target=self._run_thread, daemon=True).start()

    def _run_thread(self):
        try:
            grids = self._get_grids()
            plasma = self._get_plasma()
            n_steps = int(self.plasma["n_steps"].get())
            params = self.app.shared.get("extraction_params", {})
            mass = float(params.get("beam_mass", 1.0))
            charge = int(abs(float(params.get("beam_charge", 1))))

            # Build simulator once
            from adapters import pybemcs_adapter
            sim = pybemcs_adapter.create_simulator(mass, charge, DEFAULT_MATERIAL)
            p = pybemcs_adapter.build_domain(sim, grids, plasma)

            # Use default PY-BEMCS plasma injection (it handles current self-consistently)
            # AND inject the IBSIMU beam particles every N steps to seed the distribution
            beam = self.beam
            inject_every = 5  # inject IBSIMU beam particles every 5 steps

            # Store isBound for geometry overlay
            self._isBound = sim.isBound.copy()
            self._sim_dx = sim.dx
            self._sim_dy = sim.dy

            history = {"div": [], "min_pot": [], "t_grids": [], "remeshed": []}
            batch_size = 10  # update plots every 10 steps

            for step in range(n_steps):
                if self.stop_flag:
                    break

                # Inject IBSIMU beam particles periodically
                if beam is not None and step % inject_every == 0:
                    pybemcs_adapter.inject_beam(sim, beam)

                remeshed, min_pot, div, t_grids = sim.step(p)
                history["div"].append(div)
                history["min_pot"].append(min_pot)
                history["t_grids"].append(list(t_grids) if t_grids else [])
                history["remeshed"].append(remeshed)

                # Live plot update every batch_size steps
                if step % batch_size == 0 or step == n_steps - 1:
                    import threading
                    event = threading.Event()

                    def _update(s=step, d=div):
                        self.status_var.set(f"Step {s+1}/{n_steps}  div={d:.2f}\u00b0")
                        self._live_update(sim, history)
                        event.set()

                    self.parent.after(0, _update)
                    event.wait(timeout=2.0)  # wait for GUI to process

            # Final results
            result = {
                "history": history,
                "backflow": pybemcs_adapter.extract_cex_backflow(sim),
                "damage_map": pybemcs_adapter.get_damage_map(sim),
                "T_map": pybemcs_adapter.get_temperature_map(sim),
                "sim": sim,
                "grids": grids,
                "params": p,
            }
            self.transport_result = result
            self.app.shared["transport_results"] = result
            bf = cex_backflow.analyze(sim)
            self.parent.after(0, lambda: self._done(result, bf))
        except Exception as e:
            self.parent.after(0, lambda: self._error(str(e)))

    def _live_update(self, sim, history):
        """Update particle scatter and divergence plots mid-simulation."""
        # Particle scatter with grid geometry overlay
        self.part_ax.clear()

        # Draw grid geometry (isBound mask) as filled region
        if hasattr(self, '_isBound') and self._isBound is not None:
            ib = self._isBound
            ny, nx = ib.shape
            x_coords = np.arange(nx) * self._sim_dx
            y_coords = np.arange(ny) * self._sim_dy
            self.part_ax.contourf(x_coords, y_coords, ib.astype(float),
                                  levels=[0.5, 1.5], colors=["#555555"], alpha=0.6)
            # Also show as outline
            self.part_ax.contour(x_coords, y_coords, ib.astype(float),
                                 levels=[0.5], colors=["black"], linewidths=1.5)

        # Scatter particles
        if sim.num_p > 0:
            x = sim.p_x[:sim.num_p]; y = sim.p_y[:sim.num_p]
            cex_mask = sim.p_isCEX[:sim.num_p]
            self.part_ax.scatter(x[~cex_mask], y[~cex_mask], s=1, c="#1565C0",
                                 alpha=0.5, label=f"Primary ({int(np.sum(~cex_mask))})")
            if np.any(cex_mask):
                self.part_ax.scatter(x[cex_mask], y[cex_mask], s=2, c="#E65100",
                                     alpha=0.7, label=f"CEX ({int(np.sum(cex_mask))})")

        # Electrons
        if hasattr(sim, 'num_e') and sim.num_e > 0:
            self.part_ax.scatter(sim.e_x[:sim.num_e], sim.e_y[:sim.num_e],
                                 s=0.5, c="#4CAF50", alpha=0.3, label=f"e- ({sim.num_e})")

        self.part_ax.set_xlabel("x (mm)"); self.part_ax.set_ylabel("y (mm)")
        self.part_ax.set_title(f"Particles (n={sim.num_p})")
        self.part_ax.legend(fontsize=7, loc="upper right")
        self.part_ax.grid(True, alpha=0.3)
        self.part_fig.tight_layout(); self.part_canvas.draw()

        # Divergence history
        self.div_ax.clear()
        if history["div"]:
            self.div_ax.plot(history["div"], color="#1565C0", lw=1.5)
            self.div_ax.set_xlabel("Iteration"); self.div_ax.set_ylabel("Divergence (\u00b0)")
            self.div_ax.set_title("Beam Divergence vs Iteration")
            self.div_ax.grid(True, alpha=0.3)
        self.div_fig.tight_layout(); self.div_canvas.draw()

    def _error(self, msg):
        self.progress.stop(); self.run_btn.config(state="normal")
        self.stop_btn.config(state="disabled")
        self.status_var.set("Error"); messagebox.showerror("Error", msg)

    def _done(self, result, bf):
        self.progress.stop(); self.run_btn.config(state="normal")
        self.stop_btn.config(state="disabled")
        self.status_var.set("Done" if not self.stop_flag else "Stopped")
        self.bf_var.set(cex_backflow.report(bf))

        # Final particle scatter (reuse live update)
        self._live_update(result["sim"], result["history"])

        # Divergence history
        self.div_ax.clear()
        hist = result["history"]
        if hist["div"]:
            self.div_ax.plot(hist["div"], color="#1565C0", lw=1.5)
            self.div_ax.set_xlabel("Iteration")
            self.div_ax.set_ylabel("Divergence (\u00b0)")
            self.div_ax.set_title("Beam Divergence vs Iteration")
            self.div_ax.grid(True, alpha=0.3)
        self.div_fig.tight_layout(); self.div_canvas.draw()

        # Pass to erosion tab
        self.app.tab3.update_results(result)
