"""
Tab 3: Erosion & Thermal Results
Displays damage map, temperature map, grid temperature history.
"""
import tkinter as tk
from tkinter import ttk
import numpy as np
from gui.plot_helpers import make_canvas, placeholder


class ErosionTab:
    def __init__(self, parent, app):
        self.parent = parent
        self.app = app
        self._build_ui()

    def _build_ui(self):
        pw = ttk.PanedWindow(self.parent, orient=tk.VERTICAL)
        pw.pack(fill=tk.BOTH, expand=True)

        # Top: damage map + temperature map side by side
        topf = ttk.Frame(pw); pw.add(topf, weight=1)

        lf = ttk.LabelFrame(topf, text="Erosion Damage Map")
        lf.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.dmg_fig, self.dmg_ax, self.dmg_canvas = make_canvas(lf, figsize=(6, 3.5))
        placeholder(self.dmg_ax, "Run transport first", "Damage Map")

        rf = ttk.LabelFrame(topf, text="Temperature Map")
        rf.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.temp_fig, self.temp_ax, self.temp_canvas = make_canvas(rf, figsize=(6, 3.5))
        placeholder(self.temp_ax, "Run transport first", "Temperature (\u00b0C)")

        # Bottom: grid temperature history
        botf = ttk.LabelFrame(pw, text="Grid Temperature History")
        pw.add(botf, weight=1)
        self.hist_fig, self.hist_ax, self.hist_canvas = make_canvas(botf, figsize=(12, 3))
        placeholder(self.hist_ax, "Run transport first", "Grid Temperatures vs Iteration")

    def update_results(self, result):
        """Called by Tab 2 after transport completes."""
        # Damage map
        self.dmg_ax.clear()
        dmg = result.get("damage_map")
        if dmg is not None and dmg.size > 1:
            im = self.dmg_ax.imshow(dmg, aspect="auto", origin="lower", cmap="hot")
            self.dmg_fig.colorbar(im, ax=self.dmg_ax, label="Damage")
            self.dmg_ax.set_xlabel("x (cells)"); self.dmg_ax.set_ylabel("y (cells)")
            self.dmg_ax.set_title("Cumulative Erosion Damage")
        else:
            placeholder(self.dmg_ax, "No damage data")
        self.dmg_fig.tight_layout(); self.dmg_canvas.draw()

        # Temperature map
        self.temp_ax.clear()
        tmap = result.get("T_map")
        if tmap is not None and tmap.size > 1:
            im = self.temp_ax.imshow(tmap - 273.15, aspect="auto", origin="lower", cmap="inferno")
            self.temp_fig.colorbar(im, ax=self.temp_ax, label="\u00b0C")
            self.temp_ax.set_xlabel("x (cells)"); self.temp_ax.set_ylabel("y (cells)")
            self.temp_ax.set_title("Temperature Map (\u00b0C)")
        else:
            placeholder(self.temp_ax, "No temperature data")
        self.temp_fig.tight_layout(); self.temp_canvas.draw()

        # Grid temperature history
        self.hist_ax.clear()
        hist = result.get("history", {})
        tg_hist = hist.get("t_grids", [])
        if tg_hist and len(tg_hist[0]) > 0:
            n_grids = len(tg_hist[0])
            for gi in range(n_grids):
                temps = [t[gi] - 273.15 if gi < len(t) else 0 for t in tg_hist]
                self.hist_ax.plot(temps, label=f"Grid {gi+1}", lw=2)
            self.hist_ax.set_xlabel("Iteration")
            self.hist_ax.set_ylabel("Temperature (\u00b0C)")
            self.hist_ax.set_title("Grid Temperature vs Iteration")
            self.hist_ax.legend(fontsize=9); self.hist_ax.grid(True, alpha=0.3)
        else:
            placeholder(self.hist_ax, "No temperature history")
        self.hist_fig.tight_layout(); self.hist_canvas.draw()
