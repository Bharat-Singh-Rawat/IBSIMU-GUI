"""
Main application window with 4 tabbed workflows.
"""
import tkinter as tk
from tkinter import ttk


class IonSourceTwinApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Ion Source Digital Twin")
        self.root.geometry("1600x1000")
        self.root.minsize(1200, 700)

        # Shared state between tabs
        self.shared = {
            "ibsimu_results": None,
            "beam_handoff": None,        # beam dict for PY-BEMCS
            "extraction_params": None,
            "extraction_electrodes": None,
            "transport_results": None,
        }

        self._build_ui()

    def _build_ui(self):
        # Main notebook with 4 workflow tabs
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Import tabs (deferred to avoid circular imports)
        from gui.tab_extraction import ExtractionTab
        from gui.tab_transport import TransportTab
        from gui.tab_erosion import ErosionTab
        from gui.tab_twin import DigitalTwinTab

        self.tab1_frame = ttk.Frame(self.notebook)
        self.tab2_frame = ttk.Frame(self.notebook)
        self.tab3_frame = ttk.Frame(self.notebook)
        self.tab4_frame = ttk.Frame(self.notebook)

        self.notebook.add(self.tab1_frame, text="1. Beam Extraction (IBSIMU)")
        self.notebook.add(self.tab2_frame, text="2. Transport & CEX (PY-BEMCS)")
        self.notebook.add(self.tab3_frame, text="3. Erosion & Thermal")
        self.notebook.add(self.tab4_frame, text="4. Digital Twin")

        self.tab1 = ExtractionTab(self.tab1_frame, self)
        self.tab2 = TransportTab(self.tab2_frame, self)
        self.tab3 = ErosionTab(self.tab3_frame, self)
        self.tab4 = DigitalTwinTab(self.tab4_frame, self)

    def switch_to_tab(self, index):
        self.notebook.select(index)
