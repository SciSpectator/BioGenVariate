# main.py
import os
import tkinter as tk
from gui.main_window import GeoWorkflowGUI
import config

def main():
    """Initializes and runs the GeoExplorer application."""
    if not os.path.exists(config.RESULTS_ROOT):
        os.makedirs(config.RESULTS_ROOT)
        print(f"Created results directory at: {config.RESULTS_ROOT}")

    app = GeoWorkflowGUI()
    app.mainloop()

if __name__ == "__main__":
    main()
