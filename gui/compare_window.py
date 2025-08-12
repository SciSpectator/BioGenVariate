# gui/compare_window.py
import tkinter as tk
from tkinter import ttk, filedialog, simpledialog, messagebox
from tkinter.constants import TOP, BOTH, RIGHT, Y, LEFT, X, END, W, CENTER, NSEW
from pathlib import Path
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from scipy.stats import ranksums

from utils.statistics import analyze_gene_distribution

class CompareDistributionsWindow(tk.Toplevel):
    """
    This class creates a pop-up window for comparing gene expression distributions.
    It allows users to load their own case/control sample lists, select groups,
    and visualize their distributions against each other or against entire platforms.
    """
    def __init__(self, parent, app_ref):
        """
        Initializes the Compare Distributions window.

        Args:
            parent: The main Tkinter window that this pop-up belongs to.
            app_ref: A reference to the main application instance (GeoWorkflowGUI)
                     to access shared data like loaded GPL datasets.
        """
        super().__init__(parent)
        self.parent = parent
        self.app_ref = app_ref
        self.title("Compare Gene Expression Distributions")
        self.geometry("1200x850")

        # --- Instance Variables ---
        self.user_defined_groups = {}
        self.fig_density, self.ax_density, self.canvas_density, self.toolbar_density = None, None, None, None
        self.fig_hist, self.ax_hist, self.canvas_hist, self.toolbar_hist = None, None, None, None

        self._setup_ui()
        self.protocol("WM_DELETE_WINDOW", self._on_closing)

    def _setup_ui(self):
        """
        Creates and arranges all the widgets (buttons, lists, frames) in the window.
        """
        main_frame = ttk.Frame(self, padding="10")
        main_frame.pack(fill=BOTH, expand=True)
        main_frame.grid_columnconfigure(0, weight=2, minsize=400)
        main_frame.grid_columnconfigure(1, weight=5)
        main_frame.grid_rowconfigure(0, weight=1)

        control_pane = ttk.Frame(main_frame, padding=5)
        control_pane.grid(row=0, column=0, sticky=NSEW)
        control_pane.grid_rowconfigure(5, weight=1)
        control_pane.grid_columnconfigure(0, weight=1)

        data_loading_frame = ttk.LabelFrame(control_pane, text="1. Load Case-Control File(s)")
        data_loading_frame.grid(row=0, column=0, sticky="ew", padx=5, pady=5)
        data_loading_frame.grid_columnconfigure(0, weight=1)
        
        load_button_frame = ttk.Frame(data_loading_frame)
        load_button_frame.pack(fill=X, padx=5, pady=5)
        ttk.Button(load_button_frame, text="Load File(s)", command=self._load_case_control_file).pack(side=LEFT, fill=X, expand=True)
        ttk.Button(load_button_frame, text="Clear All", command=self._clear_user_data).pack(side=LEFT, padx=5)

        ttk.Label(control_pane, text="2. Select Groups to Compare", font=("Helvetica", 10, "bold")).grid(row=1, column=0, sticky=W, padx=5, pady=(10, 2))
        self.loaded_files_listbox = tk.Listbox(control_pane, height=6, selectmode=tk.EXTENDED)
        self.loaded_files_listbox.grid(row=2, column=0, sticky="ew", padx=5, pady=2)

        mode_frame = ttk.LabelFrame(control_pane, text="3. Choose Comparison Mode")
        mode_frame.grid(row=3, column=0, sticky="ew", padx=5, pady=10)
        self.comparison_mode = tk.StringVar(value="groups_only")
        ttk.Radiobutton(mode_frame, text="Compare Selected Case/Control Groups Only", variable=self.comparison_mode, value="groups_only").pack(anchor=W, padx=10)
        ttk.Radiobutton(mode_frame, text="Compare Groups vs. Gene(s) on Platform(s)", variable=self.comparison_mode, value="vs_gene").pack(anchor=W, padx=10)
        ttk.Radiobutton(mode_frame, text="Compare Groups vs. Entire Platform(s)", variable=self.comparison_mode, value="vs_platform").pack(anchor=W, padx=10)

        spec_frame = ttk.LabelFrame(control_pane, text="4. Specify Gene(s) and/or Platform(s)")
        spec_frame.grid(row=4, column=0, sticky="ew", padx=5, pady=5)

        ttk.Label(spec_frame, text="Gene Symbol(s) (comma-separated):").pack(anchor=W, padx=5, pady=(5,0))
        self.gene_entry = ttk.Entry(spec_frame)
        self.gene_entry.pack(fill=X, padx=5, pady=2)
        
        ttk.Label(spec_frame, text="Comparison Platform(s) (for modes 2 & 3):").pack(anchor=W, padx=5, pady=(5,0))
        self.platform_frame = ttk.Frame(spec_frame)
        self.platform_frame.pack(fill=X, padx=5, pady=2)
        self.platform_vars = {}
        if self.app_ref.gpl_datasets:
            for p in sorted(self.app_ref.gpl_datasets.keys()):
                var = tk.BooleanVar(value=False)
                ttk.Checkbutton(self.platform_frame, text=p, variable=var).pack(side=LEFT)
                self.platform_vars[p] = var

        ttk.Button(control_pane, text="Plot & Analyze Distributions", command=self._plot_and_analyze).grid(row=5, column=0, sticky="sew", padx=5, pady=15, ipady=5)
        self.status_label = ttk.Label(control_pane, text="Ready.", foreground="blue")
        self.status_label.grid(row=6, column=0, sticky=W, padx=5, pady=5)

        output_pane = ttk.Frame(main_frame)
        output_pane.grid(row=0, column=1, sticky=NSEW, padx=5)
        output_pane.grid_rowconfigure(0, weight=3); output_pane.grid_rowconfigure(1, weight=2)
        output_pane.grid_columnconfigure(0, weight=1); output_pane.grid_columnconfigure(1, weight=1)

        self.density_plot_frame = ttk.LabelFrame(output_pane, text="Overlaid Density Plots")
        self.density_plot_frame.grid(row=0, column=0, sticky=NSEW, padx=2, pady=5)
        self.hist_plot_frame = ttk.LabelFrame(output_pane, text="Overlaid Histograms & Classification")
        self.hist_plot_frame.grid(row=0, column=1, sticky=NSEW, padx=2, pady=5)
        self.stats_frame = ttk.LabelFrame(output_pane, text="Statistical Test Results (Pairwise Wilcoxon rank-sum)")
        self.stats_frame.grid(row=1, column=0, columnspan=2, sticky=NSEW, padx=2, pady=5)

    def _clear_user_data(self):
        """Resets the loaded user data and the listbox."""
        self.user_defined_groups.clear()
        self.loaded_files_listbox.delete(0, END)
        self.status_label.config(text="Cleared all loaded files.")

    def _on_closing(self):
        """Handles the window close event to prevent memory leaks from plots."""
        self._clear_figures()
        self.destroy()

    def _clear_figures(self):
        """Safely destroys any existing Matplotlib figures and canvases."""
        for fig, canvas, toolbar in [(self.fig_density, self.canvas_density, self.toolbar_density), (self.fig_hist, self.canvas_hist, self.toolbar_hist)]:
            if fig: plt.close(fig)
            if canvas: canvas.get_tk_widget().destroy()
            if toolbar: toolbar.destroy()
        self.fig_density, self.canvas_density, self.toolbar_density = None, None, None, None
        self.fig_hist, self.canvas_hist, self.toolbar_hist = None, None, None, None

    def _clear_stats_table(self):
        """Removes all widgets from the statistics frame to prepare for new results."""
        for widget in self.stats_frame.winfo_children():
            widget.destroy()

    def _ask_for_platform(self, platforms: list) -> str or None:
        """Creates a simple dialog to force the user to associate loaded samples with a platform."""
        dialog = tk.Toplevel(self)
        dialog.title("Select Platform")
        dialog.transient(self); dialog.grab_set()
        ttk.Label(dialog, text="Select the platform for the loaded samples:", wraplength=300).pack(pady=10)
        var = tk.StringVar(value=platforms[0])
        for p in platforms:
            ttk.Radiobutton(dialog, text=p, variable=var, value=p).pack(anchor=W, padx=30)
        
        result = {"platform": None}
        def on_ok():
            result["platform"] = var.get()
            dialog.destroy()
        
        ttk.Button(dialog, text="OK", command=on_ok).pack(pady=10, ipadx=20)
        self.wait_window(dialog)
        return result["platform"]

    def _load_case_control_file(self):
        """
        Opens a file dialog for the user to select one or more CSV files.
        Each file should contain sample IDs and their case/control status (1/0).
        """
        filepaths = filedialog.askopenfilenames(
            title="Select Case-Control CSV or gzipped CSV file(s)",
            filetypes=[("CSV files", "*.csv"), ("Gzipped CSV files", "*.csv.gz")]
        )
        if not filepaths: return

        for filepath in filepaths:
            try:
                df_annot = pd.read_csv(filepath, compression='gzip' if filepath.endswith('.gz') else None, low_memory=False)
                if df_annot.empty: continue

                columns = df_annot.columns.tolist()
                formatted_columns = '\n'.join(columns)
                gsm_col = simpledialog.askstring("GSM Column", f"For file '{Path(filepath).name}', enter column name for GSM IDs.\n\nAvailable:\n{formatted_columns}", parent=self)
                if not gsm_col or gsm_col not in columns: continue

                cc_col = simpledialog.askstring("Case/Control Column", f"For file '{Path(filepath).name}', enter column for labels (0=Control, 1=Case).\n\nAvailable:\n{formatted_columns}", parent=self)
                if not cc_col or cc_col not in columns: continue

                loaded_platforms = list(self.app_ref.gpl_datasets.keys())
                if not loaded_platforms:
                    messagebox.showerror("Error", "No platforms are loaded in the main app.", parent=self)
                    return
                
                platform = self._ask_for_platform(loaded_platforms)
                if not platform: continue

                df_annot.dropna(subset=[gsm_col, cc_col], inplace=True)
                df_annot[gsm_col] = df_annot[gsm_col].astype(str).str.upper()
                df_annot[cc_col] = pd.to_numeric(df_annot[cc_col], errors='coerce')
                
                if not set(df_annot[cc_col].unique()).issubset({0.0, 1.0}):
                    messagebox.showerror("Data Error", f"File '{Path(filepath).name}': The case/control column '{cc_col}' must contain only 0s and 1s.", parent=self)
                    continue

                case_gsms = df_annot[df_annot[cc_col] == 1.0][gsm_col].tolist()
                control_gsms = df_annot[df_annot[cc_col] == 0.0][gsm_col].tolist()
                
                filename = Path(filepath).name
                if case_gsms:
                    label = f"{filename} (Case)"
                    self.user_defined_groups[label] = {"platform": platform, "gsms": case_gsms}
                    self.loaded_files_listbox.insert(END, label)
                if control_gsms:
                    label = f"{filename} (Control)"
                    self.user_defined_groups[label] = {"platform": platform, "gsms": control_gsms}
                    self.loaded_files_listbox.insert(END, label)
            except Exception as e:
                messagebox.showerror("Error", f"Failed to process file {Path(filepath).name}:\n{e}", parent=self)

    def _plot_and_analyze(self):
        """
        Main logic function that collects all user inputs, fetches the relevant
        expression data, and generates the plots and statistical analyses.
        """
        self._clear_figures()
        self._clear_stats_table()
        
        data_to_plot = {}
        
        selected_indices = self.loaded_files_listbox.curselection()
        selected_groups = [self.loaded_files_listbox.get(i) for i in selected_indices]
        mode = self.comparison_mode.get()
        gene_input_str = self.gene_entry.get().strip()
        genes_to_plot = [g.strip().upper() for g in gene_input_str.split(',') if g.strip()]
        selected_platforms = [p for p, var in self.platform_vars.items() if var.get()]

        # --- Data Gathering Logic ---
        if selected_groups:
            if mode == "vs_gene" and not genes_to_plot:
                messagebox.showerror("Input Required", "Please enter gene(s) to compare your groups against platforms.", parent=self)
                return
            if mode != "groups_only" and not selected_platforms:
                messagebox.showerror("Input Required", f"Please select at least one Comparison Platform for the chosen mode.", parent=self)
                return

            # 1a. Add data for selected user-defined groups
            for label in selected_groups:
                group_info = self.user_defined_groups.get(label)
                if not group_info: continue
                platform, gsms = group_info["platform"], group_info["gsms"]
                df_plat = self.app_ref.gpl_datasets.get(platform)
                gmap_plat = self.app_ref.gpl_gene_mappings.get(platform, {})
                if df_plat is None or not gmap_plat: continue
                
                if genes_to_plot:
                    for gene in genes_to_plot:
                        gene_col = gmap_plat.get(gene)
                        if gene_col and gene_col in df_plat.columns:
                            subset_df = df_plat[df_plat['GSM'].isin(gsms)]
                            expr = pd.to_numeric(subset_df[gene_col], errors='coerce').dropna()
                            if not expr.empty: data_to_plot[f"{gene} - {label}"] = expr
                else:
                    all_gene_cols = [col for col in gmap_plat.values() if col in df_plat.columns]
                    subset_df = df_plat[df_plat['GSM'].isin(gsms)]
                    all_expr = pd.concat([pd.to_numeric(subset_df[col], errors='coerce').dropna() for col in all_gene_cols], ignore_index=True)
                    if not all_expr.empty: data_to_plot[label] = all_expr
            
            # 1b. Add data for comparison platforms based on mode
            if mode == "vs_gene":
                for platform in selected_platforms:
                    df, gmap = self.app_ref.gpl_datasets.get(platform), self.app_ref.gpl_gene_mappings.get(platform, {})
                    if df is None or not gmap: continue
                    for gene in genes_to_plot:
                        gene_col = gmap.get(gene)
                        if gene_col and gene_col in df.columns:
                            expr = pd.to_numeric(df[gene_col], errors='coerce').dropna()
                            if not expr.empty: data_to_plot[f"{gene} ({platform})"] = expr
            elif mode == "vs_platform":
                for platform in selected_platforms:
                    df, gmap = self.app_ref.gpl_datasets.get(platform), self.app_ref.gpl_gene_mappings.get(platform, {})
                    if df is None or not gmap: continue
                    all_expr = pd.concat([pd.to_numeric(df[col], errors='coerce').dropna() for col in gmap.values() if col in df.columns], ignore_index=True)
                    if not all_expr.empty: data_to_plot[f"All Genes ({platform})"] = all_expr

        else: # Case 2: No groups selected, direct platform/gene comparison
            if not genes_to_plot or not selected_platforms:
                self.status_label.config(text="No valid data to plot. Select genes and platforms.")
                return

            for platform in selected_platforms:
                df, gmap = self.app_ref.gpl_datasets.get(platform), self.app_ref.gpl_gene_mappings.get(platform, {})
                if df is None or not gmap: continue
                for gene in genes_to_plot:
                    gene_col = gmap.get(gene)
                    if gene_col and gene_col in df.columns:
                        expr = pd.to_numeric(df[gene_col], errors='coerce').dropna()
                        if not expr.empty:
                            data_to_plot[f"{gene} ({platform})"] = expr
        
        if not data_to_plot:
            self.status_label.config(text="No valid data to plot based on selections.")
            return

        self.status_label.config(text=f"Plotting {len(data_to_plot)} distributions...")
        
        # --- Plotting ---
        self.fig_density, self.ax_density = plt.subplots(figsize=(7, 5))
        for label, data in data_to_plot.items():
            sns.kdeplot(data, ax=self.ax_density, label=label, fill=True, alpha=0.4)
        self.ax_density.set_title("Overlaid Density", fontsize=12)
        self.ax_density.legend(fontsize='x-small'); self.fig_density.tight_layout()
        color_map = {text.get_text(): line.get_color() for text, line in zip(self.ax_density.get_legend().get_texts(), self.ax_density.get_lines())}
        
        self.canvas_density = FigureCanvasTkAgg(self.fig_density, master=self.density_plot_frame)
        self.canvas_density.draw()
        self.canvas_density.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
        self.toolbar_density = NavigationToolbar2Tk(self.canvas_density, self.density_plot_frame)
        self.toolbar_density.update()

        self.fig_hist, self.ax_hist = plt.subplots(figsize=(7, 5))
        if data_to_plot:
            all_values = pd.concat(data_to_plot.values())
            bins = np.linspace(all_values.min(), all_values.max(), 50)
            classifications = {label: analyze_gene_distribution(series) for label, series in data_to_plot.items()}
            for label, data in data_to_plot.items():
                self.ax_hist.hist(data, bins=bins, color=color_map.get(label), label=f"{label} ({classifications.get(label, 'N/A')})", alpha=0.6)
        
        self.ax_hist.set_title("Overlaid Histograms", fontsize=12)
        self.ax_hist.set_xlabel("Expression"); self.ax_hist.set_ylabel("Frequency")
        self.ax_hist.legend(fontsize='x-small'); self.fig_hist.tight_layout()
        
        self.canvas_hist = FigureCanvasTkAgg(self.fig_hist, master=self.hist_plot_frame)
        self.canvas_hist.draw()
        self.canvas_hist.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
        self.toolbar_hist = NavigationToolbar2Tk(self.canvas_hist, self.hist_plot_frame)
        self.toolbar_hist.update()

        self._run_stats_analysis(data_to_plot)
        self.status_label.config(text="Analysis complete.")

    def _run_stats_analysis(self, data_for_stats):
        """
        Performs a pairwise Wilcoxon rank-sum test on all combinations of the
        plotted distributions and displays the results in a table.
        """
        if len(data_for_stats) < 2: return
        self._clear_stats_table()
        
        columns = ["Distribution 1", "Distribution 2", "Wilcoxon Score (U)"]
        stats_tree = ttk.Treeview(self.stats_frame, columns=columns, show="headings")
        stats_tree.pack(side=LEFT, fill=BOTH, expand=True, padx=5, pady=5)
        
        vsb = ttk.Scrollbar(self.stats_frame, orient="vertical", command=stats_tree.yview)
        vsb.pack(side=RIGHT, fill=Y)
        stats_tree.configure(yscrollcommand=vsb.set)

        stats_tree.column("Distribution 1", width=300, anchor=W)
        stats_tree.column("Distribution 2", width=300, anchor=W)
        stats_tree.column("Wilcoxon Score (U)", width=150, anchor=CENTER)
        for col in columns:
            stats_tree.heading(col, text=col)

        for combo in itertools.combinations(data_for_stats.keys(), 2):
            label1, label2 = combo
            data1, data2 = data_for_stats[label1], data_for_stats[label2]
            try:
                if len(data1) > 1 and len(data2) > 1:
                    stat, _ = ranksums(data1, data2)
                    row_data = [label1, label2, f"{stat:.2f}"]
                else:
                    row_data = [label1, label2, "N/A"]
            except ValueError:
                row_data = [label1, label2, "Error"]
            stats_tree.insert("", "end", values=row_data)
