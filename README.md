# BioGenVariate 1.0 

<img width="2048" height="2048" alt="Gemini_Generated_Image_liw558liw558liw5" src="https://github.com/user-attachments/assets/e2061a0f-d039-4e6f-943c-5bce38670a65" />


GenVariate is a sophisticated, GUI-driven Python application designed for the in-depth analysis of gene expression data from the Gene Expression Omnibus (GEO) and custom user datasets. It provides a powerful, multi-step workflow to extract, filter, label, and analyze biological samples. The application moves beyond simple keyword searches by integrating a local AI agent, powered by the gpt-oss:20b open-source model, for automated metadata classification and features a unique, interactive tool for deep-dive analysis of specific gene expression distributions.

The suite is engineered for researchers who need to identify relevant experimental cohorts from the vast GEO database, classify samples using either AI-assisted or manual methods, and explore the descriptive metadata associated with unique gene expression patterns. It excels at analyzing user-defined regions of interest within a gene's expression distribution (such as the high or low-expression tails) and generating insightful, publication-ready visualizations and data tables.


# 🔬 Core Features

Comprehensive Experiment Extraction

Utilizes the entire GEOmetadb.sqlite database, loaded into memory for rapid, multi-threaded querying. It performs a two-stage filtering process based on user-defined keywords to search metadata across both entire experiments (GSE) and individual samples (GSM).


🧬 Multi-Species & Custom Dataset Compatibility


While GenVariate includes quick-load buttons for common human and mouse platforms from GEO, its capabilities extend far beyond. The application is designed for flexibility, allowing users to load any custom dataset from a local file. This enables the analysis and comparison of data from different microarray or sequencing platforms, other species, or even private, non-public experimental results.


🤖 Dual-Mode Sample Labeling


Offers two distinct workflows for sample classification:

- AI Agent-Powered Labeling: Leverages a local AI agent via the Ollama service, running the gpt-oss:20b open-source model. The agent analyzes sample metadata to automatically extract structured information, including Condition, Tissue, Age, and Treatment with Treatment Time.


- Manual Labeling: An interactive dialog allows for sample-by-sample classification with customizable labels, ensuring accuracy and control over the final dataset.
    

📊 Interactive Gene Distribution Explorer


A powerful visualization tool to plot and investigate expression distributions for any gene across multiple platforms. Key features include:

- Statistical Classification: Automatically analyzes and classifies the shape of each gene's distribution (e.g., Normal, Bimodal, Lognormal).

- Interactive Range Selection: Users can draw a rectangle directly on any histogram to select a specific expression range for a deep-dive analysis.
    

🎯 AI-Driven Analysis of Selected Ranges


This is the cornerstone feature of GenVariate. Selecting an expression range on a plot triggers a targeted analysis on the samples within that subset. The AI agent labeling pipeline is run specifically on these samples to provide immediate biological context to their unique expression patterns.


⚖️ Advanced Distribution Comparison Tool


A dedicated window for comparing expression profiles between different groups. It supports loading custom case/control lists, plotting them against platform backgrounds, and performing pairwise Wilcoxon rank-sum statistical tests to quantify differences.


📂 Structured & Reproducible Results


All analyses generate a comprehensive set of plots, graphs, and data tables. Results are automatically saved into a logically structured, timestamped directory. Interactive data tables provide direct web links to NCBI GEO, and all figures can be exported in multiple formats (PNG, PDF, SVG).




# Installation & Setup ⚙️




# 1.Tutorial:

Using GenVariate is a straightforward process that guides you from a broad research question to specific, actionable insights. The workflow is divided into a few key phases.


# 1.1 Initial Setup and Sample Extraction


# Step 1: Extract Relevant Experiments


This initial step is designed to find relevant studies and samples for your research. Instead of querying the GEO website directly, GenVariate performs a rapid search of a local GEOmetadb.sqlite.gz file. This file is a comprehensive database containing all the metadata from GEO, allowing for powerful and fast offline filtering.

In the main application window, locate the "Step 1: GSE Extraction" section. Configure your search using the following parameters:

    Platform Filter (optional): Enter a comma-separated list of platform IDs (e.g., GPL570,GPL96) to restrict your search to specific microarray or sequencing technologies. Leave this blank to search across all platforms.

    Filtering Tokens: This is the core of your search. Enter a comma-separated list of keywords that describe your research interest (e.g., alzheimer disease,control,brain). The tool will search for these terms in both the overall study (GSE) descriptions and the detailed individual sample (GSM) metadata.

Once you have configured your search, click "Run GSE Extraction". The process will begin in the background, and you can monitor its progress via the progress bar and the log output window.


# Step 1.5: Interactively Review and Refine Your Results


Upon completion, a new "Step 1.5: Review, Analyze, and Select GSEs" window will automatically open. This powerful interface is designed to help you quickly validate the relevance of each study (GSE) found by the tool.

The window is split into two panes:

    Left Pane (GSE List): This is a table of all studies that contained at least one of your keywords. It provides a concise summary, including the platform, the specific "Matching Tokens" found, and a count of samples that were "Direct Matches" versus those that are just associated with the study. This allows you to immediately gauge the relevance of each result. You can select one or more GSEs from this list to inspect and keep.

    Right Pane (Details and Keyword Highlighting): When you click on a GSE in the left pane, this area populates with the combined titles, descriptions, and characteristics from its samples. To make validation effortless, all of your original search tokens are highlighted, so you can instantly see them in their proper context.

Review the studies and select all the ones you wish to include in your final dataset by clicking on them in the left pane (use Ctrl-Click or Shift-Click to select multiple). When you are finished, click the "Save Kept GSEs and Continue to Step 2" button.

Saving and Finalizing Your Selection

GenVariate automatically organizes your results for you.

    Directory: A new directory is created inside the NEW_RESULTS_ROOT folder. This new folder is named using your search tokens, platform filter, and a timestamp for easy identification (e.g., alzheimer_GPL570_20251007_154500).

    File: Inside this folder, the tool saves a single, comprehensive CSV file named step1_selected_samples.csv. This file contains the detailed metadata for every sample from all the studies you selected.

Finally, the main application window will update to show the selected GSE IDs, confirming that your curated dataset is now ready for the labeling and analysis stages in Step 2.

















# Step 2: Grouping Your Samples (Case/Control Labeling):

After extracting your data, the next critical step is to classify each individual sample (GSM) into meaningful biological groups, such as 'case' vs. 'control'. BioGenVariate offers two flexible modes for this task. You can either use the data you just filtered in Step 1 or load your own dataset using the "Load External File" button.

# Option 1: Manual Labeling ( For Maximum Precision )

This mode gives you complete, sample-by-sample control over the classification process.

1. Click the "Manual Labeling" button.

2. A dialog window will appear for the first sample, displaying all of its associated metadata (title, characteristics, etc.).

3. Based on the information, choose one of the three options:

- Case (1): Assigns the sample to the 'case' group (e.g., diseased, treated).

- Control (0): Assigns the sample to the 'control' group (e.g., healthy, untreated).

- Skip (-1): Excludes the sample from the final labeled dataset.

4. The next sample will appear automatically. Repeat this process until all samples have been reviewed.

# Option 2: Automatic Labeling ( For Speed and Scalability )

This mode uses an advanced, platform-aware AI to classify samples. It's ideal for large datasets where manual labeling is impractical and long. The key feature is its ability to use different keywords for different microarray platforms.

1. Click the "Automatic Labeling" button

2. BioGenVariate will first identify all the unique platforms (e.g., GPL570, GPL1261 ) present in your dataset

3. It will then open a series of interactive dialogs, one for each platform.

4. Each dialog presents a two-pane view where you can browse the samples belonging to that specific platform to understand their context.

5. At the bottom, enter the keywords that define your groups specifically for that platform.

- Case Keywords: Words that identify the 'case' group (e.g., tumor, carcinoma, treated).

- Control Keywords: Words that identify the 'control' group (e.g., normal, healthy, untreated).

6. Click "Set Keywords & Continue to Next Platform". Once you have provided keywords for the final platform, click "Finish & Classify All Samples".

For example, your 'control' keywords for a human study on GPL570 might be healthy donor, while for a mouse study on GPL1261 they might be wild-type. BioGenVariate handles this context seamlessly.

# Saving Your Labeled Data:

Whether you chose manual or automatic labeling, the final step is saving. The tool will add a new column to your data (User_Label or Predicted_Label) containing the 0s and 1s.

The labeled datasets are saved as new CSV files inside the same results folder created during Step 1. To keep your results organized, the tool saves a separate labeled file for each platform, named like:

- step2_GPL570_auto_labeled.csv

- step2_GPL1261_manual_labeled.csv

Your data is now fully extracted, filtered, and classified, ready for downstream statistical analysis and visualization within the BioGenVariate suite.


# 2. The Gene Distribution Explorer

This tool is designed for deep-dive analysis. It allows you to visualize the expression distribution of specific genes from the loaded platforms and then perform a comprehensive, automated analysis of the metadata associated with interesting sample groups (e.g., those with very high or very low expression). This process uncovers the biological context—such as diseases, tissues, or cell types—linked to specific gene expression patterns.

# Step 1: Plot Your Gene's Distribution

First, you need to visualize the expression data for your gene(s) of interest.

1. Open the Explorer: Click the "Show Gene Distribution" button in the "Analysis Tools" section of the main window.

2.Select Platforms: In the new window that appears, check the boxes for one or more loaded platforms you wish to analyze (e.g., GPL570, GPL1261). 

3.Enter Your Gene(s): Type the official symbol(s) for the gene(s) you are interested in into the "Gene Symbol(s)" box. You can enter multiple genes separated by commas (e.g., TP53, EGFR). 

4. Plot Distributions: Click the "Plot Distributions" button. A grid of histograms will appear in the window below. Each plot shows the full expression pattern for a single gene on a single platform.

# Step 2: Interactively Analyze a Region of Interest

Now you can specify in BioGenVariate which group of samples you want to investigate.

# Primary Method: Click and Drag Selection

The most powerful way to analyze a group is to click and drag your mouse directly on any histogram. This will draw a rectangle selector. Once you release the mouse button, BioGenVariate will automatically begin analyzing all the samples that fall within the expression range you selected. This gives you precise control to investigate any part of the distribution—the tails, the center, or any specific peak.





# 3. Compare Distributions




This tool is designed for advanced comparisons. It lets you load your own lists of case and control samples and visually compare their gene expression patterns against each other or against different datasets.

# Step 1: Load Your Sample Groups

First, you need to tell BioGenVariate which samples belong to which group.

    Open the Tool: Click the "Compare Distributions" button in the main window to open the comparison interface.

    Load Your File: In the new window, click "Load File(s)" and select a CSV file that contains your sample list.

    Specify Columns: The tool will ask you for two things:

        The name of the column containing the sample IDs (GSMs).

        The name of the column containing the labels (where 1 means case and 0 means control).

    Link to a Platform: You will be asked to associate your sample list with a GPL platform that you have already loaded in the main window. Your loaded groups will now appear in the listbox.

# Step 2: Set Up the Comparison

Next, define what you want to compare.

    Select Your Groups: Click on the groups in the listbox that you want to visualize.

    Choose a Mode: Select one of the Comparison Modes. For example, you can compare your loaded groups directly against each other ("Compare Selected...Only") or compare them against the expression of a specific gene across an entire platform ("Compare Groups vs. Gene(s)...").

    Enter a Gene: If your mode requires it, type a Gene Symbol (e.g., EGFR) into the entry box.

# Step 3: Run the Analysis and View Plots

Finally, generate and view your comparison.

    Click "Plot & Analyze Distributions": The tool will gather all the necessary data and generate plots on the right side of the window.

    Interpret the Results: You will see two types of plots:

        Overlaid Density Plots: These show the overall shape and spread of the expression distributions for each group.

        Overlaid Histograms: These provide a more detailed, frequency-based view of the data.

Below the plots, a table will also show the results of a statistical test (Wilcoxon rank-sum), telling you if the differences between your groups are statistically significant.






# 4. Generating and Interpreting Results



The true power of BioGenVariate is unlocked when you analyze a region of interest from the gene expression plot. Simply clicking "Analyze LEFT Tail" or "Analyze RIGHT Tail" triggers a powerful backend cascade. The application automatically identifies statistically enriched tokens in the metadata of those specific samples, filters them for biological relevance, and clusters them using disease and tissue ontologies. This process generates a suite of visualizations, including a powerful Canonical Word Cloud where terms are colored by their ontological category (e.g., red for disease, green for tissue), providing an immediate, contextualized overview of the biological themes that define your sample set. All plots and detailed statistical tables are automatically saved to the NEW_RESULTS_ROOT directory for your records.


# Example results - tissue canonical tokens frequency results ( AC1045372.2 from GPL570 - right tail sector of distribution ) 

<img width="1781" height="1183" alt="frequencyplot_clusteredtokens_Tissue Tokens" src="https://github.com/user-attachments/assets/f38914b2-5110-4569-9b55-eded4e13a20b" />



# Example results - disease canonical tokens frequency results ( WDR96 from GPL96 - right tail sector of distribution ) 

<img width="1781" height="1183" alt="frequencyplot_clusteredtokens_Disease Tokens" src="https://github.com/user-attachments/assets/2ac801aa-dcbc-4273-890d-3e8ac95ad266" />



# Workflow Overview: 📊 

  Load GPL Datasets: Use the buttons in the UI to load the desired platform expression datasets into memory.

  Run Step 1 (GSE Extraction): Enter filter tokens (e.g., "cancer,tumor") and platform IDs (e.g., "GPL570") and click "Run GSE Extraction".

  Review & Select GSEs: A dialog will appear summarizing the extracted GEO Series (GSEs). Select the series you wish to analyze further.

  Run Step 2 (Labeling): Choose either "Automatic" or "Manual" labeling to classify the samples from the selected GSEs into case/control groups.

  Analyze in Gene Explorer: Launch the "Gene Distribution Explorer", select genes and platforms, and click "Plot Distributions". From there, you can interactively analyze tails or bins to generate word clouds, cluster graphs, and statistical reports.


  

 # Output Results File Structure: 📁

BioGenVariate/
├── main.py
├── config.py
├── requirements.txt
|
├── gui/
│   ├── __init__.py
│   ├── main_window.py
│   ├── compare_window.py
│   └── distribution_explorer.py
|
├── data_processing/
│   ├── __init__.py
│   ├── extraction.py
│   ├── labeling.py
│   └── loader.py
|
├── utils/
│   ├── __init__.py
│   ├── embedding.py
│   ├── ontology.py
│   ├── statistics.py
│   ├── patches.py
│   └── helpers.py
|
└── NEW_RESULTS_ROOT/
    └── .gitkeep

 

All results are saved to the NEW_RESULTS_ROOT/ directory, which is created automatically.

  Extraction/Labeling Results:

        NEW_RESULTS_ROOT/{query}_{platform}_{timestamp}/

            step1_extracted_gsm.csv: Raw data for all matched samples.

            step2_{...}_labeled.csv: Data with an added label column.

  Gene Distribution Analysis Results:

        NEW_RESULTS_ROOT/GeneDistributionAnalysis/{gene}_{platform}_{analysis_type}/

            Contains all plots (density, word clouds, graphs) as .png or .jpg files.

            Contains all statistics (enriched tokens, bio-specific tokens, all p-values) as .csv.gz files.
            
  

# License: 📜

This project is licensed under the MIT License - see the LICENSE.md file for details.
