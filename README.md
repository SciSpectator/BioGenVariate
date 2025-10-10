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


<img width="1024" height="158" alt="Step1" src="https://github.com/user-attachments/assets/e6ddeb44-1f98-4e1c-8570-0565d38fb533" />

_Ryc1. Step 1 visible on the GUI of GeneVariate application with Platform Filter and Filtering Tokens option_
















    
Platform Filter (optional): To restrict your search, enter a comma-separated list of platform IDs (e.g., GPL570, GPL96). This is not limited to the default platforms; any platform ID available in the geometadb file can be used. Leave this field blank to search across all platforms.


Filtering Tokens: This is the core of your search. Enter a comma-separated list of keywords that describe your research interest (e.g., alzheimer disease,control,brain). The tool will search for these terms in both the overall study (GSE) descriptions and the detailed individual sample (GSM) metadata.

Once you have configured your search, click "Run GSE Extraction". The process will begin in the background, and you can monitor its progress via the progress bar and the log output window.


# Step 1.5: Interactively Review and Refine Your Results


Upon completion, a new "Step 1.5: Review, Analyze, and Select GSEs" window will automatically open. This powerful interface is designed to help you quickly validate the relevance of each study (GSE) found by the tool.

The window is split into two panes:


<img width="1199" height="832" alt="Step1_5" src="https://github.com/user-attachments/assets/04e0c936-e0b2-42f9-a16f-283b00877d96" />

_Ryc2. Example of the review window in Step 1.5 for extracted Filtering Tokens: "Alzheimer Disease" in platform GPL570_
















Left Pane (GSE List): This is a table of all studies that contained at least one of your keywords. It provides a concise summary, including the platform, the specific "Matching Tokens" found, and a count of samples that were "Direct Matches" versus those that are just associated with the study. This allows you to immediately gauge the relevance of each result. You can select one or more GSEs from this list to inspect and keep.

Right Pane (Details and Keyword Highlighting): When you click on a GSE in the left pane, this area populates with the combined titles, descriptions, and characteristics from its samples. To make validation effortless, all of your original search tokens are highlighted, so you can instantly see them in their proper context.

Review the studies and select all the ones you wish to include in your final dataset by clicking on them in the left pane (use Ctrl-Click or Shift-Click to select multiple). When you are finished, click the "Save Kept GSEs and Continue to Step 2" button.


Saving and Finalizing Your Selection

GenVariate automatically organizes your results for you.


Directory: A new directory is created inside the NEW_RESULTS_ROOT folder. This new folder is named using your search tokens, platform filter, and a timestamp for easy identification (e.g., alzheimer_GPL570_20251007_154500).

File: Inside this folder, the tool saves a single, comprehensive CSV file named step1_selected_samples.csv. This file contains the detailed metadata for every sample from all the studies you selected.

Finally, the main application window will update to show the selected GSE IDs, confirming that your curated dataset is now ready for the labeling and analysis stages in Step 2.


# Step 2: Classifying and Labeling Your Samples

After extracting your data, the next critical step is to classify each individual sample (GSM) into meaningful biological groups. GenVariate offers two flexible and powerful modes for this task: a high-speed, AI-driven approach and a precise, expert-driven manual workflow.

# Preparing Your Data Source

Before you begin labeling, you must choose your input data. You have two options:


- Use Data from Step 1 (Default Workflow): If you have just completed Step 1.5, your curated list of samples is already loaded and ready to go.

- Load an External File: You can bypass Step 1 entirely and use your own pre-filtered dataset. This is ideal if you already have a list of samples from a previous analysis, a publication, or a private dataset.

  To do this, click the "Load External File" button in the "Step 2" panel of the main window.

- Select your .csv file. The file must contain the necessary sample metadata columns for the labeling process to be effective.
  

Once your data source is ready, choose one of the following labeling methods.

# Option 1: AI-Powered Condition Labeling 🤖

This mode uses a local AI agent (gpt-oss:20b) to analyze the rich text metadata of each sample and extract structured information automatically. It's ideal for large datasets where manual labeling is impractical.

- Click the "AI-Powered Condition Labeling" button.

  The process begins immediately in the background. The AI agent reads the full metadata for each sample (including title, source, characteristics, and protocols).
For every sample, the agent identifies and extracts key biological attributes into distinct columns:


<img width="1010" height="625" alt="Step2" src="https://github.com/user-attachments/assets/88cfe815-b20d-4bb8-9f23-beb76bd1351d" />

_Ryc3. Example of GUI preview during "AI-Powered Condition Labeling" for Filtering Tokens: "alzheimer diosease, alzheimer" in human and mouse default platforms_
































Classification is based on these criteria such as:


- Classified_Condition: The primary biological state (e.g., Alzheimer Disease, Healthy, Pancreatic Cancer).

- Classified_Tissue: The tissue of origin (e.g., Whole Blood, Brain, PBMC).

- Classified_Age: The sample's age (e.g., 35 years, Old Age, infant).

- Classified_Treatment: Any substance or procedure applied (e.g., LPS stimulation, Vehicle, None).

- Classified_Treatment_Time: The duration or time point of treatment (e.g., 24h, 6 hours, one month).


You can monitor the progress in the log window, which provides real-time updates [AI Progress] how many samples were analyzed/how many left to end the analysis, processing speed [Speed] sample/sec, and an estimated time of completion [ETA].




# Option 2: Manual Labeling ✍️

This mode gives you complete, sample-by-sample control over the classification process, allowing you to define the exact labels you need for your analysis.

- Click the "Manual Labeling" button.

A setup dialog will appear asking you to define the labels you want to assign. You can select from a list of defaults (e.g., Condition, Age, Treatment) or add your own custom labels (e.g., Genotype, Batch_Number).

Once your labels are set, a detailed dialog window will appear for the first sample. It displays all available metadata in a clear, readable format.

Based on the information, fill in the text fields for each of the labels you defined.


- Choose one of the three actions for the sample:

        OK: Saves your labels and proceeds to the next sample.

        Skip Sample: Excludes the current sample from the final labeled dataset and moves to the next one.

        Cancel All: Stops the entire manual labeling process.
  

- Saving Your Labeled Data

Whether you chose the AI-powered or manual method, GenVariate automatically saves the results.


- Directory: The labeled data is saved inside the same timestamped results folder that was created during Step 1. If you loaded an external file, a new results folder will be created.

- File Naming: The tool adds new columns to your data (e.g., Classified_Condition for AI or User_Condition for manual). To keep your results organized, it saves a separate labeled file for each platform present in your dataset, named like:
  

        step2_GPL570_ai-conditions.csv

        step2_GPL1261_manual.csv
  

Your data is now fully extracted, filtered, and classified, ready for deep-dive statistical analysis and visualization within the GenVariate suite.










# 2. Compare Distributions

This tool is designed for advanced comparisons of user-defined sample groups. It allows you to load your own lists of samples from an external file and visually and statistically compare their gene expression patterns against each other or against a larger background dataset.


# Step 1: Load and Define Your Sample Groups




First, you need a .csv file that defines your sample groups.


Open the Tool: In the main window, click the "Compare Distributions" button.

Load Your File: In the new window, click "Load File(s)" and select your .csv or .csv.gz file.

Specify Columns: A dialog will appear asking for two key pieces of information:

Sample ID Column (GSM): Select the column in your file that contains the sample IDs (e.g., GSM, Sample_ID).

Grouping/Label Column(s): Select one or more columns that define your experimental groups (e.g., classifications conducted by AI agent such as Condition , Age, Tissue, Treatment or Treatment Time).


Link to a Platform: Finally, you'll be asked to associate your file with one of the GPL platforms you've already loaded. This is essential for the tool to retrieve the correct expression data for your samples.

After completing these steps, your file will be processed, and entries for each unique group will appear in the listbox.



# Step 2: Configure the Comparison




Next, use the panels on the left to define exactly what you want to compare.


Select Your Groups: In the "Select Groups to Compare" listbox, click on the groups you want to visualize. You can select one or multiple groups.

Choose a Comparison Mode: This is the most important setting. Select one of the three modes:


_Compare Selected Groups Only: Directly compares the distributions of the groups you selected against each other._


_Compare Groups vs. Gene(s) on Platform(s): Plots your selected groups alongside the distribution of a specific gene across an entire background platform. This is useful for seeing if your group's expression is unusual compared to all other samples._

_Compare Groups vs. Entire Platform(s) (All Genes): Compares your groups against the combined distribution of all genes from a selected background platform, providing a very broad context._

    
Specify Details:


_Gene Symbol(s): If your analysis is for specific genes (required for Mode 2), enter their symbols here, separated by commas (e.g., EGFR, TNF)._


_Comparison Platform(s): If you chose Mode 2 or 3, check the box for the platform(s) you want to use as a background._


# Step 3: Analyze and Interact with the Results



Click "Plot & Analyze Distributions" to generate the results.


_Overlaid Density & Histogram Plots: These plots visualize the expression distributions. The density plot shows a smoothed shape, while the histogram shows raw sample counts._

_Statistical Test Results: A table below the plots displays the results of a pairwise Wilcoxon rank-sum test, helping you determine if the observed differences are statistically significant._


    

# Interactive Feature: Change Plot Colors

You can customize the colors of the plots for better visualization





# 3. Show Gene Distribution





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
