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

After an analysis is complete, the "Save Plots" button becomes active. Clicking it opens a dialog where you can choose which of the generated plots you want to save, select a format (png, svg, pdf), and choose a save location.









_Ryc4. Example of GUI preview of results in "Compare Distributions"_






# 3. Show Gene Distribution 📊


This  tool allows you to visualize and analyze the expression distribution of specific genes across loaded microarray platforms. Its primary purpose is to help you identify interesting subsets of samples based on their gene expression levels (e.g., high or low expressors) and then automatically characterize those samples using an integrated AI agent.


# Step 1: Plot Initial Gene Distributions


First, you'll generate plots to see how your gene(s) of interest are expressed across one or more platforms.



Open the Tool: In the main application window, click the "Show Gene Distribution" button. A new "Gene Distribution Explorer" window will open.



Select Platforms: Check the boxes for the GPL platform(s) you want to investigate. You must have already loaded these platforms in the main window.



Enter Genes: In the "Gene Symbol(s)" text box, enter the gene symbols you are interested in. You can enter multiple genes separated by commas (e.g., APP, BACE1, PSEN1).



Plot: Click the "Plot Distributions" button. The tool will generate a grid of histograms, with each plot showing the expression distribution for one gene on one platform. The title of each plot provides useful information, including the best-fit statistical distribution (e.g., Normal, Bimodal) and the total number of samples.
    


# Step 2: Interactively Select a Sample Range


Once the plots are displayed, you can isolate a group of samples within a specific expression range.


Click and Drag: Position your cursor on any histogram. Click and drag horizontally to draw a selection box over the expression range you want to analyze. The bars within your selection will change color to highlight the range.


Activate Analysis: After you make a selection, the "Analyze Selected Range" button will become active.







<img width="1269" height="875" alt="Screenshot from 2025-10-11 00-48-59" src="https://github.com/user-attachments/assets/2f84a2ee-03df-4535-ae9b-929f9b46bac5" />

_Ryc5. Example of GUI preview of gene APP distributions from selected GPL96, GPL6947, GPL7202 and GPL6885 with selected subset of distribution in GPL96 which can be analyzed by clicking "Analyze Selected Range"_
    


# Step 3: Characterize the Selected Samples with AI

This is where the tool's main power lies. By clicking the "Analyze Selected Range" button, you trigger a multi-step automated analysis of the samples you just selected.


GSE Overlay Plot: A new interactive plot window will immediately appear.


What it shows: This plot displays the overall distribution of the gene across the entire platform as a gray histogram. Overlaid on this are "rug" marks representing the individual samples from your selection, color-coded by their experiment ID (GSE).


<img width="850" height="571" alt="Screenshot from 2025-10-13 11-58-19" src="https://github.com/user-attachments/assets/f6061d09-2dc5-4414-9737-1e048e3acf08" />

_Ryc6. Example of GSE Overlay  of gene APP distributions from selected subset of distribution in GPL96 after choosing option "Analyze Selected Range" - each colored "rug" marks tagging positions across entire distribution from the platform which subset is analyzed_




Interactivity: You can click on a GSE ID in the legend to open its official page on the NCBI GEO website. You can also click the colored patch next to it to change its color on the plot. This helps you quickly identify which experiments contributed the most samples to your selected expression range.


AI Classification: In the background, the tool fetches the full metadata for every sample in your selection and sends it to the local AI agent. The AI analyzes the text descriptions (title, characteristics, protocols) to classify each sample based on categories like Condition, Tissue, Treatment, Age, and Treatment Time.


AI Analysis Plots & Table: Once the analysis is finished, several results new windows will appear:


Classified Analysis Plots: For each category the AI classified (e.g., Condition), a new plot window is generated. This plot shows the expression density of your selected samples, broken down by the AI-assigned labels (e.g., 'Alzheimer Disease' vs. 'Control'). This allows you to see if the expression level you selected corresponds to a specific biological or experimental group.






<img width="1298" height="695" alt="Screenshot from 2025-10-13 11-53-10" src="https://github.com/user-attachments/assets/d227ca46-5848-41ff-933e-579b5c9609d3" /><img width="1302" height="693" alt="Screenshot from 2025-10-11 14-16-18" src="https://github.com/user-attachments/assets/34b07f11-900a-4412-a64d-fda946fe29b5" />


_Ryc6. Example Classified Analysis Plots from selected subset of distribution in GPL96 for gene APP after choosing option "Analyze Selected Range" for classified Condition and Tissue colored differently_


Results Table: A final, detailed table lists every sample from your selection. It includes their GSM and GSE IDs, exact expression value, and all the new classifications determined by the AI. You can double-click any row to open that sample's page on the GEO website.





<img width="1799" height="851" alt="Screenshot from 2025-10-13 11-51-46" src="https://github.com/user-attachments/assets/f6336df8-f9bf-4244-bd5a-6cd5279bcc62" />

_Ryc7. Example final table from selected subset of distribution in GPL96 for gene APP after choosing option "Analyze Selected Range" with informations such as:GSM, GSE, Platform, Expression, Classified Age, Condition, Tissue, Treatment and Treatment time_


    

# Step 4: Save Your Analysis Plots

After the analysis is complete, you can save the newly generated figures for your records. Activate Save Button: The "Save Plots" button in the Gene Distribution Explorer window will become active.

Choose Plots and Format: Clicking it opens a dialog where you can: Select which of the generated plots you want to save (the GSE overlay and all the AI analysis plots).

Choose an image format (png, svg, pdf) and select a directory to save the files to.







            
  

# License: 📜

This project is licensed under the MIT License - see the LICENSE.md file for details.
