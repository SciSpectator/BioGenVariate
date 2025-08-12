# GeoExplorer 

GeoExplorer is a sophisticated, GUI-driven Python application designed for in-depth analysis of gene expression multispiecies data from the Gene Expression Omnibus (GEO). It provides a multi-step workflow to extract, filter, label, and analyze biological samples, moving beyond simple keyword searches to incorporate advanced semantic filtering, statistical analysis, and ontology-based natural language processing contextualization from sample descriptions.

The suite is engineered for researchers who need to identify relevant sample cohorts, analyze case/control groups, and explore the descriptive metadata associated with unique gene expression patterns. It excels at analyzing regions of interest in gene expression distributions (such as tails or specific bins) and generating insightful visualizations, including word clouds of biologically enriched tokens that are clustered into canonical terms based on ontology names like DOID and Uberon.


🔬 #Core Features

- Intelligent Sample Extraction: Utilizes the entire GEOmetadb.sqlite database, loaded into memory for rapid querying, and applies a two-stage semantic and lexical filtering process to find the most relevant samples.

- Bio-Aware Token Analysis: Leverages a fine-tuned PubMedBERT model to identify statistically enriched and biologically relevant keywords from sample metadata.

- Ontology-Driven Clustering: Maps relevant keywords to the DOID (Human Disease Ontology) and Uberon (multi-species anatomy ontology) to cluster findings into meaningful biological categories.

- Flexible Case-Control Labeling: Supports both fully automatic labeling based on semantic similarity and an interactive manual labeling mode for expert-driven classification.

- Interactive Gene Distribution Explorer: A powerful interface to plot expression distributions, interactively select sample tails or bins, and trigger deep-dive analysis on those subsets.

- Structured & Reproducible Results: Automatically generates a comprehensive results ( set of plots, graphs, and data tables for each analysis ) , saved in a logically structured, timestamped directory.


⚙️ # Installation & Setup

Prerequisites

    Python 3.8+

    PyTorch (with CUDA support recommended for performance)

    Sufficient RAM to load GEOmetadb and GPL expression datasets into memory (>16 GB recommended).

1. Clone the Repository

git clone https://github.com/your-username/GeoExplorer.git
cd GeoExplorer


2. Install Dependencies

It is highly recommended to use a virtual environment.

python -m venv venv
source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
pip install -r requirements.txt


The application will also download necessary NLTK data models (punkt, stopwords, wordnet) on first run.
3. Required Data Files

GeoExplorer requires several large data files to be placed in specific directories.

    GEOmetadb.sqlite.gz:

        Action: Download from GEOmetadb or another source.

        Location: Place this file in the root directory of the project.

    GPL Expression & Token Data:

        The application expects a directory structure for each platform's data. Create a main data directory (e.g., /home/user/GEOMETADB_TOKENS_PLATFORMS) and update the path in the GeoWorkflowGUI class constructor.

        Inside this directory, create a sub-folder for each platform (e.g., GPL570, GPL96).

        Expression Data: Place the gzipped CSV expression matrix for each platform in its respective folder (e.g., GPL570/gpl570_all_samples_normalized.csv.gz).

        Token Data: Place the pre-computed token file for each platform in its folder (e.g., GPL570/geometadb_tokens_GPL570.csv.gz).

    Fine-tuned Model & Keywords:

        The fine-tuned PubMedBERT model and biomedical keyword corpus are required for the AttentionSeeker functionality.

        Update the self.model_dir and self.keywords_dir paths in the GeoWorkflowGUI class to point to their locations.


🚀 Usage:

Once all dependencies are installed and data files are in place, run the main script from the project's root directory:

python geo_workflow_gui.py


The main application window will launch, and you can begin the analysis workflow.


How to Use GeoExplorer:



📊 Workflow Overview:

  Load GPL Datasets: Use the buttons in the UI to load the desired platform expression datasets into memory.

  Run Step 1 (GSE Extraction): Enter filter tokens (e.g., "cancer,tumor") and platform IDs (e.g., "GPL570") and click "Run GSE Extraction".

  Review & Select GSEs: A dialog will appear summarizing the extracted GEO Series (GSEs). Select the series you wish to analyze further.

  Run Step 2 (Labeling): Choose either "Automatic" or "Manual" labeling to classify the samples from the selected GSEs into case/control groups.

  Analyze in Gene Explorer: Launch the "Gene Distribution Explorer", select genes and platforms, and click "Plot Distributions". From there, you can interactively analyze tails or bins to generate word clouds, cluster graphs, and statistical reports.


  

📁 Output Results File Structure:

All results are saved to the NEW_RESULTS_ROOT/ directory, which is created automatically.

  Extraction/Labeling Results:

        NEW_RESULTS_ROOT/{query}_{platform}_{timestamp}/

            step1_extracted_gsm.csv: Raw data for all matched samples.

            step2_{...}_labeled.csv: Data with an added label column.

  Gene Distribution Analysis Results:

        NEW_RESULTS_ROOT/GeneDistributionAnalysis/{gene}_{platform}_{analysis_type}/

            Contains all plots (density, word clouds, graphs) as .png or .jpg files.

            Contains all statistics (enriched tokens, bio-specific tokens, all p-values) as .csv.gz files.
            

📚 Major Dependencies:

  PyTorch

  Hugging Face Transformers

  Pandas & NumPy

  SciPy & Scikit-learn

  Matplotlib & Seaborn

  BERTopic

  NetworkX
  

📜 License:

This project is licensed under the MIT License - see the LICENSE.md file for details.
