# GeoExplorer 

GeoExplorer is a sophisticated, GUI-driven Python application designed for in-depth analysis of gene expression multispiecies data from the Gene Expression Omnibus (GEO). It provides a multi-step workflow to extract, filter, label, and analyze biological samples, moving beyond simple keyword searches to incorporate advanced semantic filtering, statistical analysis, and ontology-based natural language processing contextualization from sample descriptions.

The suite is engineered for researchers who need to identify relevant sample cohorts, analyze case/control groups, and explore the descriptive metadata associated with unique gene expression patterns. It excels at analyzing regions of interest in gene expression distributions (such as tails or specific bins) and generating insightful visualizations, including word clouds of biologically enriched tokens that are clustered into canonical terms based on ontology names like DOID and Uberon.


🔬 Core Features

- Intelligent Sample Extraction: Utilizes the entire GEOmetadb.sqlite database, loaded into memory for rapid querying, and applies a two-stage semantic and lexical filtering process to find the most relevant samples.

- Bio-Aware Token Analysis: Leverages a fine-tuned PubMedBERT model to identify statistically enriched and biologically relevant keywords from sample metadata.

- Ontology-Driven Clustering: Maps relevant keywords to the DOID (Human Disease Ontology) and Uberon (multi-species anatomy ontology) to cluster findings into meaningful biological categories.

- Flexible Case-Control Labeling: Supports both fully automatic labeling based on semantic similarity and an interactive manual labeling mode for expert-driven classification.

- Interactive Gene Distribution Explorer: A powerful interface to plot expression distributions, interactively select sample tails or bins, and trigger deep-dive analysis on those subsets.

- Automated Batch Processing: Includes functionality to automate tail analysis for all genes on a platform, complete with a negative control run (shuffled metadata) to validate findings.

- Structured & Reproducible Results: Automatically generates a comprehensive set of plots, graphs, and data tables for each analysis, saved in a logically structured, timestamped directory.
