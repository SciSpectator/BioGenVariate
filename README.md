# GeoExplorer 

🔬 Core Features

Intelligent Sample Extraction: Utilizes the entire GEOmetadb.sqlite database, loaded into memory for rapid querying, and applies a two-stage semantic and lexical filtering process to find the most relevant samples.

Bio-Aware Token Analysis: Leverages a fine-tuned PubMedBERT model to identify statistically enriched and biologically relevant keywords from sample metadata.

Ontology-Driven Clustering: Maps relevant keywords to the DOID (Human Disease Ontology) and Uberon (multi-species anatomy ontology) to cluster findings into meaningful biological categories.

Flexible Case-Control Labeling: Supports both fully automatic labeling based on semantic similarity and an interactive manual labeling mode for expert-driven classification.

Interactive Gene Distribution Explorer: A powerful interface to plot expression distributions, interactively select sample tails or bins, and trigger deep-dive analysis on those subsets.

Automated Batch Processing: Includes functionality to automate tail analysis for all genes on a platform, complete with a negative control run (shuffled metadata) to validate findings.

Structured & Reproducible Results: Automatically generates a comprehensive set of plots, graphs, and data tables for each analysis, saved in a logically structured, timestamped directory.
