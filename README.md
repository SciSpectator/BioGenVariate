# GeoExplorer 
🔬 Core Features

    Intelligent Sample Extraction: Utilizes the entire GEOmetadb.sqlite database, loaded into memory for rapid querying. Samples are filtered not just by platform but through a powerful two-stage process:

        Lexical Screen: A fast, broad search using user-defined tokens and their WordNet synonyms.

        Semantic Screen: A precise filter using SciBERT embeddings to find samples semantically similar to the user's query context.

    Bio-Aware Token Analysis: Leverages a fine-tuned PubMedBERT model (AttentionSeeker class) and a curated biomedical keyword corpus to:

        Identify statistically enriched keywords in sample metadata using Z-tests.

        Filter enriched keywords for biomedical relevance based on semantic similarity, eliminating noise.

    Ontology-Driven Clustering: Maps bio-relevant keywords to canonical terms from the DOID (Human Disease Ontology) and Uberon (multi-species anatomy ontology). This allows for the clustering and visualization of tokens into meaningful biological categories (e.g., diseases, tissues).

    Flexible Case-Control Labeling:

        Automatic Mode: Classifies samples using BioBERT embeddings to measure similarity against user-defined "case" and "control" concepts.

        Manual Mode: Provides a streamlined UI for expert-driven, sample-by-sample labeling.

    Interactive Gene Distribution Explorer: A powerful pop-up interface to:

        Plot and analyze expression distributions for any gene across multiple platforms.

        Interactively select samples from histogram bins or distribution tails for deep-dive analysis.

        Generate and save publication-quality plots: density plots, word clouds, ontology-based cluster graphs, and frequency bar plots.

    Batch Processing: Automates the analysis of distribution tails for all genes on a platform, including a negative control run with shuffled metadata to validate findings.

    Robust & Self-Contained: Multi-threaded operations for a non-blocking UI, intelligent caching for ontologies and embeddings, and a clear, structured output file system for results.
