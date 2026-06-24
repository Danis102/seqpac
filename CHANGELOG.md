# Changelog

All notable changes to the `seqpac` package will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Shiny Application**: Integrated an interactive, modular Shiny application under `inst/shiny/` to run the small RNA sequence counts analysis workflow described in the vignette:
  - **Launch Page** (`module_launch.R`): A landing page introducing the workflow steps.
  - **Load / Create PAC** (`module_load.R`): Supports loading the Drosophila example PAC, uploading custom `.RData`/`.rds` files, or generating a PAC from raw FASTQ files and a pheno CSV.
  - **Filter & Normalize** (`module_filter.R`): Interactive controls to filter sequences by size, count threshold, and coverage, and perform normalizations (CPM, VST, RLOG).
  - **Annotation Explorer** (`module_annotate.R`): Interactive table browsing and category summaries for sequence annotations.
  - **Post-Filtering Analysis** (`module_analyze.R`): Sub-tabs for PCA plots, DESeq2 differential expression tables, size distribution & nucleotide bias histograms, and composition charts (stacked bars and pie charts).
  - **Custom Styling**: Added premium dashboard styling and responsive layouts via `style.css` (Inter and Outfit google fonts).
- **Package Launcher**: Added `run_seqpac_app()` in `R/run_app.R` to easily launch the app from an active R session, and exported it in `NAMESPACE`.
- **Docker Support**: Added a `Dockerfile` configuring a version-pinned R environment (`rocker/shiny-verse:latest` under emulated `linux/amd64` architecture) with all Bioconductor and CRAN dependencies preconfigured to build and run the Shiny dashboard on port `3838`.
- **Instructions**: Added `run-shiny.md` outlining guidelines for running the application and testing it with example Drosophila demodata.

### Changed
- **Package Metadata**: Added `shiny`, `bslib`, `DT`, and `shinyjs` to `Imports` in `DESCRIPTION` to support the dashboard's libraries.
