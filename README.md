# GNMF-Var

GNMF-Var is a graph-regularized non-negative matrix factorization framework for deconvolving SARS-CoV-2 variant-associated mutation signatures from wastewater sequencing data. The method decomposes quarter-specific sample-by-mutation frequency matrices into latent mutation components and sample-level mixture weights, while incorporating a spatiotemporal graph to encourage smooth and epidemiologically plausible variant dynamics across wastewater sampling sites and collection dates.

This repository contains scripts used for GNMF-Var analysis of SARS-CoV-2-positive wastewater samples from Southern Nevada. The analysis compares GNMF-Var with baseline NMF and Freyja-based lineage deconvolution to evaluate how variant-specific mutation architecture influences the timing of first detection in wastewater.

## Repository contents

- `batch1.R`: GNMF-Var analysis script for batch 1 wastewater samples.
- `batch2.R`: GNMF-Var analysis script for batch 2 wastewater samples.

## Method overview

GNMF-Var takes as input a sample-by-mutation frequency matrix and a spatiotemporal adjacency matrix connecting samples collected from the same site within one week. The model factorizes the mutation matrix into:

- sample-level mixture weights representing latent variant-associated components;
- component-specific mutation profiles;
- graph-smoothed temporal and spatial trajectories.

Each component can be annotated to SARS-CoV-2 lineages by comparing its mutation profile with curated lineage barcodes.

## Requirements

The scripts were developed in R and require the following packages:

```r
R.matlab
openxlsx
dplyr
tidyr
FAST
