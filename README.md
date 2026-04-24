# [Farming practices exerts selection pressures on the resistome of natural populations of house mice](https://www.biorxiv.org/content/10.64898/2025.12.12.693898v1.full)

Ecological Dynamics of Antibiotic Resistance Genes (ARGs) in Wild Mice Populations.\

## Authors

| Name                       | Affiliation |
|----------------------------|-------------|
| Morgane Gicquel            | IZW         |
| Aimara Planillo            | IZW         |
| Emanuel Heitlinger         | HU          |
| Sofia K Forslund-Startceva | MDC         |
| Stephanie Kramer-Schadt    | IZW         |
| Susana C. M. Ferreira      | UV          |
| Victor Hugo Jarquín-Díaz   | MDC         |

## Raw data

Available online at the European Nucleotide Archive with an accession number `PRJEB102564`\

## Overview

This repository contains data, analysis scripts, and results from a study investigating the presence, abundance, and ecological drivers of antibiotic resistance genes (ARGs) in the gut microbiome of wild mice captured across multiple farms.\

## Objectives

1)  Identify the diversity and distribution of ARGs in wild mice gut microbiomes.\
2)  Quantify the contribution of host characteristics, land cover, farm management, and climate to ARG presence and abundance.\
3)  Assess the importance of random effects such as spatial structure and farm-specific factors.\
4)  Explore trait-level patterns (mobility, localization, drug class, resistance mechanism) influencing ARG occurrence and abundance.
5)  Compare mice resistomes to livestock resistomes.\

## Data

-   Sampling units:\
    -   875 mice from 230 farms across Berlin, Brandenburg, Mecklenburg and Bavaria, Germany 
-   ARGs:\
    -   340 ARGs (presence/absence analysis)\
    -   161 ARGs (abundance analysis after filtering)\
-   Explanatory variables:\
    -   Mouse characteristics (sex, BMI)\
    -   Mouse density (number of mice trapped at the farm)\
    -   Land cover (agriculture, tree cover, imperviousness, small woody features)\
    -   Accessibility (distance to roads/paths)\
    -   Farm animal density (poultry, pig, cattle, at the municipality level)\
    -   Climatic variables (precipitation, temperature)\

## Modelling Approach

Method: Hierarchical Modelling of Species Communities (HMSC)\
Random effects: Spatial structure (latitude/longitude), Farm ID, Year\
Posterior sampling: 4 MCMC chains, 1000/1500 samples each

## Repository Structure

-   `data/raw` \# Raw datafiles
-   `data/processed` \# Preprocessed ARG datasets\
-   `scripts/` \# Analysis scripts (jSDM models, data processing, microbiome and ARG description)
-   `1_ARG_Preparing.Rmd`
    -   Aim: Preparing Variables in ARG data for jSDM
-   `2_EnvData_Extract_Explo.Rmd`
    -   Aim: Environmental data extraction for areas with mice for jSDM
    -   Outputs in supplements manuscript: *Supplementary S5*
-   `3_ARG_Exploring.Rmd`
    -   Aim: Exploring Variables in ARG data (general description of ARG variables)
-   `4.1_RunjSDM_ARG_ARO_pa.Rmd`
    -   Aim: Mice ARG jSDM presence/absence (occurrence)
    -   Outputs in main manuscript: *Figure 3 and 4*
    -   Outputs in supplements manuscript: *Supplementary S12-S13*
-   `4.2_RunjSDM_ARG_ARO_abd.Rmd`
    -   Aim: Mice ARG jSDM abundance
    -   Outputs in supplements manuscript: *Supplementary S15-S17*
-   `4.3_RunjSDM_ARG_ARO_pa_HI.Rmd`
    -   Aim: Mice ARG jSDM presence/absence (occurrence) including host genotype
    -   Outputs in supplements manuscript: *Supplementary S14*
-   `5_Map_figures.Rmd`
    -   Aim: Geographical localisation of mice from analysis
    -   Outputs in main manuscript: *Figure 1*
    -   Outputs in supplements manuscript: *Supplementary S4*
-   `6_Tax_ARG_alpha_beta_div.Rmd`
    -   Aim: Analysis of alpha and beta diversity of microbiome and ARGs
    -   Outputs in main manuscript: *Figure 2*
    -   Outputs in supplements manuscript: *Supplementary S7, S8, S10, S11*
-   `7_Mouse_Livestock_comparison.Rmd`
    -   Aim: Exploration of livestock resistomes and comparison of alpha and beta diversity with mice's resistomes
    -   Outputs in main manuscript: *Figure 5*
-   `8_ARG_genetic_context.Rmd`
    -   Aim: Screening of ARGs localised in contigs annotated as plasmids
    -   Outputs in supplements manuscript: Supplementary S9
-   `9_Compare_livestock_data.Rmd`
    -   Aim: Quantify differences in livestock density along years
    -   Outputs in supplements manuscript: *Supplementary S18*
-   `source_packages.R`
    -   Aim: Source file for libraries and functions used in the project
-   `source_themes_mg.R`
    -   Aim: Set up premade themes for ggplots jSDMs
-   `scripts/Preprocessing` \# Code snippets scripts (Raw sequencing data preprocessing, microbiome and ARG annotation)
-   `output/` \# Model output files for different analyses Model outputs can be read directly in the scripts, rather than having the model run again to save time.
-   `figures/` \# Figures for manuscript\
-   `README.md` \# Project documentation (this file)

## Description

-   Raw sequencing data has been preprocessed on the High Performance Computing Cluster from the Max Delbrück Centrum, Berlin (Max-Cluster) equipped with 5.2k Cores – 46TB RAM.\
-   Snippets of the scripts used for preprocessing are available in `scripts/Preprocessing`.\
-   The statistical analysis is presented in R-markdown documents and R-Scripts. Rmd files should be run in an R session using knitr.\
-   Scripts are numbered sequentially, but they do not need to run sequentially. The input of each is maintained at `data/`.\
