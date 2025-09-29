# ARGs in Mice Gut Microbiome

Ecological Dynamics of Antibiotic Resistance Genes (ARGs) in Wild Mice Populations  

## Overview

This repository contains data, analysis scripts, and results from a study investigating the presence, abundance, and ecological drivers of antibiotic resistance genes (ARGs) in the gut microbiome of wild mice captured across multiple farms.   
The project explores how environmental, host, and farm-related factors shape the distribution of ARGs.  


## Objectives

Identify the diversity and distribution of ARGs in wild mice gut microbiomes.  
Quantify the contribution of host characteristics, land cover, farm management, and climate to ARG presence and abundance.  
Assess the importance of random effects such as spatial structure and farm-specific factors.  
Explore trait-level patterns (mobility, localization, drug class, resistance mechanism) influencing ARG occurrence and abundance.  

## Data

Sampling units: 846 mice across multiple farms  
- ARGs:  
340 ARGs (presence/absence analysis)  
161 ARGs (abundance analysis after filtering)  
- Explanatory variables:  
Mouse characteristics (sex, BMI)  
Mouse density (number of mice trapped at the farm)  
Land cover (agriculture, tree cover, imperviousness, small woody features)  
Accessibility (distance to roads/paths)  
Farm animal density (poultry, pig, cattle, at the municipality level)  
Climatic variables (precipitation, temperature)  

## Modelling Approach

Method: Hierarchical Modelling of Species Communities (HMSC)  
Random effects: Spatial structure (latitude/longitude), Farm ID, Year  
Posterior sampling: 4 MCMC chains, 1000/1500 samples each  

## Repository Structure

- data/processed       # Preprocessed ARG datasets  
- scripts/             # Analysis scripts (HMSC models, data processing)  
- output/              # Model outputs, variance partitioning plots  
- figures/             # Visualization of figures for manuscript  
- README.md            # Project documentation (this file)  



