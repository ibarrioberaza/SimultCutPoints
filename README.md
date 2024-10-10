# SimultCutPoints
Optimal Simultaneous Categorisation of Continuous Variables in Prediction Models

## Description
This repository contains code developed based on the methodology proposed by Barrio et al. (2024) to:
1. Categorise a predictor variable in any regression model where the response variable follows an exponential family distribution.
2. Simultaneously categorise multiple continuous covariates in multivariate settings.
The methodology employs a computationally efficient approach based on a pseudo-BIC criterion to determine the optimal number of categories for each variable.

## References
Barrio, I., Roca-Pardiñas, J., Esteban, C., & Durban, M. (2024). Proposal of a general framework to categorize continuous predictor variables. arXiv preprint arXiv:2403.11983. https://doi.org/10.48550/arXiv.2403.11983

## R-scripts
There are two main files:
1. Help_Functions.R 
2. SimultCutPoints.R
3. Example.R

## Instructions
Step 1: Open R (or RStudio) and set the working directory to the folder where the files (codes and data) are located. Another option is to create a project directly from this repository.
Step 2: Open the "Example.R" file and execute the code.
Step 3: 
