# SimultCutPoints
Optimal Simultaneous Categorisation of Continuous Variables in Prediction Models

## Description
This repository contains code developed based on the methodology proposed by Barrio et al. (2026) to:
1. Categorise a predictor variable in any regression model where the response variable follows an exponential family distribution.
2. Simultaneously categorise multiple continuous covariates in multivariate settings.

The methodology employs a computationally efficient approach based on a pseudo-BIC criterion to determine the optimal number of categories for each variable.

## References
Barrio, I., Roca-Pardiñas, J., Esteban, C., & Durban, M. (2026). Proposal of a general framework to categorize continuous predictor variables. 

## R-scripts
There are three R-files:
1. Help_Functions.R: help functions required to find and optimise the values and number of cut-off points for each continuous variable.
2. SimultCutPoints.R: main function for optimal simultaneous categorisation of continuous variables in prediction models.
3. Example.R: Usage examples of the SimultCutPoints R-function.

## Instructions
Step 1: Open R (or RStudio) and set the working directory to the folder where the R-files are located. Another option is to create a project directly from this repository.

Step 2: Open the "Example.R" file and execute the code.

