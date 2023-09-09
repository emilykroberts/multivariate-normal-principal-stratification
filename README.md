
### R code for multivariate normal model for surrogate endpoint validation
---

This package can be installed using the following R code:

```r
install.packages("devtools")
```

```r
devtools::install_github("emilykroberts/multivariate-normal-principal-stratification")

library(mvnPS)
```

in development

### Resources

* [Ask a question/ Open an issue: coming soon](https://github.com/emilykroberts) (GitHub issues for bug reports, feature requests)

-----------------------------------------------------------------------------------------------------
Functionality, structure, files, and functions
-----------------------------------------------------------------------------------------------------

### Core functionality

The main functions in this repository are written to simulate data and run analysis (either of simulated or real data) for surrogate endpoint validation of two normally distributed surrogate and true outcome in a randomized clinical trial.

### Purpose of functions

First, data in the illness death format for a randomized trial can be generated using the generatedata.R file. This will provide a datasets that follows the models for each treatment arm described in the manuscripts. 
The following files contain the likelihood components or calculations for different parameters in the proposed models based on the Bayesian framework:

uniform_prior.R for the posterior computation of correlation parameters assuming a Uniform distribution, and beta_prior.R for the posterior computation assuming a Beta prior.

The run_sim.R function runs the simulations/MCMC and performs estimation and surrogacy validation based on our proposed methods. Once this file is run, .txt files are created with the results (one file per dataset over iterations of the MCMC). These results can be read in using code to create a table suitable for LaTeX or for output to a .csv file. The simulations in the manuscript are run in parallel on a high performance computer cluster. It is possible to run many iterations in parallel, though a single simulation can be run by setting the array_id variable to a fixed number and changing the number of MCMC iterations (SIM).

Based on the estimation in run_sim.R, the draws of parameters from the MCMC and other values and saved, which can be shown and written to a file using final_results.R.
      
plot_traceplots.R is a diagnostic tool to examine the traceplots of the parameter draws as desired by the user.

-----------------------------------------------------------------------------------------------------
Other information
-----------------------------------------------------------------------------------------------------

### Contributing 

If you are interested in contributing to the development please open an issue to request or contact at emily-roberts-1@uiowa.edu.

### References and other literature

Roberts, E. K., Elliott, M. R., & Taylor, J. M. (2021). Incorporating baseline covariates to validate surrogate endpoints with a constant biomarker under control arm. Statistics in Medicine.

Roberts, E. K., Elliott, M. R., & Taylor, J. M. (2022). Solutions for surrogacy validation with longitudinal outcomes for a gene therapy. Biometrics.

Conlon, A. S., Taylor, J. M., & Elliott, M. R. (2014). Surrogacy assessment using principal stratification when surrogate and outcome measures are multivariate normal. Biostatistics.
