#-------------------------
# Description of Data
#-------------------------

[1] 'CROS_COV.csv' --> This csv file has eight columns. Total number of rows are 
205 which is also the number of sites. Column 1 has the names for different sites
(given by source). Column 2 has a serial number corresponding to every site for 
indexing purpose. Columns 3-6 have: proportion of tiger habitat available (PH)
for each site,  site level encounter rates for all prey (AP) (i.e., Chital, Gaur, 
Sambar, Muntjac, Wild Pig), large prey (LP) (i.e., Chital, Gaur, Sambar), human 
intervention (HI) (via encounter rates of livestocks), respectively. Column 7 has
values for lambda-hat(ct) for 42 camera trap sites. Column 8 has values for psi-hat(os). 

[2] 'Vcr.csv' --> This file has the neighbourhood matrix 'V' only for the 42 camera
trap sites.

[3] 'Vcros.csv' --> This file has the neighbourhood matrix 'V' for all 205 sites.

#-------------------------
# R Codes
#-------------------------

[1] 'data.R' --> This R code loads the data from 'CROS_COV.csv' into R session along
with the neighbourhood matrices from 'Vcr.csv' and 'Vcros.csv'. It also loads different
R packages needed to run the other codes. The list of required R packages are - 
'matrixcalc', 'Matrix', 'mvtnorm', 'MCMCpack', 'R2cuba', 'abind', 'ggplot2'.

[2] 'utility_functions.R' --> It has user defined functions such as posterior 
distributions of (tau_e, tau_s) for constrained and unconstraned models and integrand 
for computing the marginal density. 

#-----------------------------------------
## Each of the other codes mentioned below needs 'data.R' and 'utility_functions.R' to run.
#-----------------------------------------

[3] 'marginal.R' --> This R code computes the marginal densities for different models.

[4] 'unconstrained.R' --> This R code runs MCMC to fit the unconstrained model and to 
get posterior estimates of different parameters for the unconstrained model (where spatial
confounding is not corrected).

[5] 'constrained.R' --> This R code runs MCMC to fit the constrained model and to get 
posterior estimates of different parameters for the constrained model where spatial 
confounding is corrected.

[6] 'constrained_predictive_hb.R' --> This R code runs MCMC to fit the predictive model
under constrained setup, i.e., to get posterior estimates of latent abundunces for each
of the 205 site including the sites where data lambda-hat(ct) are not obeserved. It also 
estimates the other parameters such as total abundance A(lambda), coefficient parameters
beta, precision parameters tau_e and tau_s. 

