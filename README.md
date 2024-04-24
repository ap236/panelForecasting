# Forecasting with Panel Data
Code for Pesaran, Pick and Timmermann "Forecasting with panel data: Estimation uncertainty versus parameter heterogeneity"

CPI subindices
--------------
The main file to run is CPI_main.m. You need to specify a path for Matlab to
save some files. 

Choices in this file are: 
1. saveMSFEi: set this to 1 if you want to save files that will be used later to
  produce the density plots
2. printMSFEtables: set this to 1 if you want to print the output tables
3. printQuantiles: set this to 1 if you want to print the quantile tables 

To make the density plots use CPI_desity_plots.m

There are two choices: 
1. savePlot: set this to 1 if you want to save .pdf and .fig files of the plots
  (otherwise they are produced by Matlab and you can save them manually in any
  format you like)
2. do_cdf: set this to 0 to get a density and 1 to get a cummulative density plot

House price inflation across MSAs
---------------------------------
The main file to run is HousePrice.m. The remainder is as in the CPI application.

Note that the hierarchical Bayesian model is based on random draws of the 
Gibbs sampler and results may therefore be subject to mild variations.

Andreas Pick, pick@ese.eur.nl
