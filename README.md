# Stochastic GW Background: Recoverying the SMBHB Population

This script uses the EMCEE package to recover the underlying distribution of supermassive black hole binaries (SMBHB) from the observed stochastic gravitational-wave background. This script is optimized to use the SMBHB population realizations in Klein et al. 2019, which are based on the semi-analytic galaxy formation model introduced by Barausse et al. 2012. The SMBHB populations predicted by these models are available at: 

https://people.sissa.it/~barausse/catalogs.html


This script follows two steps: 
1) Calculate the stochastic GW background for each given SMBHB population, as an "observation" of the stochastic background. 
2) Uses the stochastic GW background to recover the distribution of SMBHBs over the total-mass space. 

What it returns in paramter.pdf files at the output folder, are constraints on the total mass of the binary at each mass bin, compared against the true value from the simulations.
