Bayesian modelling of football outcomes: using the Skellam’s distribution for the goal difference 
Simulation & Estimation Workflow
============================================

1. Set Up Working Directory
---------------------------
Open R and set your working directory to the folder named “Team_17”.  
For example:
    setwd("C:/path/to/Team_17")

2. Simulation of Data
---------------------
• Run the simulation script by typing:
    source("simulation.R")
• This will generate a file called
    simulated_data.csv
  in your “Team_17” folder.  
• (You already have simulated_data.csv uploaded—run this only to verify the simulation code.)

3. Parameter Estimation
-----------------------
• Make sure simulated_data.csv is present in the “Team_17” directory.  
• Then run:
    source("estimation.R")
• This script will:
    – Perform the MCMC sampling  
    – Save the sampler outputs and sample draws  
    – Produce all diagnostic plots (ACF plots, trace plots, etc.)  
    – Compute and display Effective Sample Size (ESS) values  
• All the plots and ESS values generated are saved in the “plots and ess values” folder.  
• Note: Estimation may take several minutes, depending on your machine.

4. Required R Packages
----------------------
Before running the scripts, install these packages once:
    install.packages(c(
      "dplyr",
      "Bessel",
      "coda",
      "MASS",
      "gridExtra",
      "ggplot2"
    ))

5. Issues
--------------------
If you run into any issues, please check that:
  • Your working directory is correctly set to “Team_17”  
  • All required packages are installed  
  • simulation.R and estimation.R are in the same folder as simulated_data.csv  


