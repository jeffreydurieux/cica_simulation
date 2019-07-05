# C-ICA simulations
Code repository for C-ICA simulations
Contains:


## Simulation 1
Simulation1/Simulate_CICA.R       : R-script to generate data according to C-ICA model
Simulation1/SimulationScript.R    : R-script to run simulation study 1 

Simulation1/CreateANOVAData.R     : R-script to that combines all output from simulation 1 and creates a dataframe

Simulation1/AnalysisOnSimu...R    : Script for doing ANOVA on simulation results

Simulation1/Simulation1ANOVAdata.Rdata  : R dataframe for ANOVA


## Simulation 2

### Create FC patterns for C-ICA model

Simulation2/FCpatterns.R          : R-script used to generated FC patterns for three clusters (4 components each)

Simulation2/cluster1.Rdata        : R data for cluster 1 (4 components each)
Simulation2/cluster2.Rdata        : R data for cluster 2 (4 components each)
Simulation2/cluster3.Rdata        : R data for cluster 3 (4 components each)
Simulation2/bmaskslice.Rdata      : R data for a mask used in the analysis

### Create C-ICA data 
In this part the generated FC patterns are mixed with subject specific time courses and specific noise structures are added

Simulation2/CreateMaps/GenerateCICAdata.R : R script that simulated the data for simulation 2
