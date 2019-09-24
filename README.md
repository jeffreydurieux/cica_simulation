# C-ICA simulations

# Table of Contents
1. [Simulation 1](#Sim1)
2. [Simulation 2](#Sim2) 
    1. [Create FC patterns](#CreateMaps)
3. [Simulation 3](#Sim3)



## Folder Simulation 1 <a name="Sim1"></a>
* Simulation1/Simulate_CICA.R       : R-script to generate data according to C-ICA model
* Simulation1/SimulationScript.R    : R-script to run simulation study 1 

* Simulation1/CreateANOVAData.R     : R-script to that combines all output from simulation 1 and creates a dataframe

* Simulation1/AnalysisOnSimu...R    : Script for doing ANOVA on simulation results


## Folder Simulation 2 <a name="Sim2"></a>

* Simulation2/Simulation2Script.R           : R script for computer cluster to run simulation 2
* Simulation2/AnalysisOnSimulation2Results.R : R script for analysis

### Subfolder: CreateMaps <a name="CreateMaps"></a>

* Simulation2/FCpatterns.R          : R-script used to generated FC patterns for three clusters (4 components each)

* Simulation2/cluster1.Rdata        : R data for cluster 1 (4 components each)
* Simulation2/cluster2.Rdata        : R data for cluster 2 (4 components each)
* Simulation2/cluster3.Rdata        : R data for cluster 3 (4 components each)
* Simulation2/bmaskslice.Rdata      : R data for a mask used in the analysis

* Simulation2/CreateMaps/GenerateCICAdata.R : R script that simulated the data for simulation 2

## Folder Simulation 3 <a name="Sim3"></a>

* Simulation3/Simulation3Script.R     : R-script to do model selection (simulation 3) 
* Simulation3/SelectLowestLoss.R      : R-scriopt to select lowest loss for parallel multiple starts
* Simulation3/SequentialScreeSim3.R   : R-script for sequential scree test of sim 3




