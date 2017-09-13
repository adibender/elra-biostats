## Overview files and folders for simulations regarding evaluation of PAM modeling approach

- **`modelEvaluation`**:

    - Main folder for simulations regarding the ELRA model evaluation
    (Simulation Part A)

    - Contains scripts to setup registry and experiments
    (`setupBatchExperimentsHypocalorics.R`) and submit Jobs for
    individual simulations/simulation settings (`submitJobs.R`).
    *Problems* (functions that define how data for each simulation run is generated,
    `problems.R`) and *Algorithms* (functions that define how simulated data sets
    are evaluated, `algorithms.R`) are also defined here.

    - **`input`**: Contains scripts and data used as input in simulation studies
    (`completeData.Rds`) for example

- **`comparison`**:

    - Main folder for Simulation Part B
    - Part B, Scenario (1) is build upon http://www.ag-myresearch.com/2017_gasparrini_biomet.html and
    adapted to the setting of survival analysis
    - Part B, Scenario (2) is an extension of Scenario (1) and demonstrates
    the ability of our approach to fit a "time-varying" DLNM