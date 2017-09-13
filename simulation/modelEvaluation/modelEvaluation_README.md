## Overview simulation scripts

- `createCompleteData.R`: function that takes real data set and imputes days
for patients that did not survive until last follow up day (used later to
generate new simulated data sets), the complete data is then stored in folder
**`input`**
- `createSettingParameters.R`: Uses models created in
**`runModelBatchJobs`** to create list of lambdas (on complete data in
**`input`**). These lambdas are then used to simulate new data sets. Also
creates a list of formulas, that are later used as input (`static`) to fit
`bam` or `gam` to simulated data sets (`dynamic`).
Lambdas and formulas are stored as list objects in the **`input`** folder.
- `setupBatchExperimentHypocalorics.R`: create and set up registry for
BatchExperiments. Load complete data, lambdas and formulas from
**`input`** folder, sets up registry + adds Problems, Algorithms + Experiments.
- **`modelEvaluationELRA-files`**: created in
`setupBatchExperimentHypocalorics.R`, contains registry + all outputs
- `runSimEval.R`: extracts data from registry and saves results to
`../paper/lda/simulation/` for further processing in the paper.