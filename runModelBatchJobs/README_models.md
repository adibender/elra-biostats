### Description of steps to run the models on the data sets created in `../dataGenerationScripts/`

1.) create combinations of data and model specifications for which models should be run (see `combinations.R`).

2.) create R scripts for different specifications on different data-sets
(main analysis vs. sensitivity analysis), etc. (see `createRunModelFiles.R` and `runModelsTemplate.R`).

3.) The created scripts can be run on the cluster (all settings are already specified within the files and include setting of the cluster(s) via the
`BatchJobs` R-package)

  - server setup is specified in `BatchJobs.R`/`BatchJobsLocal.R` (both need
  to be updated to reflect your setup)
  - all scripts created in 2.) start with `submitModels`, possibly followed
  by further specifications (e.g. used data sets)

4.) When scripts described above are run a folder with a name starting with
`gamBatch` is created by the BatchJobs package and contains

  - `registry.RData`: stores information on Jobs run etc. (there are couple of
  functions to extract this information)
  - `jobs`: folder containing results of each job, if something is returned,
  the `*.out` file with echoed R-Session of the respective job
  - `results`: folder containing model objects in `.Rds` format

NOTE: These models can take up a lot of memory (upwards of 20GB) and take some time to finish (up to 6 hours, depending on the hardware)
