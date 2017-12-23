# Evaluation of the association between nutritional adequacy and survival

## How to rerun analysis:

Rerunning the analysis involves 4 steps:

1. Preprocess raw data and create data in piece-wise exponential format
2. Estimate models (main and sensitivity) + run alternative models
with cumulative measures of nutrition (added at review stage)
3. Rerun simulation analysis
    - Simulation Part B is independent of the application and could
    be run "in a vacuum"
    - For Simulation Part A the previous two steps are necessary
4. Recreate Graphs and Tables used in the publication (assumes that
all previous steps ran without errors)

To perform these steps in one, run the code below (your working directory should be set to the directory of the `rerun-analyses.R` file):
```r
source("rerun-analyses.R")
```

This will perform steps 1-4 described above.

**Remark on runtime/memory**: The complexity of the model is very high (many parameters + penalization) and the data sets are also very large (~10k subject + data splitting). Therefore, to run the code (especially simulation studies), we recommend running the code on a server or a very powerful desktop. On our servers, we were able to rerun the entire analysis within 2 days.

### Prerequisites

- For parallel computations we use `mclapply` from the **`parallel`** package
(which doesn't work on windows machines). When you execute the code on a
windows machine, `mclapply` will probably fall back to the default
`mc.cores=1` and thus code will still run, but computation time will
be increased greatly.


- For parallel processing of model fits and simulation runs of Part B we use
the **`BatchJobs`** and **`BatchExperiments`** packages (Bischl et al.
https://www.jstatsoft.org/article/view/v064i11).
For Simulation Part A we use the successor package **`batchtools`**
(https://github.com/mllg/batchtools).

- To use them it is necessary to setup your parallel execution environment (see
files `BatchJobs.R` (server) and `BatchJobsLocal.R` (local) for examples).
Setting `max.jobs=1` in `BatchJobsLocal.R` will run code sequentially, which
might take a while, especially for a full simulation rerun.
Under Linux, make sure that you have execution privileges for the scripts in
`<your R library>/BatchJobs/bin/linux-helper`.
Note: If you only want to check, whether all of the above runs as expected,
but don't want to fully replicate all simulations, reduce `n_simA` and
`n_simB` in `rerun-analyses.R`.


### Additional Notes

- Simulation Study Part A (`simulation/comparison/`) is much more general and
could be of interest for researchers interested in replicating/reusing the data structure and simulation (for example to test extensions of the method, etc.)

- Simulation Study Part B was designed to closely resemble the application
example, thus most code is hard coded (including functions in `elrapack`) and
will not be of much use for general settings.

- We currently develop an R package that facilitates working with PAMMs, including data preparation, visualization, etc.. There are also a lot of vignettes with application examples: [pamm](!http://github.com/adibender/pamm)

### Folder structure
- **`data`**: Raw data for the application example (after initial import
from SAS and minor preprocessing)
- **`dataGenerationScripts`**: Contains scripts for (further) data preprocessing. Creates folder **`dataCurrent`** and **`dataCurrentHosp`** (storing data for main and sensitivity analysis, respectively).
Run `dataImportFromSAStoCleaned.R` to process all data processing steps
- **`elrapack`**: A minimal R-package containing helper functions for data
preparation/evaluation and simulation. This package is not meant to be broadly used, but rather a convenience package for storing helper functions (will be installed locally at the beginning of the `rerun-analyses.R` script).
- **`paper`**: Contains Scripts that produce tables and figures used in the
publication.
- **`runModelBatchJobs`**: Contains scripts to rerun main and sensitivity
analyses of the application example
- **`simulation`**: Scripts to rerun simulation studies
    - **`modelEvaluation`**: Scripts to rerun *Simulation Part A*
    - **`comparison`**: Scripts to rerun *Simulation Part B*


### Session Info

Below you can find the session information of our R session:

```r
sessionInfo()

## R Under development (unstable) (2017-09-06 r73210)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Debian GNU/Linux 8 (jessie)
##
## Matrix products: default
## BLAS: /usr/lib/libblas/libblas.so.3.0
## LAPACK: /usr/lib/lapack/liblapack.so.3.0
##
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C
##
## attached base packages:
## [1] grid      parallel  stats     graphics  grDevices utils     datasets
## [8] methods   base
##
## other attached packages:
##  [1] gganimate_0.1.0.9000   reshape2_1.4.2         tables_0.8
##  [4] Hmisc_4.0-3            Formula_1.2-2          survival_2.41-3
##  [7] lattice_0.20-35        tidyr_0.7.1            gridExtra_2.3
## [10] pam_0.0.776            purrr_0.2.3            magrittr_1.5
## [13] batchtools_0.9.6       data.table_1.10.4      ggplot2_2.2.1
## [16] tsModel_0.6            dlnm_2.3.2             BatchExperiments_1.4.1
## [19] BatchJobs_1.6          BBmisc_1.11            mgcv_1.8-19
## [22] nlme_3.1-131           bindrcpp_0.2           dplyr_0.7.3
## [25] checkmate_1.8.3        lubridate_1.6.0        elrapack_0.0.2
##
## loaded via a namespace (and not attached):
##  [1] httr_1.3.1          bit64_0.9-7         splines_3.5.0
##  [4] assertthat_0.2.0    expm_0.999-2        animation_2.5
##  [7] latticeExtra_0.6-28 base64url_1.2       blob_1.1.0
## [10] progress_1.1.2      RSQLite_2.0         backports_1.1.0
## [13] glue_1.1.1          digest_0.6.12       RColorBrewer_1.1-2
## [16] colorspace_1.3-2    cowplot_0.8.0       htmltools_0.3.6
## [19] Matrix_1.2-11       plyr_1.8.4          pkgconfig_2.0.1
## [22] devtools_1.13.3     mvtnorm_1.0-6       scales_0.5.0
## [25] brew_1.0-6          git2r_0.19.0        tibble_1.3.4
## [28] htmlTable_1.9       withr_2.0.0         nnet_7.3-12
## [31] lazyeval_0.2.0      memoise_1.1.0       msm_1.6.4
## [34] fail_1.3            foreign_0.8-69      tools_3.5.0
## [37] prettyunits_1.0.2   stringr_1.2.0       sendmailR_1.2-1
## [40] munsell_0.4.3       cluster_2.0.6       compiler_3.5.0
## [43] rlang_0.1.2         rstudioapi_0.6      rappdirs_0.3.1
## [46] htmlwidgets_0.9     base64enc_0.1-3     labeling_0.3
## [49] gtable_0.2.0        DBI_0.7             curl_2.8.1
## [52] R6_2.2.2            zoo_1.8-0           knitr_1.17
## [55] bit_1.1-12          bindr_0.1           stringi_1.1.5
## [58] Rcpp_0.12.12        rpart_4.1-11        acepack_1.4.1
## [61] tidyselect_0.2.0

devtools::session_info()

## setting  value
## version  R Under development (unstable) (2017-09-06 r73210)
## system   x86_64, linux-gnu
## ui       X11
## language en_US:en
## collate  en_US.UTF-8
## tz       Europe/Berlin
## date     2017-09-13
##
## package          * version    date       source
## acepack            1.4.1      2016-10-29 CRAN (R 3.5.0)
## animation          2.5        2017-03-30 cran (@2.5)
## assertthat         0.2.0      2017-04-11 CRAN (R 3.5.0)
## backports          1.1.0      2017-05-22 CRAN (R 3.5.0)
## base             * 3.5.0      2017-09-07 local
## base64enc          0.1-3      2015-07-28 CRAN (R 3.5.0)
## base64url          1.2        2017-06-14 CRAN (R 3.5.0)
## BatchExperiments * 1.4.1      2015-03-18 CRAN (R 3.5.0)
## BatchJobs        * 1.6        2015-03-18 CRAN (R 3.5.0)
## batchtools       * 0.9.6      2017-09-06 CRAN (R 3.5.0)
## BBmisc           * 1.11       2017-03-10 CRAN (R 3.5.0)
## bindr              0.1        2016-11-13 CRAN (R 3.5.0)
## bindrcpp         * 0.2        2017-06-17 CRAN (R 3.5.0)
## bit                1.1-12     2014-04-09 CRAN (R 3.5.0)
## bit64              0.9-7      2017-05-08 CRAN (R 3.5.0)
## blob               1.1.0      2017-06-17 CRAN (R 3.5.0)
## brew               1.0-6      2011-04-13 CRAN (R 3.5.0)
## checkmate        * 1.8.3      2017-07-03 CRAN (R 3.5.0)
## cluster            2.0.6      2017-03-10 CRAN (R 3.5.0)
## colorspace         1.3-2      2016-12-14 CRAN (R 3.5.0)
## compiler           3.5.0      2017-09-07 local
## cowplot            0.8.0      2017-07-30 CRAN (R 3.5.0)
## curl               2.8.1      2017-07-21 CRAN (R 3.5.0)
## data.table       * 1.10.4     2017-02-01 CRAN (R 3.5.0)
## datasets         * 3.5.0      2017-09-07 local
## DBI                0.7        2017-06-18 CRAN (R 3.5.0)
## devtools           1.13.3     2017-08-02 CRAN (R 3.5.0)
## digest             0.6.12     2017-01-27 CRAN (R 3.5.0)
## dlnm             * 2.3.2      2017-01-16 CRAN (R 3.5.0)
## dplyr            * 0.7.3      2017-09-09 CRAN (R 3.5.0)
## elrapack         * 0.0.3      2017-09-13 local (@0.0.3)
## expm               0.999-2    2017-03-29 CRAN (R 3.5.0)
## fail               1.3        2015-10-01 CRAN (R 3.5.0)
## foreign            0.8-69     2017-06-22 CRAN (R 3.5.0)
## Formula          * 1.2-2      2017-07-10 CRAN (R 3.5.0)
## gganimate        * 0.1.0.9000 2017-09-11 Github (dgrtwo/gganimate@bf82002)
## ggplot2          * 2.2.1      2016-12-30 CRAN (R 3.5.0)
## git2r              0.19.0     2017-07-19 CRAN (R 3.5.0)
## glue               1.1.1      2017-06-21 CRAN (R 3.5.0)
## graphics         * 3.5.0      2017-09-07 local
## grDevices        * 3.5.0      2017-09-07 local
## grid             * 3.5.0      2017-09-07 local
## gridExtra        * 2.3        2017-09-09 cran (@2.3)
## gtable             0.2.0      2016-02-26 CRAN (R 3.5.0)
## Hmisc            * 4.0-3      2017-05-02 CRAN (R 3.5.0)
## htmlTable          1.9        2017-01-26 CRAN (R 3.5.0)
## htmltools          0.3.6      2017-04-28 CRAN (R 3.5.0)
## htmlwidgets        0.9        2017-07-10 CRAN (R 3.5.0)
## httr               1.3.1      2017-08-20 CRAN (R 3.5.0)
## knitr              1.17       2017-08-10 CRAN (R 3.5.0)
## labeling           0.3        2014-08-23 CRAN (R 3.5.0)
## lattice          * 0.20-35    2017-03-25 CRAN (R 3.5.0)
## latticeExtra       0.6-28     2016-02-09 CRAN (R 3.5.0)
## lazyeval           0.2.0      2016-06-12 CRAN (R 3.5.0)
## lubridate        * 1.6.0      2016-09-13 CRAN (R 3.5.0)
## magrittr         * 1.5        2014-11-22 CRAN (R 3.5.0)
## Matrix             1.2-11     2017-08-21 CRAN (R 3.5.0)
## memoise            1.1.0      2017-04-21 CRAN (R 3.5.0)
## methods          * 3.5.0      2017-09-07 local
## mgcv             * 1.8-19     2017-09-01 CRAN (R 3.5.0)
## msm                1.6.4      2016-10-04 CRAN (R 3.5.0)
## munsell            0.4.3      2016-02-13 CRAN (R 3.5.0)
## mvtnorm            1.0-6      2017-03-02 CRAN (R 3.5.0)
## nlme             * 3.1-131    2017-02-06 CRAN (R 3.5.0)
## nnet               7.3-12     2016-02-02 CRAN (R 3.5.0)
## pammtools        * 0.0.3.2    2017-09-07 Github (adibender/pammtools@2f5a6d0)
## parallel         * 3.5.0      2017-09-07 local
## pkgconfig          2.0.1      2017-03-21 CRAN (R 3.5.0)
## plyr               1.8.4      2016-06-08 CRAN (R 3.5.0)
## prettyunits        1.0.2      2015-07-13 CRAN (R 3.5.0)
## progress           1.1.2      2016-12-14 CRAN (R 3.5.0)
## purrr            * 0.2.3      2017-08-02 CRAN (R 3.5.0)
## R6                 2.2.2      2017-06-17 CRAN (R 3.5.0)
## rappdirs           0.3.1      2016-03-28 CRAN (R 3.5.0)
## RColorBrewer       1.1-2      2014-12-07 CRAN (R 3.5.0)
## Rcpp               0.12.12    2017-07-15 CRAN (R 3.5.0)
## reshape2         * 1.4.2      2016-10-22 CRAN (R 3.5.0)
## rlang              0.1.2      2017-08-09 CRAN (R 3.5.0)
## rpart              4.1-11     2017-03-13 CRAN (R 3.5.0)
## RSQLite            2.0        2017-06-19 CRAN (R 3.5.0)
## rstudioapi         0.6        2016-06-27 CRAN (R 3.5.0)
## scales             0.5.0      2017-08-24 CRAN (R 3.5.0)
## sendmailR          1.2-1      2014-09-21 CRAN (R 3.5.0)
## splines            3.5.0      2017-09-07 local
## stats            * 3.5.0      2017-09-07 local
## stringi            1.1.5      2017-04-07 CRAN (R 3.5.0)
## stringr            1.2.0      2017-02-18 CRAN (R 3.5.0)
## survival         * 2.41-3     2017-04-04 CRAN (R 3.5.0)
## tables           * 0.8        2017-01-03 CRAN (R 3.5.0)
## tibble             1.3.4      2017-08-22 CRAN (R 3.5.0)
## tidyr            * 0.7.1      2017-09-01 CRAN (R 3.5.0)
## tidyselect         0.2.0      2017-08-30 CRAN (R 3.5.0)
## tools              3.5.0      2017-09-07 local
## tsModel          * 0.6        2013-06-24 CRAN (R 3.5.0)
## utils            * 3.5.0      2017-09-07 local
## withr              2.0.0      2017-07-28 CRAN (R 3.5.0)
## zoo                1.8-0      2017-04-12 CRAN (R 3.5.0)
```
