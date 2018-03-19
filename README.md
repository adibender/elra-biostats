# Evaluation of the association between nutritional adequacy and survival
[![DOI](https://zenodo.org/badge/103416524.svg)](https://zenodo.org/badge/latestdoi/103416524)

This is the Code and Data repository for: <br>

Andreas Bender, Fabian Scheipl, Wolfgang Hartl, Andrew G Day, Helmut Küchenhoff; *Penalized estimation of complex, non-linear exposure-lag-response associations*, Biostatistics, , kxy003, https://doi.org/10.1093/biostatistics/kxy003

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

- We currently develop an R package that facilitates working with PAMMs, including data preparation, visualization, etc.. There are also a lot of vignettes with application examples: http://github.com/adibender/pammtools

### Folder structure
- **`data`**: Raw data for the application example (after initial import
from SAS and minor preprocessing)
- **`dataGenerationScripts`**: Contains scripts for (further) data preprocessing. Creates folder **`dataCurrent`** and **`dataCurrentHosp`** (storing data for main and sensitivity an analysis, respectively).
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
##  [1] tables_0.8             Hmisc_4.0-3            Formula_1.2-2
##  [4] lattice_0.20-35        pec_2.5.4              reshape2_1.4.3
##  [7] survival_2.41-3        prodlim_1.6.1          tidyr_0.7.2
## [10] gridExtra_2.3          pammtools_0.0.3.2      purrr_0.2.4
## [13] magrittr_1.5           batchtools_0.9.6       data.table_1.10.4-3
## [16] ggplot2_2.2.1          tsModel_0.6            dlnm_2.3.2
## [19] bindrcpp_0.2           dplyr_0.7.4            BatchExperiments_1.4.1
## [22] BatchJobs_1.7          BBmisc_1.11            mgcv_1.8-19
## [25] nlme_3.1-131           checkmate_1.8.3        elrapack_0.0.3
##
## loaded via a namespace (and not attached):
##  [1] bit64_0.9-7         splines_3.5.0       foreach_1.4.3
##  [4] modelr_0.1.1        assertthat_0.2.0    expm_0.999-2
##  [7] latticeExtra_0.6-28 base64url_1.2       blob_1.1.0
## [10] progress_1.1.2      timereg_1.9.1       numDeriv_2016.8-1
## [13] RSQLite_2.0         backports_1.1.1     glue_1.2.0
## [16] digest_0.6.12       RColorBrewer_1.1-2  colorspace_1.3-2
## [19] htmltools_0.3.6     cowplot_0.8.0       Matrix_1.2-11
## [22] plyr_1.8.4          psych_1.7.5         pkgconfig_2.0.1
## [25] broom_0.4.2         mvtnorm_1.0-6       scales_0.5.0
## [28] brew_1.0-6          lava_1.5            htmlTable_1.9
## [31] tibble_1.3.4        withr_2.1.0         nnet_7.3-12
## [34] lazyeval_0.2.0      mnormt_1.5-5        memoise_1.1.0
## [37] msm_1.6.5           foreign_0.8-69      tools_3.5.0
## [40] prettyunits_1.0.2   stringr_1.2.0       sendmailR_1.2-1
## [43] munsell_0.4.3       cluster_2.0.6       compiler_3.5.0
## [46] rlang_0.1.4         iterators_1.0.8     htmlwidgets_0.9
## [49] rappdirs_0.3.1      base64enc_0.1-3     labeling_0.3
## [52] gtable_0.2.0        codetools_0.2-15    DBI_0.7
## [55] R6_2.2.2            zoo_1.8-0           knitr_1.17
## [58] bit_1.1-12          bindr_0.1           stringi_1.1.6
## [61] Rcpp_0.12.14        rpart_4.1-11        acepack_1.4.1
## [64] tidyselect_0.2.3


devtools::session_info()

##Session info ------------------------------------------------------------------
## setting  value
## version  R Under development (unstable) (2017-09-06 r73210)
## system   x86_64, linux-gnu
## ui       X11
## language en_US:en
## collate  en_US.UTF-8
## tz       Europe/Berlin
## date     2017-12-24
##
##Packages ----------------------------------------------------------------------
## package          * version  date       source
## acepack            1.4.1    2016-10-29 CRAN (R 3.5.0)
## assertthat         0.2.0    2017-04-11 CRAN (R 3.5.0)
## backports          1.1.1    2017-09-25 cran (@1.1.1)
## base             * 3.5.0    2017-09-07 local
## base64enc          0.1-3    2015-07-28 CRAN (R 3.5.0)
## base64url          1.2      2017-06-14 CRAN (R 3.5.0)
## BatchExperiments * 1.4.1    2015-03-18 CRAN (R 3.5.0)
## BatchJobs        * 1.7      2017-11-28 cran (@1.7)
## batchtools       * 0.9.6    2017-09-06 CRAN (R 3.5.0)
## BBmisc           * 1.11     2017-03-10 CRAN (R 3.5.0)
## bindr              0.1      2016-11-13 CRAN (R 3.5.0)
## bindrcpp         * 0.2      2017-06-17 CRAN (R 3.5.0)
## bit                1.1-12   2014-04-09 CRAN (R 3.5.0)
## bit64              0.9-7    2017-05-08 CRAN (R 3.5.0)
## blob               1.1.0    2017-06-17 CRAN (R 3.5.0)
## brew               1.0-6    2011-04-13 CRAN (R 3.5.0)
## broom              0.4.2    2017-02-13 CRAN (R 3.5.0)
## checkmate        * 1.8.3    2017-07-03 CRAN (R 3.5.0)
## cluster            2.0.6    2017-03-10 CRAN (R 3.5.0)
## codetools          0.2-15   2016-10-05 CRAN (R 3.5.0)
## colorspace         1.3-2    2016-12-14 CRAN (R 3.5.0)
## compiler           3.5.0    2017-09-07 local
## cowplot            0.8.0    2017-07-30 CRAN (R 3.5.0)
## data.table       * 1.10.4-3 2017-10-27 cran (@1.10.4-)
## datasets         * 3.5.0    2017-09-07 local
## DBI                0.7      2017-06-18 CRAN (R 3.5.0)
## devtools           1.13.3   2017-08-02 CRAN (R 3.5.0)
## digest             0.6.12   2017-01-27 CRAN (R 3.5.0)
## dlnm             * 2.3.2    2017-01-16 CRAN (R 3.5.0)
## dplyr            * 0.7.4    2017-09-28 cran (@0.7.4)
## elrapack         * 0.0.3    2017-12-13 local (@0.0.3)
## expm               0.999-2  2017-03-29 CRAN (R 3.5.0)
## foreach            1.4.3    2015-10-13 CRAN (R 3.5.0)
## foreign            0.8-69   2017-06-22 CRAN (R 3.5.0)
## Formula          * 1.2-2    2017-07-10 CRAN (R 3.5.0)
## ggplot2          * 2.2.1    2016-12-30 CRAN (R 3.5.0)
## glue               1.2.0    2017-10-29 cran (@1.2.0)
## graphics         * 3.5.0    2017-09-07 local
## grDevices        * 3.5.0    2017-09-07 local
## grid             * 3.5.0    2017-09-07 local
## gridExtra        * 2.3      2017-09-09 cran (@2.3)
## gtable             0.2.0    2016-02-26 CRAN (R 3.5.0)
## Hmisc            * 4.0-3    2017-05-02 CRAN (R 3.5.0)
## htmlTable          1.9      2017-01-26 CRAN (R 3.5.0)
## htmltools          0.3.6    2017-04-28 CRAN (R 3.5.0)
## htmlwidgets        0.9      2017-07-10 CRAN (R 3.5.0)
## iterators          1.0.8    2015-10-13 CRAN (R 3.5.0)
## knitr              1.17     2017-08-10 CRAN (R 3.5.0)
## labeling           0.3      2014-08-23 CRAN (R 3.5.0)
## lattice          * 0.20-35  2017-03-25 CRAN (R 3.5.0)
## latticeExtra       0.6-28   2016-02-09 CRAN (R 3.5.0)
## lava               1.5      2017-03-16 CRAN (R 3.5.0)
## lazyeval           0.2.0    2016-06-12 CRAN (R 3.5.0)
## magrittr         * 1.5      2014-11-22 CRAN (R 3.5.0)
## Matrix             1.2-11   2017-08-21 CRAN (R 3.5.0)
## memoise            1.1.0    2017-04-21 CRAN (R 3.5.0)
## methods          * 3.5.0    2017-09-07 local
## mgcv             * 1.8-19   2017-09-01 CRAN (R 3.5.0)
## mnormt             1.5-5    2016-10-15 CRAN (R 3.5.0)
## modelr             0.1.1    2017-07-24 CRAN (R 3.5.0)
## msm                1.6.5    2017-12-05 cran (@1.6.5)
## munsell            0.4.3    2016-02-13 CRAN (R 3.5.0)
## mvtnorm            1.0-6    2017-03-02 CRAN (R 3.5.0)
## nlme             * 3.1-131  2017-02-06 CRAN (R 3.5.0)
## nnet               7.3-12   2016-02-02 CRAN (R 3.5.0)
## numDeriv           2016.8-1 2016-08-27 CRAN (R 3.5.0)
## pammtools        * 0.0.3.2  2017-12-10 Github (adibender/pammtools@2f5a6d0)
## parallel         * 3.5.0    2017-09-07 local
## pec              * 2.5.4    2017-08-08 CRAN (R 3.5.0)
## pkgconfig          2.0.1    2017-03-21 CRAN (R 3.5.0)
## plyr               1.8.4    2016-06-08 CRAN (R 3.5.0)
## prettyunits        1.0.2    2015-07-13 CRAN (R 3.5.0)
## prodlim          * 1.6.1    2017-03-06 CRAN (R 3.5.0)
## progress           1.1.2    2016-12-14 CRAN (R 3.5.0)
## psych              1.7.5    2017-05-03 CRAN (R 3.5.0)
## purrr            * 0.2.4    2017-10-18 cran (@0.2.4)
## R6                 2.2.2    2017-06-17 CRAN (R 3.5.0)
## rappdirs           0.3.1    2016-03-28 CRAN (R 3.5.0)
## RColorBrewer       1.1-2    2014-12-07 CRAN (R 3.5.0)
## Rcpp               0.12.14  2017-11-23 cran (@0.12.14)
## reshape2         * 1.4.3    2017-12-11 cran (@1.4.3)
## rlang              0.1.4    2017-11-05 cran (@0.1.4)
## rpart              4.1-11   2017-03-13 CRAN (R 3.5.0)
## RSQLite            2.0      2017-06-19 CRAN (R 3.5.0)
## scales             0.5.0    2017-08-24 CRAN (R 3.5.0)
## sendmailR          1.2-1    2014-09-21 CRAN (R 3.5.0)
## splines            3.5.0    2017-09-07 local
## stats            * 3.5.0    2017-09-07 local
## stringi            1.1.6    2017-11-17 cran (@1.1.6)
## stringr            1.2.0    2017-02-18 CRAN (R 3.5.0)
## survival         * 2.41-3   2017-04-04 CRAN (R 3.5.0)
## tables           * 0.8      2017-01-03 CRAN (R 3.5.0)
## tibble             1.3.4    2017-08-22 CRAN (R 3.5.0)
## tidyr            * 0.7.2    2017-10-16 cran (@0.7.2)
## tidyselect         0.2.3    2017-11-06 cran (@0.2.3)
## timereg            1.9.1    2017-05-21 CRAN (R 3.5.0)
## tools              3.5.0    2017-09-07 local
## tsModel          * 0.6      2013-06-24 CRAN (R 3.5.0)
## utils            * 3.5.0    2017-09-07 local
## withr              2.1.0    2017-11-01 cran (@2.1.0)
## zoo                1.8-0    2017-04-12 CRAN (R 3.5.0)

```
