## Supporting information for Assessing the use of environmental flows and stocking for the persistence of a flow-dependent spawner in a drying climate


### Details

This repository contains code to recreate all analyses presented in Stoessel et al. (submitted), Assessing the use of environmental flows and stocking for the persistence of a flow-dependent spawner in a drying climate.


### Usage

The entire analysis is contained within `main.R` and may take several minutes to run. Several helper scripts are provided in the `R` directory and these are sourced directly by `main.R`.

The full analysis requires several packages not available on CRAN. The `aae.pop` and `aae.pop.templates` packages define the Australian bass population model and provide methods to simulate and summarise population dynamics. The `aae.hydro` package provides methods to download and work with discharge and water temperature data from the Victorian WMIS. All `aae.x` packages are available at https://github.com/aae-stats and can be installed directly from GitHub with the `remotes` or `devtools` packages.

Several additional directories will need to be created to store model outputs: `data`, `outputs/figs`, `outputs/simulations`, and `outputs/tables`.


### Contact

For additional information, please contact Jian Yen (jian.yen [at] delwp.vic.gov.au).

