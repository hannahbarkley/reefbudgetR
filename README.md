
<!-- README.md is generated from README.Rmd. Please edit that file -->

# reefbudgetR

<!-- badges: start -->
<!-- badges: end -->

This R package provides tools for working with ReefBudget carbonate
budget data, with functions to process field-based and SfM-derived
benthic, urchin, and parrotfish census data and calculate carbonate
production and erosion.

For additional information on data analyses and methodological
approaches, see:

[Barkley, H.C., Halperin, A.A., Torres-Pulliza, D., Lamirand, M.S.,
Couch, C.S., Alagata, C.E., Weible, R.M., Smith, J.N., Oliver, T.A. and
Perry, C.T., 2025. Estimating coral reef carbonate budgets using
Structure-from-Motion photogrammetry. Coral Reefs, 44(3),
pp.937-951.](https://link.springer.com/article/10.1007/s00338-025-02660-7)

[Barkley HC, Weible RM, Halperin AA, Alagata CE, Kindinger TL,
Torres-Pulliza D, Lamirand MS, Huntington BE, Couch CS, Amir CG, et
al. 2023. Carbonate budget assessments in the U.S. Pacific Islands:
report of methods comparison results and summary of standard operating
procedures. U.S. Dept. of Commerce, NOAA Technical Memorandum
NMFS-PIFSC-154, 7979 p. doi:
10.25923/g4hg-7686.](https://repository.library.noaa.gov/view/noaa/56372/noaa_56372_DS1.pdf)

For additional metadata and downloadable data, see the [NCRMP carbonate
budgets InPort
record](https://www.fisheries.noaa.gov/inport/item/67804).

This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an ‘as is’ basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.

## Installation

You can install the development version of reefbudgetR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hannahbarkley/reefbudgetR")
```

## Example workflow for processing Hawaii data

Load benthic, urchin, and fish observation data sets:

``` r
library(reefbudgetR)

# Benthic observation data
benthic_data <- read.csv("ESD_NCRMP_CarbBudget_Benthic_HA_2024.csv", na = "", check.names = FALSE)

#Urchin observation data
urchin_data <- read.csv("ESD_NCRMP_CarbBudget_Urchin_HA_2024.csv", na = "", check.names = FALSE)

# Fish SPC data
fish_data <- read.csv("ESD_NCRMP_CarbBudget_Fish_HA_2024.csv", na = "", check.names = FALSE)
```

Process benthic data:

``` r

# Process benthic production data
prod <- process_prod(data = benthic_data)

# Summarize overall production data at site level
prod_site <- prod$summary_site

# Summarize overall production data at transect level
prod_transect <- prod$summary_transect

# Summarize production data at site level by coral group
prod_coral <- prod$summary_site_coral

# Summarize production data at site level by substrate group
prod_substrate <- prod$summary_site_substrateclass

# Return a summary of site metadata
sites_metadata <- prod$sites_metadata
```

Process urchin data:

``` r

 urch <- process_urchins(data = urchin_data)

# Summarize urchin erosion data at site level
 urch_site <- urch$site_erosion
 
 # Summarize urchin erosion data at transect level
 urch_transect <- urch$transect_erosion
 
 # Summarize urchin erosion data at transect level by taxon
 urch_taxon <- urch$transect_taxon
```

Process fish data by method type:

``` r

# Process fish erosion data for all sites
fish_site <- process_fish(data = spc_data, method_type = "Fixed", fixed_metadata = sites_metadata) 
```

Calculate net production rates :

``` r

# Calculate net production at transect site
net_site <- process_net(prod = prod_transect, urch = urch_transect, fish = fish_site, sum_by = "transect")
```
