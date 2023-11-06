
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Pacific NCRMP reefbudgetR

<!-- badges: start -->
<!-- badges: end -->

This R package provides tools for working with ReefBudget carbonate
budget data collected by NOAA's National Coral Reef Monitoring Program (NCRMP) in the U.S. Pacific Islands, with functions to process field-based and SfM-derived
benthic, urchin, and parrotfish census data and calculate carbonate
production and erosion rates. Package tools and processes are translated to R based on Indo-Pacific ReefBudget methdology and materials (Perry et al. 2018) and modified for use with Pacific NCRMP data.

For additional information on data analyses and methodological
approaches, see: Hannah C. Barkley, Rebecca M. Weible, Ariel A.
Halperin, Candace E. Alagata, Tye L. Kindinger, Damaris Torres-Pulliza,
Mia S. Lamirand, Brittany E. Huntington, Courtney S. Couch, Corinne G.
Amir, Nicole I. Besemer, Jonathan A. Charendoff, Jon Ehrenberg, Joao D.
Garriques, Andrew E. Gray, Nathan Hayes, Kurt E. Ingeman, Lori H. Luers,
Kaylyn S. McCoy, Noah V. Pomeroy, Joy N. Smith, Bernardo Vargas-Ángel,
Erica K. Towle, Jennifer C. Samson. 2023. Carbonate budget assessments
in the U.S. Pacific Islands: report of methods comparison results and
summary of standard operating procedures. U.S. Dept. of Commerce, NOAA
Technical Memorandum NMFS-PIFSC-##, p. <doi:10>… \[UPDATE when
published\].

For additional metadata and downloadable data, see:
<https://www.fisheries.noaa.gov/inport/item/67804>.

For additional information on ReefBudget, see: https://www.exeter.ac.uk/research/projects/geography/reefbudget/

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

You can install the development version of Pacific NCRMP reefbudgetR from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hannahbarkley/reefbudgetR")
```

## Example Workflow

Load benthic, urchin, and fish data sets:

``` r
library(reefbudgetR)

# Benthic observation data
benthic_data <- read.csv("ESD_CarbBudget_Benthic_OAHU_2021.csv", na = "", check.names = FALSE)

#Urchin observation data
urchin_data <- read.csv("ESD_CarbBudget_Urchins_OAHU_2021.csv", na = "", check.names = FALSE)

# Fish belt transect data
fish_data_belt <- read.csv("ESD_CarbBudget_Fixed_Belt_OAHU_2021.csv", na = "", check.names = FALSE)

# Fish fixed-site and stratified random stationary point count (SPC) data
fish_data_spc <- read.csv("ESD_CarbBudget_SPC_OAHU_2021.csv", na = "", check.names = FALSE)
```

Process benthic data by method type:

``` r
# Process Indo-Pacific ReefBudget (IPRB) benthic data
prod_iprb <- process_prod(
  data = benthic_data[benthic_data$CB_METHOD == "IPRB", ],
  method_name = "IPRB"
)

# Process NCRMP-intermediate (chords) benthic data
prod_chords <- process_prod(
  data = benthic_data[benthic_data$CB_METHOD == "Chords", ],
  method_name = "Chords"
)

# Process NCRMP-leveraged (Structure-from-Motion, SfM) benthic data
prod_sfm <- process_prod(
  data = benthic_data[benthic_data$CB_METHOD == "SfM", ],
  method_name = "SfM"
)

# Combine site-level production data
prod_site <- bind_rows(
  prod_iprb$summary_site,
  prod_chords$summary_site,
  prod_sfm$summary_site
)
```

Process urchin data by method type:

``` r
# Process Indo-Pacific ReefBudget (IPRB) urchin data
urch_iprb <- process_urchins(
  data = urchin_data[urchin_data$CB_METHOD == "IPRB", ],
  method_name = "IPRB"
)

# Process NCRMP-intermediate (chords) urchin data
urch_chords <- process_urchins(
  data = urchin_data[urchin_data$CB_METHOD == "Chords", ],
  method_name = "Chords"
)

# Process NCRMP-leveraged (Structure-from-Motion, SfM) urchin data
urch_sfm <- process_urchins(
  data = urchin_data[urchin_data$CB_METHOD == "SfM", ],
  method_name = "SfM"
)

# Combine site-level urchin erosion data
urch_site <- bind_rows(urch_iprb$site_erosion,
                       urch_chords$site_erosion,
                       urch_sfm$site_erosion)
```

Process fish data by method type:

``` r
# Process fish belt and stationary point count (SPC) data, open plot window to see strs shapefile plots
fish_site_ <- process_fish(spc_data = fish_data_spc,
                               dbase_type = "Kindinger",
                               belt_data = fish_data_belt,
                               subset_distance_m = 6000)
fish_site <- fish_site_$dat
fish_site_$assoc_site_count
```

Combine data and calculate net production rates by methodology:

``` r
# Calculate net production using IPRB methodology
net_site_iprb <- process_net(
  prod = prod_site[prod_site$CB_METHOD == "IPRB" ,],
  urch = urch_site[urch_site$CB_METHOD == "IPRB" ,],
  fish = fish_site[fish_site$CB_METHOD == "IPRB" ,],
  sum_by = "site"
)

net_site_iprb$METHOD <- "IPRB"

# Calculate net production using NCRMP-intermediate methodology
net_site_int <- process_net(
  prod = prod_site[prod_site$CB_METHOD == "Chords" ,],
  urch = urch_site[urch_site$CB_METHOD == "Chords" ,],
  fish = fish_site[fish_site$CB_METHOD == "Fixed SPC" ,],
  sum_by = "site"
)

net_site_int$METHOD <- "NCRMP-intermediate"

# Calculate net production using NCRMP-leveraged methodology
net_site_lev <- process_net(
  prod = prod_site[prod_site$CB_METHOD == "SfM" ,],
  urch = urch_site[urch_site$CB_METHOD == "SfM" ,],
  fish = fish_site[fish_site$CB_METHOD == "StRS SPC" ,],
  sum_by = "site"
)

net_site_lev$METHOD <- "NCRMP-leveraged"

# Calculate net production using NCRMP-proposed methodology
net_site_prop <- process_net(
  prod = prod_site[prod_site$CB_METHOD == "SfM" ,],
  urch = urch_site[urch_site$CB_METHOD == "Chords" ,],
  fish = fish_site[fish_site$CB_METHOD == "Fixed SPC" ,],
  sum_by = "site"
)

net_site_prop$METHOD <- "NCRMP-proposed"

# Combine net production data
net_site <- bind_rows(net_site_iprb, net_site_int, net_site_lev, net_site_prop)
```
