Process carbonate budgets data with reefbudgetR
================

``` r
# Set working directory containing data
setwd("T:/Oceanography/Carbonate Budgets/Data/")

# Load benthic, urchin, and fish data
benthic_data <- read.csv("CB_Benthic_alldata.csv", na = "", check.names = FALSE)

urchin_data <- read.csv("CB_Urchin_alldata.csv", na = "", check.names = FALSE)

fish_data <- read.csv("CB_FishBelt_alldata.csv", na = "", check.names = FALSE)
```

``` r
#### Process benthic data by method type

# Process "Indo-Pacific ReefBudget (IPRB)" methodology benthic data
prod_IPRB <- process_prod(
  data = benthic_data[benthic_data$CB_METHOD == "IPRB", ],
  transect_id = c("A1", "A2", "A3", "B1", "B2", "B3"),
  transect_length = c(10, 10, 10, 10, 10, 10),
  dbase_type = "NCRMP",
  method_name = "IPRB",
  data_type = "In water",
  full_summary = TRUE
)

# Process "Chords" methodology benthic data
prod_chords <- process_prod(
  data = benthic_data[benthic_data$CB_METHOD == "Chords", ],
  transect_id = c("A", "B", "C", "D", "E", "F"),
  transect_length = c(8, 9, 10, 10, 9, 8),
  dbase_type = "NCRMP",
  method_name = "Chords",
  data_type = "In water",
  full_summary = TRUE
)

# Process "SfM" methodology benthic data
prod_sfm <- process_prod(
  data = benthic_data[benthic_data$CB_METHOD == "SfM", ],
  transect_id = c("A", "B", "C", "D", "E", "F"),
  transect_length = c(8, 9, 10, 10, 9, 8),
  dbase_type = "NCRMP",
  method_name = "Chords",
  data_type = "SfM",
  full_summary = TRUE
)

# Combine site-level production data back together
prod_site <- bind_rows(
  prod_IPRB$summary_site,
  prod_chords$summary_site,
  prod_sfm$summary_site
)
```

``` r
#### Process urchin data by method type

# Process IPRB urchin data
urch_IPRB <- process_urchins(
  data = urchin_data[urchin_data$CB_METHOD == "IPRB", ],
  method_name = "IPRB",
  data_type = "In water",
  transect_id = c("A1", "A2", "A3", "B1", "B2", "B3"),
  transect_length = c(10, 10, 10, 10, 10, 10),
  full_summary = TRUE
)

# Process chords urchin data
urch_chords <- process_urchins(
  data = urchin_data[urchin_data$CB_METHOD == "Chords", ],
  method_name = "Chords",
  data_type = "In water",
  transect_id = c("A", "B", "C", "D", "E", "F"),
  transect_length = c(8, 9, 10, 10, 9, 8),
  full_summary = TRUE
)


# Combine site-level urchin erosion data back together
urch_site <- bind_rows(urch_IPRB$site_erosion,
                                 urch_chords$site_erosion)
```

``` r
# Proess fish data
fish_belt <- process_fish(
  data = fish_data, 
  rates_dbase = "Kindinger", 
  full_summary = TRUE)

fish_site <- fish_belt$fish_erosion_site
```

``` r
# Process net production data
net_site = process_net(
  prod = prod_site,
  urch = urch_site,
  fish = fish_site
)
```
