# Calculating formulas in Tye's Database to create Kindinger bioerosion values

dbase <- read.csv("C:/Users/Rebecca.Weible/Downloads/Kindinger et al_AppendixS2_FINAL.csv", na = "", check.names = FALSE, skip=1)
sp.fxn <- read.csv("C:/Users/Rebecca.Weible/Downloads/NCRMP_herb_fxn.grps.csv", na = "", check.names = FALSE)
GrazTypes <- read.csv("C:/Users/Rebecca.Weible/Desktop/OCC/Carbonate Budgets/Data Entry/CPM to R conversion/Input/GrazingTypesClassifications_Combo.csv", check.names = F) # call in grazing types .csv to assign grazing types for surveyed species with assigned sister sp



# Bring in bioerosion equations and calculate metrics for PROP SCARS, BITE RATES, BITE VOLUMES

format <- dbase %>%
            filter(!grepl("Acanthurus|Melichthys|Zebrasoma|Naso|Ctenochaetus|Siganus", FISH_sciname)) %>%
            filter(grepl("Scarus|Chlorurus|Calotomus|NA", FISH_genus)) %>%
            mutate("I.11-20cm" = NA,
                   "I.21-30cm" = NA,
                   "I.31-40cm" = NA,
                   "I.41-50cm" = NA,
                   "I.51-60cm" = NA,
                   "T.11-20cm" = NA,
                   "T.21-30cm" = NA,
                   "T.31-40cm" = NA,
                   "T.41-50cm" = NA,
                   "T.51-60cm" = NA) %>%
            gather(., sizebin, value, -c(RESP_var:SST_rng)) %>%
            separate(sizebin, c("phase", "sizeclass"), extra = "merge", fill = "left") %>% #pull phase letters out of sizebin column (J, I, or T)
            mutate(TL_INPUT = case_when(
              sizeclass == "11-20cm" ~ 15,
              sizeclass == "21-30cm" ~ 25,
              sizeclass == "31-40cm" ~ 35,
              sizeclass == "41-50cm" ~ 45,
              sizeclass == "51-60cm" ~ 55
            )) %>% #create column where the average size will be used for each bin for the TL to FL conversion in a later step
            select(-value) %>%
            mutate(SST_INPUT = "26") %>% #input SST of Hawaii (~79 deg F or 26 deg C)
            mutate_at(vars("SST_INPUT"), as.numeric) %>% #convert all columns that should be numeric
            mutate(FISH_sciname = replace(FISH_sciname, FISH_genus == "Chlorurus", "Chlorurus sp")) %>%
            mutate(FISH_sciname = replace(FISH_sciname, FISH_genus == "Scarus", "Scarus sp")) %>%
            mutate(FISH_sciname = replace(FISH_sciname, FISH_genus == "Calotomus", "Calotomus sp"))

#Attach functional group assignments
data.fxn <- format %>% # any fish data with scientific names
              # add fxn assignments
              left_join(sp.fxn %>% select(TAXONNAME, Herb_fxn, Herb_fxn2, starts_with("Fxn_")), by = c("FISH_sciname" = "TAXONNAME")) %>% # here, fish sci name is called "FISH_sciname" in my "data" dataframe & "TAXONNAME" in my "sp.fxn" dataframe
              # copy over Herb_fxn assignments to new column (FXN_grp) that will be the final assignments
              mutate(FXN_grp = case_when(FXN_grp == "NA" ~ Herb_fxn,
                                         TRUE ~ FXN_grp)) %>% 
              # update assignments based on the size of individual fish observed for species that shift their functional groups at a size threshold
              mutate(FXN_grp = ifelse(!is.na(Herb_fxn2) & Fxn_size_include == "Y" & TL_INPUT > Fxn_size_cmTL, Herb_fxn2, FXN_grp)) 

            # check for any species that weren't assigned a functional grouping
            data.fxn %>% filter(is.na(FXN_grp))
            
            #data.fxn.final <- data.fxn %>% 
                                # manually assign any missing groups
                                #mutate(FXN_grp = ifelse(FISH_sciname == "Kyphosus sp.", "Browser", FXN_grp))
                                  
   
  calc <- data.fxn %>%
            mutate(SST_degC = replace(SST_degC, SST_degC == "NA", 0)) %>% # replace NA values in select columns to 0 so that we can do the calculations for all species in Tye's list
            mutate_at(vars("SST_degC"), as.numeric) %>%
            mutate(METRIC = case_when((RESP_var == "bite rate" & SST_degC != "0") ~ (exp(a + (TL_INPUT * b) + (SST_INPUT * SST_degC))),
                                      RESP_var == "prop scars" ~ ((exp(a + (TL_INPUT * b)))/(1 + exp(a + (TL_INPUT * b)))),
                                      TRUE ~ (a * (TL_INPUT * b)))) %>%
            filter(RESP_var != "bite rate" | SST_degC != "0") %>% # remove combination of bite rate rows with NA bite analysis
            select(RESP_var, FXN_grp:FISH_sciname, phase, sizeclass, METRIC) %>%
            spread(., sizeclass, METRIC, fill = "T") %>%
            mutate_at(vars(c("RESP_var":"phase")), as.factor) %>% # change specific columns to factors
            mutate_if(is.character, as.numeric) %>% # make numeric values b/c need to do math in a sec
            mutate(across(where(is.numeric), round, 5)) %>%
            right_join(., GrazTypes %>% select(c(SPECIES, FISH_sciname, sister_sp)), by = "FISH_sciname")
    
    
    
    
  
# format calc output from formulas in Tye's database and replace sister species
  
formatcalc <- calc %>%
                mutate(FISH_genus = str_extract(FISH_sciname, "(\\w+)")) %>% # copy Genus from sciname and paste into FISH_genus column
                mutate(phase = replace_na(phase, "I")) %>%
                gather(., "sizeclass", "value", -c(RESP_var:phase, SPECIES:sister_sp)) %>% # gather size bins into single column
                spread(., RESP_var, value, fill = "0") %>% # spread by three metrics for calculations down the road
                select(-`<NA>`, -`bite area`) %>%
                mutate_at(vars(c("bite rate":"prop scars")), as.numeric) %>% # convert characters to numeric for next steps
                complete(., sizeclass, phase, nesting(FXN_grp, FISH_genus, FISH_sciname, SPECIES, sister_sp)) #complete phase and sizeclass for select species

# code below for filling in FXN_grp NA values with those of it's sister sp 
to_rbind <- formatcalc %>%
                filter(!is.na(FXN_grp)) %>%
                filter_at(vars(`bite rate`:`prop scars`), any_vars(. != 0)) %>%
                distinct(sister_sp, phase, sizeclass, FXN_grp) %>%
                rename(FXN_grp_replace = FXN_grp) %>% 
                left_join(formatcalc %>% filter(is.na(FXN_grp)), . , by = c("sister_sp", "phase", "sizeclass")) %>%
                select(-FXN_grp) %>% rename(FXN_grp = FXN_grp_replace) %>%
                select(sizeclass, FXN_grp, everything(.)) %>%
                mutate_at(vars(`bite rate`:`prop scars`), ~replace_na(., 0))


formatcalcfinal <- formatcalc %>% filter(!is.na(FXN_grp)) %>%
                      rbind(., to_rbind) %>%
                      mutate_at(vars(`bite rate`, `bite vol`, `prop scars`), ~replace_na(., 0)) %>%
                      arrange(FISH_sciname) %>%
                      rename(`Bite Rate` = `bite rate`, `Bite Volume` = `bite vol`, `Proportion of Scars` = `prop scars`) %>%
                      unite("sizeclass", c(phase,sizeclass)) %>% #combine phase and size range
                      # replace 0 values in bite rate/ bite volume/ bite scars with sister species specified by Grazing Types Classifications .csv file 
                      gather(., "METRIC", "VALUE", -c(sizeclass:sister_sp)) %>% #columns to row values
                      cbind(., REPLACE = .$FISH_sciname)  #create duplicate fish species column that we name "REPLACE" for later step

                
                # check for any species that weren't assigned a functional grouping
                formatcalcfinal %>% filter(is.na(FXN_grp)) %>% distinct(.$FISH_sciname)


calc.rep.sissp <- GrazTypes %>%
                      select(SPECIES, FISH_sciname, sister_sp, c(replace_biterate:replace_scars)) %>%
                      gather(., "METRIC", "REPLACE", -c(SPECIES:sister_sp)) %>%
                      mutate(METRIC = replace(METRIC, METRIC == "replace_biterate", "Bite Rate"),
                             METRIC = replace(METRIC, METRIC == "replace_volume", "Bite Volume"),
                             METRIC = replace(METRIC, METRIC == "replace_scars", "Proportion of Scars")) %>%
                      left_join(formatcalcfinal %>% select(-REPLACE), ., by = c("FISH_sciname", "SPECIES", "sister_sp", "METRIC")) %>% # merge cleaned up Grazing Types Classifications file with Tye's added on Metrics
                      left_join(., formatcalcfinal %>% select(sizeclass, METRIC, VALUE, REPLACE, FXN_grp), by = c("sizeclass", "METRIC", "REPLACE")) %>% # merge the metric values of the sister sp to replace as a separate column
                                            mutate(FXN_grp.x = case_when((FXN_grp.x != "Browser" & str_detect(REPLACE, " sp")) ~ FXN_grp.y,
                                                   TRUE ~ FXN_grp.x)) %>% # replace functional groups for species where replacement species has both a Scraper and Excavator value associated with them
                      mutate(VALUE.y = case_when(REPLACE == ""  ~ VALUE.x,
                                                 TRUE ~ VALUE.y)) %>% # replace the Metric values that are zero with the new values in VALUES.y column for select sp
                      rename(VALUE = VALUE.y,
                             FXN_grp = FXN_grp.x) %>%
                      select(-VALUE.x, -FXN_grp.y, -REPLACE, -sister_sp) %>%
                      distinct(.) %>% #remove duplicate rows
                      filter_at(vars(VALUE), any_vars(. != 0)) %>% #remove 0 values in "VALUE" column that indicates no metric values exist
                      select(SPECIES, FISH_genus, FISH_sciname, FXN_grp, everything(.)) %>%
                      spread(METRIC, VALUE, fill = "0")



# Calculate Final Erosion values from Bite Rate, Bite Volume, and Proportion of Scars

ErosionRateCalculation <- function(df, substratedensity = 1.47, percentofdayfeeding = 83.3) {
                            
                            df %>%
                              mutate_at(vars(`Bite Rate`:`Proportion of Scars`), as.numeric) %>%
                              mutate(Bites_Scars_min = `Bite Rate` * `Proportion of Scars`) %>% # Calculate Bites Leaving Scars per Minute
                              mutate(Volume_perday = Bites_Scars_min * `Bite Volume` * 60 * (12 * (percentofdayfeeding/100))) %>% # Calculate Volume Removed per Day
                              mutate(MassRemoved_perday = (Volume_perday * substratedensity)/1000) %>% # Calculate Mass removed per day (converted g to kg with /1000)
                              mutate(MassRemoved_peryear = MassRemoved_perday * 365) %>% # Calculate Mass (kg) removed per year
                              mutate(ErosionRates = MassRemoved_peryear) # Final erosion rates = mass (kg) removed per year
                          }   


test <- ErosionRateCalculation(calc.rep.sissp) %>%
          select(c(sizeclass:SPECIES), ErosionRates) %>%
          separate(., sizeclass, c("PHASE", "SIZE_CLASS"), sep = "_") %>%
          rename(TAXON_NAME = "FISH_sciname",
                 EROSION_RATE = "ErosionRates",
                 FXN_GRP = "FXN_grp") %>%
          select(-FISH_genus, -SPECIES)










  
  
  
  filter_at(vars(value), all_vars(!is.na(.))) %>% # remove rows where NA values are present in value column, because these are duplicates from Scraper/Excavator assignment
    