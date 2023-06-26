#'Calculate Kindinger erosion rate metrics: bite rate, volume, and proportion of scars
#'
#'@author Rebecca Weible
#'
#'@param kindinger_equations_csv .csv file in Kindinger et al AppendixS2 final format,
#'downloaded as is.
#'@param kindinger_functional_groups_csv .csv file in NCRMP herbivore functional groups format,
#'downloaded as is.
#'@param sister_species_csv .csv file with parrot fish sister species assignments,
#'downloaded as is.
#'@param SST_Input_Celsius average sea surface temperature in celsius from region data was 
#'collected.
#'
#'#'@import tidyverse
#'@import dplyr
#'@importFrom rlang .data
#'
#'@export calc_kindinger_erosion_metrics



calc_kindinger_erosion_metrics <- function(kindinger_equations_csv, kindinger_functional_groups_csv, sister_species_csv, SST_Input_Celsius) {
  
  # clean up and format Kindinger database equations with TL and SST constants to be used in the calculation
  format_equations <- kindinger_equations_csv %>%
                        filter(!grepl("Acanthurus|Melichthys|Zebrasoma|Naso|Ctenochaetus|Siganus", FISH_sciname)) %>% # remove non-parrots from Kindinger equations database
                        filter(grepl("Scarus|Chlorurus|Calotomus|NA", FISH_genus)) %>% # remove non-parrots genus level rows from Kindinger equations database
                        mutate("I.11-20cm" = NA,
                               "I.21-30cm" = NA,
                               "I.31-40cm" = NA,
                               "I.41-50cm" = NA,
                               "I.51-60cm" = NA,
                               "T.11-20cm" = NA,
                               "T.21-30cm" = NA,
                               "T.31-40cm" = NA,
                               "T.41-50cm" = NA,
                               "T.51-60cm" = NA) %>% #add in size and phase bins from belt survey methodology
                        gather(., sizebin, value, -c(RESP_var:SST_rng)) %>% 
                        separate(sizebin, c("phase", "sizeclass"), extra = "merge", fill = "left") %>% #pull phase letters out of size bin column (J, I, or T)
                        mutate(TL_INPUT = case_when(
                          sizeclass == "11-20cm" ~ 15,
                          sizeclass == "21-30cm" ~ 25,
                          sizeclass == "31-40cm" ~ 35,
                          sizeclass == "41-50cm" ~ 45,
                          sizeclass == "51-60cm" ~ 55
                        )) %>% #create TL_INPUT column where the average size for each bin will be substituted for TL to calculate Bite Rate, Bite Volume, and Proporotion of Scars in the next step
                        select(-value) %>%
                        mutate(SST_INPUT = SST_Input_Celsius) %>% #input SST in celcius for the location your data was collected (~79 deg F or 26 deg C)
                        mutate_at(vars("SST_INPUT"), as.numeric) %>% #convert all columns that should be numeric
                        mutate(FISH_sciname = replace(FISH_sciname, FISH_genus == "Chlorurus", "Chlorurus sp")) %>% # rename genus level to include " sp" in the scientific name column
                        mutate(FISH_sciname = replace(FISH_sciname, FISH_genus == "Scarus", "Scarus sp")) %>% # rename genus level to include " sp" in the scientific name column
                        mutate(FISH_sciname = replace(FISH_sciname, FISH_genus == "Calotomus", "Calotomus sp")) %>% # rename genus level to include " sp" in the scientific name column
                        
                        # add functional group assignments from separate .csv file
                        left_join(kindinger_functional_groups_csv %>% select(TAXONNAME, Herb_fxn, Herb_fxn2, starts_with("Fxn_")), by = c("FISH_sciname" = "TAXONNAME")) %>% # here, fish sci name is called "FISH_sciname" in Kindinger equations database & "TAXONNAME" in kndinger functional groups database              
                        mutate(FXN_grp = case_when(FXN_grp == "NA" ~ Herb_fxn,
                                                   TRUE ~ FXN_grp)) %>% # copy over Herb_fxn assignments to new column (FXN_grp) that will be the final assignments
                        # update assignments based on the size of individual fish observed for species that shift their functional groups at a size threshold
                        mutate(FXN_grp = ifelse(!is.na(Herb_fxn2) & Fxn_size_include == "Y" & TL_INPUT > Fxn_size_cmTL, Herb_fxn2, FXN_grp)) 
                        
                        
                        # check for any species that weren't assigned a functional grouping
                          #format_equations %>% filter(is.na(FXN_grp))
            
                        # manually assign any missing groups
                          #format_equations <- format_equations %>% 
                          #mutate(FXN_grp = ifelse(FISH_sciname == "Kyphosus sp.", "Browser", FXN_grp))
  
  
  
  # calculate metrics (RESP_var: Bite Rate, Bite Volume, Bite Area, Proportion of Scars) for each row in kindinger_equations
  calc_metrics <- format_equations %>%
                    mutate(SST_degC = replace(SST_degC, SST_degC == "NA", 0)) %>% # replace NA values in select columns to 0 so that we can do the calculations for all species
                    mutate_at(vars("SST_degC"), as.numeric) %>%
                    mutate(METRIC = case_when((RESP_var == "bite rate" & SST_degC != "0") ~ (exp(a + (TL_INPUT * b) + (SST_INPUT * SST_degC))),
                                              RESP_var == "prop scars" ~ ((exp(a + (TL_INPUT * b)))/(1 + exp(a + (TL_INPUT * b)))),
                                              TRUE ~ (a * (TL_INPUT^b)))) %>%
                    # clean up
                    filter(RESP_var != "bite rate" | SST_degC != "0") %>% # remove combination of bite rate rows with NA bite analysis
                    select(RESP_var, FXN_grp:FISH_sciname, phase, sizeclass, METRIC) %>%
                    distinct() %>%
                    spread(., sizeclass, METRIC, fill = "T") %>%
                    mutate_at(vars(c("RESP_var":"phase")), as.factor) %>% # change specific columns to factors
                    mutate_if(is.character, as.numeric) %>% # make numeric values b/c need to do math in a sec
                    mutate(across(where(is.numeric), round, 5)) %>%
                    right_join(., sister_species_csv %>% select(c(SPECIES, FISH_sciname, sister_sp)), by = "FISH_sciname") %>%
                    mutate(FISH_genus = str_extract(FISH_sciname, "(\\w+)")) %>% # copy Genus from sciname and paste into FISH_genus column
                    mutate(phase = replace_na(phase, "I")) %>%
                    gather(., "sizeclass", "value", -c(RESP_var:phase, SPECIES:sister_sp)) %>% # gather size bins into single column
                    spread(., RESP_var, value, fill = "0") %>% # spread by three metrics for calculations down the road
                    select(-`<NA>`, -`bite area`) %>%
                    mutate_at(vars(c("bite rate":"prop scars")), as.numeric) %>% # convert characters to numeric for next steps
                    complete(., sizeclass, phase, nesting(FXN_grp, FISH_genus, FISH_sciname, SPECIES, sister_sp)) #fill in missing phase and sizeclass for select species
                
  # Subset missing species' Functional Group assignments and fill in with those of it's sister species 
  subset_missing_fxngrp <- calc_metrics %>%
                            #subset unique functional group values for species that are assigned functional groups
                            filter(!is.na(FXN_grp)) %>% 
                            filter_at(vars(`bite rate`:`prop scars`), any_vars(. != 0)) %>%
                            distinct(sister_sp, phase, sizeclass, FXN_grp) %>% 
                            dplyr::rename(FXN_grp_replace = FXN_grp) %>% 
                            # join the unique functional groups subsetted above with the species that do not have functional group assignments
                            left_join(calc_metrics %>% filter(is.na(FXN_grp)), . , by = c("sister_sp", "phase", "sizeclass")) %>% 
                            #clean up
                            select(-FXN_grp) %>% dplyr::rename(FXN_grp = FXN_grp_replace) %>% 
                            select(sizeclass, FXN_grp, everything(.)) %>% 
                            mutate_at(vars(`bite rate`:`prop scars`), ~replace_na(., 0)) 

  fill_missing_fxngrp <- calc_metrics %>% filter(!is.na(FXN_grp)) %>% # remove old rows where FXN_grp is NA
                            rbind(., subset_missing_fxngrp) %>% # attach new rows where FXN_grp were assigned from sister species (from fill_missing_fxngrp)
                            #clean up
                            mutate_at(vars(`bite rate`, `bite vol`, `prop scars`), ~replace_na(., 0)) %>%
                            arrange(FISH_sciname) %>%
                            dplyr::rename(`Bite Rate` = `bite rate`, `Bite Volume` = `bite vol`, `Proportion of Scars` = `prop scars`) %>%
                            unite("sizeclass", c(phase,sizeclass)) %>% #combine phase and size range
                            gather(., "METRIC", "VALUE", -c(sizeclass:sister_sp)) %>% #columns to row values
                            cbind(., REPLACE = .$FISH_sciname)  #create duplicate fish species column that we name "REPLACE" for later step
                          
  
                            # check for any species that weren't assigned a functional grouping
                                #fill_missing_fxngrp %>% filter(is.na(FXN_grp)) %>% distinct(.$FISH_sciname)
  
  # Fill in missing metric values with sister species metrics (RESP_var: Bite Rate, Bite Volume, Bite Area, Proportion of Scars)
  fill_missing_metrics <- sister_species_csv %>%
                            #clean up 
                            select(SPECIES, FISH_sciname, sister_sp, c(replace_biterate:replace_scars)) %>%
                            gather(., "METRIC", "REPLACE", -c(SPECIES:sister_sp)) %>%
                            mutate(METRIC = replace(METRIC, METRIC == "replace_biterate", "Bite Rate")) %>%
                            mutate(METRIC = replace(METRIC, METRIC == "replace_volume", "Bite Volume")) %>%
                            mutate(METRIC = replace(METRIC, METRIC == "replace_scars", "Proportion of Scars")) %>%
                            #join sister species assignments for missing metrics and functional groups
                            left_join(fill_missing_fxngrp %>% select(-REPLACE), ., by = c("FISH_sciname", "SPECIES", "sister_sp", "METRIC")) %>% # merge cleaned up sister species assignments with cleaned and calculated metrics dataframe
                            left_join(., fill_missing_fxngrp %>% select(sizeclass, METRIC, VALUE, REPLACE, FXN_grp), by = c("sizeclass", "METRIC", "REPLACE")) %>% 
                            #replace functional groups for genus level metrics that have both Scraper and Excavator assignments                       
                            mutate(FXN_grp.x = case_when((FXN_grp.x != "Browser" & str_detect(REPLACE, " sp")) ~ FXN_grp.y,
                                                                                 TRUE ~ FXN_grp.x)) %>% 
                            #assign sister species metric values where they are missing
                             mutate(VALUE.y = case_when(REPLACE == ""  ~ VALUE.x,
                                                                               TRUE ~ VALUE.y)) %>%
                            #clean up
                            dplyr::rename(VALUE = VALUE.y,
                                   FXN_grp = FXN_grp.x) %>%
                            select(-VALUE.x, -FXN_grp.y, -REPLACE, -sister_sp) %>%
                            distinct(.) %>% #remove duplicate rows
                            filter_at(vars(VALUE), any_vars(. != 0)) %>% #remove 0 values in "VALUE" column that indicates no metric values exist
                            select(SPECIES, FISH_genus, FISH_sciname, FXN_grp, everything(.)) %>%
                            spread(METRIC, VALUE, fill = "0")
  
  return(fill_missing_metrics)
  
}
