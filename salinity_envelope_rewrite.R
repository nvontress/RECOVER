# Rewriting Paul's salinity envelope code so I can work with it more easily

# Clear workspace
rm(list = ls())

# Load packages
library(tidyverse)
library(janitor)
library(lubridate)
library(zoo)
library(reshape2)

# Set WD
setwd('/Users/natalie/Documents/R projects/sccf/')

# Read in data ----
cepp <- bind_rows(read_csv('RSMBN_data/CEPP_RSMBN.csv'),read_csv('RSMBN_data/CEPP_RSMBN_S308QFC.csv')) %>% clean_names() 
lowrp <- bind_rows(read_csv('RSMBN_data/LOWRP_RSMBN.csv'),read_csv('RSMBN_data/LOWP_RSMBN_S308QFC.csv')) %>% clean_names() 

stl <- read_csv('RSMBN_data/gw_flws.csv') 
  
lowrp$alt[which(lowrp$alt == 'ALT1BWR')] <- 'TSP'
cepp$alt[which(cepp$alt == 'C240')] <- 'TSP'

# Tidy data ----

# Set thresholds 
# Generate 2007 results 
cre_breaks_2007 <- c(0,450,2800,4500,1000000) # last number is arbitrarily high
sle_breaks_2007 <- c(0,350,2000,3000,1000000) # last number is arbitrarily high
labels_2007 <- c('low','opt','high','very_high')

cre_breaks_2020 <- c(0,750,2100,2600,1000000) # last number is arbitrarily high
sle_breaks_2020 <- c(0,150,1400,1700,1000000) # last number is arbitrarily high
cre_labels_2020 <- c('low','opt','stress','damaging')
sle_labels_2020 <- c('low','opt','stress','damaging')



# write function

calculate_flow_exceedences <- function(proj) {
  # proj <- lowrp
  
  proj %>%
    filter(parameter == 'FLOW') %>%
    group_by(alt,parameter) %>% 
    count() %>% 
    mutate(qpfc_check = ifelse(n == 314475,T,F)) -> qpfc_check
  
  alts <- proj %>% pull(alt) %>% unique()
  
  # Create empty objects & run loop
  proj_flows <- NULL
  for(i in 1:length(alts)){
    # Set alternative
    alt_i <- alts[i]
    
    # Calculate results
    proj %>%
      filter(parameter == 'FLOW',
             alt == alt_i) %>%
      mutate(cy = year(date), month = month(date)) %>%
      pivot_wider(names_from = c(site),values_from = data_value) %>% 
      left_join(stl)-> proj_reshaped
    
    if(pull(filter(qpfc_check,alt == alt_i),qpfc_check) == T) {
      proj_reshaped %>% rename(S80_lake = S80_QPFCSOURCE_LAKE) %>% dplyr::select(-S308_QFC) -> proj_reshaped
    } else {proj_reshaped %>% rename(S80_lake = S308_QFC) -> proj_reshaped} 
    
    
    proj_reshaped %>%
      relocate(S80_lake, .before = S80) %>%
      # SLE row sums & determine basin column
      mutate(S80_trib = S80 + TMC2EST + S48 + S49 + NSF2EST + sle_gw) %>%
      # Calculate 14-d rolling averages
      mutate(S79_14d = rollmean(S79,k = 14,fill = NA,align = 'right'), # CRE Lake + Basin
             S79_lake_14d = rollmean(S79_QPFCSOURCE_LAKE,k = 14,fill = NA,align = 'right'), # CRE Basin-specific
             S79_basin_14d = S79_14d - S79_lake_14d,
             S80_14d = rollmean(S80,k = 14,fill = NA,align = 'right'), # Don't actually think I need this one
             S80_trib_14d = rollmean(S80_trib,k = 14,fill = NA,align = 'right'), # SLE Lake + Basin
             S80_lake_14d = rollmean(S80_lake,k = 14,fill = NA,align = 'right'), # SLE Basin-specific
             S80_basin_14d = S80_trib_14d - S80_lake_14d
      ) %>%
      mutate(S79_mm = ave(S79,cy,month), # CRE Lake + Basin
             S79_lake_mm = ave(S79_QPFCSOURCE_LAKE,cy,month), # CRE Basin-specific
             S79_basin_mm = S79_mm - S79_lake_mm,
             S80_mm = ave(S80,cy,month), # Don't actually think I need this one
             S80_trib_mm = ave(S80_trib,cy,month), # SLE Lake + Basin
             S80_lake_mm = ave(S80_lake,cy,month), # SLE Basin-specific
             S80_basin_mm = S80_trib_mm - S80_lake_mm) %>% bind_rows(proj_flows) %>%
      dplyr::select(date,parameter,alt,cy,month,S79_14d,S79_basin_14d,S79_lake_14d,S80_14d,S80_trib_14d,S80_basin_14d,S80_lake_14d,S79_mm,S79_basin_mm,S79_lake_mm,S80_mm,S80_trib_mm,S80_basin_mm,S80_lake_mm) -> proj_flows
  }
  
  # Check work
  # proj_flows %>% 
  #   group_by(alt) %>% count()
  
  # Original
  # tally_exceedences_2020 <- function(flows_14d,flow_cat) {
  #   ifelse(c(NA,diff(ifelse(flows_14d == flow_cat,1,0))) == 1 |
  #            c(0,diff(cumsumbinning(ifelse(flows_14d == flow_cat,1,0),
  #                                   threshold = 14))) == 1,1,0)
  # }
  
  # Updated
  # tally_exceedences_2020 <- function(flows_14d,flow_cat) {
  #   ifelse(
  #     c(NA,diff(ifelse(flows_14d %in% flow_cat,1,0))) == 1 | 
  #       c(NA,diff(ave(ifelse(flows_14d %in% flow_cat,1,0),
  #                     data.table::rleid(c(NA,diff(ifelse(flows_14d %in% flow_cat,1,0)))),
  #                     FUN = function(x) cumsumbinning(x,threshold = 13)))) == 1, 
  #     1,0
  #   ) %>% return()
  # }
  
  proj_flows %>% 
    # Categorize flow vals
    mutate(
      ## for replicating values in original project engineering reports
      # cre_flows_mm = cut(S79_mm,breaks = cre_breaks_2007,labels_2007,include.lowest = T),
      # sle_flows_mm = cut(S80_trib_mm,breaks = sle_breaks_2007,labels_2007,include.lowest = T),
      # sle_basin_14d = ifelse(cut(S80_basin_14d,breaks = c(-1,2000,1000000),labels = c('below','exceeds'),
      #                            include.lowest = T) == 'exceeds',1,0),
      # sle_lake_14d = ifelse(cut(S80_lake_14d,breaks = c(0,2000,1000000),labels = c('below','exceeds'),
      #                           include.lowest = T) == 'exceeds',1,0)
      ## Updated results
      cre_flows_14d = cut(S79_14d,breaks = cre_breaks_2020,labels = cre_labels_2020, right = F),
      cre_lake_14d = cut(S79_lake_14d,breaks = cre_breaks_2020,labels = cre_labels_2020, right = F),
      sle_flows_14d = cut(S80_trib_14d,breaks = sle_breaks_2020,labels = sle_labels_2020, right = F),
      sle_lake_14d = cut(S80_lake_14d,breaks = sle_breaks_2020,labels = sle_labels_2020, right = F)
    )  %>% 
    # Tally exceedences
    transmute(
      alt = alt,
      # For checking my work
      cre_flows_14d = cre_flows_14d,
      cre_lake_14d = cre_lake_14d,
      sle_flows_14d = sle_flows_14d,
      sle_lake_14d = sle_lake_14d,
      # CRE
      low_tally_cre = NA,
      opt_tally_cre = NA,
      stress_tally_cre = NA,
      stress_tally_cre_lake = NA,
      stress_tally_cre_basin = NA,
      damaging_tally_cre = NA,
      damaging_tally_cre_lake = NA,
      damaging_tally_cre_basin = NA,
      # SLE
      low_tally_sle = NA,
      opt_tally_sle = NA,
      stress_tally_sle = NA,
      stress_tally_sle_lake = NA,
      stress_tally_sle_basin = NA,
      damaging_tally_sle = NA,
      damaging_tally_sle_lake = NA,
      damaging_tally_sle_basin = NA,
      # SLE
    ) -> proj_flows
  
  tallied_flows <- read_csv("\n", col_names = colnames(proj_flows),col_types = 'cffffdddddddddddddddd')
  for(i in 1:length(alts)) {
    # i <- 1
    alt_i <- alts[i]
    proj_flows_i <- filter(proj_flows,alt == alt_i)
    
    for(j in 14:nrow(proj_flows_i)) {
      # CRE tally
      proj_flows_i$low_tally_cre[j] <- ifelse(proj_flows_i$cre_flows_14d[j] == 'low' & 
                                            sum(proj_flows_i$low_tally_cre[(j-13):(j-1)],na.rm = T)==0,
                                          1,0)
      proj_flows_i$opt_tally_cre[j] <- ifelse(proj_flows_i$cre_flows_14d[j] == 'opt' & 
                                            sum(proj_flows_i$opt_tally_cre[(j-13):(j-1)],na.rm = T)==0,
                                          1,0)
      proj_flows_i$stress_tally_cre[j] <- ifelse(proj_flows_i$cre_flows_14d[j] == 'stress' & 
                                               sum(proj_flows_i$stress_tally_cre[(j-13):(j-1)],na.rm = T)==0,
                                             1,0)
      proj_flows_i$stress_tally_cre_lake[j] <- ifelse(proj_flows_i$stress_tally_cre[j] == 1 & 
                                                        proj_flows_i$cre_lake_14d[j] == 'stress',1,0)
      proj_flows_i$stress_tally_cre_basin[j] <- proj_flows_i$stress_tally_cre[j] - proj_flows_i$stress_tally_cre_lake[j]
      proj_flows_i$damaging_tally_cre[j] <- ifelse(proj_flows_i$cre_flows_14d[j] == 'damaging' & 
                                                 sum(proj_flows_i$damaging_tally_cre[(j-13):(j-1)],na.rm = T)==0,
                                               1,0)
      proj_flows_i$damaging_tally_cre_lake[j] <- ifelse(proj_flows_i$damaging_tally_cre[j] == 1 &  
                                                          proj_flows_i$cre_lake_14d[j] == 'damaging',1,0)
      proj_flows_i$damaging_tally_cre_basin[j] <- proj_flows_i$damaging_tally_cre[j] - proj_flows_i$damaging_tally_cre_lake[j]
      # SLE tally
      proj_flows_i$low_tally_sle[j] <- ifelse(proj_flows_i$sle_flows_14d[j] == 'low' & 
                                                sum(proj_flows_i$low_tally_sle[(j-13):(j-1)],na.rm = T)==0,
                                              1,0)
      proj_flows_i$opt_tally_sle[j] <- ifelse(proj_flows_i$sle_flows_14d[j] == 'opt' & 
                                                sum(proj_flows_i$opt_tally_sle[(j-13):(j-1)],na.rm = T)==0,
                                              1,0)
      proj_flows_i$stress_tally_sle[j] <- ifelse(proj_flows_i$sle_flows_14d[j] == 'stress' & 
                                                   sum(proj_flows_i$stress_tally_sle[(j-13):(j-1)],na.rm = T)==0,
                                                 1,0)
      proj_flows_i$stress_tally_sle_lake[j] <- ifelse(proj_flows_i$stress_tally_sle[j] == 1 &  
                                                        proj_flows_i$sle_lake_14d[j] == 'stress',1,0)
      proj_flows_i$stress_tally_sle_basin[j] <- proj_flows_i$stress_tally_sle[j] - proj_flows_i$stress_tally_sle_lake[j]
      proj_flows_i$damaging_tally_sle[j] <- ifelse(proj_flows_i$sle_flows_14d[j] == 'damaging' & 
                                                     sum(proj_flows_i$damaging_tally_sle[(j-13):(j-1)],na.rm = T)==0,
                                                   1,0)
      proj_flows_i$damaging_tally_sle_lake[j] <- ifelse(proj_flows_i$damaging_tally_sle[j] == 1 &  
                                                          proj_flows_i$sle_lake_14d[j] == 'damaging',1,0)
      proj_flows_i$damaging_tally_sle_basin[j] <- proj_flows_i$damaging_tally_sle[j] - proj_flows_i$damaging_tally_sle_lake[j]
    }
    bind_rows(tallied_flows,proj_flows_i) -> tallied_flows
  }
    
  tallied_flows %>%
    group_by(alt) %>%
    summarize(
      # CRE
      low_cre = sum(low_tally_cre,na.rm = T),
      opt_cre = sum(opt_tally_cre,na.rm = T),
      stress_cre = sum(stress_tally_cre,na.rm = T),
      stress_cre_lake = sum(stress_tally_cre_lake,na.rm = T),
      stress_cre_basin = sum(stress_tally_cre_basin,na.rm = T),
      damaging_cre = sum(damaging_tally_cre,na.rm = T),
      damaging_cre_lake = sum(damaging_tally_cre_lake,na.rm = T),
      damaging_cre_basin = sum(damaging_tally_cre_basin,na.rm = T),
      # SLE
      low_sle = sum(low_tally_sle,na.rm = T),
      opt_sle = sum(opt_tally_sle,na.rm = T),
      stress_sle = sum(stress_tally_sle,na.rm = T),
      stress_sle_lake = sum(stress_tally_sle,na.rm = T),
      stress_sle_basin = sum(stress_tally_sle,na.rm = T),
      damaging_sle = sum(damaging_tally_sle,na.rm = T),
      damaging_sle_lake = sum(damaging_tally_sle_lake,na.rm = T),
      damaging_sle_basin = sum(damaging_tally_sle_basin,na.rm = T)
    )%>% return()
}

calculate_flow_exceedences(lowrp) %>% as.data.frame()
beep(2)
calculate_flow_exceedences(cepp) %>% View()

# Code from trying to recalc the results presented in the original PIR/PACR ----
# mutate(
#   cre_low_months = ifelse(cre_flows_mm == 'low',1,0),
  #   cre_high_months = ifelse(cre_flows_mm == 'high',1,0),
  #   cre_very_high_months = ifelse(cre_flows_mm == 'very_high',1,0),
  #   cre_tot_high_months = ifelse(cre_flows_mm == 'high' | cre_flows_mm == 'very_high',1,0),# high + very_high
  #   sle_low_months = ifelse(sle_flows_mm == 'low',1,0),
  #   sle_high_months = ifelse(sle_flows_mm == 'high',1,0),
  #   sle_very_high_months = ifelse(sle_flows_mm == 'very_high',1,0),
  #   # sle_high_basin = ave(sle_basin_14d,alt,
  #   #                       FUN = function(x) rollapply(x, FUN = function(x) sum(x, na.rm = T), width = 14, fill = NA, align = 'right')),
  #   sle_high_lake = ifelse(ave(sle_lake_14d,alt,
  #                       FUN = function(x) c(0,diff(x))) == 1,1,0) 
  # ) %>% filter(sle_lake_14d >0 | sle_high_lake > 0) %>% dplyr::select(date, alt, sle_lake_14d, sle_high_lake) %>% 
  # group_by(alt) %>%
  # summarize(
  #   cre_low_months = sum(ifelse(month %in% c(10:12,1:7),
  #     ifelse(day(date) == 1,ave(cre_low_months,cy,month,
  #                            FUN = function(x) ceiling(mean(x))),0),0)),
  #   cre_tot_high_months = sum(ifelse(day(date) == 1,ave(cre_tot_high_months,cy,month,
  #                                                          FUN = function(x) ceiling(mean(x))),0)),
  #   cre_very_high_months = sum(ifelse(day(date) == 1,ave(cre_very_high_months,cy,month,
  #                                                        FUN = function(x) ceiling(mean(x))),0)),
  #   sle_low_months = sum(ifelse(day(date) == 1,ave(sle_low_months,cy,month,
  #                                                  FUN = function(x) ceiling(mean(x))),0)),
  #   sle_high_14d_basins = sum(sle_high_lake)
  # )
    # cre_high_months = ifelse(cre_flows_mm == 'high',1,0),
    # cre_very_high_months = ifelse(cre_flows_mm == 'very_high',1,0),
    # cre_tot_high_months = ifelse(cre_flows_mm == 'high' | cre_flows_mm == 'very_high',1,0),# high + very_high
    # sle_low_months = ifelse(sle_flows_mm == 'low',1,0),
    # sle_high_months = ifelse(sle_flows_mm == 'high',1,0),
    # sle_very_high_months = ifelse(sle_flows_mm == 'very_high',1,0),
    # sle_high_basin = ave(sle_basin_14d,alt,
    #                      FUN = function(x) rollapply(x, FUN = function(x) sum(x, na.rm = T), width = 14, by = 14,fill = NA)),
    # sle_high_lake = ave(sle_lake_14d,alt,
    #                     FUN = function(x) rollapply(x, FUN = function(x) sum(x, na.rm = T), width = 14, by = 14,fill = NA))
  
