library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(sampling)
library(survey)
library(purrr)
library(tidycensus)

pv22 <- pums_variables %>% 
  filter(year == 2022, survey == "acs1") %>%
  distinct(survey, year, var_code, var_label, data_type, level) %>%
  filter(is.na(level) | level == "person")

WAdat <- get_pums(
  variables=c("SERIALNO", "PUMA", "PWGTP", "AGEP", "ADJINC", 
              "JWMNP", "WKHP", "PERNP", "PINCP", "CIT", "SCHL", "ESR"), 
  state="WA", 
  year=2022, 
  # rep_weights="person",
  survey="acs1", 
  recode=T
) %>%
  mutate(across(c(PERNP, PINCP), ~.*as.numeric(ADJINC)),
         CIT = ifelse(CIT == "5", 0, 1),
         degree = ifelse(SCHL_label >= "Associate's degree", 1, 0)) %>%
  select(serial=SERIALNO, puma=PUMA, pweight=PWGTP, age=AGEP, transit=JWMNP, 
         workhrs=WKHP, earnings=PERNP, income=PINCP, citizen=CIT,
         degree, employment=ESR) %>%
  filter(!(employment %in% c("6", "b"))) #drop those not in the labour force

save(WAdat, file="ACS Analysis/WAdat.rda")

  
  
  
  
  
  
  
  
  



