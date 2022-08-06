### C00data_processing.R
### Compiles relevant variables from various data tables
### Working Directory = "WISDM-Precision-Medicine/"
### Author: Kushal Shah

library(tidyverse)
select <- dplyr::select
`%notin%` <- Negate(`%in%`)

### Data Table Readin ###
  PtRoster <- read.table("data/study data/PtRoster.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
  DiabScreening <- read.table("data/study data/DiabScreening.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
  DiabSocioEcon <- read.table("data/study data/DiabSocioEcon.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16", quote="\"")
  DiabLocalHbA1c <- read.table("data/study data/DiabLocalHbA1c.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
  DiabPhysExam <- read.table("data/study data/DiabPhysExam.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
  gluIndicesRCT <- read.table("data/study data/gluIndices RCT.txt", header = TRUE, sep = "|") %>%
    mutate(gluBelow70 = parse_number(as.character(gluBelow70))) %>%
    #mutate(gluAbove180 = parse_number(as.character(gluAbove180))) %>%
    mutate(gluInRange = parse_number(as.character(gluInRange)))
  STASampleResults <- read.table("data/study data/STASampleResults.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
  cpep_info <- STASampleResults %>% filter(Visit == "Randomization", ResultName == "CPEP", STAResultStatus != "Canceled")

### 203 patients from Pratley et al. ###
  patients <- (gluIndicesRCT %>% filter(period == "1) Baseline" & time == "1) Overall"))$PtID

### Compilation of Relevant Variables ###
  dat_raw_X <- select(PtRoster, RecID, PtID, SiteID, AgeAsOfEnrollDt) %>%
    filter(PtID %in% patients) %>%
    left_join(select(DiabScreening, PtID, DiagAge, Gender, PumpUse, SHNumLast12Months, DKANumLast12Months), by = "PtID") %>%
    left_join(select(DiabSocioEcon, PtID, EducationLevel, AnnualIncome, InsPrivate, InsMedicare), by = "PtID") %>%
    left_join(select(DiabLocalHbA1c %>% filter(Visit == "Screening"), PtID, HbA1cTestRes), by = "PtID") %>%
    left_join(select(gluIndicesRCT %>% filter(period == "1) Baseline" & time == "1) Overall"), PtID, gluBelow70, gluCV, gluInRange), by = "PtID") %>%
    left_join(select(cpep_info, PtID, Value), by = c("PtID"))

### Data Manipulation and Feature Creation ###
  less_bach <- c("12th grade - no diploma", "Associate Degree (AA)",
                 "High school graduate/diploma/GED", 
                 "Some college but no degree")
  bach <- c("Bachelor's Degree (BS,BA,AB)")
  greater_bach <- c("Doctorate Degree (PhD, EdD)", 
                    "Master's Degree (MA, MS, MSW, MBA, MPH)", 
                    "Professional Degree (MD, DDS, DVM, LLB, JD)")

  # Clean Feature Set
    # Categorical variables dichotomized or trichotomized
    # Final variables generated from raw variables
  dat_clean_X <- dat_raw_X %>% 
    mutate(SHNumLast12Months_binary = ifelse(SHNumLast12Months == "0", 0, 1),
           DKANumLast12Months_binary = ifelse(DKANumLast12Months == 0, 0, 1),
           DiabDuration = AgeAsOfEnrollDt - DiagAge,
           PumpUse_Binary = ifelse(PumpUse == "", 0, 1),
           EducationLevel_Tri = ifelse(EducationLevel %in% less_bach, "<Bach",
                                       ifelse(EducationLevel %in% bach, "Bach",
                                              ifelse(EducationLevel %in% greater_bach, ">Bach", NA))),
           Ins_Combined = ifelse(is.na(InsPrivate), "Medicare/other",
                                 ifelse(!is.na(InsPrivate) & is.na(InsMedicare), "Private",
                                        ifelse(!is.na(InsPrivate) & !is.na(InsMedicare), "Private and Medicare", NA))),
           gluCV_num = readr::parse_number(as.character(gluCV)),
           CPep_detected = ifelse(Value == "<0.003", "no", "yes")) %>%
    select(!c(RecID, SiteID, PumpUse, SHNumLast12Months, DKANumLast12Months, EducationLevel, AnnualIncome, 
            InsPrivate, InsMedicare, gluCV, Value))
  
### Addition of treatment and outcome information from the RCT ###
  dat_clean_XAY <- dat_clean_X %>%
    left_join(select(PtRoster, PtID, TrtGroup), by = "PtID") %>%
    left_join(select(gluIndicesRCT %>% filter(period == "2) Follow-up (26 week)" & time == "1) Overall"), PtID, gluBelow70Chg), by = "PtID") %>%
    mutate(gluBelow70Chg = -gluBelow70Chg) %>% # Changed so that "better" values (more negative gluBelow70Chg) are positive
    column_to_rownames("PtID")
  
### Omission of patients with missing data ###
  # 7 lost to follow-up, requested withdrawal from study, or discontinued intervention
  # 2 missing education data
  dat_clean_XAY_full <- na.omit(dat_clean_XAY)
  
#saveRDS(dat_clean_XAY_full, "./data/dat_clean_XAY_full.rds")
  