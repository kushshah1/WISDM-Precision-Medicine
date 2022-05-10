### dataset_creation.R
### Collects variables from different datasets and compiles them
### working directory = "Diabetes/Repository"

library(tidyverse)
select <- dplyr::select
`%notin%` <- Negate(`%in%`)

PtRoster <- read.table("data/study data/PtRoster.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
DiabScreening <- read.table("data/study data/DiabScreening.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
DiabSocioEcon <- read.table("data/study data/DiabSocioEcon.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16", quote="\"")
DiabLocalHbA1c <- read.table("data/study data/DiabLocalHbA1c.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
DiabPhysExam <- read.table("data/study data/DiabPhysExam.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
gluIndicesRCT <- read.table("data/study data/gluIndices RCT.txt", header = TRUE, sep = "|") %>%
  mutate(gluBelow70 = parse_number(as.character(gluBelow70)))

# Table 1: Baseline Characteristics
  # [PtRoster] Age
  # [DiabScreening] Diabetes duration (use DiagAge and subtract Age from PtRoster from it)
  # [DiabScreening] Age at diagnosis
  # [DiabScreening] Sex
  # [DiabScreening] Race/Ethnicity
  # [DiabSocioEcon] Annual HH income
  # [DiabSocioEcon] Highest education
  # [DiabSocioEcon] Health insurance
  # [DiabScreening$CGMUseStat] CGM use (past vs never)
  # [DiabScreening$PumpUse] Insulin pump use (both InsModPump and PumpUse have 111 users in DiabScreening whereas paper says 58+50)
  # [DiabLocalHbA1c] Screening Hba1c (filter for Visit=="Screening")
  # [????] Hba1c at randomization
  # [????] Detectable C-peptide
  # [DiabScreening$UnitsInsTotal / DiabPhysExam$Weight] Total daily insulin doses per kg (11 missing, matches paper exactly - 97+95)
  # [DiabScreening] Severe hypoglycemia event in past 12 mo
  # [DiabScreening] DKA event in past 12 mo
  # Functional activities questionnaire score
  # NIH toolbox cognition battery age-corrected fluid composite score
  # Cognition status measured by NIH toolbox cognition battery
  # Reduced hypo awareness
  # Wear hearing aids regularly
  # Near vision card last line worse than 20/40

# eTable 16: Subgroup Analysis by Baseline Characteristics
  # [already covered] Age
  # [already covered] T1D Duration
  # [already covered] Age at diagnosis
  # [gluIndicesRCT$gluBelow70 where period == "1) Baseline", time == "1) Overall"] Percent time <70
  # [gluIndicesRCT$gluCV where period == "1) Baseline", time == "1) Overall"] Coefficient of Variation
  # [?????] HbA1c (screening/randomization/different?)
  # [already covered] Gender
  # [already covered] Race/ethnicity
  # [already covered] Highest education completed
  # [already covered] NIH Toolbox age-corrected fluid composite score
  # [already covered] Cognition status measured by the NIH toolbox
  # [already covered] Hypoglycemic awareness
  # [already covered] Detectable C-peptide
  # [already covered] Severe hypo event past 12 months

# patients <- DiabScreening$PtID #213 
patients <- (gluIndicesRCT %>% filter(period == "1) Baseline" & time == "1) Overall"))$PtID # 203

# Variable Collection
dat_raw_X <- select(PtRoster, RecID, PtID, SiteID, AgeAsOfEnrollDt) %>%
  filter(PtID %in% patients) %>%
  left_join(select(DiabScreening, PtID, DiagAge, Gender, Ethnicity, Race, CGMUseStat, PumpUse, SHNumLast12Months, DKANumLast12Months), by = "PtID") %>%
  left_join(select(DiabSocioEcon, PtID, EducationLevel, AnnualIncome, InsPrivate, InsMedicare), by = "PtID") %>%
  left_join(select(DiabLocalHbA1c %>% filter(Visit == "Screening"), PtID, HbA1cTestRes), by = "PtID") %>%
  left_join(select(DiabScreening, PtID, UnitsInsTotal), by = "PtID") %>%
  left_join(select(DiabPhysExam %>% filter(Visit == "Screening"), PtID, Weight, WeightUnits), by = "PtID") %>%
  left_join(select(gluIndicesRCT %>% filter(period == "1) Baseline" & time == "1) Overall"), PtID, gluBelow70, gluCV), by = "PtID")

# Data Manipulation
less_bach <- c("12th grade - no diploma", "Associate Degree (AA)",
               "High school graduate/diploma/GED", 
               "Some college but no degree")
bach <- c("Bachelor's Degree (BS,BA,AB)")
greater_bach <- c("Doctorate Degree (PhD, EdD)", 
                  "Master's Degree (MA, MS, MSW, MBA, MPH)", 
                  "Professional Degree (MD, DDS, DVM, LLB, JD)")

dat_clean_X <- dat_raw_X %>% 
  mutate(SHNumLast12Months_binary = ifelse(SHNumLast12Months == "0", 0, 1),
         DKANumLast12Months_binary = ifelse(DKANumLast12Months == 0, 0, 1),
         Weight_Kg = ifelse(WeightUnits == "kg", Weight, Weight / 2.20462),
         DiabDuration = AgeAsOfEnrollDt - DiagAge,
         PumpUse_Binary = ifelse(PumpUse == "", 0, 1),
         InsulinDosesKg = UnitsInsTotal / Weight_Kg,
         EducationLevel_Tri = ifelse(EducationLevel %in% less_bach, "<Bach",
                                     ifelse(EducationLevel %in% bach, "Bach",
                                            ifelse(EducationLevel %in% greater_bach, ">Bach", NA))),
         Ins_Combined = ifelse(is.na(InsPrivate), "Medicare/other",
                               ifelse(!is.na(InsPrivate) & is.na(InsMedicare), "Private",
                                      ifelse(!is.na(InsPrivate) & !is.na(InsMedicare), "Private and Medicare", NA))),
         gluCV_num = readr::parse_number(as.character(gluCV))) %>%
  select(!c(RecID, SiteID, Ethnicity, Race, PumpUse, SHNumLast12Months, DKANumLast12Months, EducationLevel, AnnualIncome, InsPrivate, InsMedicare, UnitsInsTotal, Weight, WeightUnits, gluCV, Weight_Kg))

dat_clean_XAY <- dat_clean_X %>%
  left_join(select(PtRoster, PtID, TrtGroup), by = "PtID") %>%
  left_join(select(gluIndicesRCT %>% filter(period == "2) Follow-up (26 week)" & time == "1) Overall"), PtID, gluBelow70Chg), by = "PtID") %>%
  mutate(gluBelow70Chg = -gluBelow70Chg) %>% # Changed so that "better" values (more negative gluBelow70Chg) are positive
  column_to_rownames("PtID")

dat_clean_XAY_full <- na.omit(dat_clean_XAY)
#saveRDS(dat_clean_XAY_full, "./data/dat_clean_XAY_full.rds")