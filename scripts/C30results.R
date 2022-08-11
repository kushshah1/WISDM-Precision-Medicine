### C30results.R
### Compilation of results for manuscript
### Working Directory = "WISDM-Precision-Medicine/"
### Author: Kushal Shah

library(stats)
library(tidyverse)
DiabScreening <- read.table("data/study data/DiabScreening.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
DiabPhysExam <- read.table("data/study data/DiabPhysExam.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")

dat_clean_XAY_full <- readRDS("./data/dat_clean_XAY_full.rds")
dat_opt <- readRDS("./data/dat_opt.rds")

### Figure 1: Differences in treatment effect of CGM versus BGM on hypo reduction, by baseline %CV ###
  ggplot(aes(x = gluCV_num, y = 100*gluBelow70Chg, color = TrtGroup), data = dat_clean_XAY_full) + 
    geom_point() + 
    stat_smooth(method = "lm", formula = y ~ poly(x,2), se = FALSE) +
    xlab("Coefficient of Variation") +
    ylab("Decrease in % Time Under 70 mg/dL") +
    ggtitle("CGM vs. BGM Treatment Effect is Moderated by Baseline CV%")
  
  # summary(lm(gluBelow70Chg ~ gluCV_num, data = dat_clean_XAY_full))
  # summary(lm(gluBelow70Chg ~ poly(gluCV_num, 2), data = dat_clean_XAY_full))

### Table 1: Characteristics of Study Participants, Stratified by Decision Rule Subgroup ###

  # 1. Addition of extra variables for characterization of subgroups
  dat_opt <- rownames_to_column(dat_opt) %>%
    mutate(rowname = as.numeric(rowname)) %>%
    left_join(select(DiabScreening, PtID, Ethnicity, Race, UnitsInsTotal), by = c("rowname" = "PtID")) %>%
    left_join(select(DiabPhysExam %>% filter(Visit == "Screening"), PtID, Weight, WeightUnits, Height, HeightUnits), by = c("rowname" = "PtID")) %>%
    mutate(Weight_Kg = ifelse(WeightUnits == "kg", Weight, Weight / 2.20462)) %>%
    mutate(Height_M = ifelse(HeightUnits == "in", Height*0.0254, Height/100)) %>%
    mutate(BMI = Weight_Kg / (Height_M)^2) %>%
    mutate(Race_binary = ifelse(Race == "White", "White", ifelse(Race == "Unknown/not reported", NA, "Non-White"))) %>%
    mutate(Ethnicity_binary = ifelse(Ethnicity == "Not Hispanic or Latino", "Not Hispanic or Latino", 
                                     ifelse(Ethnicity == "Hispanic or Latino", "Hispanic or Latino", NA))) %>%
    mutate(InsulinDosesKg = UnitsInsTotal / Weight_Kg) %>%
    mutate(count = 1) %>%
    select(!c(Ethnicity, Race, Weight, WeightUnits, Height, HeightUnits, Weight_Kg, Height_M))

  # 2. Summarize characterizing variable averages (for continuous variables) by trt group assigned by opt decision rule
  
    # Collection of continuous variables
    cts_vars <- c("AgeAsOfEnrollDt", "DiagAge", "HbA1cTestRes", "gluBelow70", "gluInRange", 
                  "DiabDuration", "gluCV_num", "BMI", "InsulinDosesKg")
    
    dat_opt_cts <- dat_opt %>%
      select(rowname, opt, count, all_of(cts_vars))
    
    # Summary table with average/SD of selected variables
    opt_table_cts <- dat_opt_cts %>%
      group_by(opt) %>%
      summarise(n = n(),
                meanAgeAsOfEnrollDt = mean(AgeAsOfEnrollDt),
                sdAge = sd(AgeAsOfEnrollDt),
                meanDiagAge = mean(DiagAge),
                sdDiag = sd(DiagAge),
                meanHbA1cTestRes = mean(HbA1cTestRes),
                sdHbA1c = sd(HbA1cTestRes),
                meangluBelow70 = mean(gluBelow70),
                sdgluBelow70 = sd(gluBelow70),
                meangluInRange = mean(gluInRange),
                sdgluInRange = sd(gluInRange),
                meanDiabDuration = mean(DiabDuration),
                sdDiabDuration = sd(DiabDuration),
                meangluCV_num = mean(gluCV_num),
                sdgluCV = sd(gluCV_num),
                meanBMI = mean(BMI, na.rm = TRUE),
                sdBMI = sd(BMI, na.rm = TRUE),
                meanInsulinDosesKg = mean(InsulinDosesKg, na.rm = TRUE),
                sdInsulinDosesKg = sd(InsulinDosesKg, na.rm = TRUE))
  
    # Transposing and addition of column headers
    opt_table_cts_t <- t(opt_table_cts)[-c(1,2), ]
    opt_table_cts_t <- round(opt_table_cts_t, 3)
    opt_table_cts_t <- cbind(opt_table_cts_t[, 2], opt_table_cts_t[, 1])
    colnames(opt_table_cts_t) <- c("CGM (n = 173)", "BGM (n = 21)")
    
    # Addition of pvals through 2-sample t-tests for difference in mean
    p.vals <- rep(NA, 2*length(cts_vars))
    for (var in cts_vars) {
      f <- paste0(var, " ~ opt")
      t.test <- t.test(as.formula(f), dat_opt)
      p.vals[2*which(cts_vars == var) - 1] <- t.test$p.value
    }
    opt_table_cts_t <- as.data.frame(opt_table_cts_t)
    opt_table_cts_t$p.val <- p.vals
  
  
  # 3. Summarize characterizing variable proportions (for categorical variables) by trt group assigned by opt decision rule
  
    # Collection of categorical variables
    cat_vars <- c("Gender", "SHNumLast12Months_binary", "DKANumLast12Months_binary", "PumpUse_Binary", 
                  "EducationLevel_Tri", "Ins_Combined", "CPep_detected", "Race_binary", "Ethnicity_binary")
    cat_vars_split <- c("Gender", "SHNumLast12Months_binary", "DKANumLast12Months_binary", "PumpUse_Binary", 
                  "LessBach", "Bach", "GreaterBach", "MedicareOther", "Private", "PrivateMedicare", 
                  "CPep_detected", "Race_binary", "Ethnicity_binary")
    dat_opt_cat <- dat_opt %>%
      select(rowname, opt, count, all_of(cat_vars))

    # Summary table with proportions and patient count (n) for selected variables
    opt_table_cat <- dat_opt_cat %>%
      group_by(opt) %>%
      summarise(n = n(),
                propMale = sum(count[Gender == "M"]) / n,
                nMale = sum(count[Gender == "M"]),
                propSHNumLast12Months_binary = sum(count[SHNumLast12Months_binary == 1]) / n,
                nSHN = sum(SHNumLast12Months_binary),
                propDKANumLast12Months_binary = sum(count[DKANumLast12Months_binary == 1]) / n,
                nDKA = sum(DKANumLast12Months_binary),
                propPumpUse_Binary = sum(count[PumpUse_Binary == 1]) / n,
                nPumpUse = sum(PumpUse_Binary),
                propLessBach = sum(count[EducationLevel_Tri == "<Bach"]) / n,
                nLessBach = sum(count[EducationLevel_Tri == "<Bach"]),
                propBach = sum(count[EducationLevel_Tri == "Bach"]) / n,
                nBach = sum(count[EducationLevel_Tri == "Bach"]),
                propGreaterBach = sum(count[EducationLevel_Tri == ">Bach"]) / n,
                nGreaterBach = sum(count[EducationLevel_Tri == ">Bach"]),
                propMedicareOther = sum(count[Ins_Combined == "Medicare/other"]) / n,
                nMedicareOther = sum(count[Ins_Combined == "Medicare/other"]),
                propPrivate = sum(count[Ins_Combined == "Private"]) / n,
                nPrivate = sum(count[Ins_Combined == "Private"]),
                propPrivateMedicare = sum(count[Ins_Combined == "Private and Medicare"]) / n,
                nPrivateMedicare = sum(count[Ins_Combined == "Private and Medicare"]),
                propCPepDetected = sum(count[CPep_detected == "yes"]) / n,
                nCPepDetected = sum(count[CPep_detected == "yes"]),
                propWhite = sum(count[Race_binary == "White"], na.rm = TRUE) / n, # will be overridden by mutate()
                nWhite = sum(count[Race_binary == "White"], na.rm = TRUE),
                propNotHispanic = sum(count[Ethnicity_binary == "Not Hispanic or Latino"], na.rm = TRUE) / n, # will be overridden by mutate()
                nNotHispanic = sum(count[Ethnicity_binary == "Not Hispanic or Latino"], na.rm = TRUE)) %>%
      mutate(propWhite = ifelse(opt == 2, nWhite/(n-2), nWhite/n),
             propNotHispanic = ifelse(opt == 2, nNotHispanic / (n-1), nNotHispanic/n))
  
    # Transposing and addition of column headers
    opt_table_cat_t <- t(opt_table_cat)[-c(1,2), ]
    opt_table_cat_t <- round(opt_table_cat_t, 3)
    opt_table_cat_t <- cbind(opt_table_cat_t[, 2], opt_table_cat_t[, 1])
    colnames(opt_table_cat_t) <- c("CGM (n = 173)", "BGM (n = 21)")

    # Addition of pvals through 2-proportion t-tests
    p.vals <- rep(NA, 2*length(cat_vars_split))
    for (var in cat_vars_split) {
      prop.test <- prop.test(x = as.numeric(opt_table_cat_t[2*which(cat_vars_split == var), ]), n = c(173, 21))
      p.vals[2*which(cat_vars_split == var) - 1] <- prop.test$p.value
    }
    opt_table_cat_t <- as.data.frame(opt_table_cat_t)
    opt_table_cat_t$p.val <- p.vals
  
### Supplementary Table S1 ###

  decision_list_CV <- readRDS("./data/decision_list.rds")
  policy_tree_CV <- readRDS("./data/policy_tree_CV.rds")
  all_CV <- rbind(policy_tree_CV, decision_list_CV)

