library(ggplot2)
DiabScreening <- read.table("data/study data/DiabScreening.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
DiabPhysExam <- read.table("data/study data/DiabPhysExam.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")

dat_clean_XAY_full <- readRDS("./data/dat_clean_XAY_full.rds")
dat_with_opt_trt <- readRDS("./data/dat_with_opt_trt_194.rds")

ggplot(aes(x = gluCV_num, y = 100*gluBelow70Chg, color = TrtGroup), data = dat_clean_XAY_full) + 
  geom_point() + 
  #geom_smooth(method = "lm") +
  #geom_smooth() +
  stat_smooth(method = "lm", formula = y ~ poly(x,2), se = FALSE) +
  #stat_smooth(method = "gam", formula = y ~ s(x, k = 3), se = FALSE) +
  xlab("Coefficient of Variation") +
  ylab("Decrease in % Time Under 70 mg/dL") +
  ggtitle("CGM vs. BGM Treatment Effect is Moderated by Baseline CV%")

summary(lm(gluBelow70Chg ~ gluCV_num, data = dat_clean_XAY_full))
summary(lm(gluBelow70Chg ~ poly(gluCV_num, 2), data = dat_clean_XAY_full))

dat_with_opt_trt <- rownames_to_column(dat_with_opt_trt) %>%
  mutate(rowname = as.numeric(rowname)) %>%
  left_join(select(DiabScreening, PtID, Ethnicity, Race, UnitsInsTotal), by = c("rowname" = "PtID")) %>%
  left_join(select(DiabPhysExam %>% filter(Visit == "Screening"), PtID, Weight, WeightUnits, Height, HeightUnits), by = c("rowname" = "PtID")) %>%
  mutate(Weight_Kg = ifelse(WeightUnits == "kg", Weight, Weight / 2.20462)) %>%
  mutate(Height_M = ifelse(HeightUnits == "in", Height*0.0254, Height/100)) %>%
  mutate(BMI = Weight_Kg / (Height_M)^2) %>%
  mutate(Race_binary = ifelse(Race == "White", "White", ifelse(Race == "Unknown/not reported", NA, "Non-White"))) %>%
  mutate(Ethnicity_binary = ifelse(Ethnicity == "Not Hispanic or Latino", "Not Hispanic or Latino", ifelse(Ethnicity == "Hispanic or Latino", "Hispanic or Latino", NA))) %>%
  mutate(InsulinDosesKg = UnitsInsTotal / Weight_Kg) %>%
  mutate(count = 1) %>%
  select(!c(Ethnicity, Race, Weight, WeightUnits, Height, HeightUnits, Weight_Kg, Height_M))

opt_table <- dat_with_opt_trt %>%
  group_by(opt) %>%
  summarize(n = n(),
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
            meanSHNumLast12Months_binary = mean(SHNumLast12Months_binary),
            nSHN = sum(SHNumLast12Months_binary),
            meanDKANumLast12Months_binary = mean(DKANumLast12Months_binary),
            nDKA = sum(DKANumLast12Months_binary),
            meanDiabDuration = mean(DiabDuration),
            sdDiab = sd(DiabDuration),
            meanPumpUse_Binary = mean(PumpUse_Binary),
            nPumpUse = sum(PumpUse_Binary),
            #InsulinDosesKg = mean(InsulinDosesKg),
            meangluCV_num = mean(gluCV_num),
            sdgluCV = sd(gluCV_num),
            propMale = sum(count[Gender == "M"]) / n,
            nMale = sum(count[Gender == "M"]),
            propCPepDetected = sum(count[CPep_detected == "yes"]) / n,
            nCPepDetected = sum(count[CPep_detected == "yes"]),
            propMedicareOther = sum(count[Ins_Combined == "Medicare/other"]) / n,
            nMedicareOther = sum(count[Ins_Combined == "Medicare/other"]),
            propPrivate = sum(count[Ins_Combined == "Private"]) / n,
            nPrivate = sum(count[Ins_Combined == "Private"]),
            propPrivateMedicare = sum(count[Ins_Combined == "Private and Medicare"]) / n,
            nPrivateMedicare = sum(count[Ins_Combined == "Private and Medicare"]),
            propLessBach = sum(count[EducationLevel_Tri == "<Bach"]) / n,
            nLessBach = sum(count[EducationLevel_Tri == "<Bach"]),
            propBach = sum(count[EducationLevel_Tri == "Bach"]) / n,
            nBach = sum(count[EducationLevel_Tri == "Bach"]),
            propGreaterBach = sum(count[EducationLevel_Tri == ">Bach"]) / n,
            nGreaterBach = sum(count[EducationLevel_Tri == ">Bach"]),
            meanBMI = mean(BMI),
            sdBMI = sd(BMI),
            propWhite = sum(count[Race_binary == "White"], na.rm = TRUE) / n,
            nWhite = sum(count[Race_binary == "White"], na.rm = TRUE),
            propNotHispanic = sum(count[Ethnicity_binary == "Not Hispanic or Latino"], na.rm = TRUE) / n,
            nNotHispanic = sum(count[Ethnicity_binary == "Not Hispanic or Latino"], na.rm = TRUE),
            propNotHispanicAndWhite = sum(count[Ethnicity_binary == "Not Hispanic or Latino" & Race_binary == "White"], na.rm = TRUE) / n,
            nNotHispanicAndWhite = sum(count[Ethnicity_binary == "Not Hispanic or Latino" & Race_binary == "White"], na.rm = TRUE),
            meanInsulinDosesKg = mean(InsulinDosesKg, na.rm = TRUE),
            sdInsulinDosesKg = sd(InsulinDosesKg, na.rm = TRUE)
            )

opt_table_t <- t(opt_table)[-c(1,2), ]
opt_table_t <- round(opt_table_t, 3)
#colnames(opt_table_t) <- c("BGM (n = 19)", "CGM (n = 165)")
colnames(opt_table_t) <- c("BGM (n = 21)", "CGM (n = 173)")
opt_table_t <- cbind(opt_table_t[, 2], opt_table_t[, 1])
colnames(opt_table_t) <- c("CGM", "BGM")
#write.table(opt_table_t, file = "./data/opt_table_t.txt", sep = ",", quote = FALSE, row.names = T)


# T-Tests

t.test(gluCV_num ~ opt, dat_with_opt_trt)
t.test(AgeAsOfEnrollDt ~ opt, dat_with_opt_trt)
t.test(DiabDuration ~ opt, dat_with_opt_trt)
t.test(DiagAge ~ opt, dat_with_opt_trt)
t.test(HbA1cTestRes ~ opt, dat_with_opt_trt)
t.test(InsulinDosesKg ~ opt, dat_with_opt_trt)
t.test(BMI ~ opt, dat_with_opt_trt)
t.test(gluBelow70 ~ opt, dat_with_opt_trt)
t.test(gluInRange ~ opt, dat_with_opt_trt)

prop.test(x = c(83, 10), n = c(173, 21)) # Gender
prop.test(x = c(157, 21), n = c(173, 21)) # Non-Hispanic
prop.test(x = c(162, 21), n = c(173, 21)) # White
prop.test(x = c(70, 5), n = c(173, 21)) # <Bach
prop.test(x = c(55, 7), n = c(173, 21)) # Bach
prop.test(x = c(48, 9), n = c(173, 21)) # >Bach
prop.test(x = c(49, 3), n = c(173, 21)) # Private
prop.test(x = c(58, 9), n = c(173, 21)) # Private/Medicare
prop.test(x = c(66, 9), n = c(173, 21)) # Medicare/other
prop.test(x = c(93, 9), n = c(173, 21)) # Insulin Pump Use
prop.test(x = c(37, 9), n = c(173, 21)) # C-Pep
prop.test(x = c(27, 1), n = c(173, 21)) # SHN
prop.test(x = c(8, 0), n = c(173, 21)) # DKA


