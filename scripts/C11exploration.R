library(ggplot2)
DiabScreening <- read.table("data/study data/DiabScreening.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
DiabPhysExam <- read.table("data/study data/DiabPhysExam.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")

dat_clean_XAY_full <- readRDS("./data/dat_clean_XAY_full.rds")
dat_with_opt_trt <- readRDS("./data/dat_with_opt_trt.rds")

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
  left_join(select(DiabScreening, PtID, Ethnicity, Race), by = c("rowname" = "PtID")) %>%
  left_join(select(DiabPhysExam %>% filter(Visit == "Screening"), PtID, Weight, WeightUnits, Height, HeightUnits), by = c("rowname" = "PtID")) %>%
  mutate(Weight_Kg = ifelse(WeightUnits == "kg", Weight, Weight / 2.20462)) %>%
  mutate(Height_M = ifelse(HeightUnits == "in", Height*0.0254, Height/100)) %>%
  mutate(BMI = Weight_Kg / (Height_M)^2) %>%
  mutate(Race_binary = ifelse(Race == "White", "White", ifelse(Race == "Unknown/not reported", NA, "Non-White"))) %>%
  mutate(Ethnicity_binary = ifelse(Ethnicity == "Not Hispanic or Latino", "Not Hispanic or Latino", ifelse(Ethnicity == "Hispanic or Latino", "Hispanic or Latino", NA))) %>%
  mutate(count = 1) %>%
  select(!c(Ethnicity, Race, Weight, WeightUnits, Height, HeightUnits, Weight_Kg, Height_M))

opt_table <- dat_with_opt_trt %>%
  group_by(opt) %>%
  summarize(n = n(),
            AgeAsOfEnrollDt = mean(AgeAsOfEnrollDt),
            DiagAge = mean(DiagAge),
            HbA1cTestRes = mean(HbA1cTestRes),
            gluBelow70 = mean(gluBelow70),
            gluInRange = mean(gluInRange),
            SHNumLast12Months_binary = mean(SHNumLast12Months_binary),
            DKANumLast12Months_binary = mean(DKANumLast12Months_binary),
            DiabDuration = mean(DiabDuration),
            PumpUse_Binary = mean(PumpUse_Binary),
            InsulinDosesKg = mean(InsulinDosesKg),
            gluCV_num = mean(gluCV_num),
            propMale = sum(count[Gender == "M"]) / n,
            propCPepDetected = sum(count[CPep_detected == "yes"]) / n,
            propMedicareOther = sum(count[Ins_Combined == "Medicare/other"]) / n,
            propPrivate = sum(count[Ins_Combined == "Private"]) / n,
            propPrivateMedicare = sum(count[Ins_Combined == "Private and Medicare"]) / n,
            propLessBach = sum(count[EducationLevel_Tri == "<Bach"]) / n,
            propBach = sum(count[EducationLevel_Tri == "Bach"]) / n,
            propGreaterBach = sum(count[EducationLevel_Tri == ">Bach"]) / n,
            BMI = mean(BMI),
            propWhite = sum(count[Race_binary == "White"], na.rm = TRUE) / n,
            propNotHispanic = sum(count[Ethnicity_binary == "Not Hispanic or Latino"], na.rm = TRUE) / n
            )

opt_table_t <- t(opt_table)[-c(1,2), ]
opt_table_t <- round(opt_table_t, 2)
colnames(opt_table_t) <- c("BGM (n = 19)", "CGM (n = 165)")
write.table(opt_table_t, file = "./data/opt_table_t.txt", sep = ",", quote = FALSE, row.names = T)
