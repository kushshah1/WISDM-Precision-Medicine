### readin.R
### Reads in all different data tables (file necessary because the readin process is slightly different for different tables, as seen below)
### working directory = "Diabetes/Repository"

library(tidyverse)
# library(lubridate) # has functions like year(), etc. if needed
select <- dplyr::select

### Patient Data ###

# 16 patients withdraw pre-randomization
# 9 patients withdraw post-randomization
  PtFinalStat <- read.table("data/PtFinalStat.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")

# 219 rows, 16 of them with RandDt=NA (those who withdrew pre-randomization)
  PtRoster <- read.table("data/PtRoster.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")

  ScreenCompVisit <- read.table("data/ScreenCompVisit.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
  Visits <- read.table("data/Visits.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16") 

# Baseline Characteristics
  DiabScreening <- read.table("data/DiabScreening.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
  DiabSocioEcon <- read.table("data/DiabSocioEcon.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16", quote="\"") # Quote default here treats both ' and " as character delimiters. Removed ' from that list, since there are words like Master's and Bachelor's
  DiabPhysExam <- read.table("data/DiabPhysExam.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")

### CGM Data ###
  DeviceCGM <- read.table("data/DeviceCGM.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
  DeviceCGM <- DeviceCGM %>% 
  mutate(date = as.POSIXct(DeviceDtTm, tz = "UTC"))
  # All CGM data, >10 million records
  cgmAnalysisRCT <- read.table("data/cgmAnalysis RCT.txt", header = TRUE, sep = "|")
  cgmAnalysisRCT <- cgmAnalysisRCT %>% 
    mutate(date = as.POSIXct(DeviceDtTm, format = "%d%b%Y:%H:%M:%S", tz = "UTC"))
  cgmAnalysisEXT <- read.table("data/cgmAnalysis Ext.txt", header = TRUE, sep = "|")
  cgmAnalysisEXT <- cgmAnalysisEXT %>% 
    mutate(date = as.POSIXct(DeviceDtTm, format = "%d%b%Y:%H:%M:%S", tz = "UTC"))
  
### Outcomes ###
  gluIndicesRCT <- read.table("data/gluIndices RCT.txt", header = TRUE, sep = "|") %>%
    mutate(gluBelow70 = parse_number(as.character(gluBelow70)))
  
### HbA1c ###
  DiabLocalHbA1c <- read.table("data/DiabLocalHbA1c.txt", header = TRUE, sep = "|", fileEncoding = "UTF-16")
  