library(readxl)
library(ggplot2)

experimentInfo <- read.csv("C:/Users/guypw/OneDrive/Documents/PhD/UpgradeReport/data/BrainBankClinicalData.csv")

experimentInfo <- as.data.frame(experimentInfo)

experimentInfo[experimentInfo$Status == "Ctrl", "Status"] <- "Control"

experimentInfo$caseControlGender <- 1

# Gender Split Between Case and Controls
ggplot(data=experimentInfo, aes(x=Status, y=caseControlGender, fill=Sex)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("") + ylab("Count of Patients") +
  guides(fill=guide_legend(title="Sex"))

# Age Distribution between Case and Controls
ggplot(data=experimentInfo, aes(x=Status, y=Age, fill=Sex)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("Age at Death (Years)") +
  guides(fill=guide_legend(title="Sex"))

ggplot(data=experimentInfo[!is.na(experimentInfo$Form) & !experimentInfo$Form == "",], aes(x=Form, y=caseControlGender, fill=Sex)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("") + ylab("Count of Patients") +
  guides(fill=guide_legend(title="Sex"))
