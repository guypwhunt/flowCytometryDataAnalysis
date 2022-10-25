library(readxl)
library(ggplot2)

experimentInfo <- read.csv("C:/Users/guypw/OneDrive/Documents/PhD/UpgradeReport/data/AnswerAlsClinicalData.csv")

experimentInfo <- as.data.frame(experimentInfo)

experimentInfo[experimentInfo$Subject.Group != "ALS", "Subject.Group"] <- "Control"
experimentInfo[experimentInfo$Subject.Group == "ALS", "Subject.Group"] <- "Case"

experimentInfo <- experimentInfo[!is.na(experimentInfo$Age.At.Death) & experimentInfo$Subject.Group == "Case" | experimentInfo$Subject.Group == "Control",]

experimentInfo <- experimentInfo[experimentInfo$Sex != "",]

experimentInfo$caseControlGender <- 1

# Gender Split Between Case and Controls
ggplot(data=experimentInfo, aes(x=Subject.Group, y=caseControlGender, fill=Sex)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("") + ylab("Count of Patients") +
  guides(fill=guide_legend(title="Sex"))

# Age Distribution between Case and Controls
ggplot(data=experimentInfo, aes(x=Subject.Group, y=Age.At.Death, fill=Sex)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("Age at Death (Years)") +
  guides(fill=guide_legend(title="Sex"))
