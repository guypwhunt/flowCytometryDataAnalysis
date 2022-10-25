library(readxl)
library(ggplot2)

experimentInfo <- read_excel("data/metadata/clinicalData.xlsx")

experimentInfo <- as.data.frame(experimentInfo)

visitOneFlowDF <- experimentInfo[experimentInfo$visit == 1 &
                                   experimentInfo$experiment == "flowCytometry",]

visitOneFlowDF$caseControlGender <- 1


# Gender Split Between Case and Controls
ggplot(data=visitOneFlowDF, aes(x=caseControl, y=caseControlGender, fill=gender)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("") + ylab("Count of Patients") +
  guides(fill=guide_legend(title="Sex"))

# Age Distribution between Case and Controls
ggplot(data=visitOneFlowDF, aes(x=caseControl, y=ageAtVisit, fill=gender)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("Age at First Visit (Years)") +
  guides(fill=guide_legend(title="Sex"))

visitOneFlowDFCase <- visitOneFlowDF[visitOneFlowDF$caseControl == "Case",]

# Gender Split Between Fast and Slow Progressets
ggplot(data=visitOneFlowDFCase, aes(x=fastSlow, y=caseControlGender, fill=gender)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("") + ylab("Count of Patients") +
  guides(fill=guide_legend(title="Sex"))

# Age Distribution between Case and Controls
ggplot(data=visitOneFlowDFCase, aes(x=fastSlow, y=ageAtVisit, fill=gender)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("Age at First Visit (Years)") +
  guides(fill=guide_legend(title="Sex"))

ggplot(data=visitOneFlowDFCase, aes(x=fastSlow, y=ageAtOnset, fill=gender)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("Age at Disease Onset (Years)") +
  guides(fill=guide_legend(title="Sex"))

# alsfrsR
ggplot(data=visitOneFlowDFCase, aes(x=fastSlow, y=alsfrsR, fill=gender)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("ALSFRS-R Score") +
  guides(fill=guide_legend(title="Sex"))

# diseaseDurationInYears
ggplot(data=visitOneFlowDFCase, aes(x=fastSlow, y=diseaseDurationInYears, fill=gender)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("Disease Duration (Years)") +
  guides(fill=guide_legend(title="Sex"))

# diseaseDurationInYears
ggplot(data=visitOneFlowDFCase, aes(x=fastSlow, y=diagnosticDelayInYears, fill=gender)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("Diagnostic Delay (Years)") +
  guides(fill=guide_legend(title="Sex"))



# Gender Split Between Bulbar and Limb Onset
ggplot(data=visitOneFlowDFCase, aes(x=BulbarLimb, y=caseControlGender, fill=gender)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("") + ylab("Count of Patients") +
  guides(fill=guide_legend(title="Sex"))

# Age Distribution between Bulbar and Limb
ggplot(data=visitOneFlowDFCase, aes(x=BulbarLimb, y=ageAtVisit, fill=gender)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("Age at First Visit (Years)") +
  guides(fill=guide_legend(title="Sex"))

ggplot(data=visitOneFlowDFCase, aes(x=BulbarLimb, y=ageAtOnset, fill=gender)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("Age at Disease Onset (Years)") +
  guides(fill=guide_legend(title="Sex"))

# alsfrsR
ggplot(data=visitOneFlowDFCase, aes(x=BulbarLimb, y=alsfrsR, fill=gender)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("ALSFRS-R Score") +
  guides(fill=guide_legend(title="Sex"))

# diseaseDurationInYears
ggplot(data=visitOneFlowDFCase, aes(x=BulbarLimb, y=diseaseDurationInYears, fill=gender)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("Disease Duration (Years)") +
  guides(fill=guide_legend(title="Sex"))

# diseaseDurationInYears
ggplot(data=visitOneFlowDFCase, aes(x=BulbarLimb, y=diagnosticDelayInYears, fill=gender)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("") + ylab("Diagnostic Delay (Years)") +
  guides(fill=guide_legend(title="Sex"))

# Bulbar limb Split by Progressions
ggplot(data=visitOneFlowDFCase, aes(x=fastSlow, y=caseControlGender, fill=BulbarLimb)) +
  geom_bar(stat="identity") +
  theme_minimal() +
  xlab("") + ylab("Count of Patients") +
  guides(fill=guide_legend(title="Site of Onset"))

experimentInfoDf <- experimentInfo[experimentInfo$visit > 1 &
                                     experimentInfo$experiment == "flowCytometry",]

experimentInfoDf$visit <- as.factor(experimentInfoDf$visit)

ggplot(data=experimentInfoDf,
       aes(x=visit, y=timeFromVisit1InYears, fill=visit)) +
  geom_boxplot() +
  theme_minimal() +
  xlab("Visit") + ylab("Time From First Visit (Years)") +
  guides(fill=guide_legend(title="Visit"))

