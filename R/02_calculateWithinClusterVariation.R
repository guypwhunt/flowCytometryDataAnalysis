try(source("R/01_functions.R"))
loadlibraries()

directoryName <- "bCells"
columnNames <- c("IgD...PerCP.Cy5.5.A", "CD24...BV605.A", "CD27...BV650.A")
clusterName <- "meta_clusters_flowsom"
numberOfClusters <- seq(3,4)

for (number in numberOfClusters) {
  calculateVarianceWithinClusters(directoryName, columnNames, clusterName, number)
}

listOfVariances <- c()

for (number in numberOfClusters) {
  x <- readRDS(paste0("data/", directoryName, "/clusteringOutput/", clusterName, number,"flowSomDf.rds"))
  x <- x^2
  x <- sum(x)
  x <- list(x)
  names(x) <- number
  listOfVariances <- append(listOfVariances, x)
}



#---
#Date: 24/06/2022
#Author: Thomas P Spargo <thomas.spargo@kcl.ac.uk>
#
#Purpose: Compute elbow plots for clustering based on within-cluster variance for clustering solutions of 1:k clusters.
#---
#Instructions

#This operates based on trailing command line arguments:
# Format in bash script will be
# Rscript ./path/to/elbow_simplified.R \
#         "./path/to/dataset/with/clustering/variables/file.csv" \                  #this string will become args[1]
#         "./path/to/file/detailing/assignments/for/each/clustering/solution.csv"   #this string will become args[2]

#The first column of both files should be the row numbers, which will be used as IDs to match the data to clustering positions
#The first row should be a header in both files

#args[1] dataset column 2:ncols() each should be a variable used in the clustering (it should not include any cluster asssingments)
#args[2] dataset column 2:ncols() each should be cluster assignments between 1:k classes

#The elbow plot figure will be saved into the working directory, the x axis labels will be defined by 2:ncols() column names for args[2]

#---
#Begin script
library(ggplot2)
library(dplyr)

#Extract command line arguments for datasets to use
args <- commandArgs(trailingOnly = TRUE)

#Read the CSV files
data <- read.csv(args[1],header=TRUE)
clust <- read.csv(args[2],header=TRUE)


#Define function
#clust is a vector with cluster assignments for the k_i-cluster solution
#data is the dataset used in clustering
#clust_ID is a vector containing row IDs associated with the clustering vector
calcVar <- function(clust,data,clust_ID){

  clust_wIDs<- data.frame(clust_ID,clust) #Recombine the clust_ID and clust vectors into a dataframe
  colnames(clust_wIDs) <- c("ID","clust") #Rename to appropriate names

  colnames(data)[1] <- "ID"               #Rename first column of main dataset to be "ID", since this will be matched within left-join


  #Join the cluster assignments to the main dataset by the ID index
  df <- data %>%
    left_join(clust_wIDs, by="ID")

  #Compute total variance across all clusters and predictor variables
  totvar <- df %>%
    select(-ID) %>%                                           #Drop the ID column
    group_by(clust) %>%                                       #group by clusters
    dplyr::summarise(across(everything(), ~ var(.))) %>%      #Identify variance for each variable by cluster
    ungroup() %>%                                             #Remove grouping by cluster
    select(-clust) %>%                                        #drop cluster identifying column
    sum()                                                     #sum all variances

  return(totvar)                                              #Return total within-cluster variances for this clustering solution

}

#Apply calcVar across columns 2:ncol() of clust, return total within-cluster variances
clustVar <- sapply(clust[-1],calcVar,data=data,clust_ID=clust[1])

#Generate factor variable indicating the x-axis labels from column names from clustering data
#Factor ensures correct ordering
xlevs <- factor(colnames(clust[-1]),levels=colnames(clust[-1]))

#Plot the elbow plot across variances from clustVar
elbow_plot <- ggplot(data.frame(clustVar),aes(x=xlevs, y=clustVar,group=1))+
  geom_point()+
  geom_line()+
  xlab("Cluster fit")+
  ylab("Sum of within-cluster variance")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5, hjust=1)) #Orient x axis label to allow for longer names, hjust=1 is right aligned, vjust is center aligned


#Generate file name for saving elbow plot, based on input dataset name
dataname <- args[1] %>%
  gsub(".*/","",.) %>% #sub out the directory tree
  sub("\\..*?$","",.)  #sub out the appending filetype, by pattern matching the remaining full stop (this could cut more than just file extention if . is used in name)

#Save the figure into the working directory
ggsave(plot=elbow_plot,
       filename = paste0("Elbow_plot",dataname,"_",Sys.Date(),".pdf"),
       units="mm",width=200,height=150,
       device=cairo_pdf)







