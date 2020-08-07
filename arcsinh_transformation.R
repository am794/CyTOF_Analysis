#############################
# CyTOF data transformation #
#############################

source("/Users/amandava/Desktop/amandava/HIPC_IOF/transformation_functions.R")
source("/Users/amandava/Desktop/amandava/FCS_trans/FCSTrans_FCS3.1_R3.0.r")
require("flowCore")
library("gdata")
library(ggplot2)
library(colorspace)
library(gridExtra)
install.packages("RColorBrewer")
library("RColorBrewer")
library("MASS")

setwd("/Volumes/Samsung_T5/CYTOF_Stanford/FlowRepository_FR-FCM-ZYAJ_files_Holden_21HealthyCyTOF/")
file_list <- list.files(path = "/Volumes/Samsung_T5/CYTOF_Stanford/FlowRepository_FR-FCM-ZYAJ_files_Holden_21HealthyCyTOF/",pattern = ".fcs")
fcs_raw <- read.flowSet(file_list,transformation = FALSE, truncate_max_range = FALSE)

for(i in 1:21){
  fcs_files <- read.FCS(file_list[[i]])
  channels_1 <- c(2,9)
  channels_2 <- c(3,4,10:41,43,44)
  channels_3 <- c(1,5:8,42,45:49)  ##Remove
  cofactor <- 2
  expr <- exprs(fcs_raw[[i]])
  markers <- gsub(pattern = ".*_", replacement = "", x = as.vector(fcs_raw[[i]]@parameters@data$desc))
  colnames(expr)[channels_2] <- markers[channels_2]
  colnames(expr)[channels_1] <- c("Event_length","Bead")
  exprs_trans <- as.data.frame(cbind(asinh(expr[,channels_2]/cofactor),expr[,channels_1]))

  ##Linear transformation
  a=4096
  normalize <- function(data){
    return(data*4095/(max(data)-min(data)))
  }
  trans_data <- as.data.frame(apply(exprs_trans,2,normalize))
  write.table(trans_data,file = paste0("Transformed_",file_list[i],".txt"),quote=FALSE,sep = "\t",col.names=TRUE,row.names=FALSE)
}

  #Kernel density estimation and adding contours to the plot
  k <- 11; my.cols <- rev(brewer.pal(k, "RdYlBu"))
  #trans_data <- as.data.frame(cbind(apply(exprs_trans[,1:34],2,normalize),exprs_trans[,35:38]))
  trans_data <- as.data.frame(apply(exprs_trans,2,normalize))
  plot(trans_data[,14],trans_data[,24], pch=".", col="grey", main="Bead vs DNA1",xlim=c(0,4096),ylim=c(0,4096))
  z <- kde2d(trans_data[,14], trans_data[,24], n=50)
  contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
