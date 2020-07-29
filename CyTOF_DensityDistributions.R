############################
# CyTOF data distributions #
############################
source("/Users/amandava/Desktop/amandava/HIPC_IOF/transformation_functions.R")
source("/Users/amandava/Desktop/amandava/FCSTrans_FCS3.1_R3.0.r")
require("flowCore")
library("gdata")
library(ggplot2)
library(colorspace)
library(gridExtra)
install.packages("RColorBrewer")
library("RColorBrewer")
library("MASS")

setwd("/Users/amandava/Desktop/Pertussis/FCS/")
file_list <- list.files(path = "/Users/amandava/Desktop/Pertussis/FCS/")
fcs_raw <- read.flowSet(file_list,transformation = FALSE, truncate_max_range = FALSE)

for(i in 1:50){
  fcs_files <- read.FCS(file_list[[i]])
  cols_keep <- c(1:4,10:42,44,46,47,49,51)
  channels_1 <- c(3,4,10:42,44,46,47,49,51)
  channels_2 <- c(1,2)
  cofactor <- 2
  expr <- exprs(fcs_raw[[i]])
  markers <- gsub(pattern = ".*_", replacement = "", x = as.vector(fcs_raw[[1]]@parameters@data$desc))
  #colnames(expr)[which(!is.na(markers))] <- markers[which(!is.na(markers))]
  colnames(expr)[channels_1] <- markers[channels_1]
  #pregating_channels <- c("Event_length","DNA-1","DNA-2","Viability")
  pregating_channels <- c(2,46,47,49)
  lineage_channels <- c(3,4,10:42,44,51)
  instrument_channels <- c("Time","Event_length","Center","Offset","Width","Residual")
  exprs_trans <- cbind(asinh(expr[,c(pregating_channels,lineage_channels)]/cofactor),expr[,instrument_channels])
  a=4096
  normalize=function(data){
    return(
      a*(data-min(data))/(max(data)-min(data))
    )
}
  #trans_data <- as.data.frame(cbind(normalize(exprs_trans[,1:41]),exprs_trans[,42:47]))
  trans_data <- as.data.frame(cbind(normalize(exprs_trans[,1:41])))
  write.table(trans_data,file = paste0("Transformed_",file_list[i],".txt"),quote=FALSE,sep = "\t",col.names=TRUE,row.names=FALSE)
}


k <- 11; my.cols <- rev(brewer.pal(k, "RdYlBu"))
quartz()
plot(trans_data[,c(2,3)], pch=".", col="grey", main="DNA1 vs DNA2",xlim=c(0,4000),ylim=c(0,4000))
z <- kde2d(trans_data[,2], trans_data[,3], n=50)
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
quartz()
plot(exprs_trans[,c(21,38)], pch=".", col="grey", main="CD45 vs Viability")
z <- kde2d(exprs_trans[,21], exprs_trans[,38], n=50)
contour(z, drawlabels=FALSE, nlevels=k, col=my.cols, add=TRUE)
