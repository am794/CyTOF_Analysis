###########Pertussis median intensity values#############
setwd("/Users/amandava/Desktop/Pertussis/New_Donors/Gated/")
files_list <- list.files(path="/Users/amandava/Desktop/Pertussis/New_Donors/Gated/",pattern = "Transformed")
#function(file){
table <- c()
for(file in files_list){
  flock_table <- read.table(paste0("/Users/amandava/Desktop/Pertussis/New_Donors/Gated/",file,"/flock_results_all.txt"),header=T,sep="\t")
  #colnames(flock_table)[41:83] <- c("")
  ls <- c(c())
  for (i in colnames(flock_table)[41:83]) {
    rows <- which(flock_table[,i] == "0")
    for(j in names(flock_table)[1:40]){
      #median_func <- function(x){median(x)}
      #j=names(flock_table)[20]
      assign(paste0(i,"_",j),median(flock_table[rows,j]))
      ls[[i]] <- as.data.frame(cbind(ls[[i]],get(noquote(paste0(i,"_",j))))) 
      #colnames(ls[[i]])[match(j,names(flock_table)[1:40])]=paste0(i,"_",j)
      colnames(ls[[i]])[match(j,names(flock_table)[1:40])]=paste0(j)
    }
  }
  table <- cbind(table,unlist(ls)) 
}
colnames(table) <- files_list
write.xlsx(table,"/Users/amandava/Desktop/Pertussis/New_Donors/Gated/median_Intensity.xlsx",sep="\t",quote=FALSE,row.names=TRUE)

