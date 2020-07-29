#############Tregs###########
# Parent is CD4+ T cells
#############################

files <- list.files("/Users/amandava/Desktop/Pertussis/New_Donors/Gated/Gated/",pattern = ".fcs")
i = 100
cd4_for_dafi2 <-  read.table(paste0("/Users/amandava/Desktop/Pertussis/New_Donors/Gated/Gated/",files[i],"/flock_results_all.txt"),header=TRUE)
config_table <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2968/config/2968.config",sep="\t")
parent_for_cd4 <- paste0("pop",config_for_dafi2[21,8])  
cd4_for_tregs <- cd4_for_dafi2[which(cd4_for_dafi2[,parent_for_cd4]=="0"),1:40]
tab <- cd4_for_tregs

##Feature selection step
  n=31
  x <- config_table[n,2]
  y <- config_table[n,3]
  min_x <- (config_table[n,4]*4096/200)
  max_x <- (config_table[n,5]*4096/200)
  min_y <- (config_table[n,6]*4096/200)
  max_y <- (config_table[n,7]*4096/200)
  #tab <- as.data.frame(exprs(popFlowFrameList[[config_table[n,8]]]))
  
  #Wilk's Generalized variance
  generalised_variance <- matrix(nrow=40,ncol=40)
  tot_var <- matrix(nrow=40,ncol=40)
  for(i in 1:40){
    for(j in 1:40){
      if(i != j){
        generalised_variance[i,j] <- (var(tab[,i])*var(tab[,j]))-(cov(tab[,i],tab[,j]))^2 ##det(cov(tab[,c(i,j)])) or eigen(cov(tab[,c(i,j)]))
      }else{
        generalised_variance[i,j] <- 0
      }
    }
  }
  
  #Pick the top 100 feature combinations with highest dispersion(around 13%)
  #p=round(((ncol(tab)-1)*(ncol(tab)-2))*15/(2*100)) #Number of marker combinations to consider everytime
  gv <- as.data.frame(generalised_variance)
  row.names(gv) <- colnames(gv) <- colnames(tab)[1:40]
  gv[lower.tri(gv,diag=T)] <- 0
  m_gv <- melt(t(gv))
  m_gv_order <- m_gv[order(abs(m_gv$value),decreasing=TRUE),][1:100,]
  
  #FLOCK based density biased downsampling
  write.table(tab[,c(x,y)],file="/Users/amandava/Desktop/FLOCK_sampling/filtered_file.txt",quote=FALSE,sep="\t",row.names = F)
  system(paste0("mv /Users/amandava/Desktop/FLOCK_sampling/filtered_file.txt /Users/amandava/Desktop/FLOCK_sampling/Input/"))
  system(paste0("/Users/amandava/Desktop/FLOCK_sampling/FLOCK_sampling /Users/amandava/Desktop/FLOCK_sampling/Input 50"))
  flock_samp <- read.table("/Users/amandava/Desktop/FLOCK_sampling/Input_sampled/filtered_file_sampled.txt",header=TRUE)
  flock_samp_mapping <- read.table("/Users/amandava/Desktop/FLOCK_sampling/Input_sampled/filtered_file_mapping",header=TRUE)
  row.names(tab) <- seq(1,dim(tab)[1],by=1)
  events_map <- which(flock_samp_mapping$CompressedEventID != 0)
  mapped_events <- flock_samp_mapping[events_map,]
  final_mapped <- mapped_events[order(mapped_events$CompressedEventID),]
  
  #FR test for number of runs and p value
  #FR test for number of runs and p value
  FR_runs <- c()
  FR_pvals <- c()
  nn_edges <- c()
  for(i in 1:100){
    dim_1 <- match(m_gv_order[i,1],names(tab))
    dim_2 <- match(m_gv_order[i,2],names(tab))
    features_data <- cbind(flock_samp,tab[final_mapped$OriginalEventID,dim_1],tab[final_mapped$OriginalEventID,dim_2])
    colnames(features_data)[3:4] <- c(paste0(colnames(tab)[dim_1]),paste0(colnames(tab)[dim_2]))
    features_data_x <- match(colnames(tab)[x],names(features_data))
    features_data_y <- match(colnames(tab)[y],names(features_data))
    features_data$c_x <- features_data$c_y <- features_data$class <- 0
    features_data$c_x[which(features_data[,features_data_x] < min_x | features_data[,features_data_x] > max_x)] <- 1
    features_data$c_x[which(features_data[,features_data_x] >= min_x & features_data[,features_data_x] <= max_x)] <- 2  
    features_data$c_y[which(features_data[,features_data_y] < min_y | features_data[,features_data_y] > max_y)] <- 1
    features_data$c_y[which(features_data[,features_data_y] >= min_y & features_data[,features_data_y] <= max_y)] <- 2 
    features_data$class[which(features_data$c_x == "2" & features_data$c_y == "2")] <- 2
    features_data$class[which(features_data$c_x == "1" | features_data$c_y == "1")] <- 1
    col1 <- match(m_gv_order[i,1],names(features_data))
    col2 <- match(m_gv_order[i,2],names(features_data))
    a1 <- features_data[which(features_data$class == "1"),c(col1,col2)]
    a2 <- features_data[which(features_data$class == "2"),c(col1,col2)]
    FR_runs[i] <- getFR(a1,a2)$runs
    FR_pvals[i] <- getFR(a1,a2)$pNorm
    nn_edges[i] <- getNN(a1,a2)
  }
  ordered_combs <- m_gv_order[order(FR_runs),1:2]
  list <- m_gv_order[order(FR_runs),1:2][1:6,]
  dims <- union(union(match(list[,1],names(tab)),c(x,y)),match(list[,2],names(tab)))
  #return(ordered_combs)
  #return(dims)
}

for (i in 1:length(files)) {
  gated <- read.table(paste0("/Users/amandava/Desktop/Pertussis/New_Donors/Gated/Gated/",files[i],"/flock_results_all.txt"),header=TRUE)  
  for_tregs <- gated[which(cd4_for_dafi2[,parent_for_cd4]=="0"),dims]
  write.table(for_tregs,paste0("/Users/amandava/Desktop/Pertussis/New_Donors/Gated_tregs/tregs_",files[i],".txt"),sep="\t",quote=FALSE,row.names = FALSE)
  }
