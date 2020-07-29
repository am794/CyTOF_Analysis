#############
# MFI plots #
#############

library("openxlsx")
library("reshape2")
library("ggplot2")
MFI <- read.xlsx("/Users/amandava/Desktop/Pertussis/New_Donors/Gated/MFI_pertussis_an3.xlsx")

##Boxplots of MFI
for(i in 1:56){
  if(i>4){
    quartz()
    print(ggplot(data=MFI,aes(x=factor(visit),y=MFI[,i],fill=factor(cohort)))+scale_fill_brewer(palette="Blues")+
            geom_boxplot(aes(fill=cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(cohort)),shape=21,position=position_jitterdodge(),size=3)+
            stat_summary(fun.y=median, geom="smooth", aes(group=factor(cohort),colour=factor(cohort)),position=position_dodge(0.2),lwd=0.5)+
            xlab("Visit")+ ylab("MFI") + ggtitle(paste0("Boxplot of MFI for ", colnames(MFI)[i])))
  }
}
#GBP2_CD4Temra <- ddply(MFI,~visit+cohort,summarize,median=median(GBP2_CD4Temra))
#assign(paste("median_",colnames(MFI)[7],sep=""),ddply(MFI,~visit+cohort,summarize,median=median(IgM_Bcells)))


###Fold change
median <- ddply(MFI,.(visit,cohort),colwise(median,na.rm=TRUE))
data <- c()
fc_median <- c()
order_df <- df[order(median[,2]),]
for(i in 1:10){
  if(i<6){
    fc_median <- order_df[i,5:56]/order_df[1,5:56]
    }else{
      fc_median <- order_df[i,5:56]/order_df[6,5:56]
    }
    data <- rbind(data,fc_median)
  }
fc_median_data <- cbind(order_df[,1:4],data)

for(i in 1:56){
  if(i>4){
    quartz()
    print(ggplot(data=fc_median_data,aes(x=factor(visit),y=fc_median_data[,i],fill=factor(cohort)))+scale_fill_brewer(palette="Blues")+
            geom_line(aes(group=factor(cohort),colour=factor(cohort)),lwd=1)+
            geom_point(aes(colour=factor(cohort)),shape=21,size=3)+
            xlab("Visit")+ ylab("Foldchange(Median) ") + ggtitle(paste0("Fold change of median of MFI for ", colnames(MFI)[i] ," with Visit 1 as baseline")))
  }
}

write.xlsx(median[,-3],file="/Users/amandava/Desktop/Pertussis/New_Donors/MFI_Foldchange/median.xlsx",row.names = TRUE,col.names = TRUE)
write.xlsx(fc_median_data[,-3],file="/Users/amandava/Desktop/Pertussis/New_Donors/MFI_Foldchange/Foldchange_median.xlsx",row.names = FALSE,col.names = TRUE)


##############Generating MFI tables#############
setwd("/Users/amandava/Desktop/Pertussis/New_Donors/Gated/")
files <- list.files(pattern = "Transformed")
len <- length(list.files(pattern="Transformed"))
preprocessed <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2968/Preprocessed_TXT/Transformed_2968-WC-1_02_normalized.fcs.txt",header=T)
config <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2968/config/2968.config")
for(i in 1:len){
  mfi <- read.table(paste0(files[i],"/DUnSup_pop_MFI.txt"))[,c(-1)]
  colnames(mfi) <- colnames(preprocessed)
  rownames(mfi) <- config[,12]
}

system(paste0("cut -f 20,33,5,39 /Users/amandava/Desktop/FLOCK_sampling/Transformed_2691_WC-1_01_normalized.fcs.txt > /Users/amandava/Desktop/FLOCK_sampling/b.txt"))

system(paste0("for i in */; do awk '{print $18}' ORS='\t' $i/DUnSup_pop_MFI.txt | awk '{print $4 "\t" $5 "\t" $6 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $21 "\t" $22 "\t" $23 "\t" $24 "\t" $25 "\t" $26 "\t" $27 "\t" $28 "\t" $29 "\t" $30 "\t" $31 "\t" $37 "\t" $38 "\t" $40 "\t" $42 "\t" $43}' >> a1.txt; done; "))
