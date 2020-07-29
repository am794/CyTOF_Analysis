#########################
# Base population plots #
#########################

library("openxlsx")
library("reshape2")
library("ggplot2")
pertussis_events <- read.xlsx("/Users/amandava/Desktop/Pertussis/New_Donors/Events/Analysis_3/Events.xlsx",sheet=2)
#config <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2968/config/2968.config")
#parents <- c(0,1,2,3,4,4,4,4,4,9,9,9,7,7,7,7,15,15,15,15,15,15,21,21,21,21,22,22,22,22,21,4,32,32,34,34,36,36,38,38,37,37,16)
parents <- c()
events <- as.data.frame(t(pertussis_events[,3:45]))
data <- c()
for(i in 1:length(parents)){
  if(i>1){
    fun1 <- function(x){(x*100/events[parents[i],])}
    dat <-  fun1(events[i,])
    data <- rbind(data,dat)
  }
}
data_fin <- cbind(pertussis_events[,1:3],t(data))
write.xlsx(data_fin,file="/Users/amandava/Desktop/Pertussis/New_Donors/Proportions_Analysis2.xlsx",row.names = TRUE,col.names = TRUE)

write.xlsx(final_duplexes,file="/Users/amandava/Desktop/Pertussis/New_Donors/Events/Analysis_3/duplexes.xlsx",row.names = TRUE,col.names = TRUE)

#WRT live cells
events <- as.data.frame(t(pertussis_events[,c(8:10,44,46)]))
rownames(events) <- c("Live_Cells","CD3_CD14_duplexes","CD19_CD3_duplexes","mDC_CD3_duplexes","pDC_CD3_duplexes")
data <- c()
for(i in 1:5){
  if(i>1){
    fun1 <- function(x){(x*100/events[1,])}
    dat <-  fun1(events[i,])
    data <- rbind(data,dat)
  }
}
final_duplexes <- cbind(pertussis_events[,1:4],t(data))
quartz()
ggplot(data=final_duplexes,aes(x=factor(Visit),y=CD3_CD14_duplexes,fill=factor(Cohort)))+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_jitter(aes(colour=factor(Cohort)))+
  scale_fill_brewer(palette="Blues")+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of CD3CD14 duplexes with live cells as parent ")

quartz()
ggplot(data=final_duplexes,aes(x=factor(Visit),y=CD3_CD14_duplexes,fill=factor(Cohort)))+
  scale_fill_brewer(palette="Blues")+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Median of percentages of CD3CD14 duplexes with live cells as parent")

quartz()
ggplot(data=final_duplexes,aes(x=factor(Visit),y=CD19_CD3_duplexes,fill=factor(Cohort)))+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape=NA)+geom_jitter(aes(colour=factor(Cohort)))+
  scale_fill_brewer(palette="Blues")+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of CD19CD3 duplexes with live cells as parent ")

quartz()
ggplot(data=final_duplexes,aes(x=factor(Visit),y=CD19_CD3_duplexes,fill=factor(Cohort)))+
  scale_fill_brewer(palette="Blues")+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Median of percentages of CD19CD3 duplexes with live cells as parent")

quartz()
ggplot(data=final_duplexes,aes(x=factor(Visit),y=mDC_CD3_duplexes,fill=factor(Cohort)))+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape=NA)+geom_jitter(aes(colour=factor(Cohort)))+
  scale_fill_brewer(palette="Blues")+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of mDC_CD3 duplexes with live cells as parent ")

quartz()
ggplot(data=final_duplexes,aes(x=factor(Visit),y=mDC_CD3_duplexes,fill=factor(Cohort)))+
  scale_fill_brewer(palette="Blues")+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Median of percentages of mDCCD3 duplexes with live cells as parent")

quartz()
ggplot(data=final_duplexes,aes(x=factor(Visit),y=pDC_CD3_duplexes,fill=factor(Cohort)))+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape=NA)+geom_jitter(aes(colour=factor(Cohort)))+
  scale_fill_brewer(palette="Blues")+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of pDC_CD3 duplexes with live cells as parent ")

quartz()
ggplot(data=final_duplexes,aes(x=factor(Visit),y=pDC_CD3_duplexes,fill=factor(Cohort)))+
  scale_fill_brewer(palette="Blues")+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Median of percentages of pDC_CD3 duplexes with live cells as parent")

CD3CD14 <- ddply(final_duplexes,~Visit+Cohort,summarize,median=median(CD3_CD14_duplexes),mean=mean(CD3_CD14_duplexes))

###Fold change
func1 <- function(df){
  fc_median <- c()
  fc_mean <- c()
  order_df <- df[order(df[,2]),]
  for(i in 1:10){
    if(i<6){
      fc_median[i] <- order_df[i,3]/order_df[1,3]
      fc_mean[i] <- order_df[i,4]/order_df[1,4]
    }else{
      fc_median[i] <- order_df[i,3]/order_df[6,3]
      fc_mean[i] <- order_df[i,4]/order_df[6,4]
    }
  }
  final_fc <- cbind(order_df,fc_median,fc_mean)
  return(final_fc)
}
CD3CD14_fc <- func1(CD3CD14)
CD19CD3_fc <- func1(CD19CD3)
pDCCD3_fc <- func1(pDCCD3)
mDCCD3_fc <- func1(mDCCD3)



CD19CD3 <- ddply(final_duplexes,~Visit+Cohort,summarize,median=median(CD19_CD3_duplexes),mean=mean(CD19_CD3_duplexes))
mDCCD3 <- ddply(final_duplexes,~Visit+Cohort,summarize,median=median(mDC_CD3_duplexes),mean=mean(mDC_CD3_duplexes))
pDCCD3 <- ddply(final_duplexes,~Visit+Cohort,summarize,median=median(pDC_CD3_duplexes),mean=mean(pDC_CD3_duplexes))
summary_duplexes <- rbind(CD3CD14,CD19CD3,mDCCD3,pDCCD3)
population <- rep(c("CD3CD14","CD19CD3","mDCCD3","pDCCD3"),c(10,10,10,10))
fin <- cbind(population,summary_duplexes)

quartz()
ggplot(data=CD3CD14_fc,aes(x=factor(Visit),y=fc_median,fill=factor(Cohort)))+geom_line(aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ylab("Fold Change(Median)")+ ggtitle("Fold Change(Median) for CD3CD14 duplexes with Visit 1 as baseline")
+scale_color_manual(values=c("#CC6666", "#9999CC"))

quartz()
ggplot(data=CD19CD3_fc,aes(x=factor(Visit),y=fc_median,fill=factor(Cohort)))+geom_line(aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ylab("Fold Change(Median)")+ ggtitle("Fold Change(Median) for CD19CD3 duplexes with Visit 1 as baseline")

quartz()
ggplot(data=pDCCD3_fc,aes(x=factor(Visit),y=fc_median,fill=factor(Cohort)))+geom_line(aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ylab("Fold Change(Median)")+ ggtitle("Fold Change(Median) for pDC_CD3 duplexes with Visit 1 as baseline")

quartz()
ggplot(data=mDCCD3_fc,aes(x=factor(Visit),y=fc_median,fill=factor(Cohort)))+geom_line(aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ylab("Fold Change(Median)")+ ggtitle("Fold Change(Median) for mDC_CD3 duplexes with Visit 1 as baseline")

#Mean
quartz()
ggplot(data=CD3CD14_fc,aes(x=factor(Visit),y=fc_mean,fill=factor(Cohort)))+geom_line(aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ylab("Fold Change(Mean)")+ ggtitle("Fold Change(Mean) for CD3CD14 duplexes with Visit 1 as baseline")
+scale_color_manual(values=c("#CC6666", "#9999CC"))

quartz()
ggplot(data=CD19CD3_fc,aes(x=factor(Visit),y=fc_mean,fill=factor(Cohort)))+geom_line(aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ylab("Fold Change(Mean)")+ ggtitle("Fold Change(Mean) for CD19CD3 duplexes with Visit 1 as baseline")

quartz()
ggplot(data=pDCCD3_fc,aes(x=factor(Visit),y=fc_mean,fill=factor(Cohort)))+geom_line(aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ylab("Fold Change(Mean)")+ ggtitle("Fold Change(Mean) for pDC_CD3 duplexes with Visit 1 as baseline")

quartz()
ggplot(data=mDCCD3_fc,aes(x=factor(Visit),y=fc_mean,fill=factor(Cohort)))+geom_line(aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ylab("Fold Change(Mean)")+ ggtitle("Fold Change(Mean) for mDC_CD3 duplexes with Visit 1 as baseline")

##Line plots
quartz()
ggplot(data=final_duplexes,aes(x=factor(Visit),y=CD3_CD14_duplexes,fill=factor(Subject)))+geom_line()

##Mean
quartz()
ggplot(data=final_duplexes,aes(x=factor(Visit),y=CD3_CD14_duplexes,fill=factor(Cohort)))+
  scale_fill_brewer(palette="Blues")+
  stat_summary(fun.y=mean, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Mean of percentages of CD3CD14 duplexes with live cells as parent")

quartz()
ggplot(data=final_duplexes,aes(x=factor(Visit),y=CD19_CD3_duplexes,fill=factor(Cohort)))+
  scale_fill_brewer(palette="Blues")+
  stat_summary(fun.y=mean, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Mean of percentages of CD19CD3 duplexes with live cells as parent")

quartz()
ggplot(data=final_duplexes,aes(x=factor(Visit),y=mDC_CD3_duplexes,fill=factor(Cohort)))+
  scale_fill_brewer(palette="Blues")+
  stat_summary(fun.y=mean, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Mean of percentages of mDCCD3 duplexes with live cells as parent")

quartz()
ggplot(data=final_duplexes,aes(x=factor(Visit),y=pDC_CD3_duplexes,fill=factor(Cohort)))+
  scale_fill_brewer(palette="Blues")+
  stat_summary(fun.y=mean, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),lwd=1)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Mean of percentages of pDC_CD3 duplexes with live cells as parent")

########### All percentages
events <- as.data.frame(t(pertussis_events[,c(8:47)]))
data <- c()
for(i in 1:40){
  if(i>1){
    fun1 <- function(x){(x*100/events[1,])}
    dat <-  fun1(events[i,])
    data <- rbind(data,dat)
  }
}
final_percentages <- cbind(pertussis_events[,1:4],t(data))
write.xlsx(final_percentages,file="/Users/amandava/Desktop/Pertussis/New_Donors/Events/Analysis_3/All_percentages_Pertussis.xlsx",row.names = TRUE,col.names = TRUE)

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`CD3-CD19+`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of B cells with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`CD3+Tcells`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of T cells with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`CD4+Tcells`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of CD4+ T cells with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`CD8+Tcells`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of CD8+ T cells with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`CD45CD14+`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Monocytes with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`Tregs`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Tregs with live cells as parent ")



###CD4 Memory
quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`TcmCD4`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Tcm CD4+ with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`TemCD4`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Tem CD4+ with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`TemraCD4`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Temra CD4+ with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`NaiveCD4`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Naive CD4 with live cells as parent ")

#CD8 Memory
quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`TcmCD8`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Tcm CD8+ with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`TemCD8`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Tem CD8+ with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`TemraCD8`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Temra CD8+ with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`NaiveCD8`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Naive CD8 with live cells as parent ")


#Other populatins basophils, mDC, pDC, classical, nonclassical, intermediate monocytes
quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`Classical_Monocytes`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Classical Monocytes with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`Non-Classical_Monocytes`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Non-Classical Monocytes with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`Intermediate_Monocytes`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Intermediate Monocytes with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`Basophils`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of Basophils with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`CD123+CD1c-`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of all pDCs with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`CD123+CD1c+`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of all mDCs with live cells as parent ")

quartz()
ggplot(data=final_percentages,aes(x=factor(Visit),y=final_percentages$`Nkcells`,fill=factor(Cohort)))+scale_fill_brewer(palette="Blues")+
  geom_boxplot(aes(fill=Cohort),width=0.3,outlier.shape = NA)+geom_point(aes(colour=factor(Cohort)),shape=21,position=position_jitterdodge(),size=3)+
  stat_summary(fun.y=median, geom="smooth", aes(group=factor(Cohort),colour=factor(Cohort)),position=position_dodge(0.2),lwd=0.5)+
  xlab("Visit")+ ylab("Percentage") + ggtitle("Boxplot of percentages of NK cells with live cells as parent ")


mfi_all_marker_cp <- read.xlsx("/Users/amandava/Desktop/Pertussis/New_Donors/Gated/MFI_All_markers_CP.xlsx",sheet=1)
