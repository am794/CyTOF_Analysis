########################
# Pertussis analysis 2 #
########################

library("openxlsx")
library("reshape2")
library("ggplot2")
pertussis_events <- read.xlsx("/Users/amandava/Desktop/Pertussis/New_Donors/Events_Percentages_Analysis2.xlsx",rowNames = TRUE)
config <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2968/config/2968.config")
parents <- c(0,1,2,3,4,4,4,4,4,9,9,9,7,7,7,7,15,15,15,15,15,15,21,21,21,21,22,22,22,22,21,4,32,32,34,34,36,36,38,38,37,37,16)
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


#################################################
# Plotting cluster centroids for DNA-1 vs DNA-2 #
#################################################

centroids_1 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2968/Gated_TXT/Transformed_2968-WC-1_02_normalized.fcs/centroids1.txt",header=TRUE)
centroids_2 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2968/Gated_TXT/Transformed_2968-WC-2_02_normalized.fcs/centroids1.txt",header=TRUE)
centroids_3 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2968/Gated_TXT/Transformed_2968-WC-3_01_normalized.fcs/centroids1.txt",header=TRUE)
centroids_4 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2968/Gated_TXT/Transformed_2968-WC-4_01_normalized.fcs/centroids1.txt",header=TRUE)
centroids_5 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2968/Gated_TXT/Transformed_2968-WC-5_01_normalized.fcs/centroids1.txt",header=TRUE)

quartz()
plot(centroids_1$DNA.1,centroids_1$DNA.2,xlim=c(0,4096),ylim=c(0,4096),pch=".",col="blue")

quartz()
ggplot(data=centroids_1,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2968(ap) visit 1")+
  geom_rect(aes(xmin=1331.2,xmax=2355.2,ymin=1638.4,ymax=2293.76),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")

quartz()
ggplot(data=centroids_2,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2968(ap) visit 2")+
  geom_rect(aes(xmin=1331.2,xmax=2355.2,ymin=1638.4,ymax=2293.76),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")

quartz()
ggplot(data=centroids_3,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2968(ap) visit 3")+
  geom_rect(aes(xmin=1331.2,xmax=2355.2,ymin=1638.4,ymax=2293.76),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")

quartz()
ggplot(data=centroids_4,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2968(ap) visit 4")+
  geom_rect(aes(xmin=1331.2,xmax=2355.2,ymin=1638.4,ymax=2293.76),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")

quartz()
ggplot(data=centroids_5,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2968(ap) visit 5")+
  geom_rect(aes(xmin=1331.2,xmax=2355.2,ymin=1638.4,ymax=2293.76),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")


centroids_6 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2687/Gated_TXT/Transformed_2687-WC-1_01_normalized.fcs/centroids1.txt",header=TRUE)
centroids_7 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2687/Gated_TXT/Transformed_2687-WC-2_01_normalized.fcs/centroids1.txt",header=TRUE)
centroids_8 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2687/Gated_TXT/Transformed_2687-WC-3_01_normalized.fcs/centroids1.txt",header=TRUE)
centroids_9 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2687/Gated_TXT/Transformed_2687-WC-4_01_normalized.fcs/centroids1.txt",header=TRUE)
centroids_10 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2687/Gated_TXT/Transformed_2687-WC-5_01_normalized.fcs/centroids1.txt",header=TRUE)

quartz()
ggplot(data=centroids_6,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2687(wp) visit 1")+
  geom_rect(aes(xmin=1331.2,xmax=2396.16,ymin=1638.4,ymax=2457.6),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")

quartz()
ggplot(data=centroids_7,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2687(wp) visit 2")+
  geom_rect(aes(xmin=1331.2,xmax=2396.16,ymin=1638.4,ymax=2457.6),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")

quartz()
ggplot(data=centroids_8,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2687(wp) visit 3")+
  geom_rect(aes(xmin=1331.2,xmax=2396.16,ymin=1638.4,ymax=2457.6),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")

quartz()
ggplot(data=centroids_9,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2687(wp) visit 4")+
  geom_rect(aes(xmin=1331.2,xmax=2396.16,ymin=1638.4,ymax=2457.6),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")

quartz()
ggplot(data=centroids_10,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2687(wp) visit 5")+
  geom_rect(aes(xmin=1331.2,xmax=2396.16,ymin=1638.4,ymax=2457.6),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")

centroids_11 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2383/Gated_TXT/Transformed_2383-WC-1_02_normalized.fcs/centroids1.txt",header=TRUE)
centroids_12 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2383/Gated_TXT/Transformed_2383-WC-2_01_normalized.fcs/centroids1.txt",header=TRUE)
centroids_13 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2383/Gated_TXT/Transformed_2383-WC-3_01_normalized.fcs/centroids1.txt",header=TRUE)
centroids_14 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2383/Gated_TXT/Transformed_2383-WC-4_01_normalized.fcs/centroids1.txt",header=TRUE)
centroids_15 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2383/Gated_TXT/Transformed_2383-WC-5_01_normalized.fcs/centroids1.txt",header=TRUE)

quartz()
ggplot(data=centroids_11,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2383(wp) visit 1")+
  geom_rect(aes(xmin=1331.2,xmax=2396.16,ymin=1638.4,ymax=2334.72),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")
quartz()
ggplot(data=centroids_12,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2383(wp) visit 2")+
  geom_rect(aes(xmin=1331.2,xmax=2396.16,ymin=1638.4,ymax=2334.72),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")
quartz()
ggplot(data=centroids_13,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2383(wp) visit 3")+
  geom_rect(aes(xmin=1331.2,xmax=2396.16,ymin=1638.4,ymax=2334.72),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")
quartz()
ggplot(data=centroids_14,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2383(wp) visit 4")+
  geom_rect(aes(xmin=1331.2,xmax=2396.16,ymin=1638.4,ymax=2334.72),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")
quartz()
ggplot(data=centroids_15,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2383(wp) visit 5")+
  geom_rect(aes(xmin=1331.2,xmax=2396.16,ymin=1638.4,ymax=2334.72),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")


centroids_16 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2912/Gated_TXT/Transformed_2912_WC-1_01_normalized.fcs/centroids1.txt",header=TRUE)
centroids_17 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2912/Gated_TXT/Transformed_2912_WC-2_01_normalized.fcs/centroids1.txt",header=TRUE)
centroids_18 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2912/Gated_TXT/Transformed_2912_WC-3_01_normalized.fcs/centroids1.txt",header=TRUE)
centroids_19 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2912/Gated_TXT/Transformed_2912_WC-4_01_normalized.fcs/centroids1.txt",header=TRUE)
centroids_20 <- read.table("/Users/amandava/Desktop/Pertussis/New_Donors/2912/Gated_TXT/Transformed_2912_WC-5_01_normalized.fcs/centroids1.txt",header=TRUE)
quartz()
ggplot(data=centroids_16,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2912(ap) visit 1")+
  geom_rect(aes(xmin=1372.16,xmax=2457.6,ymin=1699.84,ymax=2457.6),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")

quartz()
ggplot(data=centroids_17,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2912(ap) visit 2")+
  geom_rect(aes(xmin=1372.16,xmax=2457.6,ymin=1699.84,ymax=2457.6),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")

quartz()
ggplot(data=centroids_18,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2912(ap) visit 3")+
  geom_rect(aes(xmin=1372.16,xmax=2457.6,ymin=1699.84,ymax=2457.6),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")

quartz()
ggplot(data=centroids_19,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2912(ap) visit 4")+
  geom_rect(aes(xmin=1372.16,xmax=2457.6,ymin=1699.84,ymax=2457.6),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")

quartz()
ggplot(data=centroids_20,aes(x=DNA.1,y=DNA.2))+geom_point(shape=".",colour="blue")+xlim(0,4096)+ylim(0,4096)+
  xlab("DNA.1")+ylab("DNA.2")+ggtitle("Cluster centroids of DNA.1 and DNA.2 for subject 2912(ap) visit 5")+
  geom_rect(aes(xmin=1372.16,xmax=2457.6,ymin=1699.84,ymax=2457.6),colour="black",fill=NA,size=0.1)+
  annotate("text",x=2600,y=2000,label="DAFi gate")



