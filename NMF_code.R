library(NMF)
library(plyr)
library(data.table)
library(survival)
library(survminer)
library(tidyverse)

MES_events_total<-read.csv("spladder_MES.csv")
A3_events_total<-read.csv("spladder_A3.csv")
A5_events_total<-read.csv("spladder_A5.csv")
ES_events_total<-read.csv("spladder_ES.csv")
IR_events_total<-read.csv("spladder_IR.csv")

MES_events<-na.omit(ME_events_total)
A3_events<-na.omit(A3_events_total)
A5_events<-na.omit(A5_events_total)
ES_events<-na.omit(ES_events_total)
IR_events<-na.omit(IR_events_total)

MES_events_SD<-transform(MES_data, SD=apply(MES_events,1, sd, na.rm = TRUE))
A3_events_SD<-transform(A3_data, SD=apply(A3_events,1, sd, na.rm = TRUE))
A5_events_SD<-transform(A5_data, SD=apply(A5_events,1, sd, na.rm = TRUE))
ES_events_SD<-transform(ES_data, SD=apply(ES_events,1, sd, na.rm = TRUE))
IR_events_SD<-transform(IR_data, SD=apply(IR_events,1, sd, na.rm = TRUE))

MES_0.07_nmf<- subset(MES_events_SD, SD >0.07)
A3_0.06_nmf<- subset(A3_events_SD, SD >0.06)
A5_0.055_nmf<- subset(A5_events_SD, SD >0.055)
ES_0.075_nmf<- subset(ES_events_SD, SD >0.075)
IR_0.1_nmf<- subset(IR_events_SD, SD >0.1)

MES_0.07_nmf <- subset(MES_0.07_nmf, select = -c(SD))
A3_0.06_nmf <- subset(A3_0.06_nmf, select = -c(SD))
A5_0.055_nmf <- subset(A5_0.055_nmf, select = -c(SD))
ES_0.075_nmf <- subset(IR_0.075_nmf, select = -c(SD))
IR_0.1_nmf <- subset(IR_0.1_nmf, select = -c(SD))

estim.r.MES<-nmf(MES_0.07_nmf,2:4,nrun=100,seed=123456,method="lee")
estim.r.A3<-nmf(A3_0.06_nmf,2:4,nrun=100,seed=123456,method="lee")
estim.r.A5<-nmf(A5_0.055_nmf,2:4,nrun=100,seed=123456,method="lee")
estim.r.ES<-nmf(ES_0.075_nmf,2:4,nrun=100,seed=123456,method="lee")
estim.r.IR<-nmf(IR_0.1_nmf,2:4,nrun=100,seed=123456,method="lee")

MES_optimal_k_index<-which.max(estim.r.MES.lee[["measures"]][["cophenetic"]])
MES_optimal_k<-estim.r.MES.lee[["measures"]][["rank"]][MES_optimal_k_index]
A3_optimal_k_index<-which.max(estim.r.A3.lee[["measures"]][["cophenetic"]])
A3_optimal_k<-estim.r.A3.lee[["measures"]][["rank"]][A3_optimal_k_index]
A5_optimal_k_index<-which.max(estim.r.A5.lee[["measures"]][["cophenetic"]])
A5_optimal_k<-estim.r.A5.lee[["measures"]][["rank"]][A5_optimal_k_index]
ES_optimal_k_index<-which.max(estim.r.ES.lee[["measures"]][["cophenetic"]])
ES_optimal_k<-estim.r.ES.lee[["measures"]][["rank"]][ES_optimal_k_index]
IR_optimal_k_index<-which.max(estim.r.IR.lee[["measures"]][["cophenetic"]])
IR_optimal_k<-estim.r.IR.lee[["measures"]][["rank"]][IR_optimal_k_index]

estim.r.MES_optimal<-nmf(MES_0.07_nmf,MES_optimal_k,nrun=100,seed=123456,method="lee")
estim.r.A3_optimal<-nmf(A3_0.06_nmf,A3_optimal_k,nrun=100,seed=123456,method="lee")
estim.r.A5_optimal<-nmf(A5_0.055_nmf,A5_optimal_k,nrun=100,seed=123456,method="lee")
estim.r.ES_optimal<-nmf(ES_0.075_nmf,ES_optimal_k,nrun=100,seed=123456,method="lee")
estim.r.IR_optimal<-nmf(IR_0.1_nmf,IR_optimal_k,nrun=100,seed=123456,method="lee")

MES_clusters<-data.frame(predict(estim.r.MES_optimal))
A3_clusters<-data.frame(predict(estim.r.A3_optimal))
A5_clusters<-data.frame(predict(estim.r.A5_optimal))
ES_clusters<-data.frame(predict(estim.r.ES_optimal))
IR_clusters<-data.frame(predict(estim.r.IR_optimal))

MES_clusters<-data.frame(predict(estim.r.ME.lee.r2))
A3_clusters<-data.frame(predict(estim.r.A3.lee.r2))
A5_clusters<-data.frame(predict(estim.r.A5.r2))
ES_clusters<-data.frame(predict(estim.r.ES.0.075.lee.r2))
IR_clusters<-data.frame(predict(estim.r.IR.lee.r2))


clusters<- tibble::lst(MES_clusters,A3_clusters,A5_clusters,ES_clusters,IR_clusters)

for (i in 1:length(clusters)){
  rownames(clusters[[i]])<-gsub("\\.","-",rownames(clusters[[i]]))
  colnames(clusters[[i]])[1]<-names(clusters)[i]
  clusters[[i]]<-setDT(as.data.frame((clusters[[i]])), keep.rownames = TRUE)[]
  colnames(clusters[[i]])[1]<-"id"
}

clusters[[6]]<-pdac_76_survival
test<-Reduce(function(x, y) merge(x, y, by="id", all=TRUE), clusters)

res<- survfit(Surv(as.numeric(PFI.time),as.numeric(PFI)) ~ IR_1 ,data = pdac_185_survival)

