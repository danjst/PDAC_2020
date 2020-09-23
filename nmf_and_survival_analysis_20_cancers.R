nmf_part<-function(x,pval_index){
  estim<-nmf(x,2,nrun=150,seed=123456,method="lee")
  #pval_index<-(ncol(tcga_IR_pvals)+1)
  clusters<-data.frame(predict(estim,what="consensus"))
  g=clusters[,1]
  clus_1<-sum(g=="1")
  clus_2<-sum(g=="2")
  tcga_IR_pvals[4,pval_index]<<-clus_1
  tcga_IR_pvals[5,pval_index]<<-clus_2
  tcga_IR_pvals[6,pval_index]<<-length(extractFeatures(estim)[[1]])
  tcga_IR_pvals[7,pval_index]<<-length(extractFeatures(estim)[[2]])
  
  tcga_IR_pvals[8,pval_index]<<-nrow(x)
  
  tcga_IR_pvals[9,pval_index]<<-ncol(x)
  clusters<-data.frame(setDT(as.data.frame((clusters)), keep.rownames = TRUE)[])
  colnames(clusters)<-c("Name","cluster")
  clusters$Name<-substr(clusters$Name,start=1,stop=12)
  clusters$Name<-gsub("\\.","-",clusters$Name)
  survival_table<<-merge(clusters,all_survival_data,c("Name"))
  res.os<- survfit(Surv(as.numeric(OS.time),as.numeric(OS)) ~ cluster ,data = survival_table)
  res.pfi<- survfit(Surv(as.numeric(PFI.time),as.numeric(PFI)) ~ cluster ,data = survival_table)
  res.dss<- survfit(Surv(as.numeric(DSS.time),as.numeric(DSS)) ~ cluster ,data = survival_table)
  os.p<-surv_pvalue(res.os)$pval
  pfi.p<-surv_pvalue(res.pfi)$pval
  dss.p<-surv_pvalue(res.dss)$pval
  all_pval<-c(os.p,pfi.p,dss.p)
  tcga_IR_pvals[1:3,pval_index]<<-all_pval
  return(estim)
  
}

surv_table<-function(x){
  
  clusters<-data.frame(predict(x,what="consensus"))
  clusters<-data.frame(setDT(as.data.frame((clusters)), keep.rownames = TRUE)[])
  colnames(clusters)<-c("Name","cluster")
  clusters$Name<-substr(clusters$Name,start=1,stop=12)
  clusters$Name<-gsub("\\.","-",clusters$Name)
  survival_table<<-merge(clusters,all_survival_data,c("Name"))
  res.os<- survfit(Surv(as.numeric(OS.time),as.numeric(OS)) ~ cluster ,data = survival_table)
  res.pfi<- survfit(Surv(as.numeric(PFI.time),as.numeric(PFI)) ~ cluster ,data = survival_table)
  res.dss<- survfit(Surv(as.numeric(DSS.time),as.numeric(DSS)) ~ cluster ,data = survival_table)
  os.p<-surv_pvalue(res.os)$pval
  pfi.p<-surv_pvalue(res.pfi)$pval
  dss.p<-surv_pvalue(res.dss)$pval
  all_pval<-c(os.p,pfi.p,dss.p)
  return(all_pval)
}


fix_everything<-function(x,y=0.1){
  x<-x[grepl("-01A|-01B",x)]
matches<-gsub("-",".",x)

  matches <- unique (grep(paste(matches,collapse="|"),colnames(tcga), value=TRUE))
  matches<-matches[!duplicated(substring(matches, first=0,last = 12))]
  my_data<-tcga[,matches]
  rownames(my_data)<-tcga$event_id
  my_data<-na.omit(my_data)
  my_data<-transform(my_data, SD=apply(my_data,1, sd, na.rm = TRUE))
  my_data<- subset(my_data, SD >y)
  my_data<- subset(my_data, select = -c(SD))
  return(my_data)
}
#how to run fix_everything:
#gbm_fixed<-fix_everything(gbm7,y=0.1)

kirc_samples<-c(colnames(tcga)[1:6],colnames(kirc_fixed))
all_kirc_events<-tcga[ ,colnames(tcga) %in% kirc_samples, ]
colnames(all_kirc_events)[7:170]<-substr(colnames(all_kirc_events)[7:170],start=1,stop=12)
colnames(all_kirc_events)[7:170]<-gsub("\\.","-",colnames(all_kirc_events)[7:170])

# first filter out events with more than 10% NA
kirc_less_10_percent_NA <- all_kirc_events[rowSums(is.na(all_kirc_events)) < 16, ]

###Next filtering for less than 0.05 abs NA difference between 2 clusters
#make kirc_cluster1 and kirc_cluster2 first
kirc_10na_2<-kirc_less_10_percent_NA[,kirc_cluster2]

kirc_2_NAs <- apply(kirc_10na_2, 1, function(x) sum(is.na(x)))
kirc_2_NAs_percent<-kirc_2_NAs/107
kirc_1_minus_2_percent<-kirc_1_NAs_percent-kirc_2_NAs_percent
kirc_1_minus_2_percent<-abs(kirc_1_minus_2_percent)

kirc_over_0.05<-c()
for (i in 1:nrow(kirc_less_10_percent_NA)){
  if (kirc_1_minus_2_percent[[i]]>0.05){kirc_over_0.05<-c(kirc_over_0.05,names(kirc_1_minus_2_percent[i]))}
}
kirc_less_10_percent_NA_also_0.05 <- kirc_less_10_percent_NA[!row.names(kirc_less_10_percent_NA)%in%kirc_over_0.05,]
###doing the 0.1 mean psi difference 
kirc_10na_0.05_1<-kirc_less_10_percent_NA_also_0.05[,kirc_cluster1]
kirc_10na_0.05_2<-kirc_less_10_percent_NA_also_0.05[,kirc_cluster2]

kirc_rowmeans_1<-rowMeans(kirc_10na_0.05_1,na.rm = T)
kirc_rowmeans_2<-rowMeans(kirc_10na_0.05_2,na.rm = T)

kirc_rowmeans_difference<-kirc_rowmeans_1-kirc_rowmeans_2
kirc_rowmeans_difference_greater_0.1<-kirc_rowmeans_difference[abs(kirc_rowmeans_difference)>0.1]
kirc_final_event_rownames<-names(kirc_rowmeans_difference_greater_0.1)

kirc_final_event_rownames<-names(kirc_rowmeans_difference_greater_0.1)
kirc_final_ds_events <- kirc_less_10_percent_NA_also_0.05[row.names(kirc_less_10_percent_NA_also_0.05)%in%kirc_final_event_rownames,]      
kirc_final_event_meanpsi_values<-as.numeric(paste(unlist(kirc_rowmeans_difference_greater_0.1)))

kirc_final_ds_events<-cbind(kirc_final_ds_events,kirc_final_event_meanpsi_values)

kirc_final_ds_events_x<-cbind(kirc_final_ds_events[1:6],kirc_final_event_meanpsi_values)
kirc_final_ds_events_x$gene_no_version<-substr(kirc_final_ds_events_x$gene_old_ensembl,start=1,stop=15)

kirc_final_ds_events_x<-join_all(list(kirc_final_ds_events_x,gencodev22_bed_2020), by = c('gene_no_version'), type = "left", match = "all")   

printedx<-function(x){
  x<-as.character(x)
  return(x)
}

make_surv_table<-function(estim){
  clusters<-data.frame(predict(estim,what="consensus"))
  clusters<-data.frame(setDT(as.data.frame((clusters)), keep.rownames = TRUE)[])
  colnames(clusters)<-c("Name","cluster")
  clusters$Name<-substr(clusters$Name,start=1,stop=12)
  clusters$Name<-gsub("\\.","-",clusters$Name)
  survival_table<-merge(clusters,all_survival_data,c("Name"))
  return(survival_table)
}
