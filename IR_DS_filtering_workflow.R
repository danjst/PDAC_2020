

# first filter out events with more than 10% NA
pdac_less_10_percent_NA <- all_pdac_events[rowSums(is.na(all_pdac_events)) < 0.1 * ncol(all_pdac_events), ]

###Next filtering for less than 0.05 abs NA difference between 2 clusters
#make pdac_cluster1 and pdac_cluster2 first
pdac_10na_2<-pdac_less_10_percent_NA[,pdac_cluster2]
pdac_10na_1<-pdac_less_10_percent_NA[,pdac_cluster1]

pdac_1_NAs <- apply(pdac_10na_1, 1, function(x) sum(is.na(x)))
pdac_1_NAs_percent<-pdac_2_NAs/40
pdac_2_NAs <- apply(pdac_10na_2, 1, function(x) sum(is.na(x)))
pdac_2_NAs_percent<-pdac_2_NAs/36
pdac_1_minus_2_percent<-pdac_1_NAs_percent-pdac_2_NAs_percent
pdac_1_minus_2_percent<-abs(pdac_1_minus_2_percent)

pdac_over_0.05<-c()
for (i in 1:nrow(pdac_less_10_percent_NA)){
  if (pdac_1_minus_2_percent[[i]]>0.05){pdac_over_0.05<-c(pdac_over_0.05,names(pdac_1_minus_2_percent[i]))}
}
pdac_less_10_percent_NA_also_0.05 <- pdac_less_10_percent_NA[!row.names(pdac_less_10_percent_NA)%in%pdac_over_0.05,]

#join with DS list
#re-adjust p-values here
#check for significance here, filter padj<0.05
###doing the 0.1 mean psi difference 
pdac_10na_0.05_1<-pdac_less_10_percent_NA_also_0.05[,pdac_cluster1]
pdac_10na_0.05_2<-pdac_less_10_percent_NA_also_0.05[,pdac_cluster2]

pdac_rowmeans_1<-rowMeans(pdac_10na_0.05_1,na.rm = T)
pdac_rowmeans_2<-rowMeans(pdac_10na_0.05_2,na.rm = T)

pdac_rowmeans_difference<-pdac_rowmeans_1-pdac_rowmeans_2
pdac_rowmeans_difference_greater_0.1<-pdac_rowmeans_difference[abs(pdac_rowmeans_difference)>0.1]
pdac_final_event_rownames<-names(pdac_rowmeans_difference_greater_0.1)

pdac_final_ds_events <- pdac_less_10_percent_NA_also_0.05[row.names(pdac_less_10_percent_NA_also_0.05)%in%pdac_final_event_rownames,]
pdac_final_event_meanpsi_values<-as.numeric(paste(unlist(pdac_rowmeans_difference_greater_0.1)))

pdac_final_ds_events<-cbind(pdac_final_ds_events,pdac_final_event_meanpsi_values)

pdac_final_ds_events_x<-cbind(pdac_final_ds_events[1:6],pdac_final_event_meanpsi_values)
pdac_final_ds_events_x$gene_no_version<-substr(pdac_final_ds_events_x$gene_old_ensembl,start=1,stop=15)

pdac_final_ds_events_x<-join_all(list(pdac_final_ds_events_x,gencodev22_bed_2020), by = c('gene_no_version'), type = "left", match = "all")
