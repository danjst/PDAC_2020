library(lineup)
prad_rbp_correalation<-corbetw2mat(prad_sig_rbp_norm_counts_t,prad_top_events_psi_467_t,what="all")

prad_rbp_correalation_pvalues<- data.frame(prad_rbp_correalation) # this and line below is to make dataframe with same colnames nad rownames
prad_rbp_correalation_pvalues[]<-0 

for (i in 1:ncol(prad_sig_rbp_norm_counts_t)){
  
  for (j in 1:ncol(prad_top_events_psi_467_t)){
    x<-cor.test(prad_sig_rbp_norm_counts_t[,i],prad_top_events_psi_467_t[,j])$p.value
    prad_rbp_correalation_pvalues[i,j]<-x
  }
}

prad_rbp_correalation_pvalues_interm<-as.vector(t(prad_rbp_correalation_pvalues)) #there are two things going on here: 1) t implicitly converts a data.frame to a matrix, 2) a matrix is just a special vector with dim attribute and as.vector or c removes it 

pval_adj<-p.adjust(prad_rbp_correalation_pvalues_interm,method="BH")
prad_correlation_padj<-matrix(pval_adj,ncol=545,byrow=TRUE)

prad_all_cor_adjusted<-list()

index<-1
for (x in 1:nrow(prad_rbp_correalation)){
  
  for (i in 1:ncol(prad_rbp_correalation)){
    if (abs(prad_rbp_correalation[x,i])>=0.7 & prad_correlation_padj[x,i]<0.05){
      prad_all_cor_adjusted[[index]]<-c(rownames(prad_rbp_correalation)[x],colnames(prad_rbp_correalation)[i],prad_rbp_correalation[x,i],prad_rbp_correalation_pvalues[x,i],prad_correlation_padj[x,i])
      index=index+1
      
    }
  }
}

 prad_all_cor_rbp<-sapply(prad_all_cor_adjusted, "[[", 1)
 prad_all_cor_event<-sapply(prad_all_cor_adjusted, "[[", 2)
 prad_all_cor_val<-sapply(prad_all_cor_adjusted, "[[", 3)
 prad_all_cor_pval<-sapply(prad_all_cor_adjusted,"[[",4)
 prad_all_cor_padj<-sapply(prad_all_cor_adjusted, "[[", 5)
prad_cor_results<-cbind(prad_all_cor_rbp,prad_all_cor_event,prad_all_cor_val,prad_all_cor_pval,prad_all_cor_padj)
prad_cor_results<-data.frame(prad_cor_results)
colnames(prad_cor_results)<-c("gene","event_id","correlation_value","cor_pval","cor_padj")
prad_cor_results<-join_all(list(prad_cor_results,prad_rbp), by = c('gene'), type = "left", match = "all")
colnames(prad_cor_results)[1]<-"RBP_ensembl"
colnames(prad_cor_results)[12]<-"RBP_gene_name"
#join all with event information (coordinates, mean psi, etc)
write.csv(prad_cor_results,"prad_cor_results.csv")
