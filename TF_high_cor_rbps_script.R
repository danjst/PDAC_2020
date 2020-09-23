TF_chea3_IR1<-data.frame(TF_chea3_IR1)


the_list<-list()

for (i in 1:nrow(TF_chea3_IR1)){
sig_rbp_targets<-strsplit(TF_chea3_IR1[i,18],",")[[1]]
matches<-c()
for (j in 1:length(sig_rbp_targets)){
  
  if (sig_rbp_targets[j] %in% high_cor_rbps){
    matches<-c(matches,sig_rbp_targets[j])
  }
  
}
the_list[[i]]<-matches


}
names(the_list)<-TF_chea3_IR1$gene_name

#for (is.na(TF_chea3_IR1$number.of.sig.RBP.target.matches))
num_highcor<-c()
for (i in 1:length(the_list)){
  num_highcor<-c(num_highcor,length(the_list[[i]]))
}

for (i in 1:length(the_list)){
  if (length(the_list[[i]])>=2){
    the_list[[i]]<-paste(the_list[[i]],collapse=",")
  }
}

for (i in 1:length(the_list)){
  if (length(is.na(the_list[[i]]))==0){
    the_list[[i]]<-0
    
  }
}

for (i in 1:length(the_list)){
  TF_chea3_IR1[i,19]<-the_list[[i]]
}
colnames(TF_chea3_IR1)<-"high_cor_rbps"

TF_chea3_IR1<-cbind(TF_chea3_IR1,num_highcor)

TF_chea3_IR1$high_correlation_rbps[TF_chea3_IR1$high_correlation_rbps=="0"]<-NA





