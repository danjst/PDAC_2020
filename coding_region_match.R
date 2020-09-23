#Coding Region Coordinate Matching

coding <- data.frame(matrix(ncol = 2, nrow = 100))
colnames(coding)<-c("event_id","transcript")
for (i in 1:nrow(the_4852_trial)){
  for (j in 1:nrow(coding_Bed)){
    if (the_4852_trial$contig[i]==coding_Bed$chromosome[j] & the_4852_trial$strand[i]==coding_Bed$strand_bed[j] & the_4852_trial$intron_start[i]>=coding_Bed$coding_1[j] & the_4852_trial$intron_end[i]<=coding_Bed$coding_2[j]){
      coding[i,1]<-the_4852_trial$event_id[i]
      coding[i,2]<-coding_Bed$transcript[j]
      break
    }
  }
}


final<-fuzzy_left_join(the_4852_trial, coding_Bed,
                by = c("contig" = "chromosome", "strand" = "strand_bed","intron_start"="coding_1","intron_end"="coding_2"), 
                match_fun = list(`==`, `==`,`>=`,`<=`))


unknown_4852_search2<-unknown_4852 %>% 
  mutate(dummy=TRUE) %>%
  left_join(X3_prime_bed %>% mutate(dummy=TRUE)) %>%
  filter(contig==chromosome, strand==strand_bed,intron_start>=three_start,intron_end<=three_end) %>% group_by(event_id) %>%slice(1) %>%
  select(-dummy)

DS_262_pseudo<-DS_262_possible_pseudo%>%
  mutate(dummy=TRUE) %>%
  left_join(pseudogene_bed %>% mutate(dummy=TRUE)) %>%
  filter(contig==chromosome, strand==strand_bed,intron_start>=full_1,intron_end<=full_2) %>% group_by(event_id) %>%slice(1) %>%
  select(-dummy)

DS_262_unknown_match<-DS_262_unknown%>% 
  mutate(dummy=TRUE) %>%
  left_join(coding_Bed %>% mutate(dummy=TRUE)) %>%
  filter(contig==chromosome, strand==strand_bed,intron_start>=full_1,intron_end<=full_2) %>% group_by(event_id) %>%slice(1) %>%
  select(-dummy)

final3<-the_4852_trial %>%
  mutate(dummy=TRUE) %>% left_join(coding_Bed %>% mutate(dummy=TRUE))




#############################
#UTR Matching

for (row in 1:nrow(utr_match_262x)){
    if (utr_match_262x[row,]$strand_bed=="+"){
        if (utr_match_262x[row,]$intron_start>=utr_match_262x[row,]$full_1 & utr_match_262x[row,]$intron_end<=utr_match_262x[row,]$coding_1){
            utr_match_262x[row,]$utr="5' UTR"}
        else {utr_match_262x[row,]$utr="3' UTR"}
    }
    if (utr_match_262x[row,]$strand_bed=="-"){
        if (utr_match_262x[row,]$intron_start>=utr_match_262x[row,]$full_1 & utr_match_262x[row,]$intron_end<=utr_match_262x[row,]$coding_1){
            utr_match_262x[row,]$utr="3' UTR"}
        else {utr_match_262x[row,]$utr="5' UTR"}
    }
}

