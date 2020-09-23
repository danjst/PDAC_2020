setwd("/Users/danieljasontan/metrics_files/")

counts <- read.csv("/Users/danieljasontan/metrics_files/Spladder_ES_0.075.csv", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
#here its ok to have features as rows, samples as columns, the below function (combine_groupings) will flip it
#sdev_cutoff <- read.csv("/Users/danieljasontan/metrics_files/TE_psi_0.1_events.csv",row.names = 1, stringsAsFactors = FALSE)[,1]
#counts <- subset(counts, rownames(counts) %in% sdev_cutoff)
#skip this: counts <- as.matrix(counts[76:ncol(counts)])
groupings <- read.csv("/Users/danieljasontan/metrics_files/ES_lee_clusters.csv", header = FALSE, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
groupings <- groupings[colnames(counts),,drop=FALSE]
counts <- cbind(t(counts),groupings)

combine_groupings <- function(counts, groupings){
  groupings_ordered <- groupings[colnames(counts),,drop=FALSE]
  return(cbind(t(counts),groupings))
}

ss <- function(x) {
  # calculate sum of squares per one feature
  return(sum((x-mean(x))^2))
}

rmsstd_per_cluster <- function(x) {
  # calculate rmsstd for one cluster
  pooled_ss <- apply(x, 2, ss)
  return(sqrt(sum(pooled_ss)/(nrow(x)*(nrow(x)-1))))
}

rmsstd <- function(x, groupings = ncol(x)) {
  # calculates "rmsstd" across all clusters (sum of rmsstd by cluster)
  # before this function, features are columns, samples are rows
  # groupings is the column of x with the groups
  rmsstd_individual <- c()
  for (each in unique(x[,ncol(x)])){
    cluster <- subset(x, x[,ncol(x)] == each)
    rmsstd_individual <- c(rmsstd_individual, rmsstd_per_cluster(cluster[1:ncol(cluster)-1]))
  }
  return(sum(rmsstd_individual))
}

ss_within <- function(x, groupings = ncol(x)) {
  # calculates the sum of square within across all samples
  ss_within_individual <- c()
  for (each in unique(x[,ncol(x)])){
    cluster <- subset(x, x[,ncol(x)] == each)
    ss_within_individual <- c(ss_within_individual, sum(apply(cluster[1:ncol(cluster)-1],2,ss)))
  }
  return(sum(ss_within_individual))
}

rs <- function(x, groupings = ncol(x)) {
  # calculates rs index across all samples
  ss_within_all_clusters <- ss_within(x)
  ss_tot <- sum(apply(x[1:ncol(x)-1],2,ss))
  return((ss_tot - ss_within_all_clusters) / ss_tot)
}

library(clusterCrit)

crit <- intCriteria(as.matrix(counts[,1:ncol(counts)-1]),counts[,ncol(counts)],c("Det_Ratio","SD_Scat", "SD_Dis"))
det_ind <- crit$det_ratio
nclust <- length(unique(counts[,ncol(counts)]))
SD_index <- nclust*crit$sd_scat + crit$sd_dis

internal_valid <- c(rmsstd(counts), rs(counts), det_ind, SD_index)
names(internal_valid) <- c("RMSSTD","RS","Det Ind", "SD Ind")
internal_valid
