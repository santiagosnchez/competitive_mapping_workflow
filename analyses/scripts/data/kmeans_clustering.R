library(ggplot2)
library(cowplot)
library(lemon)
library(ggridges)
library(tidyr)
library(dplyr)
library(gridExtra)
library(RColorBrewer)
library(ggsci)
library(parallel)
library(scales)
theme_set(theme_cowplot())

theme_set(theme_cowplot())
source("scripts/data/enrichment_fun.R")

# read data
rlogdata = read.csv("tables/MDS_plot/rlog-transformed_exprdata.csv", row.names=1)
rlogdata.means = data.frame(CbrF=rowMeans(rlogdata[,1:3]),
                            CbrM=rowMeans(rlogdata[,4:6]),
                            HF1F=rowMeans(rlogdata[,7:9]),
                            HF1M=rowMeans(rlogdata[,10:12]),
                            CniF=rowMeans(rlogdata[,13:15]),
                            CniM=rowMeans(rlogdata[,16:18]))
scaledata = t(scale(t(rlogdata.means)))

# clustering
# https://2-bitbio.com/2017/10/clustering-rnaseq-data-using-k-means.html
# set.seed(1)
kClust <- kmeans(scaledata, centers=15)
kClusters <- kClust$cluster

# centroids
clust.centroid = function(i, dat, clusters) {
  ind = (clusters == i)
  colMeans(dat[ind,])
}
kClustcentroids <- as.data.frame(sapply(levels(factor(kClusters)), clust.centroid, scaledata, kClusters))

colnames(kClustcentroids) = paste("M",1:15,sep="")
kClustcentroids$samples = rownames(kClustcentroids)
df.kClusterCentroids = pivot_longer(kClustcentroids, cols=-samples, names_to="models", values_to="expression")
df.kClusterCentroids$sex = c("female","male")[ factor(substr(df.kClusterCentroids$samples, 4, 4)) ]
df.kClusterCentroids$species = substr(df.kClusterCentroids$samples, 1, 3)

ggplot(df.kClusterCentroids, aes(x=species, y=expression, group=sex, color=sex)) +
    geom_line() +
    scale_x_discrete(limits=c("Cbr","HF1","Cni")) +
    facet_wrap(~models)
