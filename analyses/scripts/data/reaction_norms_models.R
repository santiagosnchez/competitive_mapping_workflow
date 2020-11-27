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

# read data
rlogdata = read.csv("tables/MDS_plot/rlog-transformed_exprdata.csv", row.names=1)
models = read.csv("tables/species_by_sex/df.species_sex.csv", row.names=1)
models.levels = sort(unique(as.character(models$models)))[c(8,5,13,7,4,12,3,11,6,2,10,1,9)]

# get numeric names of models
models_num = as.character(models$models2)

# clean and merge data
rlogdata = rlogdata[,c(1:6,13:18)]
scaledata = t(scale(t(rlogdata)))
scaledata = scaledata[rownames(models),]
scaledata = data.frame(CbrF=rowMeans(scaledata[,1:3]), CbrM=rowMeans(scaledata[,4:6]), CniF=rowMeans(scaledata[,7:9]), CniM=rowMeans(scaledata[,10:12]))

# add data
scaledata$models = models$models
scaledata$genes = rownames(models)

# melt dataframe
df.models = pivot_longer(scaledata, c(-models, -gene), names_to="samples", values_to="normalized expression")
df.models$species = substr(df.models$samples, 1, 3)
df.models$species = c("C. briggsae","C. nigoni")[ factor(df.models$species) ]
df.models$sex = substr(df.models$samples, 4, 4)
df.models$sex = c("female","male")[ factor(df.models$sex) ]
df.models$sex[ grep("sex-neutral", df.models$models) ] = "sex-neutral"
df.models$models_num = models_num[ df.models$models ]

# get centroids
df.models.cent = as.data.frame(df.models %>% group_by(models, samples) %>% summarize(expr=mean(expression), upper=quantile(expression, .95), lower=quantile(expression, .05)) %>% select(models, samples, expression = expr, upper, lower))
df.models.cent$species = substr(df.models.cent$samples, 1, 3)
df.models.cent$sex = substr(df.models.cent$samples, 4, 4)
df.models.cent$gene = "centroid"
df.models.cent$species = c("C. briggsae","C. nigoni")[ factor(df.models.cent$species) ]
df.models.cent$sex = c("female","male")[ factor(df.models.cent$sex) ]

# colors
cols.sex = c("#80CDC1","#018571","#BEBEBE")

# plot
ggplot(df.models, aes(x=species, y=expression, color=sex, group=interaction(gene, sex))) +
    geom_line(alpha=0.1) +
    geom_boxplot(alpha=0.5) +
    geom_line(data=df.models.cent, color="black", linetype=2, size=1.5) +
    scale_color_manual(values=cols.sex) +
    scale_x_discrete(expand=c(0.1,0.1)) +
    facet_rep_wrap(~models)
