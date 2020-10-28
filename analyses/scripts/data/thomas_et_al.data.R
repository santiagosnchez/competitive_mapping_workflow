# Analysis of Thomas et al. 2012 Curr Biol RNAseq data C. briggsae

library(limma)
library(edgeR)
library(DESeq2)
library(MASS)
library(ggplot2)
library(cowplot)
library(lemon)
library(tidyr)
library(dplyr)
library(gridExtra)
library(RColorBrewer)

theme_set(theme_cowplot())

# read raw count data
counts.species = read.csv("../counts/species_counts.txt", row.name=1, sep="\t")
counts.thomas = read.csv("../counts/thomas_Cbr.counts.txt", row.name=1, sep="\t")
rownames(counts.thomas) = rownames(counts.species)
counts = counts.species[,1:6]
colnames(counts) = c(paste0("Cbr_F",1:3), paste0("Cbr_M",1:3))

# group treatments
groups = gsub("[1-3]$","", colnames(counts))
sex = substr(groups,5,5)
analysis = rep("our_data",6)

# coldata
coldata = data.frame(samples=colnames(counts), groups=groups, analysis=analysis, sex=sex)
coldata.thomas = coldata
coldata.thomas$samples = colnames(counts.thomas)
coldata.thomas$analysis = rep("thomas",6)

# prepare DSseq2
dds = DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~sex)
dds = estimateSizeFactors(dds)
filter_in = rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 3
# Thomas et al.
dds.thomas = DESeqDataSetFromMatrix(countData = counts.thomas, colData = coldata.thomas, design = ~sex)
dds.thomas = estimateSizeFactors(dds.thomas)
filter_in.thomas = rowSums( counts(dds.thomas, normalized=TRUE) >= 10 ) >= 3
filter_in_both = apply(cbind(filter_in, filter_in.thomas), 1, all)
# filter data both
counts.both = cbind(counts, counts.thomas)
counts.both = counts.both[filter_in_both,]
# separate
dds = dds[filter_in_both,]
dds.thomas = dds.thomas[filter_in_both,]
dds.rlog = assays(rlog(dds))
dds.rlog.thomas = assays(rlog(dds.thomas))

# females
keep = c(1:3,7:9)
dds.females = DESeqDataSetFromMatrix(countData = counts.both[,keep], colData = rbind(coldata,coldata.thomas)[keep,], design = ~analysis)
dds.females = DESeq(dds.females)
dds.females.res = results(dds.females)
# males
keep = c(4:6,10:12)
dds.males = DESeqDataSetFromMatrix(countData = counts.both[,keep], colData = rbind(coldata,coldata.thomas)[keep,], design = ~analysis)
dds.males = DESeq(dds.males)
dds.males.res = results(dds.males)
# females vs males thomas
keep = c(7:12)
dds.sex.thomas = DESeqDataSetFromMatrix(countData = counts.both[,keep], colData = rbind(coldata,coldata.thomas)[keep,], design = ~sex)
dds.sex.thomas = DESeq(dds.sex.thomas)
dds.sex.thomas.res = results(dds.sex.thomas)
# females vs males our data
keep = c(1:6)
dds.sex.our = DESeqDataSetFromMatrix(countData = counts.both[,keep], colData = rbind(coldata,coldata.thomas)[keep,], design = ~sex)
dds.sex.our = DESeq(dds.sex.our)
dds.sex.our.res = results(dds.sex.our)

# identify hermaphroditic genes
herm = dds.sex.our.res$padj < 0.05 &
    dds.sex.thomas.res$padj < 0.05 & dds.males.res$padj > 0.05 &
    dds.females.res$padj < 0.05 & dds.sex.our.res$log2FoldChange > 0 &
    dds.sex.thomas.res$log2FoldChange > 0

df.herm = na.omit(data.frame(logFC.sex.our=dds.sex.our.res$log2FoldChange,
                    logFC.sex.thomas=dds.sex.thomas.res$log2FoldChange,
                    herm=c("male-female(fog2)","hermaphroditic")[ factor(herm) ]))
write.csv(df.herm, file="tables/Thomas_et_al_DESeq_logFC.csv", row.names=F)

ggplot(df.herm, aes(logFC.sex.our, logFC.sex.thomas, color=herm)) +
    geom_vline(xintercept=0, linetype=2) +
    geom_hline(yintercept=0, linetype=2) +
    geom_point(alpha=0.3) +
    labs(x="log2-FC male:hermaphrodite (our data)", y="log2-FC male:fog2-female (Thomas et al.)") +
    scale_color_manual(values=c("darkblue","darkgreen"), name="", guide=guide_legend(override.aes=list(alpha=1))) +
    background_grid(major="xy", minor="xy")

# save figure
ggsave(file="figures/suppl_Thomas_et_al_hermaphroditic.png")
ggsave(file="figures/suppl_Thomas_et_al_hermaphroditic.pdf", device="pdf", useDingbats=F)

# save data
herm_genes = na.omit(rownames(counts.both)[ herm ])
cat(herm_genes, sep="\n", file="tables/Thomas_et_al_herm_genes.txt")
