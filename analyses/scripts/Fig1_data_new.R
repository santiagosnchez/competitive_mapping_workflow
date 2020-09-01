# generates files and data for Fig. 1

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
source("scripts/enrichment_fun.R")

# read raw count data
counts.species = read.csv("../counts/species_counts.txt", row.name=1, sep="\t")
counts.ase = read.csv("../counts/ase_counts.txt", row.name=1, sep="\t")
counts.species = counts.species[ rownames(counts.ase), ]

# merge ase counts and data
counts = cbind(counts.species[,1:6], counts.ase[,seq(1,12,2)] + counts.ase[,seq(2,12,2)], counts.species[,7:12])
colnames(counts) = c(paste0("Cbr_F",1:3), paste0("Cbr_M",1:3), paste0("HF1_F",1:3), paste0("HF1_M",1:3), paste0("Cni_F",1:3), paste0("Cni_M",1:3))

# group treatments
groups = gsub("[1-3]$","", colnames(counts))
sex = substr(groups,5,5)
species = substr(groups,1,3)

# get factor data
coldata = data.frame(samples=colnames(counts), groups=groups, species=species, sex=sex)

# prepare DSseq2
dds = DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~groups)
dds = estimateSizeFactors(dds)
filter_in = rowSums( counts(dds, normalized=TRUE) >= 10 ) >= 3
dds = dds[filter_in,]
dds_rlog = rlog(counts(dds))
rownames(dds_rlog) = rownames(dds)
write.csv(dds_rlog, "tables/rlog-transformed_exprdata.csv", quote=F)

# expression distance data
mds = plotMDS(dds_rlog)
df.mds = data.frame(dim1=mds$x, dim2=mds$y, species, sex)
write.csv(df.mds, "tables/expression_dist.mds_coords.csv")

# per species, F1 sex-biased differential expression
# C. briggsae
keep = species == "Cbr"
dds.briggsae = DESeqDataSetFromMatrix(countData = counts(dds)[,keep], colData = coldata[keep,], design = ~ sex)
dds.briggsae = DESeq(dds.briggsae)
dds.briggsae.res = results(dds.briggsae)
write.csv(as.data.frame(dds.briggsae.res), file="tables/deseq_res.briggsae.sex.csv", quote=F)
dds.briggsae.sum = as.data.frame(dds.briggsae.res) %>% mutate(Cbr = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>% select(Cbr)

# C. nigoni
keep = species == "Cni"
dds.nigoni = DESeqDataSetFromMatrix(countData = counts(dds)[,keep], colData = coldata[keep,], design = ~ sex)
dds.nigoni = DESeq(dds.nigoni)
dds.nigoni.res = results(dds.nigoni)
write.csv(as.data.frame(dds.nigoni.res), file="tables/deseq_res.nigoni.sex.csv", quote=F)
dds.nigoni.sum = as.data.frame(dds.nigoni.res) %>% mutate(Cni = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>% select(Cni)

# F1 hybrid
keep = species == "HF1"
dds.hybrids = DESeqDataSetFromMatrix(countData = counts(dds)[,keep], colData = coldata[keep,], design = ~ sex)
dds.hybrids = DESeq(dds.hybrids)
dds.hybrids.res = results(dds.hybrids)
write.csv(as.data.frame(dds.hybrids.res), file="tables/deseq_res.hybrids.sex.csv", quote=F)
dds.hybrids.sum = as.data.frame(dds.hybrids.res) %>% mutate(F1 = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>% select(F1)

# consolidate results
df.sex = cbind(dds.briggsae.sum, dds.nigoni.sum, dds.hybrids.sum)
df.sex[ df.sex == 0 ] = "sex-neutral"
df.sex[ df.sex == 1 ] = "male"
df.sex[ df.sex == -1 ] = "female"
df.sex[ is.na(df.sex) ] = "no expression"
df.sex$chromosome = sub("\\..*","", rownames(dds))
df.sex$gene = rownames(dds)
df.sex2 = pivot_longer(df.sex, cols=c(-chromosome,-gene), names_to="species", values_to="sex")

# generate counts and summaries
df.sex.count.spp = as.data.frame(df.sex2 %>% group_by(species, sex) %>% summarize(count=n()) %>% mutate(perc = round((count / sum(count))*100,0)))
write.csv(df.sex.count.spp, file="tables/df.counts.perc.DE.sex.csv", row.names=F)

# per chromosome
df.sex.count.chr = as.data.frame(df.sex2 %>% filter(sex != "no expression") %>% group_by(species, sex, chromosome) %>% summarize(count=n()))

# enrichments
df.sex.count.chr$p.value = NA
df.sex.count.chr$enrichment = NA
x = matrix((df.sex.count.chr %>% filter(species == "Cbr"))$count, ncol=3, nrow=6)
df.sex.count.chr$p.value[1:18] = as.vector(enrichment(x))
df.sex.count.chr$enrichment[1:18] = as.vector(enrichment(x, odds.ratio=T))
x = matrix((df.sex.count.chr %>% filter(species == "Cni"))$count, ncol=3, nrow=6)
df.sex.count.chr$p.value[19:36] = as.vector(enrichment(x))
df.sex.count.chr$enrichment[19:36] = as.vector(enrichment(x, odds.ratio=T))
x = matrix((df.sex.count.chr %>% filter(species == "F1"))$count, ncol=3, nrow=6)
df.sex.count.chr$p.value[37:54] = as.vector(enrichment(x))
df.sex.count.chr$enrichment[37:54] = as.vector(enrichment(x, odds.ratio=T))
df.sex.count.chr$chr.pseudo = (1:6)[ factor(df.sex.count.chr$chromosome) ]
write.csv(df.sex.count.chr, file="tables/DE_per_sex_chr_enrichment.csv", row.names=F, quote=F)

# expression divergence
# females
# between species
keep = (species == "Cbr" | species == "Cni") & sex == "F"
dds.briggsae.vs.nigoni.F = DESeqDataSetFromMatrix(countData = counts(dds)[,keep], colData = coldata[keep,], design = ~ species)
dds.briggsae.vs.nigoni.F = DESeq(dds.briggsae.vs.nigoni.F)
dds.briggsae.vs.nigoni.res.F = results(dds.briggsae.vs.nigoni.F)
dds.briggsae.vs.nigoni.sum.F = as.data.frame(dds.briggsae.vs.nigoni.res.F) %>%
    mutate(CniF = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
    select(CniF)
dds.briggsae.vs.nigoni.res.F = cbind(as.data.frame(dds.briggsae.vs.nigoni.res.F),dds.briggsae.vs.nigoni.sum.F)
write.csv(dds.briggsae.vs.nigoni.res.F, file="tables/deseq_res.briggsae.vs.nigoni.sex.F.csv", quote=F)

# hybrids vs briggsae
keep = (species == "Cbr" | species == "HF1") & sex == "F"
dds.briggsae.vs.hybrids.F = DESeqDataSetFromMatrix(countData = counts(dds)[,keep], colData = coldata[keep,], design = ~ species)
dds.briggsae.vs.hybrids.F = DESeq(dds.briggsae.vs.hybrids.F)
dds.briggsae.vs.hybrids.res.F = results(dds.briggsae.vs.hybrids.F)
dds.briggsae.vs.hybrids.sum.F = as.data.frame(dds.briggsae.vs.hybrids.res.F) %>%
    mutate(HF1F = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
    select(HF1F)
dds.briggsae.vs.hybrids.res.F = cbind(as.data.frame(dds.briggsae.vs.hybrids.res.F),dds.briggsae.vs.hybrids.sum.F)
write.csv(dds.briggsae.vs.hybrids.res.F, file="tables/deseq_res.briggsae.vs.hybrids.sex.F.csv", quote=F)

# hybrids vs nigoni
keep = (species == "Cni" | species == "HF1") & sex == "F"
dds.nigoni.vs.hybrids.F = DESeqDataSetFromMatrix(countData = counts(dds)[,keep], colData = coldata[keep,], design = ~ species)
dds.nigoni.vs.hybrids.F = DESeq(dds.nigoni.vs.hybrids.F)
dds.nigoni.vs.hybrids.res.F = results(dds.nigoni.vs.hybrids.F)
dds.nigoni.vs.hybrids.sum.F = as.data.frame(dds.nigoni.vs.hybrids.res.F) %>%
    mutate(HF1F = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
    select(HF1F)
dds.nigoni.vs.hybrids.res.F = cbind(as.data.frame(dds.nigoni.vs.hybrids.res.F),dds.nigoni.vs.hybrids.sum.F)
write.csv(dds.nigoni.vs.hybrids.res.F, file="tables/deseq_res.nigoni.vs.hybrids.sex.F.csv", quote=F)

# males
# between species
keep = (species == "Cbr" | species == "Cni") & sex == "M"
dds.briggsae.vs.nigoni.M = DESeqDataSetFromMatrix(countData = counts(dds)[,keep], colData = coldata[keep,], design = ~ species)
dds.briggsae.vs.nigoni.M = DESeq(dds.briggsae.vs.nigoni.M)
dds.briggsae.vs.nigoni.res.M = results(dds.briggsae.vs.nigoni.M)
dds.briggsae.vs.nigoni.sum.M = as.data.frame(dds.briggsae.vs.nigoni.res.M) %>%
    mutate(CniM = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
    select(CniM)
dds.briggsae.vs.nigoni.res.M = cbind(as.data.frame(dds.briggsae.vs.nigoni.res.M),dds.briggsae.vs.nigoni.sum.M)
write.csv(dds.briggsae.vs.nigoni.res.M, file="tables/deseq_res.briggsae.vs.nigoni.sex.M.csv", quote=F)

# hybrids vs briggsae
keep = (species == "Cbr" | species == "HF1") & sex == "M"
dds.briggsae.vs.hybrids.M = DESeqDataSetFromMatrix(countData = counts(dds)[,keep], colData = coldata[keep,], design = ~ species)
dds.briggsae.vs.hybrids.M = DESeq(dds.briggsae.vs.hybrids.M)
dds.briggsae.vs.hybrids.res.M = results(dds.briggsae.vs.hybrids.M)
dds.briggsae.vs.hybrids.sum.M = as.data.frame(dds.briggsae.vs.hybrids.res.M) %>%
    mutate(HF1M = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
    select(HF1M)
dds.briggsae.vs.hybrids.res.M = cbind(as.data.frame(dds.briggsae.vs.hybrids.res.M),dds.briggsae.vs.hybrids.sum.M)
write.csv(dds.briggsae.vs.hybrids.res.M, file="tables/deseq_res.briggsae.vs.hybrids.sex.M.csv", quote=F)

# hybrids vs nigoni
keep = (species == "Cni" | species == "HF1") & sex == "M"
dds.nigoni.vs.hybrids.M = DESeqDataSetFromMatrix(countData = counts(dds)[,keep], colData = coldata[keep,], design = ~ species)
dds.nigoni.vs.hybrids.M = DESeq(dds.nigoni.vs.hybrids.M)
dds.nigoni.vs.hybrids.res.M = results(dds.nigoni.vs.hybrids.M)
dds.nigoni.vs.hybrids.sum.M = as.data.frame(dds.nigoni.vs.hybrids.res.M) %>%
    mutate(HF1M = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>%
    select(HF1M)
dds.nigoni.vs.hybrids.res.M = cbind(as.data.frame(dds.nigoni.vs.hybrids.res.M),dds.nigoni.vs.hybrids.sum.M)
write.csv(dds.nigoni.vs.hybrids.res.M, file="tables/deseq_res.nigoni.vs.hybrids.sex.M.csv", quote=F)

# consolidate
# expression divergence per sex
df.exprdiv = rbind(data.frame(logFC=dds.briggsae.vs.nigoni.res.M$log2FoldChange,
                              DE=c("sig","no sig","sig")[ factor(dds.briggsae.vs.nigoni.sum.M[,1]) ],
                              species=c("briggsae","nigoni")[ factor(dds.briggsae.vs.nigoni.res.M$log2FoldChange > 0) ],
                              sex="male",
                              p.value=dds.briggsae.vs.nigoni.res.M$padj,
                              genes=rownames(dds.briggsae.vs.nigoni.res.M), stringsAsFactors=F),
                  data.frame(logFC=dds.briggsae.vs.nigoni.res.F$log2FoldChange,
                              DE=c("sig","no sig","sig")[ factor(dds.briggsae.vs.nigoni.sum.F[,1]) ],
                              species=c("briggsae","nigoni")[ factor(dds.briggsae.vs.nigoni.res.F$log2FoldChange > 0) ],
                              sex="female",
                              p.value=dds.briggsae.vs.nigoni.res.F$padj,
                              genes=rownames(dds.briggsae.vs.nigoni.res.F), stringsAsFactors=F))
df.exprdiv$chromosome = sub("\\..*","",df.exprdiv$genes)
df.exprdiv = na.omit(df.exprdiv)
df.exprdiv$species_DE = paste(df.exprdiv$species, df.exprdiv$DE)

write.csv(df.exprdiv, file="tables/df.expr_div.per_sex.csv", row.names=F)

# DE expression enrichment per sex
df.exprdiv2 = df.exprdiv
df.exprdiv2$species_DE[ grep("no sig", df.exprdiv2$species_DE) ] = "no sig"
df.DE.exprdiv.enrich = as.data.frame(df.exprdiv2 %>% group_by(sex, species_DE, chromosome) %>% summarise(count=n()))
df.DE.exprdiv.enrich$p.value = NA
df.DE.exprdiv.enrich$enrichment = NA
x = matrix((df.DE.exprdiv.enrich %>% filter(sex == "male"))$count, ncol=3, nrow=6)
df.DE.exprdiv.enrich$p.value[1:18] = as.vector(enrichment(x))
df.DE.exprdiv.enrich$enrichment[1:18] = as.vector(enrichment(x, odds.ratio=T))
x = matrix((df.DE.exprdiv.enrich %>% filter(sex == "female"))$count, ncol=3, nrow=6)
df.DE.exprdiv.enrich$p.value[19:36] = as.vector(enrichment(x))
df.DE.exprdiv.enrich$enrichment[19:36] = as.vector(enrichment(x, odds.ratio=T))
df.DE.exprdiv.enrich$chr.pseudo = (1:6)[ factor(df.DE.exprdiv.enrich$chromosome) ]
write.csv(df.DE.exprdiv.enrich, file="tables/df.enrichment.exprdiv.DE.sex.csv", row.names=F)

# allele-specific expression cis divergence
# prep data
counts.ase.filter = counts.ase[filter_in, c(grep("Cbr", colnames(counts.ase)), grep("Cni", colnames(counts.ase)))]
keep = species == "Cbr" | species == "Cni"
coldata.ase = coldata[keep,]
coldata.ase$samples = colnames(counts.ase)

# females
keep = coldata.ase$sex == "F"
dds.ase.F = DESeqDataSetFromMatrix(countData = counts.ase.filter[,keep], colData = coldata.ase[keep,], design = ~ species)
dds.ase.F = DESeq(dds.ase.F)
dds.ase.res.F = results(dds.ase.F)
dds.ase.sum.F = as.data.frame(dds.ase.res.F) %>% mutate(CniF = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>% select(CniF)
dds.ase.res.F = cbind(as.data.frame(dds.ase.res.F),dds.ase.sum.F)
write.csv(dds.ase.res.F, file="tables/deseq_res.ase.females.csv", quote=F)


# males
keep = coldata.ase$sex == "M"
auto = rownames(counts.ase.filter)[ grep('^X', rownames(counts.ase.filter), invert=T)]
dds.ase.M = DESeqDataSetFromMatrix(countData = counts.ase.filter[auto,keep], colData = coldata.ase[keep,], design = ~ species)
dds.ase.M = DESeq(dds.ase.M)
dds.ase.res.M = results(dds.ase.M)
dds.ase.sum.M = as.data.frame(dds.ase.res.M) %>% mutate(CniM = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>% select(CniM)
dds.ase.res.M = cbind(as.data.frame(dds.ase.res.M),dds.ase.sum.M)
write.csv(dds.ase.res.M, file="tables/deseq_res.ase.males.csv", quote=F)

# consolidate ASE
df.ase = rbind(data.frame(logFC=dds.ase.res.M$log2FoldChange,
                              DE=c("sig","no sig","sig")[ factor(dds.ase.sum.M[,1]) ],
                              species=c("briggsae","nigoni")[ factor(dds.ase.res.M$log2FoldChange > 0) ],
                              sex="male",
                              p.value=dds.ase.res.M$padj,
                              genes=rownames(dds.ase.res.M), stringsAsFactors=F),
                  data.frame(logFC=dds.ase.res.F$log2FoldChange,
                              DE=c("sig","no sig","sig")[ factor(dds.ase.sum.F[,1]) ],
                              species=c("briggsae","nigoni")[ factor(dds.ase.res.F$log2FoldChange > 0) ],
                              sex="female",
                              p.value=dds.ase.res.F$padj,
                              genes=rownames(dds.ase.res.F), stringsAsFactors=F))
df.ase$chromosome = sub("\\..*","",df.ase$genes)
df.ase = na.omit(df.ase)
df.ase$species_DE = paste(df.ase$species, df.ase$DE)

write.csv(df.ase, file="tables/df.ase.per_sex.csv", row.names=F)

# ASE enrichment per sex
df.ase2 = df.ase
df.ase2$species_DE[ grep("no sig", df.ase2$species_DE) ] = "no sig"
df.DE.ase.enrich = as.data.frame(df.ase2 %>% group_by(sex, species_DE, chromosome) %>% summarise(count=n()))
df.DE.ase.enrich$p.value = NA
df.DE.ase.enrich$enrichment = NA
x = matrix((df.DE.ase.enrich %>% filter(sex == "female"))$count, ncol=3, nrow=6)
df.DE.ase.enrich$p.value[1:18] = as.vector(enrichment(x))
df.DE.ase.enrich$enrichment[1:18] = as.vector(enrichment(x, odds.ratio=T))
x = matrix((df.DE.ase.enrich %>% filter(sex == "male"))$count, ncol=3, nrow=5)
df.DE.ase.enrich$p.value[19:33] = as.vector(enrichment(x))
df.DE.ase.enrich$enrichment[19:33] = as.vector(enrichment(x, odds.ratio=T))
df.DE.ase.enrich$chr.pseudo = (1:6)[ factor(df.DE.ase.enrich$chromosome) ]
write.csv(df.DE.ase.enrich, file="tables/df.enrichment.ase.DE.sex.csv", row.names=F)
