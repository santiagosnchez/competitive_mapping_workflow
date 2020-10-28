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

# hermaphroditic genes
herm_genes = scan("tables/Thomas_et_al_herm_genes.txt", what=character(), sep="\n")

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

################
### Figure 1 ###
################

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
df.sex$sex_class = NA

df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("sex-neutral","sex-neutral","sex-neutral"))), 1, all) ] = "sex-neutral-conserved"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("female","female","female"))), 1, all) ] = "female-conserved"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("sex-neutral","female","female"))), 1, all) ] = "female-F1-nigoni"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("female","sex-neutral","female"))), 1, all) ] = "female-F1-briggsae"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("sex-neutral","sex-neutral","female"))), 1, all) ] = "female-hybrid"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("female","female","sex-neutral"))), 1, all) ] = "female-sex-neutral"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("male","male","male"))), 1, all) ] = "male-conserved"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("male","sex-neutral","male"))), 1, all) ] = "male-F1-briggsae"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("sex-neutral","male","male"))), 1, all) ] = "male-F1-nigoni"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("sex-neutral","sex-neutral","male"))), 1, all) ] = "male-hybrid"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("male","male","sex-neutral"))), 1, all) ] = "male-sex-neutral"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("female","male","female"))), 1, all) ] = "antagonistic-female-briggsae"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("male","female","female"))), 1, all) ] = "antagonistic-female-nigoni"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("male","female","male"))), 1, all) ] = "antagonistic-male-briggsae"
df.sex$sex_class[ apply(t(apply(df.sex[,1:3], 1, function(x) x  == c("female","male","male"))), 1, all) ] = "antagonistic-male-nigoni"



df.sex2 = pivot_longer(df.sex, cols=c(-chromosome,-gene,-sex_class), names_to="species", values_to="sex")
#df.sex2$sex[ as.character(df.sex2$gene) %in% herm_genes & df.sex2$species == "Cbr" ] = "hermaphrodite"

# generate counts and summaries
df.sex.count.spp = as.data.frame(df.sex2 %>%
            group_by(species, sex) %>%
            summarize(count=n()) %>%
            mutate(perc = round((count / sum(count))*100,0)))
write.csv(df.sex.count.spp, file="tables/df.counts.perc.DE.sex.csv", row.names=F)

# per chromosome
df.sex.count.chr = as.data.frame(df.sex2 %>% filter(sex != "no expression" | sex != "hermaphrodite") %>% group_by(species, sex, chromosome) %>% summarize(count=n()))

# enrichments
df.sex.count.chr$p.value = NA
df.sex.count.chr$enrichment = NA
x = matrix((df.sex.count.chr %>% filter(species == "Cbr"))$count, ncol=4, nrow=6)
df.sex.count.chr$p.value[1:24] = as.vector(enrichment(x))
df.sex.count.chr$enrichment[1:24] = as.vector(enrichment(x, odds.ratio=T))
x = matrix((df.sex.count.chr %>% filter(species == "Cni"))$count, ncol=3, nrow=6)
df.sex.count.chr$p.value[25:42] = as.vector(enrichment(x))
df.sex.count.chr$enrichment[25:42] = as.vector(enrichment(x, odds.ratio=T))
x = matrix((df.sex.count.chr %>% filter(species == "F1"))$count, ncol=3, nrow=6)
df.sex.count.chr$p.value[43:60] = as.vector(enrichment(x))
df.sex.count.chr$enrichment[43:60] = as.vector(enrichment(x, odds.ratio=T))
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
df.exprdiv$sex[ as.character(df.exprdiv$genes) %in% herm_genes & df.exprdiv$sex == "female" ] = "hermaphrodite"

write.csv(df.exprdiv, file="tables/df.expr_div.per_sex.csv", row.names=F)

# DE expression enrichment per sex
df.exprdiv2 = df.exprdiv %>% filter(sex != "hermaphrodite")
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

# get trans effects

getTransEffects = function(x){
	x = as.data.frame(x)
	x$T = rep(c("P","H"), each=6)
	x$S = rep(rep(c("br","ni"), each=3), 2)
	fit = lm(x ~ T/S, data=x)
	fit2 = car::linearHypothesis(fit, c("TH:Sni = TP:Sni"))
	p.value = fit2$`Pr(>F)`[2]
	list(x, fit, fit2, p.value)
}

# data set
keep = coldata.ase$sex == "F"
df.trans.female = cbind(cpm(counts(dds.briggsae.vs.nigoni.F), log=T), cpm(counts.ase.filter[,keep], log=T))
keep = coldata.ase$sex == "M"
df.trans.male = cbind(cpm(counts(dds.briggsae.vs.nigoni.M)[auto,], log=T), cpm(counts.ase.filter[auto,keep], log=T))
# tests Wald
excl = c(12365, 13145, 13146, 13246)
df.trans.female2 = df.trans.female[-excl,]
trans_effects.female = mclapply(split(df.trans.female2, 1:dim(df.trans.female2)[1]), getTransEffects, mc.cores=4)
trans_effects.male = mclapply(split(df.trans.male, 1:dim(df.trans.male)[1]), getTransEffects, mc.cores = 4)
# extract pvalues and FDR analysis
trans_p.values.fdr.female = p.adjust(sapply(trans_effects.female, function(x) x[[4]] ), method="BH")
trans_p.values.fdr.male = p.adjust(sapply(trans_effects.male, function(x) x[[4]] ), method="BH")
tmp = rep(NA, dim(df.trans.female)[1])
names(tmp) = rownames(df.trans.female)
names(trans_p.values.fdr.female) = rownames(df.trans.female2)
tmp[ names(trans_p.values.fdr.female) ] = trans_p.values.fdr.female
trans_p.values.fdr.female = tmp


# consolidate ASE
nas = rep(NA, dim(dds.ase.res.F)[1]-dim(dds.ase.res.M)[1])
df.ase = rbind(data.frame(logFC.sp=dds.briggsae.vs.nigoni.res.M$log2FoldChange,
                          logFC.ase=c(dds.ase.res.M$log2FoldChange,nas),
                          DE=c(c("sig","no sig","sig")[ factor(dds.ase.sum.M[,1]) ],nas),
                          species=c(c("briggsae","nigoni")[ factor(dds.ase.res.M$log2FoldChange > 0) ],nas),
                          sex="male",
                          p.value.sp=dds.briggsae.vs.nigoni.res.M$padj,
                          p.value.ase=c(dds.ase.res.M$padj,nas),
                          genes=rownames(dds.ase.res.F),
                          trans_effect=c(trans_p.values.fdr.male,nas), stringsAsFactors=F),
               data.frame(logFC.sp=dds.briggsae.vs.nigoni.res.F$log2FoldChange,
                          logFC.ase=dds.ase.res.F$log2FoldChange,
                          DE=c("sig","no sig","sig")[ factor(dds.ase.sum.F[,1]) ],
                          species=c("briggsae","nigoni")[ factor(dds.ase.res.F$log2FoldChange > 0) ],
                          sex="female",
                          p.value.sp=dds.briggsae.vs.nigoni.res.F$padj,
                          p.value.ase=dds.ase.res.F$padj,
                          genes=rownames(dds.ase.res.F),
                          trans_effect=trans_p.values.fdr.female, stringsAsFactors=F))
df.ase$chromosome = sub("\\..*","",df.ase$genes)
df.ase = na.omit(df.ase)
df.ase$species_DE = paste(df.ase$species, df.ase$DE)
df.ase$sex = as.character(df.ase$sex)
df.ase$sex[ df.ase$sex == "female" & df.ase$genes %in% herm_genes ] = "hermaphrodite"

write.csv(df.ase, file="tables/df.ase.per_sex.csv", row.names=F)

# ASE enrichment per sex
df.ase2 = df.ase
df.ase2$species_DE[ grep("no sig", df.ase2$species_DE) ] = "no sig"
df.DE.ase.enrich = as.data.frame(df.ase2 %>% filter(sex != "hermaphrodite") %>% group_by(sex, species_DE, chromosome) %>% summarise(count=n()))
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

################
### Figure 2 ###
################


##############################
### expression inheritance ###
##############################

dt.briggsae.vs.nigoni.female = read.csv("tables/deseq_res.briggsae.vs.nigoni.sex.F.csv", row.names=1)
dt.briggsae.vs.nigoni.male = read.csv("tables/deseq_res.briggsae.vs.nigoni.sex.M.csv", row.names=1)
dt.briggsae.vs.hybrids.female = read.csv("tables/deseq_res.briggsae.vs.hybrids.sex.F.csv", row.names=1)
dt.briggsae.vs.hybrids.male = read.csv("tables/deseq_res.briggsae.vs.hybrids.sex.M.csv", row.names=1)
dt.nigoni.vs.hybrids.female = read.csv("tables/deseq_res.nigoni.vs.hybrids.sex.F.csv", row.names=1)
dt.nigoni.vs.hybrids.male = read.csv("tables/deseq_res.nigoni.vs.hybrids.sex.M.csv", row.names=1)

# consolidate
# expression divergence per sex between parents and F1

models.inheritance.male = paste(dt.briggsae.vs.nigoni.male[,"CniM"], dt.briggsae.vs.hybrids.male[,"HF1M"], dt.nigoni.vs.hybrids.male[,"HF1M"])
models.inheritance.female = paste(dt.briggsae.vs.nigoni.female[,"CniF"], dt.briggsae.vs.hybrids.female[,"HF1F"], dt.nigoni.vs.hybrids.female[,"HF1F"])

# classify genes based on F1 expression relative to parent species
# order "Cbr.vs.Cni Cbr.vs.F1 Cni.vs.F1"
# 0 = no DE
# -1 = C. briggsae (in Cbr vs Cni or Cbr vs F1); C. nigoni (in Cni vs F1)
# 1 = C. nigoni (in Cbr cs Cni); F1 (in Cni vs F1 or Cni vs F1)
get_inheritance_class = function(x){
	class_inheritance = vector("character", length=length(x))
	class_inheritance[ x == "0 0 0" ] = "no change"
	class_inheritance[ x == "0 1 1" | x == "-1 1 1" | x == "1 1 1" ] = "overdominant"
	class_inheritance[ x == "0 -1 -1" | x == "-1 -1 -1" | x == "1 -1 -1" ] = "underdominant"
	class_inheritance[ x == "-1 0 1" | x == "1 0 -1" ] = "C. briggsae dominant"
	class_inheritance[ x == "1 1 0" | x == "-1 -1 0" ] = "C. nigoni dominant"
	class_inheritance[ x == "1 1 -1" | x == "-1 -1 1" ] = "additive"
	class_inheritance[ class_inheritance == "" ] = "ambiguous"
	return(class_inheritance)
}

class_inheritance.female = get_inheritance_class(models.inheritance.female)
class_inheritance.male = get_inheritance_class(models.inheritance.male)
class.levels = c("no change","ambiguous","additive","C. briggsae dominant", "C. nigoni dominant","overdominant","underdominant")

# hermaphroditic genes
herm_genes = scan("tables/Thomas_et_al_herm_genes.txt", what=character(), sep="\n")

df.inherit = data.frame(female=class_inheritance.female, male=class_inheritance.male, gene=rownames(dt.briggsae.vs.nigoni.female))
df.inherit2 = pivot_longer(df.inherit, -gene, names_to="sex", values_to="class") %>% arrange(sex)
df.inherit2[ as.character(df.inherit2$gene) %in% herm_genes & as.character(df.inherit2$sex) == "female", "sex" ] = "(F) hermaphrodite"
df.inherit2[ as.character(df.inherit2$gene) %in% herm_genes & as.character(df.inherit2$sex) == "male", "sex" ] = "(M) hermaphrodite"
df.inherit2$class = factor(df.inherit2$class, levels=class.levels)
df.inherit2$chromosome = sub("\\..*","", df.inherit2$gene)
df.inherit2$logFC.F1.vs.briggsae = c(dt.briggsae.vs.hybrids.female[,"log2FoldChange"],dt.briggsae.vs.hybrids.male[,"log2FoldChange"])
df.inherit2$logFC.F1.vs.nigoni = c(dt.nigoni.vs.hybrids.female[,"log2FoldChange"],dt.nigoni.vs.hybrids.male[,"log2FoldChange"])
df.inherit2$pvalue.logFC.F1.vs.briggsae = c(dt.briggsae.vs.hybrids.female[,"padj"],dt.briggsae.vs.hybrids.male[,"padj"])
df.inherit2$pvalue.logFC.F1.vs.nigoni = c(dt.nigoni.vs.hybrids.female[,"padj"],dt.nigoni.vs.hybrids.male[,"padj"])

# expression distance F1 from zero
df.inherit2$F1_dist_from_zero = pointDistance(as.matrix(df.inherit2[,c("logFC.F1.vs.briggsae","logFC.F1.vs.nigoni")]),
											matrix(0, ncol=2, nrow=dim(df.inherit2)[1]),
											lonlat=F, allpairs=F)

df.inherit2 = cbind(df.inherit2, bn.bf.nf=c(models.inheritance.female, models.inheritance.male))

write.csv(df.inherit2, file="tables/expr_inheritance_logFC_F1_parents.csv", row.names=F)

# enrichment and counts per chromosome
df.inherit3 = (df.inherit2 %>% group_by(sex, class, chromosome) %>% count())
x = matrix((df.inherit3 %>% filter(sex == "(F) hermaphrodite"))$n, ncol=7, nrow=6)
df.inherit3[1:42,"enrichment"] = as.vector(enrichment(x, odds.ratio=T))
df.inherit3[1:42,"p.value"] = as.vector(enrichment(x, odds.ratio=F))
x = matrix((df.inherit3 %>% filter(sex == "(M) hermaphrodite"))$n, ncol=7, nrow=6)
df.inherit3[43:84,"enrichment"] = as.vector(enrichment(x, odds.ratio=T))
df.inherit3[43:84,"p.value"] = as.vector(enrichment(x, odds.ratio=F))
x = matrix((df.inherit3 %>% filter(sex == "female"))$n, ncol=7, nrow=6)
df.inherit3[85:126,"enrichment"] = as.vector(enrichment(x, odds.ratio=T))
df.inherit3[85:126,"p.value"] = as.vector(enrichment(x, odds.ratio=F))
x = matrix((df.inherit3 %>% filter(sex == "male"))$n, ncol=7, nrow=6)
df.inherit3[127:168,"enrichment"] = as.vector(enrichment(x, odds.ratio=T))
df.inherit3[127:168,"p.value"] = as.vector(enrichment(x, odds.ratio=F))
df.inherit3$chr.pseudo = (1:6)[ factor(df.inherit3$chromosome)]
df.inherit3$sig = ""
df.inherit3$sig[ df.inherit3$p.value < 0.05 & abs(log2(df.inherit3$enrichment)) > 0.5 ] = "*"
write.csv(df.inherit3, file="tables/expr_inheritance_chr_enrich.csv", row.names=F)

##################################
### allele-specific expression ###
##################################

# classify allele specific expression
# x = "parents hybrids_ASE trans_effects"
get_regulation_type = function(x){
    regulation_type = rep(NA, dim(x)[1])
    regulation_type[ x[,"p.value.sp"] > 0.05 & x[,"p.value.ase"] > 0.05 ] = "conserved"
    regulation_type[ is.na(regulation_type) & x[,"trans_effect"] > 0.05 & x[,"p.value.ase"] < 0.05 & x[,"p.value.sp"] < 0.05 ] = "cis-only"
    regulation_type[ is.na(regulation_type) & x[,"trans_effect"] < 0.05 & x[,"p.value.ase"] > 0.05 & x[,"p.value.sp"] < 0.05 ] = "trans-only"
    regulation_type[ is.na(regulation_type) & x[,"trans_effect"] < 0.05 & x[,"p.value.ase"] < 0.05 & x[,"p.value.sp"] > 0.05 ] = "cis-trans (compensatory)"
    regulation_type[ is.na(regulation_type) & x[,"trans_effect"] < 0.05 & ((x[,"logFC.ase"] > 0 & x[,"logFC.sp"] < 0) | (x[,"logFC.ase"] < 0 & x[,"logFC.sp"] > 0)) ] = "cis x trans (compensatory)"
    regulation_type[ is.na(regulation_type) & x[,"trans_effect"] < 0.05 & ((x[,"logFC.ase"] > 0 & x[,"logFC.sp"] > 0) | (x[,"logFC.ase"] < 0 & x[,"logFC.sp"] < 0)) ] = "cis + trans (enhancing)"
    regulation_type[ is.na(regulation_type) ] = "ambiguous"
	return(regulation_type)
}
# read data
ase.div = read.csv("tables/df.ase.per_sex.csv", stringsAsFactors=F) # allele-specific expression (logFC)
ase.div$type = get_regulation_type(ase.div)
ase.div[ as.character(ase.div$genes) %in% herm_genes & as.character(ase.div$sex) == "female", "sex" ] = "(F) hermaphrodite"
ase.div[ as.character(ase.div$genes) %in% herm_genes & as.character(ase.div$sex) == "male", "sex" ] = "(M) hermaphrodite"
write.csv(ase.div, file="tables/df.cis_trans.regulation_type.auto.csv", row.names=F)

# counts and enrichments

ase.div.enrich = as.data.frame(ase.div %>% group_by(sex, type, chromosome) %>% count())
x = matrix((ase.div.enrich %>% filter(sex == "(F) hermaphrodite"))$n, ncol=7, nrow=6)
ase.div.enrich$enrichment = NA
ase.div.enrich$enrichment[1:42] = as.vector(enrichment(x, odds.ratio=T))
ase.div.enrich$p.value = NA
ase.div.enrich$p.value[1:42] = as.vector(enrichment(x, odds.ratio=F))
x = matrix((ase.div.enrich %>% filter(sex == "(M) hermaphrodite"))$n, ncol=7, nrow=5)
ase.div.enrich$enrichment[43:77] = as.vector(enrichment(x, odds.ratio=T))
ase.div.enrich$p.value[43:77] = as.vector(enrichment(x, odds.ratio=F))
x = matrix((ase.div.enrich %>% filter(sex == "female"))$n, ncol=7, nrow=6)
ase.div.enrich$enrichment[78:119] = as.vector(enrichment(x, odds.ratio=T))
ase.div.enrich$p.value[78:119] = as.vector(enrichment(x, odds.ratio=F))
ase.div.enrich$chr.pseudo = (1:6)[ factor(ase.div.enrich$chromosome)]
x = matrix((ase.div.enrich %>% filter(sex == "male"))$n, ncol=7, nrow=5)
ase.div.enrich$enrichment[120:154] = as.vector(enrichment(x, odds.ratio=T))
ase.div.enrich$p.value[120:154] = as.vector(enrichment(x, odds.ratio=F))
ase.div.enrich$chr.pseudo = (1:6)[ factor(ase.div.enrich$chromosome)]

ase.div.enrich$sig = ""
ase.div.enrich$sig[ ase.div.enrich$p.value < 0.05 & abs(log2(ase.div.enrich$enrichment)) > 0.5 ] = "*"
write.csv(ase.div.enrich, file="tables/df.cis_trans.counts.chr_enrich.csv", row.names=F)

# merge expr inheritance with regulation type

tmp1 = df.inherit2
tmp2 = ase.div
rownames(tmp1) = paste(tmp1$gene,tmp1$sex)
rownames(tmp2) = paste(tmp2$genes,tmp2$sex)
df.inherit_and_ase = merge(tmp1, tmp2, by=0, all=T)
df.inherit_and_ase = df.inherit_and_ase %>% dplyr::select(-Row.names, -sex.y, -chromosome.y, -genes) %>% rename(sex="sex.x", chromosome="chromosome.x")
write.csv(df.inherit_and_ase, file="tables/df.cis_trans.inherit.merge.csv", row.names=F)

################
### Figure 3 ###
################

# matrix fig 3

class.vs.type.mat = na.omit(df.inherit_and_ase) %>%
    filter(!(class == "no change" | class == "ambiguous" | type == "conserved" | type == "ambiguous")) %>%
    type.convert() %>%
    group_by(sex, class, type, chromosome, .drop=F) %>%
    count()
write.csv(class.vs.type.mat, file="tables/df.cis_trans.inherit.mat.csv", row.names=F)



# sex by species interactions
keep = species == "Cbr" | species == "Cni"
dds.briggsae.vs.nigoni.sex = DESeqDataSetFromMatrix(countData = counts(dds)[,keep], colData = coldata[keep,], design = ~ species * sex)
dds.briggsae.vs.nigoni.sex = DESeq(dds.briggsae.vs.nigoni.sex)
condition_names = resultsNames(dds.briggsae.vs.nigoni.sex)[2:4]
dds.briggsae.vs.nigoni.sex.res.species = results(dds.briggsae.vs.nigoni.sex, name=condition_names[1])
dds.briggsae.vs.nigoni.sex.res.sex = results(dds.briggsae.vs.nigoni.sex, name=condition_names[2])
dds.briggsae.vs.nigoni.sex.res.interaction = results(dds.briggsae.vs.nigoni.sex, name=condition_names[3])
df.species_sex = data.frame(logFC_species=dds.briggsae.vs.nigoni.sex.res.species$log2FoldChange, pvalue_species=dds.briggsae.vs.nigoni.sex.res.species$padj,
                            logFC_sex=dds.briggsae.vs.nigoni.sex.res.sex$log2FoldChange, pvalue_sex=dds.briggsae.vs.nigoni.sex.res.sex$padj,
                            logFC_interaction=dds.briggsae.vs.nigoni.sex.res.interaction$log2FoldChange, pvalue_interaction=dds.briggsae.vs.nigoni.sex.res.interaction$padj)
df.species_sex.sum = df.species_sex %>%
    mutate(species = ifelse(pvalue_species > 0.05, 0, ifelse(logFC_species > 0, 1, -1)),
           sex = ifelse(pvalue_sex > 0.05, 0, ifelse(logFC_sex > 0, 1, -1)),
           interaction = ifelse(pvalue_interaction > 0.05, 0, ifelse(logFC_interaction > 0, 1, -1))) %>%
           select(species, sex, interaction)
    select(CniM)

# paste together

models.species_sex = paste(df.species_sex.sum[,2], df.species_sex.sum[,3], df.species_sex.sum[,4])

# classify genes based on F1 expression relative to parent species
# first column:
# 0 = no DE
# -1 = C. briggsae
# 1 = C. nigoni
# second column:
# 0 = no DE
# 1 = males
# -1 = females
# third column:
# no interaction = 0
# interaction = 1 = -1
get_species_sex_class = function(x){
	class_species_sex = vector("character", length=length(x))
	class_species_sex[ x == "0 0 0" ] = "conserved sex-neutral no_interaction"
	class_species_sex[ x == "0 1 0" ] = "conserved male-biased no_interaction"
	class_species_sex[ x == "0 -1 0" ] = "conserved female-biased no_interaction"
	class_species_sex[ x == "1 0 0" ] = "nigoni sex-neutral no_interaction"
	class_species_sex[ x == "-1 0 0" ] = "briggsae sex-neutral no_interaction"
	class_species_sex[ x == "1 1 0" ] = "nigoni male-biased no_interaction"
	class_species_sex[ x == "-1 1 0" ] = "briggsae male-biased no_interaction"
	class_species_sex[ x == "1 -1 0" ] = "nigoni female-biased no_interaction"
	class_species_sex[ x == "-1 -1 0" ] = "briggsae female-biased no_interaction"
	class_species_sex[ x == "0 0 -1" | x == "0 0 1" ] = "conserved sex-neutral interaction"
	class_species_sex[ x == "-1 -1 -1" | x == "-1 -1 1" ] = "briggsae female-biased interaction"
	class_species_sex[ x == "-1 1 -1" | x == "-1 1 1" ] = "briggsae male-biased interaction"
	class_species_sex[ x == "1 -1 -1" | x == "1 -1 1" ] = "nigoni female-biased interaction"
	class_species_sex[ x == "1 1 -1" | x == "1 1 1" ] = "nigoni male-biased interaction"
	class_species_sex[ class_species_sex == "" ] = "ambiguous"
	return(class_species_sex)
}

species_sex_class = get_species_sex_class(models.species_sex)
species_sex_class[ species_sex_class == "conserved sex-neutral interaction" ] = "ambiguous"

models.levels = sort(unique(species_sex_class))[-1][c(8,5,13,7,4,12,3,11,6,2,10,1,9)]

df.species_sex_class = as.data.frame(t(as.data.frame(strsplit(species_sex_class, " "))), stringsAsFactor=F)
colnames(df.species_sex_class) = c("species","sex","interaction")
rownames(df.species_sex_class) = rownames(dds.briggsae.vs.nigoni.sex.res.interaction)

df.species_sex = cbind(df.species_sex, models=species_sex_class, df.species_sex_class)
df.species_sex = df.species_sex[ df.species_sex[,"models"] != "ambiguous", ]
df.species_sex$models = factor(df.species_sex$models, models.levels)

dnds = read.csv("tables/dnds_propcons_position_domain.csv", row.names=1)
x = gsub("\\.t[1-9]","", rownames(dds.briggsae.vs.nigoni.sex.res.species))
dnds = dnds[x,]
rownames(dnds) = rownames(dds.briggsae.vs.nigoni.sex.res.species)
df.species_sex$Ka = dnds[rownames(df.species_sex),"dn"]
df.species_sex$Ks = dnds[rownames(df.species_sex),"ds"]
df.species_sex$Ks_ENCc = dnds[rownames(df.species_sex),"ds_ENCc"]
df.species_sex$domain = dnds[rownames(df.species_sex),"domain"]
models2 = c("C-N-N","B-N-N","N-N-N","C-M-N","B-M-N","N-M-N","B-M-I","N-M-I","C-F-N","B-F-N","N-F-N","B-F-I","N-F-I")[df.species_sex$models]
df.species_sex$models2 = models2

df.aes.inherit.merged = read.csv("tables/df.cis_trans.inherit.merge.csv")

df.species_sex$class_male = df.species_sex$type_male = df.species_sex$sex2 = NA
tmp = filter(df.aes.inherit.merged, sex == "male" | sex == "(M) hermaphrodite")
tmp = tmp[ tmp$gene %in% rownames(df.species_sex),]
df.species_sex[ as.character(tmp$gene), "type_male"] = as.character(tmp$type)
df.species_sex[ as.character(tmp$gene), "class_male"] = as.character(tmp$class)
df.species_sex[ as.character(tmp$gene), "sex2"] = as.character(tmp$sex)
df.species_sex$sex2 = gsub("(M) hermaphrodite", "hermaphrodite", df.species_sex$sex2)
df.species_sex$sex2 = gsub("male", "male-female", df.species_sex$sex2)

df.species_sex$class_female = df.species_sex$type_female = NA
tmp = filter(df.aes.inherit.merged, sex == "female" | sex == "(F) hermaphrodite")
tmp = tmp[ tmp$gene %in% rownames(df.species_sex),]
df.species_sex[ as.character(tmp$gene), "type_female"] = as.character(tmp$type)
df.species_sex[ as.character(tmp$gene), "class_female"] = as.character(tmp$class)
df.species_sex$chromosome = gsub("\\..*","", rownames(df.species_sex))

write.csv(df.species_sex, file="tables/df.species_sex.csv")

# heatmap of species * sex models
df.models = unique(df.species_sex[, 8:10])
df.models$models = unique(df.species_sex[, 7])
df.models$models2 = unique(df.species_sex[, 14])
df.models$counts = table(df.species_sex$models)[ df.models$models ]
rownames(df.models) = 1:dim(df.models)[1]
df.models2 = pivot_longer(df.models, c(-models, -models2, -counts), names_to="key", values_to="value")
df.models2$sex = sapply(strsplit(as.character(df.models2$models), " "), function(x) x[2] )

write.csv(df.models2, file="tables/species_sex_models_heatmap.csv", row.names=F)

# counts and propotions
df.species_sex_class.expr_class_type = data.frame(sex=rep(c("female","male"), each=dim(df.species_sex)[1]),
												  class=c(as.character(df.species_sex$class_female), as.character(df.species_sex$class_male)),
												  type=c(as.character(df.species_sex$type_female), as.character(df.species_sex$type_male)),
												  models=rep(as.character(df.species_sex$models), 2),
                                                  models2=rep(as.character(df.species_sex$models2), 2),
												  gene=rep(rownames(df.species_sex), 2))
df.species_sex_class.expr_class_type$chromosome = sub("\\..*","", df.species_sex_class.expr_class_type$gene)

df.species_sex_class.expr_class.counts = as.data.frame(df.species_sex_class.expr_class_type %>% group_by(sex, models, models2, class) %>% summarize(n = n()) %>% na.omit() %>% mutate(prop=n/sum(n)))
df.species_sex_class.expr_type.counts = as.data.frame(df.species_sex_class.expr_class_type %>% group_by(sex, models, models2, type) %>% summarize(n = n()) %>% na.omit() %>% mutate(prop=n/sum(n)))

write.csv(df.species_sex_class.expr_type.counts, file="tables/df.species_sex_class.expr_type.counts.csv", row.names=F)
write.csv(df.species_sex_class.expr_class.counts, file="tables/df.species_sex_class.expr_class.counts.csv", row.names=F)

# scaled data for reaction norms

dds.briggsae.vs.nigoni.sex.rlog = rlog(counts(dds.briggsae.vs.nigoni.sex), blind=F)
dds.briggsae.vs.nigoni.sex.rlog = data.frame(CbrF=rowMeans(dds.briggsae.vs.nigoni.sex.rlog[,1:3]),
                                            CbrM=rowMeans(dds.briggsae.vs.nigoni.sex.rlog[,4:6]),
                                            CniF=rowMeans(dds.briggsae.vs.nigoni.sex.rlog[,7:9]),
                                            CniM=rowMeans(dds.briggsae.vs.nigoni.sex.rlog[,10:12]))
scaledata = t(scale(t(dds.briggsae.vs.nigoni.sex.rlog)))
rownames(scaledata) = rownames(dds.briggsae.vs.nigoni.sex.res.species)

# vector with models
models = df.species_sex$models2
names(models) = rownames(df.species_sex)

# subsample
scaledata = as.data.frame(scaledata[names(models),])

# add data
scaledata$models = models
scaledata$gene = names(models)

# melt dataframe
df.models = pivot_longer(scaledata, c(-models, -gene), names_to="samples", values_to="expression")
df.models$species = substr(df.models$samples, 1, 3)
df.models$species = c("C. briggsae","C. nigoni")[ factor(df.models$species) ]
df.models$sex = substr(df.models$samples, 4, 4)
df.models$sex = c("female","male")[ factor(df.models$sex) ]
df.models$sex[ grep("sex-neutral", df.models$models) ] = "sex-neutral"
write.csv(df.models, file="tables/df.full_species_sex_reaction_norms.csv", row.names=F)

# get centroids
df.models.cent = as.data.frame(df.models %>% group_by(models, samples) %>% summarize(expr=mean(expression), upper=quantile(expression, .95), lower=quantile(expression, .05)) %>% select(models, samples, expression = expr, upper, lower))
df.models.cent$species = substr(df.models.cent$samples, 1, 3)
df.models.cent$sex = substr(df.models.cent$samples, 4, 4)
df.models.cent$gene = "centroid"
df.models.cent$species = c("C. briggsae","C. nigoni")[ factor(df.models.cent$species) ]
df.models.cent$sex = c("female","male")[ factor(df.models.cent$sex) ]
write.csv(df.models.cent, file="tables/df.centroides_and_CI.species_sex_reaction_norms.csv", row.names=F)


# cis divergence and molecular evolution
dnds = read.csv("tables/dnds_propcons_position_domain.csv", row.names=1)
tmp1 = filter(cis_trans.data, sex == "male" | sex == "(M) hermaphrodite")
tmp1 = filter(cis_trans.data, sex == "female" | sex == "(F) hermaphrodite")
tmp1 = filter(cis_trans.data, sex == "male" | sex == "(M) hermaphrodite")
tmp2 = filter(cis_trans.data, sex == "female" | sex == "(F) hermaphrodite")
tmp1.1 = gsub("\\.t[1-9]","", tmp1$genes)
tmp2.1 = gsub("\\.t[1-9]","", tmp2$genes)
dnds.1 = dnds[tmp1.1,]
dnds.2 = dnds[tmp2.1,]
cis_trans.data_molevol = cbind(cis_trans.data, rbind(dnds.1, dnds.2))
write.csv(cis_trans.data_molevol, file="tables/df.cis_divergence.molevol.csv", row.names=F)

# trans-only = "yellow"
# cis-only = "cyan blue"
# cis+trans = "green"
# cis x trans = "pale pink"
# cis-trans = "dark pink"

# spermatogenesis genes

genes = rownames(df.species_sex)
spermatogenesis_genes = scan("tables/c_briggsae_sperm_genes_Ma_et_al_2014.txt", what=character(), sep="\n")
tmp1 = sapply(strsplit(genes, "__"), function(x) sub("\\.t.*", "", x[1]) )
df.sperm_herm = df.species_sex %>% dplyr::select(models2)
df.sperm_herm$herm = df.sperm_herm$sperm = "not"
df.sperm_herm[ genes[tmp1 %in% spermatogenesis_genes], "sperm" ] = "spermatogenesis"
df.sperm_herm[ genes[genes %in% herm_genes], "herm" ] = "hermaphrodite"
cat(genes[tmp1 %in% spermatogenesis_genes], file="tables/spermatogenesis_orthologs_Ma_et_al.txt", sep="\n")

# enrichment
tmp1 = df.sperm_herm %>%
    group_by(models2, sperm, .drop=F) %>%
    summarise(n=n()) %>%
    filter(sperm=="spermatogenesis") %>% rename(group="sperm")
tmp1$enrichment = enrichment(matrix((df.sperm_herm %>% group_by(models2, sperm) %>% count())$n, ncol=2, nrow=13, byrow=T), odds.ratio=T)[,2]
tmp1$p.value = enrichment(matrix((df.sperm_herm %>% group_by(models2, sperm) %>% count())$n, ncol=2, nrow=13, byrow=T), odds.ratio=F)[,2]

tmp2 = df.sperm_herm %>%
    mutate(herm = factor(herm)) %>%
    group_by(models2, herm, .drop=FALSE) %>%
    count() %>%
    filter(herm=="hermaphrodite") %>% rename(group="herm")
tmp2$enrichment = enrichment(matrix((df.sperm_herm %>% mutate(herm = factor(herm)) %>% group_by(models2, herm, .drop=FALSE) %>% summarize(n = n()) %>% mutate(n = n+1))$n, ncol=2, nrow=13, byrow=T), odds.ratio=T)[,1]
tmp2$p.value = enrichment(matrix((df.sperm_herm %>% mutate(herm = factor(herm)) %>% group_by(models2, herm, .drop=FALSE) %>% summarize(n = n()) %>% mutate(n = n+1))$n, ncol=2, nrow=13, byrow=T), odds.ratio=F)[,1]

df.sperm_herm.enrich = bind_rows(tmp1, tmp2)
write.csv(df.sperm_herm.enrich, file="tables/df.sperm_herm.enrich.csv", row.names=F)
cat("# enrichment analysis in \"hermaphrodite\" group is based on +1 pseudocounts\n", file="tables/df.sperm_herm.enrich.csv", append=T)
