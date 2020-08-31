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
keep = (species == "Cbr" | species == "Cni") & sex == "F"
dds.briggsae.vs.nigoni = DESeqDataSetFromMatrix(countData = counts(dds)[,keep], colData = coldata[keep,], design = ~ species)
dds.briggsae.vs.nigoni = DESeq(dds.briggsae.vs.nigoni)
dds.briggsae.vs.nigoni.res = results(dds.briggsae.vs.nigoni)
write.csv(as.data.frame(dds.briggsae.vs.nigoni.res), file="tables/deseq_res.briggsae.vs.nigoni.sex.csv", quote=F)
dds.briggsae.vs.nigoni.sum = as.data.frame(dds.briggsae.vs.nigoni.res) %>% mutate(CniF = ifelse(padj > 0.05, 0, ifelse(log2FoldChange > 0, 1, -1))) %>% select(Cni)

# males
e.briggsae.vs.nigoni.male = e$counts[,c(4:6,16:18)]
e.briggsae.vs.nigoni.male = DGEList(e.briggsae.vs.nigoni.male)
e.briggsae.vs.nigoni.male$samples$group = factor(groups[c(4:6,16:18)])
e.briggsae.vs.nigoni.male = calcNormFactors(e.briggsae.vs.nigoni.male, method = "TMM")
e.briggsae.vs.nigoni.male = estimateCommonDisp(e.briggsae.vs.nigoni.male)
e.briggsae.vs.nigoni.male = estimateGLMTrendedDisp(e.briggsae.vs.nigoni.male)
e.briggsae.vs.nigoni.male = estimateTagwiseDisp(e.briggsae.vs.nigoni.male)
v.briggsae.vs.nigoni.male = voom(e.briggsae.vs.nigoni.male, design.species, plot=F)
vfit.briggsae.vs.nigoni.male = lmFit(v.briggsae.vs.nigoni.male, design.species)
vfit.briggsae.vs.nigoni.male = eBayes(vfit.briggsae.vs.nigoni.male)
dt.briggsae.vs.nigoni.male = decideTests(vfit.briggsae.vs.nigoni.male, p.value = 0.05)
tt.briggsae.vs.nigoni.male = topTable(vfit.briggsae.vs.nigoni.male, n = Inf)
tt.briggsae.vs.nigoni.male = tt.briggsae.vs.nigoni.male[ rownames(dt.briggsae.vs.nigoni.male), ]
write.csv(tt.briggsae.vs.nigoni.male, file="tt.briggsae.vs.nigoni.males.csv", quote=F)

# consolidate
# expression divergence per sex
df.exprdiv = rbind(data.frame(logFC=tt.briggsae.vs.nigoni.male[,1],
                              DE=c("sig","no sig","sig")[ factor(dt.briggsae.vs.nigoni.male[,2]) ],
                              species=c("nigoni","briggsae")[ factor(tt.briggsae.vs.nigoni.male[,1] < 0) ],
                              sex="male",
                              p.value=tt.briggsae.vs.nigoni.male[,"adj.P.Val"],
                              genes=rownames(tt.briggsae.vs.nigoni.male)),
                  data.frame(logFC=tt.briggsae.vs.nigoni.female[,1],
                              DE=c("sig","no sig","sig")[ factor(dt.briggsae.vs.nigoni.female[,2]) ],
                              species=c("nigoni","briggsae")[ factor(tt.briggsae.vs.nigoni.female[,1] < 0) ],
                              sex="female",
                              p.value=tt.briggsae.vs.nigoni.female[,"adj.P.Val"],
                              genes=rownames(tt.briggsae.vs.nigoni.female)))
df.exprdiv$chromosome = sub("\\..*","",df.exprdiv$genes)
df.exprdiv$species_DE = paste(df.exprdiv$species, df.exprdiv$DE)
write.csv(df.exprdiv, file="df.expr_div.per_sex.csv", row.names=F)

# DE expression enrichment per sex
df.exprdiv2 = df.exprdiv
df.exprdiv2$species_DE[ grep("no sig", df.exprdiv2$species_DE) ] = "no sig"
df.DE.enrich = as.data.frame(df.exprdiv2 %>% group_by(sex, species_DE, chromosome) %>% count())
df.DE.enrich$p.value = NA
df.DE.enrich$enrichment = NA
x = matrix((df.DE.enrich %>% filter(sex == "male"))$n, ncol=3, nrow=6)
df.DE.enrich$p.value[1:18] = as.vector(enrichment(x))
df.DE.enrich$enrichment[1:18] = as.vector(enrichment(x, odds.ratio=T))
x = matrix((df.DE.enrich %>% filter(sex == "female"))$n, ncol=3, nrow=6)
df.DE.enrich$p.value[19:36] = as.vector(enrichment(x))
df.DE.enrich$enrichment[19:36] = as.vector(enrichment(x, odds.ratio=T))
df.DE.enrich$chr.pseudo = (1:6)[ factor(df.DE.enrich$chromosome) ]
write.csv(df.DE.enrich, file="df.enrichment.DE.sex.csv", row.names=F)
