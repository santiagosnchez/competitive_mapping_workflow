# generates files and data for Fig. 1

library(limma)
library(Glimma)
library(edgeR)
library(MASS)
library(ggplot2)
library(cowplot)
library(lemon)
library(tidyr)
library(dplyr)
library(gridExtra)
library(RColorBrewer)

theme_set(theme_cowplot())
source("enrichment_fun.R")

# read raw count data
counts.species = read.csv("../counts/species_species.txt", row.name=1, sep="\t")
counts.ase = read.csv("../counts/ase_counts.txt", row.name=1, sep="\t")
counts.species = counts.species[ rownames(counts.ase), ]

# merge ase counts and data
counts = cbind(counts.species[,1:6], counts.ase[,seq(1,12,2)] + counts.ase[,seq(2,12,2)], counts.species[,7:12])
colnames(counts) = c(paste0("Cbr_F",1:3), paste0("Cbr_M",1:3), paste0("HF1_F",1:3), paste0("HF1_M",1:3), paste0("Cni_F",1:3), paste0("Cni_M",1:3))

# group treatments
groups = gsub("[1-3]$","", colnames(counts))
sex = substr(groups,5,5)
species = substr(groups,1,3)

# log-transformed expression data
e = DGEList(counts)
e$samples$group = factor(groups)
e = calcNormFactors(e, method = "TMM")
keep.exprs = filterByExpr(e, as.factor(groups))
e = e[keep.exprs,, keep.lib.sizes=FALSE]
e = estimateCommonDisp(e)
e = estimateGLMTrendedDisp(e)
e = estimateTagwiseDisp(e)
v = voom(e, plot=T)
write.csv(v$E, "voom_log-transformed_exprdata.csv")

# expression distance data
mds = plotMDS(v$E)
df.mds = data.frame(dim1=mds$x, dim2=mds$y, species, sex)
write.csv(df.mds, "expression_dist.mds_coords.csv")

# per species, F1 sex-biased differential expression
# C. briggsae
e.briggsae = DGEList(e$counts[,1:6])
e.briggsae$samples$group = factor(sex[1:6])
e.briggsae = calcNormFactors(e.briggsae, method = "TMM")
e.briggsae = estimateCommonDisp(e.briggsae)
e.briggsae = estimateGLMTrendedDisp(e.briggsae)
e.briggsae = estimateTagwiseDisp(e.briggsae)
design.sex = model.matrix(~droplevels(factor(sex[1:6])))
colnames(design.sex)[2] = "male"
v.briggsae = voom(e.briggsae, design.sex, plot=F)
vfit.briggsae = lmFit(v.briggsae, design.sex)
vfit.briggsae = eBayes(vfit.briggsae)
dt.briggsae = decideTests(vfit.briggsae, p.value = 0.05)
tt.briggsae = topTable(vfit.briggsae, n = Inf)
tt.briggsae = tt.briggsae[rownames(dt.briggsae),]
write.csv(tt.briggsae, file="tt.briggsae.sex.csv", quote=F)

# C. nigoni
e.nigoni = DGEList(e$counts[,13:18])
e.nigoni$samples$group = factor(sex[1:6])
e.nigoni = calcNormFactors(e.nigoni, method = "TMM")
e.nigoni = estimateCommonDisp(e.nigoni)
e.nigoni = estimateGLMTrendedDisp(e.nigoni)
e.nigoni = estimateTagwiseDisp(e.nigoni)
v.nigoni = voom(e.nigoni, design.sex, plot=F)
vfit.nigoni = lmFit(v.nigoni, design.sex)
vfit.nigoni = eBayes(vfit.nigoni)
dt.nigoni = decideTests(vfit.nigoni, p.value = 0.05)
tt.nigoni = topTable(vfit.nigoni, n = Inf)
tt.nigoni = tt.nigoni[rownames(dt.nigoni),]
write.csv(tt.nigoni, file="tt.nigoni.sex.csv", quote=F)

# F1 hybrid
e.hybrids = DGEList(e$counts[,7:12])
e.hybrids$samples$group = factor(sex[1:6])
e.hybrids = calcNormFactors(e.hybrids, method = "TMM")
e.hybrids = estimateCommonDisp(e.hybrids)
e.hybrids = estimateGLMTrendedDisp(e.hybrids)
e.hybrids = estimateTagwiseDisp(e.hybrids)
v.hybrids = voom(e.hybrids, design.sex, plot=F)
vfit.hybrids = lmFit(v.hybrids, design.sex)
vfit.hybrids = eBayes(vfit.hybrids)
dt.hybrids = decideTests(vfit.hybrids, p.value = 0.05)
tt.hybrids = topTable(vfit.hybrids, n = Inf)
tt.hybrids = tt.hybrids[rownames(dt.hybrids),]
write.csv(tt.hybrids, file="tt.hybrids.sex.csv", quote=F)

# consolidate results
df.sex = data.frame(Cbr=dt.briggsae[,2], F1=dt.hybrids[,2], Cni=dt.nigoni[,2])
df.sex[ df.sex == 0 ] = "sex-neutral"
df.sex[ df.sex == 1 ] = "male"
df.sex[ df.sex == -1 ] = "female"
df.sex$chromosome = sub("\\..*","", rownames(df.sex))
df.sex2 = gather(df.sex, species, sex, -chromosome)

# generate counts and summaries
df.sex.count.spp = as.data.frame(df.sex2 %>% group_by(species, sex) %>% summarize(n=n()) %>% mutate(perc = round((n / sum(n))*100,0)))
df.sex.count.chr = as.data.frame(df.sex2 %>% group_by(species, sex, chromosome) %>% summarize(n=n()))
write.csv(df.sex.count.spp, file="df.counts.perc.DE.sex.csv", row.names=F)

# enrichments
df.sex.count.chr$p.value = NA
df.sex.count.chr$enrichment = NA
x = matrix((df.sex.count.chr %>% filter(species == "Cbr"))$n, ncol=3, nrow=6)
df.sex.count.chr$p.value[1:18] = as.vector(enrichment(x))
df.sex.count.chr$enrichment[1:18] = as.vector(enrichment(x, odds.ratio=T))
x = matrix((df.sex.count.chr %>% filter(species == "Cni"))$n, ncol=3, nrow=6)
df.sex.count.chr$p.value[19:36] = as.vector(enrichment(x))
df.sex.count.chr$enrichment[19:36] = as.vector(enrichment(x, odds.ratio=T))
x = matrix((df.sex.count.chr %>% filter(species == "F1"))$n, ncol=3, nrow=6)
df.sex.count.chr$p.value[37:54] = as.vector(enrichment(x))
df.sex.count.chr$enrichment[37:54] = as.vector(enrichment(x, odds.ratio=T))
df.sex.count.chr$chr.pseudo = (1:6)[ factor(df.sex.count.chr$chromosome) ]
write.csv(df.sex.count.chr, file="DE_per_sex_chr_enrichment.csv", row.names=F)

# expression divergence
# females
e.briggsae.vs.nigoni.female = e$counts[,c(1:3,13:15)]
e.briggsae.vs.nigoni.female = DGEList(e.briggsae.vs.nigoni.female)
e.briggsae.vs.nigoni.female$samples$group = factor(groups[c(1:3,13:15)])
e.briggsae.vs.nigoni.female = calcNormFactors(e.briggsae.vs.nigoni.female, method = "TMM")
e.briggsae.vs.nigoni.female = estimateCommonDisp(e.briggsae.vs.nigoni.female)
e.briggsae.vs.nigoni.female = estimateGLMTrendedDisp(e.briggsae.vs.nigoni.female)
e.briggsae.vs.nigoni.female = estimateTagwiseDisp(e.briggsae.vs.nigoni.female)
design.species = model.matrix(~droplevels(factor(species[c(1:3,13:15)])))
colnames(design.species)[2] = "nigoni"
v.briggsae.vs.nigoni.female = voom(e.briggsae.vs.nigoni.female, design.species, plot=F)
vfit.briggsae.vs.nigoni.female = lmFit(v.briggsae.vs.nigoni.female, design.species)
vfit.briggsae.vs.nigoni.female = eBayes(vfit.briggsae.vs.nigoni.female)
dt.briggsae.vs.nigoni.female = decideTests(vfit.briggsae.vs.nigoni.female, p.value = 0.05)
tt.briggsae.vs.nigoni.female = topTable(vfit.briggsae.vs.nigoni.female, n = Inf)
tt.briggsae.vs.nigoni.female = tt.briggsae.vs.nigoni.female[ rownames(dt.briggsae.vs.nigoni.female), ]
write.csv(tt.briggsae.vs.nigoni.female, file="tt.briggsae.vs.nigoni.females.csv", quote=F)

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
