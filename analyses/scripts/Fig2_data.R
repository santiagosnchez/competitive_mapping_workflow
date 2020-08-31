# generates data/analyses for Fig. 2

library(limma)
library(Glimma)
library(edgeR)
library(MASS)
library(parallel)
library(raster)
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

####################################################
### expression divergence between species and F1 ###
####################################################

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
#write.csv(tt.briggsae.vs.nigoni.female, file="tt.briggsae.vs.nigoni.females.csv", quote=F)

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
#write.csv(tt.briggsae.vs.nigoni.male, file="tt.briggsae.vs.nigoni.males.csv", quote=F)

## C. briggsae vs F1
# females
e.briggsae.vs.hybrids.female = e$counts[,c(1:3,7:9)]
e.briggsae.vs.hybrids.female = DGEList(e.briggsae.vs.hybrids.female)
e.briggsae.vs.hybrids.female$samples$group = factor(sub("[1-3]$","",rownames(e.briggsae.vs.hybrids.female$samples)))
e.briggsae.vs.hybrids.female = calcNormFactors(e.briggsae.vs.hybrids.female, method = "TMM")
e.briggsae.vs.hybrids.female = estimateCommonDisp(e.briggsae.vs.hybrids.female)
e.briggsae.vs.hybrids.female = estimateGLMTrendedDisp(e.briggsae.vs.hybrids.female)
e.briggsae.vs.hybrids.female = estimateTagwiseDisp(e.briggsae.vs.hybrids.female)
design.briggsae.vs.hybrids.female = model.matrix(~droplevels(factor(species[c(1:3,7:9)])))
colnames(design.briggsae.vs.hybrids.female)[2] = "F1"
v.briggsae.vs.hybrids.female = voom(e.briggsae.vs.hybrids.female, design.briggsae.vs.hybrids.female, plot=F)
vfit.briggsae.vs.hybrids.female = lmFit(v.briggsae.vs.hybrids.female, design.briggsae.vs.hybrids.female)
vfit.briggsae.vs.hybrids.female = eBayes(vfit.briggsae.vs.hybrids.female)
dt.briggsae.vs.hybrids.female = decideTests(vfit.briggsae.vs.hybrids.female, p.value = 0.05)
tt.briggsae.vs.hybrids.female = topTable(vfit.briggsae.vs.hybrids.female, n = Inf)
tt.briggsae.vs.hybrids.female = tt.briggsae.vs.hybrids.female[ rownames(dt.briggsae.vs.hybrids.female), ]
write.csv(tt.briggsae.vs.hybrids.female, file="tt.briggsae.vs.hybrids.females.csv", quote=F)

# males
e.briggsae.vs.hybrids.male = e$counts[,c(4:6,10:12)]
e.briggsae.vs.hybrids.male = DGEList(e.briggsae.vs.hybrids.male)
e.briggsae.vs.hybrids.male$samples$group = factor(sub("[1-3]$","",rownames(e.briggsae.vs.hybrids.male$samples)))
e.briggsae.vs.hybrids.male = calcNormFactors(e.briggsae.vs.hybrids.male, method = "TMM")
e.briggsae.vs.hybrids.male = estimateCommonDisp(e.briggsae.vs.hybrids.male)
e.briggsae.vs.hybrids.male = estimateGLMTrendedDisp(e.briggsae.vs.hybrids.male)
e.briggsae.vs.hybrids.male = estimateTagwiseDisp(e.briggsae.vs.hybrids.male)
design.briggsae.vs.hybrids.male = model.matrix(~droplevels(factor(species[c(4:6,10:12)])))
colnames(design.briggsae.vs.hybrids.male)[2] = "F1"
v.briggsae.vs.hybrids.male = voom(e.briggsae.vs.hybrids.male, design.briggsae.vs.hybrids.male, plot=F)
vfit.briggsae.vs.hybrids.male = lmFit(v.briggsae.vs.hybrids.male, design.briggsae.vs.hybrids.male)
vfit.briggsae.vs.hybrids.male = eBayes(vfit.briggsae.vs.hybrids.male)
dt.briggsae.vs.hybrids.male = decideTests(vfit.briggsae.vs.hybrids.male, p.value = 0.05)
tt.briggsae.vs.hybrids.male = topTable(vfit.briggsae.vs.hybrids.male, n = Inf)
tt.briggsae.vs.hybrids.male = tt.briggsae.vs.hybrids.male[ rownames(dt.briggsae.vs.hybrids.male), ]
write.csv(tt.briggsae.vs.hybrids.male, file="tt.briggsae.vs.hybrids.males.csv", quote=F)

## C. nigoni vs F1
# females
e.nigoni.vs.hybrids.female = e$counts[,c(13:15,7:9)]
e.nigoni.vs.hybrids.female = DGEList(e.nigoni.vs.hybrids.female)
e.nigoni.vs.hybrids.female$samples$group = factor(sub("[1-3]$","",rownames(e.briggsae.vs.hybrids.female$samples)))
e.nigoni.vs.hybrids.female = calcNormFactors(e.nigoni.vs.hybrids.female, method = "TMM")
e.nigoni.vs.hybrids.female = estimateCommonDisp(e.nigoni.vs.hybrids.female)
e.nigoni.vs.hybrids.female = estimateGLMTrendedDisp(e.nigoni.vs.hybrids.female)
e.nigoni.vs.hybrids.female = estimateTagwiseDisp(e.nigoni.vs.hybrids.female)
design.nigoni.vs.hybrids.female = model.matrix(~droplevels(factor(species[c(13:15,7:9)])))
colnames(design.nigoni.vs.hybrids.female)[2] = "F1"
v.nigoni.vs.hybrids.female = voom(e.nigoni.vs.hybrids.female, design.nigoni.vs.hybrids.female, plot=F)
vfit.nigoni.vs.hybrids.female = lmFit(v.nigoni.vs.hybrids.female, design.nigoni.vs.hybrids.female)
vfit.nigoni.vs.hybrids.female = eBayes(vfit.nigoni.vs.hybrids.female)
dt.nigoni.vs.hybrids.female = decideTests(vfit.nigoni.vs.hybrids.female, p.value = 0.05)
tt.nigoni.vs.hybrids.female = topTable(vfit.nigoni.vs.hybrids.female, n = Inf)
tt.nigoni.vs.hybrids.female = tt.nigoni.vs.hybrids.female[ rownames(dt.nigoni.vs.hybrids.female), ]
write.csv(tt.nigoni.vs.hybrids.female, file="tt.nigoni.vs.hybrids.females.csv", quote=F)

# males
e.nigoni.vs.hybrids.male = e$counts[,c(16:18,10:12)]
e.nigoni.vs.hybrids.male = DGEList(e.nigoni.vs.hybrids.male)
e.nigoni.vs.hybrids.male$samples$group = factor(sub("[1-3]$","",rownames(e.briggsae.vs.hybrids.male$samples)))
e.nigoni.vs.hybrids.male = calcNormFactors(e.nigoni.vs.hybrids.male, method = "TMM")
e.nigoni.vs.hybrids.male = estimateCommonDisp(e.nigoni.vs.hybrids.male)
e.nigoni.vs.hybrids.male = estimateGLMTrendedDisp(e.nigoni.vs.hybrids.male)
e.nigoni.vs.hybrids.male = estimateTagwiseDisp(e.nigoni.vs.hybrids.male)
design.nigoni.vs.hybrids.male = model.matrix(~droplevels(factor(species[c(16:18,10:12)])))
colnames(design.nigoni.vs.hybrids.male)[2] = "F1"
v.nigoni.vs.hybrids.male = voom(e.nigoni.vs.hybrids.male, design.nigoni.vs.hybrids.male, plot=F)
vfit.nigoni.vs.hybrids.male = lmFit(v.nigoni.vs.hybrids.male, design.nigoni.vs.hybrids.male)
vfit.nigoni.vs.hybrids.male = eBayes(vfit.nigoni.vs.hybrids.male)
dt.nigoni.vs.hybrids.male = decideTests(vfit.nigoni.vs.hybrids.male, p.value = 0.05)
tt.nigoni.vs.hybrids.male = topTable(vfit.nigoni.vs.hybrids.male, n = Inf)
tt.nigoni.vs.hybrids.male = tt.nigoni.vs.hybrids.male[ rownames(dt.nigoni.vs.hybrids.male), ]
write.csv(tt.nigoni.vs.hybrids.male, file="tt.nigoni.vs.hybrids.males.csv", quote=F)

# consolidate
# expression divergence per sex between parents and F1

models.inheritance.male = paste(dt.briggsae.vs.nigoni.male[,2], dt.briggsae.vs.hybrids.male[,2], dt.nigoni.vs.hybrids.male[,2])
models.inheritance.female = paste(dt.briggsae.vs.nigoni.female[,2], dt.briggsae.vs.hybrids.female[,2], dt.nigoni.vs.hybrids.female[,2])

# classify genes based on F1 expression relative to parent species
# order "Cbr.vs.Cni Cbr.vs.F1 Cni.vs.F1"
# 0 = no DE
# -1 = C. briggsae
# 1 = C. nigoni
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

df.inherit = data.frame(female=class_inheritance.female, male=class_inheritance.male, gene=rownames(dt.briggsae.vs.nigoni.female))
df.inherit2 = gather(df.inherit, sex, class, -gene)
df.inherit2$class = factor(df.inherit2$class, levels=class.levels)
df.inherit2$chromosome = sub("\\..*","", rownames(dt.briggsae.vs.hybrids.male))
df.inherit2$logFC.F1.vs.briggsae = c(tt.briggsae.vs.hybrids.female[,1],tt.briggsae.vs.hybrids.male[,1])
df.inherit2$logFC.F1.vs.nigoni = c(tt.nigoni.vs.hybrids.female[,1],tt.nigoni.vs.hybrids.male[,1])

# expression distance F1 from zero
df.inherit2$F1_dist_from_zero = pointDistance(as.matrix(df.inherit2[,c("logFC.F1.vs.briggsae","logFC.F1.vs.nigoni")]),
																							matrix(0, ncol=2, nrow=dim(df.inherit2)[1]),
																							lonlat=F, allpairs=F)
write.csv(df.inherit2, file="expr_inheritance_logFC_F1_parents.csv", row.names=F)

# enrichment and counts per chromosome
df.inherit3 = (df.inherit2 %>% group_by(sex, class, chromosome) %>% count())
x = matrix((df.inherit3 %>% filter(sex == "female"))$n, ncol=7, nrow=6)
df.inherit3[1:42,"enrichment"] = as.vector(enrichment(x, odds.ratio=T))
df.inherit3[1:42,"p.value"] = as.vector(enrichment(x, odds.ratio=F))
x = matrix((df.inherit3 %>% filter(sex == "male"))$n, ncol=7, nrow=6)
df.inherit3[43:84,"enrichment"] = as.vector(enrichment(x, odds.ratio=T))
df.inherit3[43:84,"p.value"] = as.vector(enrichment(x, odds.ratio=F))
df.inherit3$chr.pseudo = (1:6)[ factor(df.inherit3$chromosome)]
df.inherit3$sig = ""
df.inherit3$sig[ df.inherit3$p.value < 0.05 & abs(log2(df.inherit3$enrichment)) > 0.5 ] = "*"
write.csv(df.inherit3, file="expr_inheritance_chr_enrich.csv", row.names=F)

##################################
### allele-specific expression ###
##################################

# male
xchr = rownames(counts.ase)[grep("^X\\.", rownames(counts.ase))]
counts.ase = counts.ase[keep.exprs, c(seq(1,12,2), seq(2,12,2))]
e.AS.male = DGEList(counts.AS[! rownames(counts.AS) %in% xchr ,c(4:6,10:12)])
e.AS.male$samples$group = factor(sub("[1-3]_","_",rownames(e.AS.male$samples)))
e.AS.male = calcNormFactors(e.AS.male, method = "TMM")
e.AS.male = estimateCommonDisp(e.AS.male)
e.AS.male = estimateGLMTrendedDisp(e.AS.male)
e.AS.male = estimateTagwiseDisp(e.AS.male)
design.AS.male = model.matrix(~droplevels(factor(species[c(4:6,10:12)])))
colnames(design.AS.male)[2] = "nigoni"
v.AS.male = voom(e.AS.male, design.AS.male, plot=F)
vfit.AS.male = lmFit(v.AS.male, design.AS.male)
vfit.AS.male = eBayes(vfit.AS.male)
dt.AS.male = decideTests(vfit.AS.male, p.value = 0.05)
tt.AS.male = topTable(vfit.AS.male, n = Inf)
tt.AS.male = tt.AS.male[ rownames(dt.AS.male), ]
write.csv(tt.AS.male, file="tt.AS.males.csv", quote=F)

# female
e.AS.female = DGEList(counts.ase[,c(1:3,7:9)])
e.AS.female$samples$group = factor(sub("[1-3]_","_",rownames(e.AS.female$samples)))
e.AS.female = calcNormFactors(e.AS.female, method = "TMM")
e.AS.female = estimateCommonDisp(e.AS.female)
e.AS.female = estimateGLMTrendedDisp(e.AS.female)
e.AS.female = estimateTagwiseDisp(e.AS.female)
design.AS.female = model.matrix(~droplevels(factor(species[c(1:3,7:9)])))
colnames(design.AS.female)[2] = "nigoni"
v.AS.female = voom(e.AS.female, design.AS.female, plot=F)
vfit.AS.female = lmFit(v.AS.female, design.AS.female)
vfit.AS.female = eBayes(vfit.AS.female)
dt.AS.female = decideTests(vfit.AS.female, p.value = 0.05)
tt.AS.female = topTable(vfit.AS.female, n = Inf)
tt.AS.female = tt.AS.female[ rownames(dt.AS.female), ]
write.csv(tt.AS.female, file="tt.AS.females.csv", quote=F)

# get trans effects: is sp1/sp2 in parents different than sp1/sp2 in allele-specific expression

# function that uses Wald-type posthoc hypothesis testing of linear regression
# returns model object pvalue and DF
getTransEffects = function(x){
	x = as.data.frame(x)
	x$T = rep(c("P","H"), each=6)
	x$S = rep(rep(c("br","ni"), each=3), 2)
	fit = lm(x ~ T/S, data=x)
	fit2 = car::linearHypothesis(fit, c("TH:Sni = TP:Sni"))
	p.value = fit2$`Pr(>F)`[2]
	list(x, fit, fit2, p.value)
}

# female
df.trans = cbind(v.briggsae.vs.nigoni.female$E, v.AS.female$E)
trans_effects_tests = mclapply(split(df.trans, 1:dim(df.trans)[1]), getTransEffects, mc.cores = 4)
trans_p.values = sapply(trans_effects_tests, function(x) x[[4]] )
trans_p.values.fdr.female = p.adjust(trans_p.values, method="BH")
# male
df.trans = cbind(v.briggsae.vs.nigoni.male$E[rownames(v.AS.male$E),], v.AS.male$E)
trans_effects_tests = mclapply(split(df.trans, 1:dim(df.trans)[1]), getTransEffects, mc.cores = 4)
trans_p.values = sapply(trans_effects_tests, function(x) x[[4]] )
trans_p.values.fdr.male = p.adjust(trans_p.values, method="BH")

# tmp.xchr.male = df.inherit2 %>% filter(sex == "male" & substr(gene,1,1) == "X")
# xchr = as.character(tmp.xchr.male$gene)
# tmp.xchr.male$type = as.character(tmp.xchr.male$class)
# tmp.xchr.male$type[ tmp.xchr.male$type == "C. nigoni dominant" ] = "cis only"
# tmp.xchr.male$type[ tmp.xchr.male$type == "C. briggsae dominant" ] = "trans only"
# tmp.xchr.male$type[ models.inheritance.male[ which(rownames(dt.briggsae.vs.nigoni.male) %in% xchr) ] == "0 1 1" ] = "cis-trans (compensatory)"
# tmp.xchr.male$type[ models.inheritance.male[ which(rownames(dt.briggsae.vs.nigoni.male) %in% xchr) ] == "0 -1 -1" ] = "cis-trans (compensatory)"
# tmp.xchr.male$type[ models.inheritance.male[ which(rownames(dt.briggsae.vs.nigoni.male) %in% xchr) ] == "1 -1 -1" ] = "cis-trans (enhancing)"
# tmp.xchr.male$type[ models.inheritance.male[ which(rownames(dt.briggsae.vs.nigoni.male) %in% xchr) ] == "-1 -1 -1" ] = "cis-trans (enhancing)"
# tmp.xchr.male$type[ models.inheritance.male[ which(rownames(dt.briggsae.vs.nigoni.male) %in% xchr) ] == "1 1 1" ] = "cis-trans (enhancing)"
# tmp.xchr.male$type[ models.inheritance.male[ which(rownames(dt.briggsae.vs.nigoni.male) %in% xchr) ] == "-1 1 1" ] = "cis-trans (enhancing)"
# tmp.xchr.male$type[ tmp.xchr.male$type == "additive" ] = "cis-trans (enhancing)"
# tmp.xchr.male$type[ tmp.xchr.male$type == "no change" ] = "conserved"

auto = rownames(dt.briggsae.vs.nigoni.male)[ !rownames(dt.briggsae.vs.nigoni.male) %in% xchr ]
# consolidate
df.AS.auto = data.frame(sex=rep(c("female","male"), each=length(auto)),
									species_AS=c(c("Cbr","no DE","Cni")[factor(dt.AS.female[auto,2])], c("Cbr","no DE","Cni")[factor(dt.AS.male[auto,2])]),
									species_parents=c(c("Cbr","no DE","Cni")[factor(dt.briggsae.vs.nigoni.female[auto,2])], c("Cbr","no DE","Cni")[factor(dt.briggsae.vs.nigoni.male[auto,2])]),
									logFC_F1_cis=c(tt.AS.female[auto,1],tt.AS.male[auto,1]),
									logFC_parents=c(tt.briggsae.vs.nigoni.female[auto,1],tt.briggsae.vs.nigoni.male[auto,1]),
									trans_effects = c("non sig","sig")[ as.factor(c(trans_p.values.fdr.female[!rownames(dt.briggsae.vs.nigoni.male) %in% xchr],trans_p.values.fdr.male) < 0.05) ],
									p.value=c(trans_p.values.fdr.female[ !rownames(dt.briggsae.vs.nigoni.male) %in% xchr ],trans_p.values.fdr.male),
									genes=c(auto,auto),
									type=NA)

# classify allele specific expression
# x = "parents hybrids_ASE trans_effects"
get_regulation_type = function(x){
	regulation_type = vector("character", length=length(x))
	regulation_type[ x == "Cbr Cbr non sig" ] = "cis only"
	regulation_type[ x == "Cni Cni non sig" ] = "cis only"
	regulation_type[ x == "Cbr Cni non sig" ] = "cis only"
	regulation_type[ x == "Cni Cbr non sig" ] = "cis only"
	regulation_type[ x == "Cni no DE sig" ] = "trans only"
	regulation_type[ x == "Cbr no DE sig" ] = "trans only"
	regulation_type[ x == "Cbr no DE sig" ] = "trans only"
	regulation_type[ x == "Cni no DE sig" ] = "trans only"
	regulation_type[ x == "Cni Cni sig" ] = "cis-trans (enhancing)" # cis + trans
	regulation_type[ x == "Cbr Cbr sig" ] = "cis-trans (enhancing)" # cis + trans
	regulation_type[ x == "Cbr Cbr sig" ] = "cis-trans (enhancing)" # cis + trans
	regulation_type[ x == "Cni Cni sig" ] = "cis-trans (enhancing)" # cis + trans
	regulation_type[ x == "no DE no DE non sig" ] = "conserved"
	regulation_type[ x == "Cni Cbr sig" ] = "cis-trans (compensatory)" # cis x trans
	regulation_type[ x == "Cbr Cni sig" ] = "cis-trans (compensatory)" # cis x trans
	regulation_type[ x == "Cni Cbr sig" ] = "cis-trans (compensatory)" # cis x trans
	regulation_type[ x == "Cbr Cni sig" ] = "cis-trans (compensatory)" # cis x trans
	regulation_type[ x == "no DE Cbr sig" ] = "cis-trans (compensatory)"
	regulation_type[ x == "no DE Cni sig" ] = "cis-trans (compensatory)"
	regulation_type[ x == "no DE Cni sig" ] = "cis-trans (compensatory)"
	regulation_type[ x == "no DE Cbr sig" ] = "cis-trans (compensatory)"
	regulation_type[ regulation_type == "" ] = "ambiguous"
	return(regulation_type)
}

tmp = paste(df.AS.auto[!is.na(df.AS.auto$species_AS),"species_parents"],df.AS.auto[!is.na(df.AS.auto$species_AS),"species_AS"],df.AS.auto[!is.na(df.AS.auto$species_AS),"trans_effects"])
type.levels = c("conserved","ambiguous","trans only","cis only","cis-trans (enhancing)","cis-trans (compensatory)")
df.AS.auto$type = factor(get_regulation_type(tmp), levels=type.levels)
df.AS.auto$chromosome = sub("\\..*","", df.AS.auto$gene)
write.csv(df.AS.auto, file="df.cis_trans.regulation_type.auto.csv", row.names=F)

# counts

df.AS.auto2 = as.data.frame(df.AS.auto %>% group_by(sex, type, chromosome) %>% count())
x = matrix((df.AS.auto2 %>% filter(sex == "female"))$n, ncol=6, nrow=6)
df.AS.auto2$enrichment = NA
df.AS.auto2$enrichment[1:36] = as.vector(enrichment(x, odds.ratio=T))
df.AS.auto2$p.value = NA
df.AS.auto2$p.value[1:36] = as.vector(enrichment(x, odds.ratio=F))
x = matrix((df.AS.auto2 %>% filter(sex == "male"))$n, ncol=6, nrow=6)
df.AS.auto2$enrichment[37:72] = as.vector(enrichment(x, odds.ratio=T))
df.AS.auto2$p.value[37:72] = as.vector(enrichment(x, odds.ratio=F))
df.AS.auto2$chr.pseudo = (1:6)[ factor(df.AS.auto2$chromosome)]

df.AS.auto2$sig = ""
df.AS.auto2$sig[ df.AS.auto2$p.value < 0.05 & abs(log2(df.AS.auto2$enrichment)) > 0.5 ] = "*"
write.csv(df.AS.auto2, file="df.cis_trans.counts.chr_enrich.csv", row.names=F)
