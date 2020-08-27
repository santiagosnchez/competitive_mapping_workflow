# DE and stats
library(limma)
library(edgeR)
library(MASS)
library(car)
# data wrangling
library(tidyr)
library(dplyr)
library(tibble)
# plotting
library(gridExtra)
library(RColorBrewer)
library(ggsci)
library(ggplot2)
library(cowplot)
library(lemon)
library(ggridges)
# multithread processing
library(parallel)

theme_set(theme_cowplot())

# compare counts
counts_species = read.table("real_counts_species.txt", head=T, sep="\t", row.names=1)
counts_fc_hf1 = read.table("real_counts_HF1.txt", head=T, sep="\t", row.names=1)
counts_cm_hf1 = read.table("compmap_counts_HF1.txt", head=T, sep="\t", row.names=1)
counts_cm_hf1 = counts_cm_hf1[, colnames(counts_fc_hf1)]
counts_cm_hf1 = counts_cm_hf1[ rownames(counts_fc_hf1),]
colnames(counts_fc_hf1) = sub("$", "_featureCounts", colnames(counts_fc_hf1))
colnames(counts_cm_hf1) = sub("$", "_CompMap", colnames(counts_cm_hf1))
compare_counts = cbind(counts_fc_hf1, counts_cm_hf1, gene=rownames(counts_fc_hf1))

# melt data set
df.compare_counts = compare_counts %>%
    pivot_longer(cols=-gene, names_to="samples", values_to="counts") %>%
    extract(samples, c("sample","species","method"), "(HF1_[1-3])_(Cbr|Cni)_(CompMap|featureCounts)") %>%
    pivot_wider(names_from=method, values_from=counts)

# plot with ggplot
ggplot(df.compare_counts, aes(log2(featureCounts), log2(CompMap), color=sample)) +
    geom_point() +
    geom_smooth(method="lm", se=F, color="black", linetype=2) +
    facet_rep_grid(sample~species) +
    labs(x="True log2 counts (featureCounts)",y="Competitive mapping log2 counts (CompMap)") +
    background_grid() +
    xlim(4,12) + ylim(4,12)

# counts edgeR objects: DGEList
# species
e.SP = DGEList(counts_species)
groups.SP = sub("_.*","",colnames(counts_species))
e.SP$samples$group = factor(groups.SP)
e.SP = calcNormFactors(e.SP, method = "TMM")
e.SP = estimateCommonDisp(e.SP)
e.SP = estimateGLMTrendedDisp(e.SP)
e.SP = estimateTagwiseDisp(e.SP)
design.SP = model.matrix(~groups.SP)
colnames(design.SP)[2] = "Cni"
v.SP = voom(e.SP, design.SP, plot=F)
vfit.SP = lmFit(v.SP, design.SP)
vfit.SP = eBayes(vfit.SP)
dt.SP = decideTests(vfit.SP, p.value = 0.05)
tt.SP = topTable(vfit.SP, n = Inf)
tt.SP = tt.SP[ rownames(dt.SP), ]

# alleles real
e.AS.real = DGEList(counts_fc_hf1)
e.AS.real$samples$group = factor(groups.SP)
e.AS.real = calcNormFactors(e.AS.real, method = "TMM")
e.AS.real = estimateCommonDisp(e.AS.real)
e.AS.real = estimateGLMTrendedDisp(e.AS.real)
e.AS.real = estimateTagwiseDisp(e.AS.real)
v.AS.real = voom(e.AS.real, design.SP, plot=F)
vfit.AS.real = lmFit(v.AS.real, design.SP)
vfit.AS.real = eBayes(vfit.AS.real)
dt.AS.real = decideTests(vfit.AS.real, p.value = 0.05)
tt.AS.real = topTable(vfit.AS.real, n = Inf)
tt.AS.real = tt.AS.real[ rownames(dt.AS.real), ]

# alleles CompMap
e.AS.compmap = DGEList(counts_cm_hf1)
e.AS.compmap$samples$group = factor(groups.SP)
e.AS.compmap = calcNormFactors(e.AS.compmap, method = "TMM")
e.AS.compmap = estimateCommonDisp(e.AS.compmap)
e.AS.compmap = estimateGLMTrendedDisp(e.AS.compmap)
e.AS.compmap = estimateTagwiseDisp(e.AS.compmap)
v.AS.compmap = voom(e.AS.compmap, design.SP, plot=F)
vfit.AS.compmap = lmFit(v.AS.compmap, design.SP)
vfit.AS.compmap = eBayes(vfit.AS.compmap)
dt.AS.compmap = decideTests(vfit.AS.compmap, p.value = 0.05)
tt.AS.compmap = topTable(vfit.AS.compmap, n = Inf)
tt.AS.compmap = tt.AS.compmap[ rownames(dt.AS.compmap), ]

# read real fc values
simulated_real_fc = read.csv("../sim_fold_change_comb.csv", row.names=1)

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

trans.real = cbind(v.SP$E, v.AS.real$E)
trans.compmap = cbind(v.SP$E, v.AS.compmap$E)
trans_effects_tests = mclapply(split(trans.real, 1:dim(trans.real)[1]), getTransEffects, mc.cores = 4)
trans_p.values = sapply(trans_effects_tests, function(x) x[[4]] )
trans_p.values.fdr = p.adjust(trans_p.values, method="BH")
trans_effects.real = trans_p.values.fdr
trans_effects_tests = mclapply(split(trans.compmap, 1:dim(trans.compmap)[1]), getTransEffects, mc.cores = 4)
trans_p.values = sapply(trans_effects_tests, function(x) x[[4]] )
trans_p.values.fdr = p.adjust(trans_p.values, method="BH")
trans_effects.compmap = trans_p.values.fdr

# combine results
combined_DE = data.frame(species.logFC=tt.SP$logFC, species.pvalue=tt.SP$adj.P.Val,
    alleles.logFC.real=tt.AS.real$logFC, alleles.pvalue.real=tt.AS.real$adj.P.Val, trans_effects.real=trans_effects.real,
    alleles.logFC.compmap=tt.AS.compmap$logFC, alleles.pvalue.compmap=tt.AS.compmap$adj.P.Val, trans_effects.compmap=trans_effects.compmap,
    gene=rownames(tt.SP))

# classify true counts
class_real = rep(NA, 1000)
class_real[ combined_DE[,"species.pvalue"] > 0.05 & combined_DE[,"alleles.pvalue.real"] > 0.05 ] = "conserved"
class_real[ is.na(class_real) & combined_DE[,"trans_effects.real"] > 0.05 & combined_DE[,"alleles.pvalue.real"] < 0.05 & combined_DE[,"species.pvalue"] < 0.05 ] = "cis only"
class_real[ is.na(class_real) & combined_DE[,"trans_effects.real"] < 0.05 & combined_DE[,"alleles.pvalue.real"] > 0.05 & combined_DE[,"species.pvalue"] < 0.05 ] = "trans only"
class_real[ is.na(class_real) & combined_DE[,"trans_effects.real"] < 0.05 & combined_DE[,"alleles.pvalue.real"] < 0.05 & combined_DE[,"species.pvalue"] > 0.05 ] = "cis-trans (compensatory)"
class_real[ is.na(class_real) & combined_DE[,"trans_effects.real"] < 0.05 & ((combined_DE[,"alleles.logFC.real"] > 0 & combined_DE[,"species.logFC"] < 0) | (combined_DE[,"alleles.logFC.real"] < 0 & combined_DE[,"species.logFC"] > 0)) ] = "cis-trans x (compensatory)"
class_real[ is.na(class_real) & combined_DE[,"trans_effects.real"] < 0.05 & ((combined_DE[,"alleles.logFC.real"] > 0 & combined_DE[,"species.logFC"] > 0) | (combined_DE[,"alleles.logFC.real"] < 0 & combined_DE[,"species.logFC"] < 0)) ] = "cis-trans (enhancing)"
class_real[ is.na(class_real) ] = "ambiguous"
# classify compmap counts
class_compmap = rep(NA, 1000)
class_compmap[ combined_DE[,"species.pvalue"] > 0.05 & combined_DE[,"alleles.pvalue.compmap"] > 0.05 ] = "conserved"
class_compmap[ is.na(class_compmap) & combined_DE[,"trans_effects.compmap"] > 0.05 & combined_DE[,"alleles.pvalue.compmap"] < 0.05 & combined_DE[,"species.pvalue"] < 0.05 ] = "cis only"
class_compmap[ is.na(class_compmap) & combined_DE[,"trans_effects.compmap"] < 0.05 & combined_DE[,"alleles.pvalue.compmap"] > 0.05 & combined_DE[,"species.pvalue"] < 0.05 ] = "trans only"
class_compmap[ is.na(class_compmap) & combined_DE[,"trans_effects.compmap"] < 0.05 & combined_DE[,"alleles.pvalue.compmap"] < 0.05 & combined_DE[,"species.pvalue"] > 0.05 ] = "cis-trans (compensatory)"
class_compmap[ is.na(class_compmap) & combined_DE[,"trans_effects.compmap"] < 0.05 & ((combined_DE[,"alleles.logFC.compmap"] > 0 & combined_DE[,"species.logFC"] < 0) | (combined_DE[,"alleles.logFC.compmap"] < 0 & combined_DE[,"species.logFC"] > 0)) ] = "cis-trans x (compensatory)"
class_compmap[ is.na(class_compmap) & combined_DE[,"trans_effects.compmap"] < 0.05 & ((combined_DE[,"alleles.logFC.compmap"] > 0 & combined_DE[,"species.logFC"] > 0) | (combined_DE[,"alleles.logFC.compmap"] < 0 & combined_DE[,"species.logFC"] < 0)) ] = "cis-trans (enhancing)"
class_compmap[ is.na(class_compmap) ] = "ambiguous"
combined_DE = cbind(combined_DE, class_real, class_compmap)

# classify by fold-change
# cm_class_fc = rep(NA, 999)
# cm_class_fc[ abs(combined[,"species_fc_cm"]) < 1 & abs(combined[,"allele_fc_cm"]) < 1 ] = "conserved"
# cm_class_fc[ is.na(cm_class_fc) & abs(combined[,"species_fc_cm"] - combined[,"allele_fc_cm"]) < 1 & abs(combined[,"allele_fc_cm"]) > 1 & abs(combined[,"species_fc_cm"]) > 1 ] = "cis only"
# cm_class_fc[ is.na(cm_class_fc) & abs(combined[,"species_fc_cm"] - combined[,"allele_fc_cm"]) > 1 & abs(combined[,"allele_fc_cm"]) < 1 & abs(combined[,"species_fc_cm"]) > 1 ] = "trans only"
# cm_class_fc[ is.na(cm_class_fc) & abs(combined[,"species_fc_cm"] - combined[,"allele_fc_cm"]) > 1 & abs(combined[,"allele_fc_cm"]) > 1 & abs(combined[,"species_fc_cm"]) < 1 ] = "cis-trans (compensatory)"
# cm_class_fc[ is.na(cm_class_fc) & abs(combined[,"species_fc_cm"] - combined[,"allele_fc_cm"]) > 1 & ((combined[,"allele_fc_cm"] > 1 & combined[,"species_fc_cm"] < -1) | (combined[,"allele_fc_cm"] < -1 & combined[,"species_fc_cm"] > 1)) ] = "cis-trans x (compensatory)"
# cm_class_fc[ is.na(cm_class_fc) & abs(combined[,"species_fc_cm"] - combined[,"allele_fc_cm"]) > 1 & ((combined[,"allele_fc_cm"] > 1 & combined[,"species_fc_cm"] > 1) | (combined[,"allele_fc_cm"] < -1 & combined[,"species_fc_cm"] < -1)) ] = "cis-trans (enhancing)"
# cm_class_fc[ is.na(cm_class_fc) ] = "ambiguous"
# combined = cbind(combined, cm_class_fc)

# data wrangling

# counts
combined_DE %>%
    select(class.real, class.compmap) %>%
    pivot_longer(cols=everything(), names_to="method", values_to="class") %>%
    extract(method, "method", "\\.(.*)") %>%
    group_by(method, class) %>%
    count() -> combined_DE_counts

# tmp with class
tmp = combined_DE %>%
    select(class.real, class.compmap) %>%
    pivot_longer(cols=everything(), names_to="method", values_to="class") %>%
    select(class)

combined_DE %>%
    select(species.logFC, alleles.logFC.real, alleles.logFC.compmap, class.real, class.compmap) %>%
    pivot_longer(cols=c(-species.logFC, -class.real, -class.compmap), names_to="method", values_to="alleles.logFC") %>%
    extract(method, "method", "\\.(real|compmap)$") %>%
    add_column(tmp) %>%
    select(-class.real, -class.compmap) %>%
    select(species.logFC, alleles.logFC, class, method) -> combined_DE_logFC

# combined_DE %>%
#     pivot_longer(cols=c(-species.logFC, -species.pvalue, -gene, -class.real, -class.compmap, -trans_effects.pvalue.compmap, -trans_effects.pvalue.real), names_to="key", values_to="value") %>%
#     extract(key, c("effect","pvalue","method"), "(species|alleles)\\.(pvalue|logFC)\\.(compmap|real)") %>%
#     pivot_wider(names_from=pvalue, values_from=value) %>%
#     pivot_wider(names_from=effect, values_from=c(logFC, pvalue)) %>%
#     pivot_longer(cols=c(-species.logFC, -species.pvalue, -trans_effects.pvalue.real, -trans_effects.pvalue.compmap, -gene, -method, -logFC_alleles, -pvalue_alleles), values_to="class") %>%
#     select(-name) %>%
#     unique() -> df.combined_DE

brewer.pal(10, "Paired") -> cols
cols = cols[-c(7,8)]
cols.class = c("black","grey",cols[ c(4,6,2,8,7) ])
cols.type = c("black","grey",pal_jama()(5))
cols.type[3] = "#09BC8A"

p1 = ggplot(combined_DE_logFC, aes(x=species.logFC, y=alleles.logFC, color=class)) +
    geom_point(alpha=0.5) +
    facet_rep_wrap(~method) +
    background_grid() +
    labs(x="Species (Cni/Cbr) log2-FC", y="Alleles (Cni/Cbr) log2-FC") +
    scale_color_manual(values=cols.type)
p2 = ggplot(combined_DE_counts, aes(x=class, y=n, fill=class)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=n), nudge_y=20) +
    facet_rep_wrap(~method) +
    background_grid(major="y") +
    labs(x="", y="number of genes") +
    scale_fill_manual(values=cols.type) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

plot_grid(p1, p2, nrow=2, rel_heights=c(0.4,0.6))
