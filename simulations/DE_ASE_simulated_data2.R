# DE and stats
library(DESeq2)
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
    pivot_wider(names_from=method, values_from=counts) %>%
    mutate(ratio = log2(CompMap) / log2(featureCounts))

# plot with ggplot
biplots = ggplot(df.compare_counts, aes(log2(featureCounts), log2(CompMap), color=sample)) +
    geom_point(show.legend=F) +
    geom_smooth(method="lm", se=F, color="black", linetype=2) +
    facet_rep_grid(sample~species) +
    labs(x="True log2 counts (featureCounts)",y="Competitive mapping log2 counts (CompMap)") +
    background_grid() +
    xlim(4,12) + ylim(4,12)

histogram = ggplot(df.compare_counts, aes(x=log2(CompMap)/log2(featureCounts), fill=sample)) +
    geom_histogram(show.legend=F) +
    background_grid() +
    facet_rep_grid(sample~species) +
    xlim(0.75,1.25)

plot_grid(biplots, histogram, labels=c("A","B"), ncol=2)
ggsave("log2_count_comp_true_vs_compmap.png")
ggsave("log2_count_comp_true_vs_compmap.pdf", device="pdf", useDingbats=F)

# gene lengths
br.fasta = readDNAStringSet("../briggsae.cds.small.fa")
ni.fasta = readDNAStringSet("../nigoni.cds.small.fa")
br.gene_len = width(br.fasta)
ni.gene_len = width(ni.fasta)
names(br.gene_len) = paste(names(br.fasta), names(ni.fasta), sep="__")
names(ni.gene_len) = paste(names(br.fasta), names(ni.fasta), sep="__")

# divergence
div = read.csv("../ortho_div.txt", head=F)
divergence = div[,3]
names(divergence) = paste(div[,1], div[,2], sep="__")

# merge with df
df.compare_counts$Cbr.gene_length = br.gene_len[ df.compare_counts$gene ]
df.compare_counts$Cni.gene_length = ni.gene_len[ df.compare_counts$gene ]
df.compare_counts$divergence = divergence[ df.compare_counts$gene ]

ggplot(df.compare_counts, aes(ratio, Cbr.gene_length)) + geom_point()


# DE with DESeq2
# species
e.SP = DGEList(counts_species)
groups.SP = sub("_.*","",colnames(counts_species))
e.SP$samples$group = factor(groups.SP)
dds.SP = DESeqDataSetFromMatrix(countData = e.SP$counts, colData = e.SP$samples, design = ~ group)
dds.SP = DESeq(dds.SP)
dds.SP.res = results(dds.SP)
#dds.SP = rlog(dds.SP, blind=FALSE)

# ASE
# real
e.AS.real = DGEList(counts_fc_hf1)
e.AS.real$samples$group = factor(groups.SP)
dds.AS.real = DESeqDataSetFromMatrix(countData = e.AS.real$counts, colData = e.AS.real$samples, design = ~ group)
dds.AS.real = DESeq(dds.AS.real)
dds.AS.real.res = results(dds.AS.real)
#dds.AS.real = rlog(dds.AS.real, blind=FALSE)

# compmap
e.AS.compmap = DGEList(counts_cm_hf1)
e.AS.compmap$samples$group = factor(groups.SP)
dds.AS.compmap = DESeqDataSetFromMatrix(countData = e.AS.compmap$counts, colData = e.AS.compmap$samples, design = ~ group)
dds.AS.compmap = DESeq(dds.AS.compmap)
dds.AS.compmap.res = results(dds.AS.compmap)
#dds.AS.compmap = rlog(dds.AS.compmap, blind=FALSE)

# trans real
# e.AS.trans.real = DGEList(cbind(counts_species,counts_fc_hf1))
# e.AS.trans.real$samples$group = factor(rep(groups.SP, 2))
# e.AS.trans.real$samples$treat = factor(rep(c("P","H"), each=6))
# dds.AS.trans.real = DESeqDataSetFromMatrix(countData = e.AS.trans.real$counts, colData = e.AS.trans.real$samples, design = ~ treat/group)
# dds.AS.trans.real = DESeq(dds.AS.trans.real)
# dds.AS.trans.real.res = results(dds.AS.trans.real, contrast=c("treat","P","H"))
# #dds.AS.trans.real = rlog(dds.AS.trans.real, blind=FALSE)
#
# # trans compmap
# e.AS.trans.compmap = DGEList(cbind(counts_species,counts_cm_hf1))
# e.AS.trans.compmap$samples$group = factor(rep(groups.SP, 2))
# e.AS.trans.compmap$samples$treat = factor(rep(c("P","H"), each=6))
# dds.AS.trans.compmap = DESeqDataSetFromMatrix(countData = e.AS.trans.compmap$counts, colData = e.AS.trans.compmap$samples, design = ~ treat/group)
# dds.AS.trans.compmap = DESeq(dds.AS.trans.compmap)
# dds.AS.trans.compmap.res = results(dds.AS.trans.compmap, contrast=c("treat","P","H"))
#dds.AS.trans.compmap = rlog(dds.AS.trans.compmap, blind=FALSE)

getTransEffects = function(x){
	x = as.data.frame(x)
	x$T = rep(c("P","H"), each=6)
	x$S = rep(rep(c("br","ni"), each=3), 2)
	fit = lm(x ~ T/S, data=x)
	fit2 = car::linearHypothesis(fit, c("TH:Sni = TP:Sni"))
	p.value = fit2$`Pr(>F)`[2]
	list(x, fit, fit2, p.value)
}

trans.real = cbind(cpm(e.SP, log=T), cpm(e.AS.real, log=T))
trans.compmap = cbind(cpm(e.SP, log=T), cpm(e.AS.compmap, log=T))
trans_effects_tests = mclapply(split(trans.real, 1:dim(trans.real)[1]), getTransEffects, mc.cores = 4)
trans_p.values = sapply(trans_effects_tests, function(x) x[[4]] )
trans_p.values.fdr = p.adjust(trans_p.values, method="BH")
trans_effects.real = trans_p.values.fdr
trans_effects_tests = mclapply(split(trans.compmap, 1:dim(trans.compmap)[1]), getTransEffects, mc.cores = 4)
trans_p.values = sapply(trans_effects_tests, function(x) x[[4]] )
trans_p.values.fdr = p.adjust(trans_p.values, method="BH")
trans_effects.compmap = trans_p.values.fdr

# read real fc values
simulated_real_fc = read.csv("../sim_fold_change_comb.csv", row.names=1)
simulated_real_fc = simulated_real_fc[ sub("__.*","", rownames(dds.AS.compmap.res)), ]

# combine results
combined_DE_dds = data.frame(species.logFC=dds.SP.res$log2FoldChange, species.pvalue=dds.SP.res$padj,
    alleles.logFC.real=dds.AS.real.res$log2FoldChange, alleles.pvalue.real=dds.AS.real.res$padj, trans_effects.pvalue.real=trans_effects.real,
    alleles.logFC.compmap=dds.AS.compmap.res$log2FoldChange, alleles.pvalue.compmap=dds.AS.compmap.res$padj, trans_effects.pvalue.compmap=trans_effects.compmap,
    gene=rownames(dds.AS.real.res))

cislogFC_real_comp = ggplot(combined_DE_dds, aes(x=alleles.logFC.real, y=alleles.logFC.compmap)) +
    geom_point(alpha=0.2) +
    geom_smooth(method="lm", se=TRUE, color="grey") +
    labs(x="cis logFC (real)", y="cis logFC (CompMap)") +
    facet_wrap(~"cis logFC divergence") +
    background_grid()
translogFC_real_comp = ggplot(combined_DE_dds, aes(x=species.logFC-alleles.logFC.real, y=species.logFC-alleles.logFC.compmap)) +
    geom_point(alpha=0.2) +
    geom_smooth(method="lm", se=TRUE, color="grey") +
    labs(x="trans logFC (real)", y="trans logFC (CompMap)") +
    facet_wrap(~"trans logFC divergence") +
    background_grid()

cis_pval = ggplot(combined_DE_dds, aes(x=-log10(alleles.pvalue.real), y=-log10(alleles.pvalue.compmap))) +
    geom_point(alpha=0.2) +
    geom_smooth(method="lm", se=TRUE, color="grey") +
    labs(x="-log10 p-value (real)", y="-log10 p-value (CompMap)") +
    facet_wrap(~"cis divergence p-values") +
    background_grid()
trans_pval = ggplot(combined_DE_dds, aes(x=-log10(trans_effects.pvalue.real), y=-log10(trans_effects.pvalue.compmap))) +
    geom_point(alpha=0.2) +
    geom_smooth(method="lm", se=TRUE, color="grey") +
    labs(x="-log10 p-value (real)", y="-log10 p-value (CompMap)") +
    facet_wrap(~"trans divergence p-values") +
    background_grid()

volcano_real = ggplot(combined_DE_dds, aes(x=alleles.logFC.real, y=-log10(alleles.pvalue.real))) +
    geom_point(alpha=0.2) +
    geom_hline(yintercept=-log10(0.05), linetype=2) +
    labs(x="cis logFC (real)", y="-log10 p-value (real)") +
    facet_wrap(~"volcano plot (real)") +
    background_grid()
volcano_comp = ggplot(combined_DE_dds, aes(x=alleles.logFC.compmap, y=-log10(alleles.pvalue.compmap))) +
    geom_point(alpha=0.2) +
    geom_hline(yintercept=-log10(0.05), linetype=2) +
    labs(x="cis logFC (CompMap)", y="-log10 p-value (CompMap)") +
    facet_wrap(~"volcano plot (CompMap)") +
    background_grid()

plot_grid(cislogFC_real_comp, translogFC_real_comp, cis_pval, trans_pval, volcano_real, volcano_comp, ncol=2)
ggsave("cis_and_trans_logFC_and_pvalues.png")
ggsave("cis_and_trans_logFC_and_pvalues.pdf", device="pdf", useDingbats=F)

# classify true counts
class.real = rep(NA, 1000)
class.real[ combined_DE_dds[,"species.pvalue"] > 0.05 & combined_DE_dds[,"alleles.pvalue.real"] > 0.05 ] = "conserved"
class.real[ is.na(class.real) & combined_DE_dds[,"trans_effects.pvalue.real"] > 0.05 & combined_DE_dds[,"alleles.pvalue.real"] < 0.05 & combined_DE_dds[,"species.pvalue"] < 0.05 ] = "cis only"
class.real[ is.na(class.real) & combined_DE_dds[,"trans_effects.pvalue.real"] < 0.05 & combined_DE_dds[,"alleles.pvalue.real"] > 0.05 & combined_DE_dds[,"species.pvalue"] < 0.05 ] = "trans only"
class.real[ is.na(class.real) & combined_DE_dds[,"trans_effects.pvalue.real"] < 0.05 & combined_DE_dds[,"alleles.pvalue.real"] < 0.05 & combined_DE_dds[,"species.pvalue"] > 0.05 ] = "cis-trans (compensatory)"
class.real[ is.na(class.real) & combined_DE_dds[,"trans_effects.pvalue.real"] < 0.05 & ((combined_DE_dds[,"alleles.logFC.real"] > 0 & combined_DE_dds[,"species.logFC"] < 0) | (combined_DE_dds[,"alleles.logFC.real"] < 0 & combined_DE_dds[,"species.logFC"] > 0)) ] = "cis-trans x (compensatory)"
class.real[ is.na(class.real) & combined_DE_dds[,"trans_effects.pvalue.real"] < 0.05 & ((combined_DE_dds[,"alleles.logFC.real"] > 0 & combined_DE_dds[,"species.logFC"] > 0) | (combined_DE_dds[,"alleles.logFC.real"] < 0 & combined_DE_dds[,"species.logFC"] < 0)) ] = "cis-trans (enhancing)"
class.real[ is.na(class.real) ] = "ambiguous"
# classify compmap counts
class.compmap = rep(NA, 1000)
class.compmap[ combined_DE_dds[,"species.pvalue"] > 0.05 & combined_DE_dds[,"alleles.pvalue.compmap"] > 0.05 ] = "conserved"
class.compmap[ is.na(class.compmap) & combined_DE_dds[,"trans_effects.pvalue.compmap"] > 0.05 & combined_DE_dds[,"alleles.pvalue.compmap"] < 0.05 & combined_DE_dds[,"species.pvalue"] < 0.05 ] = "cis only"
class.compmap[ is.na(class.compmap) & combined_DE_dds[,"trans_effects.pvalue.compmap"] < 0.05 & combined_DE_dds[,"alleles.pvalue.compmap"] > 0.05 & combined_DE_dds[,"species.pvalue"] < 0.05 ] = "trans only"
class.compmap[ is.na(class.compmap) & combined_DE_dds[,"trans_effects.pvalue.compmap"] < 0.05 & combined_DE_dds[,"alleles.pvalue.compmap"] < 0.05 & combined_DE_dds[,"species.pvalue"] > 0.05 ] = "cis-trans (compensatory)"
class.compmap[ is.na(class.compmap) & combined_DE_dds[,"trans_effects.pvalue.compmap"] < 0.05 & ((combined_DE_dds[,"alleles.logFC.compmap"] > 0 & combined_DE_dds[,"species.logFC"] < 0) | (combined_DE_dds[,"alleles.logFC.compmap"] < 0 & combined_DE_dds[,"species.logFC"] > 0)) ] = "cis-trans x (compensatory)"
class.compmap[ is.na(class.compmap) & combined_DE_dds[,"trans_effects.pvalue.compmap"] < 0.05 & ((combined_DE_dds[,"alleles.logFC.compmap"] > 0 & combined_DE_dds[,"species.logFC"] > 0) | (combined_DE_dds[,"alleles.logFC.compmap"] < 0 & combined_DE_dds[,"species.logFC"] < 0)) ] = "cis-trans (enhancing)"
class.compmap[ is.na(class.compmap) ] = "ambiguous"
combined_DE_dds = cbind(combined_DE_dds, class.real=as.character(class.real), class.compmap=as.character(class.compmap))

# data wrangling

# counts
combined_DE_dds %>%
    select(class.real, class.compmap) %>%
    pivot_longer(cols=everything(), names_to="met", values_to="class") %>%
    extract(met, "met", "\\.(.*)") %>%
    group_by(met, class) %>% summarize(counts=n()) -> combined_DE_counts

# tmp with class
tmp = combined_DE_dds %>%
    select(class.real, class.compmap) %>%
    pivot_longer(cols=everything(), names_to="method", values_to="class") %>%
    select(class)

combined_DE_dds %>%
    select(species.logFC, alleles.logFC.real, alleles.logFC.compmap, class.real, class.compmap) %>%
    pivot_longer(cols=c(-species.logFC, -class.real, -class.compmap), names_to="method", values_to="alleles.logFC") %>%
    extract(method, "method", "\\.(real|compmap)$") %>%
    add_column(tmp) %>%
    select(-class.real, -class.compmap) %>%
    select(species.logFC, alleles.logFC, class, method) -> combined_DE_logFC

tmp = combined_DE_dds %>%
    select(class.real) %>%
    summarize(match=sum(class.compmap == class.real), total=n()) %>%
    mutate(error.match=abs(match-total)/total)

combined_DE_dds %>%
    select(class.real, class.compmap) %>%
    group_by(class.real) %>%
    summarize(match=sum(class.compmap == class.real), total=n()) %>%
    mutate(error.match=abs(match-total)/total) %>%
    mutate(error.type="type 1") %>%
    rename(class.real = "class") %>%
    mutate_if(is.factor, as.character) -> combined_DE_class_error_type1

combined_DE_dds %>%
    select(class.real, class.compmap) %>%
    group_by(class.compmap) %>%
    summarize(match=sum(class.compmap == class.real), total=n()) %>%
    mutate(error.match=abs(match-total)/total) %>%
    mutate(error.type="type 2") %>%
    rename(class.compmap = "class") %>%
    mutate_if(is.factor, as.character) -> combined_DE_class_error_type2

combined_DE_class_error = rbind(combined_DE_class_error_type1, combined_DE_class_error_type2)

ggplot(combined_DE_class_error, aes(x=class, y=error.match, color=error.type)) +
    geom_point(position=position_dodge(0.5), size=4) +
    geom_linerange(aes(ymax=error.match, ymin=0), position=position_dodge(0.5), size=2) +
    scale_y_continuous(labels = scales::percent, limits=c(0,0.1), expand=c(0,0)) +
    scale_color_manual(values=c("darkblue","darkgreen"), name="error type") +
    labs(x="", y="error rate") +
    background_grid(major="y") +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

ggsave("cis_and_trans_error_rate.png")
ggsave("cis_and_trans_error_rate.pdf", device="pdf", useDingbats=F)

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
p2 = ggplot(combined_DE_counts, aes(x=class, y=counts, fill=class)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=counts), nudge_y=20) +
    facet_rep_wrap(~met) +
    background_grid(major="y") +
    labs(x="", y="number of genes") +
    scale_fill_manual(values=cols.type) +
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

plot_grid(p1, p2, nrow=2, rel_heights=c(0.4,0.6))
ggsave("cis_trans_allele.vs.species_logFC.png")
ggsave("cis_trans_allele.vs.species_logFC.pdf", device="pdf", useDingbats=F)
