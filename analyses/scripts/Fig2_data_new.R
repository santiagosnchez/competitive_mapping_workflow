# generates data/analyses for Fig. 2

library(raster)
library(tidyr)
library(dplyr)

theme_set(theme_cowplot())
source("scripts/enrichment_fun.R")

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

# expression distance F1 from zero
df.inherit2$F1_dist_from_zero = pointDistance(as.matrix(df.inherit2[,c("logFC.F1.vs.briggsae","logFC.F1.vs.nigoni")]),
																							matrix(0, ncol=2, nrow=dim(df.inherit2)[1]),
																							lonlat=F, allpairs=F)
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

# matrix fig 3

class.vs.type.mat = na.omit(df.inherit_and_ase) %>%
    filter(!(class == "no change" | class == "ambiguous" | type == "conserved" | type == "ambiguous")) %>%
    type.convert() %>%
    group_by(sex, class, type, chromosome, .drop=F) %>%
    count()
write.csv(class.vs.type.mat, file="tables/df.cis_trans.inherit.mat.csv", row.names=F)
