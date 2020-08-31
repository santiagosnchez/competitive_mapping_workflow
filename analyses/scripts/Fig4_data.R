library(limma)
library(Glimma)
library(edgeR)
library(MASS)
library(car)
library(ggplot2)
library(cowplot)
library(lemon)
library(ggridges)
library(tidyr)
library(dplyr)
library(stringr)
library(gridExtra)
library(RColorBrewer)
library(ggsci)
library(parallel)
library(WGCNA)
library(raster) # distance between points
library(eulerr) # for venn diagrams
theme_set(theme_cowplot())

# counts

counts = read.csv("counts_both_AS_final.csv", row.name=1)
counts.addH = cbind(counts[,1:6], counts[,7:12] + counts[,19:24], counts[,13:18])
colnames(counts.addH) = sub("_c[B9]","", colnames(counts.addH))
groups = gsub("1|2|3|c","", colnames(counts.addH))
species = sapply(groups, function(x) substr(x,1,1) )
species = sub("9","Cni",species)
species = sub("B","Cbr",species)
species = sub("H","F1",species)
sex = sapply(groups, function(x) substr(x,2,2 ) )

# log-transformed expression data

e = DGEList(counts.addH)
e$samples$group = factor(groups)
e = calcNormFactors(e, method = "TMM")
keep.exprs = filterByExpr(e, as.factor(groups))
e = e[keep.exprs,, keep.lib.sizes=FALSE]
e = estimateCommonDisp(e)
e = estimateGLMTrendedDisp(e)
e = estimateTagwiseDisp(e)

###########################################
### expression divergence species * sex ###
###########################################

## C. briggsae vs C. nigoni + sex models

e.briggsae.vs.nigoni.sex = e$counts[,c(1:6,13:18)]
e.briggsae.vs.nigoni.sex = DGEList(e.briggsae.vs.nigoni.sex)
e.briggsae.vs.nigoni.sex$samples$group = factor(sub("1|2|3","",rownames(e.briggsae.vs.nigoni.sex$samples)))
e.briggsae.vs.nigoni.sex = calcNormFactors(e.briggsae.vs.nigoni.sex, method = "TMM")
e.briggsae.vs.nigoni.sex = estimateCommonDisp(e.briggsae.vs.nigoni.sex)
e.briggsae.vs.nigoni.sex = estimateGLMTrendedDisp(e.briggsae.vs.nigoni.sex)
e.briggsae.vs.nigoni.sex = estimateTagwiseDisp(e.briggsae.vs.nigoni.sex)
design.briggsae.vs.nigoni.sex = model.matrix(~species[c(1:6,13:18)] * sex[c(1:6,13:18)])
colnames(design.briggsae.vs.nigoni.sex)[2:4] = c("nigoni","male","nigoni:male")
v.briggsae.vs.nigoni.sex = voom(e.briggsae.vs.nigoni.sex, design.briggsae.vs.nigoni.sex, plot=F)
vfit.briggsae.vs.nigoni.sex = lmFit(v.briggsae.vs.nigoni.sex, design.briggsae.vs.nigoni.sex)
vfit.briggsae.vs.nigoni.sex = eBayes(vfit.briggsae.vs.nigoni.sex)
dt.briggsae.vs.nigoni.sex = decideTests(vfit.briggsae.vs.nigoni.sex, p.value = 0.05)
tt.briggsae.vs.nigoni.sex = topTable(vfit.briggsae.vs.nigoni.sex, n = Inf)
tt.briggsae.vs.nigoni.sex = tt.briggsae.vs.nigoni.sex[ rownames(dt.briggsae.vs.nigoni.sex), ]
write.csv(tt.briggsae.vs.nigoni.sex, file="tt.briggsae.vs.nigoni.sex.csv", quote=F)

# consolidate
# expression divergence per sex between parents and F1

models.species_sex = paste(dt.briggsae.vs.nigoni.sex[,2], dt.briggsae.vs.nigoni.sex[,3], dt.briggsae.vs.nigoni.sex[,4])

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
rownames(df.species_sex_class) = rownames(dt.briggsae.vs.nigoni.sex)
df.species_sex_class$models = paste(df.species_sex_class[,1], df.species_sex_class[,2], df.species_sex_class[,3])
df.species_sex_class = df.species_sex_class[ df.species_sex_class[,1] != "ambiguous", ]
df.species_sex_class$models = factor(df.species_sex_class$models, models.levels)
df.species_sex_class$logFC_species = tt.briggsae.vs.nigoni.sex[rownames(df.species_sex_class),1]
df.species_sex_class$logFC_sex = tt.briggsae.vs.nigoni.sex[rownames(df.species_sex_class),2]

# add dnds data
dnds = read.csv("dnds_propcons_position_domain.csv", row.names=1)
df.species_sex_class$Ka = dnds[ rownames(df.species_sex_class),"dn"]
df.species_sex_class$Ks = dnds[ rownames(df.species_sex_class),"ds"]
df.species_sex_class$Ks_ENCc = dnds[ rownames(df.species_sex_class),"ds_ENCc"]

# add cis-trans and expr inheritance

expr_inherit = read.csv("expr_inheritance_logFC_F1_parents.csv")
cis_trans = read.csv("df.cis_trans.regulation_type.csv")
class.levels = c("no change","ambiguous","additive","C. briggsae dominant", "C. nigoni dominant","overdominant","underdominant")
type.levels = c("conserved","ambiguous","trans only","cis only","cis-trans (enhancing)","cis-trans (compensatory)")

x = subset(expr_inherit, sex == "female")
rownames(x) = as.character(x$gene)
df.species_sex_class$class_female = factor(x[ rownames(df.species_sex_class), "class" ], class.levels)
x = subset(expr_inherit, sex == "male")
rownames(x) = as.character(x$gene)
df.species_sex_class$class_male = factor(x[ rownames(df.species_sex_class), "class" ], class.levels)

x = subset(cis_trans, sex == "female")
rownames(x) = as.character(x$genes)
df.species_sex_class$type_female = factor(x[ rownames(df.species_sex_class), "type" ], type.levels)
x = subset(cis_trans, sex == "male")
rownames(x) = as.character(x$genes)
df.species_sex_class$type_male = factor(x[ rownames(df.species_sex_class), "type" ], type.levels)
models.levels = sort(unique(as.character(df.species_sex_class$models)))[c(8,5,13,7,4,12,3,11,6,2,10,1,9)]
df.species_sex_class$models = factor(df.species_sex_class$models, models.levels)
df.species_sex_class$models_num = paste0("M",1:length(models.levels))[ df.species_sex_class$models ]

write.csv(df.species_sex_class, file="df.species_sex_class.csv")

df.species_sex_class.expr_class_type = data.frame(sex=rep(c("female","male"), each=dim(df.species_sex_class)[1]),
												  class=factor(c(as.character(df.species_sex_class$class_female), as.character(df.species_sex_class$class_male)), class.levels),
												  type=factor(c(as.character(df.species_sex_class$type_female), as.character(df.species_sex_class$type_male)), type.levels),
												  models=rep(as.character(df.species_sex_class$models), 2),
												  gene=rep(rownames(df.species_sex_class), 2))
df.species_sex_class.expr_class_type$chromosome = sub("\\..*","", df.species_sex_class.expr_class_type$gene)
df.species_sex_class.expr_class_type$models = factor(df.species_sex_class.expr_class_type$models, models.levels)
df.species_sex_class.expr_class_type$models_num = paste0("M",1:length(models.levels))[ df.species_sex_class.expr_class_type$models ]

df.species_sex_class.expr_class.counts = as.data.frame(df.species_sex_class.expr_class_type %>% group_by(sex, models, models_num, class) %>% count())
df.species_sex_class.expr_type.counts = as.data.frame(df.species_sex_class.expr_class_type %>% group_by(sex, models, models_num, type) %>% count())
df.species_sex_class.expr_class.counts$class = factor(df.species_sex_class.expr_class.counts$class, class.levels)
df.species_sex_class.expr_type.counts$type = factor(df.species_sex_class.expr_type.counts$type, type.levels)

write.csv(df.species_sex_class.expr_type.counts, file="df.species_sex_class.expr_type.counts.csv", row.names=F)
write.csv(df.species_sex_class.expr_class.counts, file="df.species_sex_class.expr_class.counts.csv", row.names=F)


# heatmap of species * sex models
df.models = unique(df.species_sex_class[, 1:3])
rownames(df.models) = 1:dim(df.models)[1]
df.models$models = paste(df.models[,1], df.models[,2], df.models[,3])
df.models2 = gather(df.models, key, value, -models)
df.models2$counts = table(df.species_sex_class$models)[ df.models2$models ]
df.models2$sex = sapply(strsplit(as.character(df.models2$models), " "), function(x) x[2] )
df.models2$models = factor(df.models2$models, models.levels)
df.models2$models_num = paste0("M",1:length(models.levels))[ df.models2$models ]

write.csv(df.models2, file="species_sex_models_heatmap.csv", row.names=F)
