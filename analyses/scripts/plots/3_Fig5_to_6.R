# libraries
library(ggplot2)
library(cowplot)
library(lemon)
library(ggridges)
library(ggstance)
library(ggsci)
library(ggforce)
library(tidyr)
library(dplyr)
library(gtable)
library(tibble)
library(ggplotify)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(viridis)
library(scales)
library(venneuler)
theme_set(theme_cowplot())

source("scripts/enrichment_fun.R")

# colors

# brewer.pal(10, "Paired") -> cols
# cols = cols[-c(7,8)]
# cols.spp = cols[c(6,2,8)]
# cols.sex = c("#80CDC1","#018571","#BEBEBE")
# cols.heatmap = c(cols.spp[1],"grey",cols.sex[1],"black",cols.sex[2],cols.spp[2],"grey",cols.sex[3])
# cols.class = c("black","grey",cols[ c(4,6,2,8,7) ])
# cols.type = c("black","grey",pal_jama()(5))
# cols.type[3] = "#09BC8A"
# cols.type[7]= "#C49991"
# cols.type = cols.type[c(1,2,4,3,5,6,7)]
# cols.type[4] = "#476A6F"
# cols.heatmap1 = cols.heatmap[c(1,2,2)]

#cols.inherit_old = rev(c("#A68481","#4F3B39","#6F7862","#B379AF","#4C5E91"))
cols.inherit = c("#4C5E91","#B379AF","#3F784C","#C17817","#533E2D")
cols.class = c("black","grey",cols.inherit)
# brewer.pal(10, "Paired") -> cols
# cols = cols[-c(7,8)]
# cols.class = c("black","grey",cols[ c(4,6,2,8,7) ])
cols.type1 = c("#D4C41C","#39B8E3","#3B9133","#A3627D","#DBA9BE")
cols.type = c("black","grey",cols.type1)
# cols sex
cols.sex = c("#80CDC1","#018571","#A68481","#BEBEBE")
# heatmap
cols.heatmap = c(cols.inherit[1],"grey",cols.sex[1],"black",cols.sex[3],cols.inherit[2],"grey",cols.sex[4])

########################
# Supplementary Fig. X #
########################

# reaction norms

df.models.cent = read.csv("tables/species_by_sex/df.centroides_and_CI.species_sex_reaction_norms.csv")

models.levels = c("C-N-N","B-N-N","N-N-N","C-M-N","B-M-N","B-M-I","N-M-N","N-M-I","C-F-N","B-F-N","B-F-I","N-F-N","N-F-I")
df.models.cent$models = factor(df.models.cent$models, models.levels)

cons = ggplot(filter(df.models.cent, models %in% c("C-N-N","B-N-N","N-N-N")),
        aes(x=species, y=expression, group=sex)) +
    geom_ribbon(aes(ymax=upper, ymin=lower), fill="grey", alpha=0.5) +
    geom_line(aes(color=sex), size=1.5) +
    facet_rep_grid("sex-neutral"~models) +
    scale_y_continuous(breaks=c(-1,0,1)) +
    scale_x_discrete(expand=c(0.1,0.1)) +
    labs(x="",y="") +
    #background_grid(major="y") +
    theme(strip.background = element_rect(fill="grey90", color=NA)) +
    theme(axis.text.x=element_text(face="italic", size=10)) +
    scale_color_manual(values=cols.sex[c(1,3)])

male_biased = ggplot(filter(df.models.cent, models %in% c("C-M-N","B-M-N","B-M-I","N-M-N","N-M-I")),
        aes(x=species, y=expression, group=sex)) +
    geom_ribbon(aes(ymax=upper, ymin=lower), fill="grey", alpha=0.5) +
    geom_line(aes(color=sex), size=1.5, show.legend=F) +
    facet_rep_grid("male-biased"~models) +
    scale_y_continuous(breaks=c(-1,0,1)) +
    scale_x_discrete(expand=c(0.1,0.1)) +
    labs(x="",y="normalized expression") +
    #background_grid(major="y") +
    theme(strip.background = element_rect(fill="grey90", color=NA)) +
    theme(axis.text.x=element_text(face="italic", size=10)) +
    scale_color_manual(values=cols.sex[c(1,3)])

female_biased = ggplot(filter(df.models.cent, models %in% c("C-F-N","B-F-N","B-F-I","N-F-N","N-F-I")),
        aes(x=species, y=expression, group=sex)) +
    geom_ribbon(aes(ymax=upper, ymin=lower), fill="grey", alpha=0.5) +
    geom_line(aes(color=sex), size=1.5, show.legend=F) +
    facet_rep_grid("female-biased"~models) +
    scale_y_continuous(breaks=c(-1,0,1)) +
    scale_x_discrete(expand=c(0.1,0.1)) +
    labs(x="",y="") +
    #background_grid(major="y") +
    theme(strip.background = element_rect(fill="grey90", color=NA)) +
    theme(axis.text.x=element_text(face="italic", size=10)) +
    scale_color_manual(values=cols.sex[c(1,3)])

row1 = plot_grid(cons + theme(plot.margin = unit(c(t = 0, r = 9, b = 0, l = 0), "cm")),
        male_biased + theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "cm")),
        female_biased + theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "cm")),
        nrow=3)

row1
ggsave("figures/Fig5a_species_sex_reaction_norms.pdf", device="pdf", useDingbats=F)
ggsave("figures/Fig5a_species_sex_reaction_norms.png")

# read class count and expression data

############
# figure 5 #
############

df.species_sex_class = read.csv("tables/species_by_sex/df.species_sex.csv", row.names=1)
df.species_sex_class.expr_type.counts = read.csv("tables/species_by_sex/df.species_sex_class.expr_type.counts.csv")
df.species_sex_class.expr_class.counts = read.csv("tables/species_by_sex/df.species_sex_class.expr_class.counts.csv")
df.models2 = read.csv("tables/species_by_sex/species_sex_models_heatmap.csv")

df.species_sex_class.expr_class.counts$class = factor(df.species_sex_class.expr_class.counts$class, class.levels)
df.species_sex_class.expr_type.counts$type = factor(df.species_sex_class.expr_type.counts$type, type.levels)
df.species_sex_class.expr_class.counts$models2 = factor(df.species_sex_class.expr_class.counts$models2, models.levels)
df.species_sex_class.expr_type.counts$models2 = factor(df.species_sex_class.expr_type.counts$models2, models.levels)
df.species_sex_class.expr_class.counts$sex2 = sapply(strsplit(as.character(df.species_sex_class.expr_class.counts$models), " "), function(x) x[2] )
df.species_sex_class.expr_type.counts$sex2 = sapply(strsplit(as.character(df.species_sex_class.expr_type.counts$models), " "), function(x) x[2] )
df.species_sex_class.expr_class.counts$sex2 = factor(df.species_sex_class.expr_class.counts$sex2, c("sex-neutral","male-biased","female-biased"))
df.species_sex_class.expr_type.counts$sex2 = factor(df.species_sex_class.expr_type.counts$sex2, c("sex-neutral","male-biased","female-biased"))
#df.species_sex_class$models = factor(df.species_sex_class$models, models.levels)
df.species_sex_class$models2 = factor(df.species_sex_class$models2, models.levels)
df.species_sex_class$sex = factor(df.species_sex_class$sex, c("sex-neutral","male-biased","female-biased"))
df.models2$models2 = factor(df.models2$models2, models.levels)
df.models2$sex = factor(df.models2$sex, c("sex-neutral","male-biased","female-biased"))
models.counts = unique(as.character(df.models2$counts))[ order(unique(df.models2$models)) ]
df.models2$counts = factor(as.character(df.models2$counts), models.counts)

# plots

models_grid = ggplot(df.models2, aes(x=key, y=counts, fill=value)) +
	geom_tile(color="black", show.legend=F) +
	scale_fill_manual(values=scales:::alpha(cols.heatmap, 0.7)) +
	scale_x_discrete(limits=c("species","sex","interaction"), position="top") +
	#scale_y_discrete(labels=unique(as.character(df.models2$counts))[ order(unique(df.models2$models)) ]) +
	xlab("") + ylab("") +
    facet_rep_wrap(~sex, strip.position="right", drop=T, scales="free_y", ncol=1) +
    theme(strip.background = element_rect(fill="grey90", color=NA)) +
	theme(
		axis.line=element_blank(),
		axis.ticks=element_blank(),
		axis.text.x=element_text(angle=90, hjust=0, vjust=0.5))

expr_div = ggplot(df.species_sex_class, aes(y=models2, x=abs(logFC_species))) +
	geom_density_ridges(aes(fill=sex), scale=2, panel_scaling=F, color=NA, alpha=0.7, show.legend=F) +
	geom_boxploth(alpha=0.5, outlier.shape=NA, width=0.5) +
	scale_fill_manual(values=rev(cols.sex[c(1,3,4)])) +
	scale_x_continuous(position="top", limits=c(0,10), labels=c("0","","5","","10")) +
	#scale_y_discrete(labels=models.levels2) +
	xlab(expression(paste("log"[2]," expression divergence"))) +
	ylab("") +
    facet_rep_wrap(~sex, strip.position="right", drop=T, scales="free_y", ncol=1) +
	background_grid(major="x", size.major = 0.2, colour.major = "grey75", minor="x", size.minor = 0.2, colour.minor = "grey75") +
	theme(strip.background=element_blank(), strip.text=element_blank())

aa_div = ggplot(df.species_sex_class, aes(y=models2, x=Ka/Ks_ENCc)) +
	geom_density_ridges(aes(fill=sex), scale=1, panel_scaling=F, color=NA, alpha=0.7, show.legend=F) +
	geom_boxploth(alpha=0.5, outlier.shape=NA, width=0.5) +
	scale_fill_manual(values=rev(cols.sex[c(1,3,4)])) +
	scale_x_continuous(position="top", limits=c(0,0.5), breaks=c(0,0.2,0.4,0.6), labels=c("0","0.2","0.4","0.6")) +
	#scale_y_discrete(limits=rev(models.levels)) +
	xlab(expression(paste(italic("K")[a],"/",italic("K")[s],"\'"))) +
	ylab("") +
    facet_rep_wrap(~sex, strip.position="right", drop=T, scales="free_y", ncol=1) +
	background_grid(major="x", size.major = 0.2, colour.major = "grey75", minor="x", size.minor = 0.2, colour.minor = "grey75") +
    theme(strip.background=element_blank(), strip.text=element_blank()) +
	theme(axis.text.y=element_blank())

prop_class = ggplot(df.species_sex_class.expr_class.counts, aes(x=models2, y=n, fill=class)) +
	geom_bar(stat="identity", position="fill", show.legend=F, alpha=0.7) +
	scale_y_continuous(position="right", labels=c("0","","0.5","","1")) +
	#scale_x_discrete(limits=rev(models.levels)) +
	scale_fill_manual(values=cols.class) +
	ylab("proportion of genes") + xlab("") +
	coord_flip() +
    facet_rep_grid(sex2~sex, drop=T, scales="free_y") +
	background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
	theme(strip.background.x = element_rect(fill="grey90", color=NA)) +
    theme(strip.background.y=element_blank(), strip.text.y=element_blank()) +
	theme(axis.text.y=element_blank())

prop_type = ggplot(df.species_sex_class.expr_type.counts, aes(x=models2, y=n, fill=type)) +
	geom_bar(stat="identity", position="fill", show.legend=F, alpha=0.7) +
	scale_y_continuous(position="right", labels=c("0","","0.5","","1")) +
	#scale_x_discrete(limits=rev(models.levels)) +
	scale_fill_manual(values=cols.type) +
	facet_rep_grid(sex2~sex, drop=T, scales="free_y") +
	ylab("proportion of genes") + xlab("") +
	coord_flip() +
	background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
	theme(strip.background.x = element_rect(fill="grey90", color=NA)) +
    theme(strip.background.y=element_blank(), strip.text.y=element_blank()) +
	theme(axis.text.y=element_blank())

row2 = plot_grid(models_grid + theme(plot.margin = unit(c(t=-0.3, r=-0.1, b=0, l=-0.2), "cm")),
 					# expr_div + theme(plot.margin = unit(c(t=0.2, r=0, b=0, l=-0.3), "cm")),
					# aa_div + theme(plot.margin = unit(c(t=0.15, r=0, b=0, l=-0.1), "cm")),
					prop_class + theme(plot.margin = unit(c(t=0.3, r=0, b=0, l=-0.2), "cm")),
					prop_type + theme(plot.margin = unit(c(t=0.3, r=0, b=0, l=-0.2), "cm")),
					rel_widths=c(0.15,0.4,0.4), nrow=1)

# row2 = plot_grid(models_grid + theme(plot.margin = unit(c(t=-0.3, r=-0.1, b=0, l=-0.2), "cm")),
#  					expr_div + theme(plot.margin = unit(c(t=0.2, r=0, b=0, l=-0.3), "cm")),
# 					aa_div + theme(plot.margin = unit(c(t=0.15, r=0, b=0, l=-0.1), "cm")),
# 					prop_class + theme(plot.margin = unit(c(t=0.3, r=0, b=0, l=-0.2), "cm")),
# 					prop_type + theme(plot.margin = unit(c(t=0.3, r=0, b=0, l=-0.2), "cm")),
# 					rel_widths=c(0.15,0.2,0.18,0.4,0.4), nrow=1)

row2
ggsave("figures/Fig5b_species_sex_expr_div_dnds_counts.pdf", device="pdf", useDingbats=F)
ggsave("figures/Fig5b_species_sex_expr_div_dnds_counts.png")

tmp1 =  as.data.frame(df.species_sex_class %>%
    group_by(sex, models2, chromosome) %>%
    summarize(n=n()) %>%
    mutate(region=c("auto","X")[ factor(chromosome == "X")]) %>%
    group_by(sex, models2, region) %>%
    summarise(n=sum(n)))
x = enrichment(t(matrix(tmp1$n, ncol=2, nrow=13, byrow=T)), odds.ratio=T)
y = enrichment(t(matrix(tmp1$n, ncol=2, nrow=13, byrow=T)), odds.ratio=F)
tmp1 = as.data.frame(tmp1 %>%
    add_column(enrichment=as.vector(x)) %>%
    add_column(p.value=as.vector(y)) %>%
    group_by(sex, models2) %>%
    summarize(region, n, enrichment, p.value, rel_enrich=enrichment/sum(enrichment)))
tmp2 =  as.data.frame(df.species_sex_class %>%
    filter(!(domain == "tipA" | domain == "tipB")) %>%
    mutate(region = gsub("A$|B$","", as.character(domain))) %>%
    group_by(sex, models2, region) %>%
    summarize(n=n()))
x = enrichment(t(matrix(tmp2$n, ncol=2, nrow=13, byrow=T)), odds.ratio=T)
y = enrichment(t(matrix(tmp2$n, ncol=2, nrow=13, byrow=T)), odds.ratio=F)
tmp2 = as.data.frame(tmp2 %>%
    add_column(enrichment=as.vector(x)) %>%
    add_column(p.value=as.vector(y)) %>%
    group_by(sex, models2) %>%
    summarize(region, n, enrichment, p.value, rel_enrich=enrichment/sum(enrichment)))

df.enrich_arms_X2 = rbind(tmp1, tmp2)
df.enrich_arms_X2$models2 = factor(df.enrich_arms_X2$models2, models.levels2)
df.enrich_arms_X2$sex = factor(df.enrich_arms_X2$sex, c("sex-neutral","male-biased","female-biased"))

arms_X = ggplot(df.enrich_arms_X2 %>% filter(region == "arm" | region == "X"), aes(x=log2(enrichment), y=models2, fill=sex)) +
    geom_bar(stat="identity", show.legend=F, alpha=0.7) +
    scale_x_continuous(expand=c(0,0), position="top") +
    facet_rep_grid(sex~region, scales="free_y") +
    scale_fill_manual(values=rev(cols.sex[c(1,3,4)])) +
    xlab(expression(paste("enrichment [log"[2]," odds ratio]",sep=""))) +
    ylab("") +
    theme(axis.text.y=element_blank()) +
    background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
    theme(strip.background = element_rect(fill="grey90", color=NA))

row3 = plot_grid(expr_div + theme(plot.margin = unit(c(t=0.2, r=0, b=0, l=-0.3), "cm")),
					aa_div + theme(plot.margin = unit(c(t=0.15, r=0, b=0, l=-0.1), "cm")),
                    arms_X + theme(plot.margin = unit(c(t=0.3, r=0, b=0, l=-0.2), "cm")),
					rel_widths=c(0.25,0.25,0.5), nrow=1)
row3

plot_grid(row2, row3, ncol=1)
ggsave("figures/Fig5b2_species_sex_expr_div_dnds_counts_pre.pdf", device="pdf", useDingbats=F)
#ggsave("figures/Fig5b2_species_sex_expr_div_dnds_counts.png")

############
# figure 6 #
############

# spermatogenesis / hermaphrodite genes

herm_genes = scan("tables/Thomas_et_al_herm_genes.txt", what=character(), sep="\n")
spermatogenesis_genes = scan("tables/spermatogenesis_orthologs_Ma_et_al.txt", what=character(), sep="\n")
BMI_genes = rownames(df.species_sex_class)[ df.species_sex_class$models2 == "B-M-I"]

HS = sum(herm_genes %in% spermatogenesis_genes)
BMIH = sum(BMI_genes %in% herm_genes)
BMIS = sum(BMI_genes %in% spermatogenesis_genes)
BMIHS = sum(BMI_genes %in% spermatogenesis_genes & BMI_genes %in% herm_genes)
BMI = length(BMI_genes) - (BMIH-BMIHS + BMIS-BMIHS + BMIHS)
H = length(herm_genes) - (BMIH-BMIHS + HS-BMIHS + BMIHS)
S = length(spermatogenesis_genes) - (BMIS-BMIHS + HS-BMIHS + BMIHS)
venn = venneuler(c(H=H, S=S, "H&S"=HS, "BMI&H"=BMIH, "BMI&S"=BMIS, "BMI&H&S"=BMIHS))
df.venn = data.frame(venn$centers, diameters=venn$diameters, r=venn$diameters/2, label=c("hermaphrodite","spermatogenesis","B-M-I"), number=c(H,S,BMI))
df.venn = rbind(df.venn, data.frame(x=0.5, y=0.5, diameters=NA, r=NA, label=NA, number=HS))

x = locator()
labels = as.data.frame(x)
labels$text = c("S","H","BMI","BMI&S","S&H","H&BMI","BMI&H&S")

herm_logFC = read.csv("tables/Thomas_et_al_DESeq_logFC.csv")
herm_logFC$herm = as.character(herm_logFC$herm)
herm_logFC$herm = sub("hermaphroditic","hermaphrodite", herm_logFC$herm)
sperm_herm_models.enrich = read.csv("tables/df.sperm_herm.enrich.csv", comment.char="#")
sperm_herm_models.enrich$sex = c("female-biased","male-biased","sex-neutral")[ factor(sapply(strsplit(as.character(sperm_herm_models.enrich$models2), "-"), function(x) x[2] )) ]
sperm_herm_models.enrich$sig = ""
sperm_herm_models.enrich$sig[ sperm_herm_models.enrich$p.value < 0.05 & log2(sperm_herm_models.enrich$enrichment) > 0.5 ] = "*"
sperm_herm_models.enrich$models2 = factor(sperm_herm_models.enrich$models2, models.levels)
sperm_herm_models.enrich$sex = factor(sperm_herm_models.enrich$sex, c("sex-neutral","male-biased","female-biased"))

# load expression inheritance data
class.levels = c("no change","ambiguous","C. briggsae dominant", "C. nigoni dominant","additive","overdominant","underdominant")
class.legend = c("no change","ambiguous",expression(paste(italic("C. briggsae"), " dominant")),
                expression(paste(italic("C. nigoni"), " dominant")),"additive","overdominant","underdominant")
sex.levels = c("female","male","(F) hermaphrodite","hermaphrodite")
df.inherit2 = read.csv("tables/expression_inheritance/expr_inheritance_logFC_F1_parents.csv")
df.inherit3 = read.csv("tables/expression_inheritance/expr_inheritance_chr_enrich.csv")
df.inherit2$class = factor(df.inherit2$class, levels=class.levels)
df.inherit3$class = factor(df.inherit3$class, levels=class.levels)
df.inherit2$sex = factor(df.inherit2$sex, levels=sex.levels)
df.inherit3$sex = factor(df.inherit3$sex, levels=sex.levels)

# load cis-trans regulatory data
type.levels = c("conserved","ambiguous","trans-only","cis-only", "cis + trans (enhancing)", "cis x trans (compensatory)", "cis-trans (compensatory)")
type.legend = c("conserved","ambiguous",expression(paste(italic("trans"), "-only")), expression(paste(italic("cis"),"-only")),
                expression(paste(italic("cis + trans"), " (enhancing)")), expression(paste(italic("cis x trans"), " (compensatory)")),
                expression(paste(italic("cis-trans"), " (compensatory)")))
df.AS = read.csv("tables/cis_trans/df.cis_trans.regulation_type.auto.csv")
df.AS2 = read.csv("tables/cis_trans/df.cis_trans.counts.chr_enrich.csv")
df.AS$type = factor(df.AS$type, levels=type.levels)
df.AS2$type = factor(df.AS2$type, levels=type.levels)
df.AS$sex = factor(df.AS$sex, levels=sex.levels)
df.AS2$sex = factor(df.AS2$sex, levels=sex.levels)

logFC_biplot = ggplot(herm_logFC, aes(x=logFC.sex.our, y=logFC.sex.thomas, color=herm)) +
    geom_point(alpha=0.3) +
    scale_x_continuous(limits=c(-10,10)) +
    scale_y_continuous(limits=c(-10,13)) +
    labs(x="logFC male:female (this study)", y="logFC male:female (Thomas et al.)") +
    geom_hline(yintercept=0, linetype=2) +
    geom_vline(xintercept=0, linetype=2) +
    scale_color_manual(values=cols.sex[c(2,4)], name="", guide=guide_legend(override.aes=list(alpha=1))) +
    background_grid(major="xy") +
    theme(legend.position="bottom", axis.title=element_text(size=12))
logFC_biplot

herm_sperm_enrich = ggplot(sperm_herm_models.enrich, aes(x=log2(enrichment), y=models2, fill=sex)) +
    geom_bar(stat="identity", alpha=0.7) +
    scale_fill_manual(values=rev(cols.sex[-2]), name="") +
    labs(y="", x=expression(paste("enrichment (log"[2], " odds ratio)"))) +
    geom_text(aes(x=log2(enrichment)/2, label=sig), size=5, position=position_dodge(width=0.5), show.legend=F, color="white") +
    background_grid(major="x") +
    theme(legend.position="bottom") +
    theme(strip.background = element_rect(fill="grey90", color=NA)) +
    facet_rep_grid(sex~group, scales="free_y")

venn_diag = ggplot(df.venn) +
    geom_circle(aes(x0=x, y0=y, r=r, fill=label), alpha=0.5, color=NA, show.legend=T) +
    scale_fill_manual(values=cols.sex[c(3,2,4)], name="") +
    ylim(0,1.3) +
    #geom_text(data=labels, aes(x=x, y=y, label=text)) +
    theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        axis.ticks=element_blank()) +
    theme(legend.position=c(0.2,0.8))

top = plot_grid(logFC_biplot, herm_sperm_enrich, venn_diag, nrow=1, rel_widths=c(0.3,0.4,0.3))

############################
## expression inheritance ##
############################

# biplot of expression difference between parents and F1
# colored by expression inheritance class
pA = ggplot(df.inherit2 %>% filter(grepl("\\(F\\) hermaphrodite", sex)) %>%
            #mutate(sex = sub("\\(M\\) hermaphrodite","male", sex)) %>%
            mutate(sex = sub("\\(F\\) ","", sex)),
                aes(x=logFC.F1.vs.briggsae, y=logFC.F1.vs.nigoni, color=class)) +
	geom_point(show.legend=F, alpha=0.3) +
	scale_color_manual(values=cols.class) +
	geom_hline(yintercept=0, linetype=2) +
	geom_vline(xintercept=0, linetype=2) +
	xlab(expression(paste("log"[2],"(F1/",italic("Cbri"),")", sep=''))) + # bold("expr. diff."),
	ylab(expression(paste("log"[2],"(F1/",italic("Cnig"),")", sep=''))) + # bold("expr. diff."),
	background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
	#facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	theme(strip.background = element_rect(fill="grey90", color=NA)) +
	theme(strip.background = element_blank(), strip.text=element_blank(),
				axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))
#p1

# F1 expression distance from origin (0,0)
pB = ggplot(df.inherit2 %>% filter(grepl("\\(F\\) hermaphrodite", sex)) %>%
            #mutate(sex = sub("\\(M\\) hermaphrodite","male", sex)) %>%
            mutate(sex = sub("\\(F\\) ","", sex)),
                aes(y=class, x=F1_dist_from_zero, fill=class)) +
	geom_boxploth(alpha=0.7, outlier.shape=NA, width=0.5, show.legend=F) +
	geom_density_ridges(alpha=0.7, color=NA, show.legend=F) +
	scale_y_discrete(limits=rev(class.levels)) +
	scale_fill_manual(values=cols.class) +
	#facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	xlab("F1 expression distance") +
	ylab("") +
	background_grid(major="x", minor="x", size.major = 0.2, colour.major = "grey75", size.minor = 0.2, colour.minor = "grey75") +
	xlim(0,7) +
	theme(strip.background = element_blank(), strip.text=element_blank(),
	 			axis.text.y=element_blank(),
	 			axis.ticks.y=element_blank())

# stacked barplot with gene counts for each class per chromosome
pC = ggplot(df.inherit3 %>% filter(grepl("\\(F\\) hermaphrodite", sex)) %>%
            #mutate(sex = sub("\\(M\\) hermaphrodite","male", sex)) %>%
            mutate(sex = sub("\\(F\\) ","", sex)),
                aes(x=chromosome, y=n, fill=class)) +
	geom_bar(stat="identity", show.legend=F, alpha=0.7) +
    #scale_y_continuous(breaks=seq(0,3000,500)) +
	scale_fill_manual(values=cols.class) +
	ylab("number of genes") +
	#facet_rep_wrap(~sex, ncol=1, strip.position="right", scales="free") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	theme(strip.background = element_rect(fill="grey90", color=NA)) +
	theme(strip.background = element_blank(), strip.text=element_blank())
#p2

# for pD
# grey shade for each chromosome
rectsD = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-3,2,-3), y2=c(2,-3,2))

# gene enrichment of each expression inheritance class
pD = ggplot(df.inherit3 %>% filter(grepl("\\(F\\) hermaphrodite", sex)) %>%
            #mutate(sex = sub("\\(M\\) hermaphrodite","male", sex)) %>%
            mutate(sex = sub("\\(F\\) ","", sex)),
                aes(x=chr.pseudo, y=log2(enrichment), fill=class)) +
 	geom_rect(data=rectsD, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
	geom_bar(stat="identity", position=position_dodge(), alpha=0.7, show.legend=T) +
	scale_y_continuous(breaks=seq(-3,2,1)) +
	scale_x_continuous(breaks=1:6, labels=df.inherit3$chromosome[1:6]) +
	scale_fill_manual(values=cols.class, labels=class.legend) +
	#facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	ylab(expression(paste("enrichment [log"[2]," odds ratio]",sep=""))) +
	xlab("chromosome") +
	geom_text(aes(y=log2(enrichment)/2, label=sig), size=5, position=position_dodge(width=0.9), show.legend=F, color="white") +
	theme(strip.background = element_rect(fill="grey90", color=NA),
				legend.title=element_blank(),
				legend.text=element_text(size=10),
                legend.text.align = 0,
				axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

# extract legends to separate plot objects
legend.class = as_grob(get_legend(pD))
pD = pD + theme(legend.position = "none")

# layout plots into a grid
row1 = plot_grid(pA + ggtitle("A") + theme(plot.margin = unit(c(0, 0, 0.2, 0), "cm")),
								pB + ggtitle("B") + theme(plot.margin = unit(c(0, 0.5, 0.2, 0), "cm")),
								pC + ggtitle("C") + theme(plot.margin = unit(c(0, 0.5, 0.2, 0), "cm")),
								pD + ggtitle("D") + theme(plot.margin = unit(c(0, 0.5, 0.2, 0), "cm")),
								legend.class, #+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
								nrow=1, align="both", rel_widths=c(0.4,0.4,0.35,0.7,0.3))

##########################################
## cis trans expression regulation type ##
##########################################

# biplot of allele-specific expression in F1 (y-axis) vs
# expression difference between parents (cis + trans)
pE = ggplot(df.AS %>% filter(grepl("\\(F\\) hermaphrodite", sex)) %>%
            mutate(sex = sub("\\(M\\) hermaphrodite","male", sex)) %>%
            mutate(sex = sub("\\(F\\) ","", sex)),
                aes(logFC.sp, logFC.ase, color=type)) +
	geom_point(alpha=0.3, show.legend=F) +
	scale_color_manual(values=cols.type) +
	geom_hline(yintercept=0, linetype=2) +
	geom_vline(xintercept=0, linetype=2) +
	ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cnig"),"/",italic("Cbri"),")]", sep=''))) +
	xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cnig"),"/",italic("Cbri"),")]", sep=''))) +
	background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
	#facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	theme(strip.background = element_rect(fill="grey90", color=NA)) +
	theme(strip.background = element_blank(), strip.text=element_blank(),
				axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

# boxplots of expression divergence between parent species
# for each regulation type
pF = ggplot(df.AS %>% filter(grepl("\\(F\\) hermaphrodite", sex)) %>%
            mutate(sex = sub("\\(M\\) hermaphrodite","male", sex)) %>%
            mutate(sex = sub("\\(F\\) ","", sex)),
                aes(y=type, x=abs(logFC.sp), fill=type)) +
	geom_boxploth(alpha=0.7, outlier.shape=NA, width=0.5, show.legend=F) +
	geom_density_ridges(alpha=0.7, color=NA, show.legend=F) +
	scale_y_discrete(limits=rev(type.levels)) +
	scale_fill_manual(values=cols.type) +
	#facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	xlab(expression(paste("|log"[2],"-FC expression divergence|", sep=""))) +
	ylab("") +
	background_grid(major="x", minor="x", size.major = 0.2, colour.major = "grey75", size.minor = 0.2, colour.minor = "grey75") +
	xlim(0,7) +
	theme(strip.background = element_blank(), strip.text=element_blank(),
	 			axis.text.y=element_blank(),
	 			axis.ticks.y=element_blank())

# stacked barplot of counts of each regulation type per chromosome
pG = ggplot(df.AS %>% filter(grepl("\\(F\\) hermaphrodite", sex)) %>%
            mutate(sex = sub("\\(M\\) hermaphrodite","male", sex)) %>%
            mutate(sex = sub("\\(F\\) ","", sex)),
                aes(x=chromosome, fill=type)) +
	geom_bar(show.legend=F, alpha=0.7) +
	scale_fill_manual(values=cols.type) +
    #scale_y_continuous(breaks=seq(0,3000,500)) +
	ylab("number of genes") +
	#facet_rep_wrap(~sex, ncol=1, strip.position="right", scales="free_y") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	theme(strip.background = element_rect(fill="grey90", color=NA)) +
	theme(strip.background = element_blank(), strip.text=element_blank())

# for pH
# grey shade for each chromosome
rectsH = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-1.2,1.2,-1.2), y2=c(1.2,-1.2,1.2))

# gene enrichment of cis and trans regulation types in each chromosome
pH = ggplot(df.AS2 %>% filter(grepl("\\(F\\) hermaphrodite", sex)) %>%
            mutate(sex = sub("\\(M\\) hermaphrodite","male", sex)) %>%
            mutate(sex = sub("\\(F\\) ","", sex)),
                aes(x=chr.pseudo, y=log2(enrichment), fill=type)) +
 	geom_rect(data=rectsH, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
	geom_bar(stat="identity", position=position_dodge(), alpha=0.7, show.legend=T) +
	scale_y_continuous(breaks=seq(-1,1,0.5), labels=c("-1","","0","","1")) +
	scale_x_continuous(breaks=1:6, labels=df.AS2$chromosome[1:6]) +
	scale_fill_manual(values=cols.type, labels=type.legend) +
	#facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	ylab(expression(paste("enrichment [log"[2]," odds ratio]",sep=""))) +
	xlab("chromosome") +
	geom_text(aes(y=log2(enrichment)/2, label=sig), size=5, position=position_dodge(width=0.9), show.legend=F, color="white") +
	theme(strip.background = element_rect(fill="grey90", color=NA),
				legend.title=element_blank(),
				legend.text=element_text(size=10),
                legend.text.align = 0,
				axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

# extract legends to separate plot objects
legend.type = as_grob(get_legend(pH))
pH = pH + theme(legend.position = "none")

# layout plots into a grid
row2 = plot_grid(pE + ggtitle("E") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
								pF + ggtitle("F") + theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")),
								pG + ggtitle("G") + theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")),
								pH + ggtitle("H") + theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")),
								legend.type, # + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
								nrow=1, align="both", rel_widths=c(0.4,0.4,0.35,0.7,0.3))
bottom = plot_grid(row1, row2, nrow=2)

plot_grid(top, bottom, nrow=2)
ggsave("figures/Fig6_herm_spermatogenesis_inherit_ase_pre.pdf", device="pdf", useDingbats=F)

# her-1 WBGene00038586, V.g17207, Cnig_chr_V.g20275
df.inherit2 %>% filter(grepl("Cnig_chr_V.g20275", gene))
df.AS %>% filter(grepl("Cnig_chr_V.g20275", genes))
df.species_sex_class %>% filter(grepl("Cnig_chr_V.g20275", rownames(df.species_sex_class)))

# tra-1 WBGene00033987, III.g8535.t1__Cnig_chr_III.g9780.t2
df.inherit2 %>% filter(grepl("Cnig_chr_III.g9780", gene))
df.AS %>% filter(grepl("Cnig_chr_III.g9780", genes))
df.species_sex_class %>% filter(grepl("Cnig_chr_III.g9780", rownames(df.species_sex_class)))

# tra-2 WBGene00032357, II.g4896.t1__Cnig_chr_II.g5233.t3
df.inherit2 %>% filter(grepl("Cnig_chr_II.g5233", gene))
df.AS %>% filter(grepl("Cnig_chr_II.g5233", genes))
df.species_sex_class %>% filter(grepl("Cnig_chr_II.g5233", rownames(df.species_sex_class)))

# tra-3 WBGene00040309, IV.g10811.t1__Cnig_chr_IV.g12476.t2
df.inherit2 %>% filter(grepl("Cnig_chr_IV.g12476", gene))
df.AS %>% filter(grepl("Cnig_chr_IV.g12476", genes))
df.species_sex_class %>% filter(grepl("Cnig_chr_IV.g12476", rownames(df.species_sex_class)))

# fog-1 WBGene00035332, I.g1335__Cnig_chr_I.g907
# not recovered as a 1-1 ortholog
df.inherit2 %>% filter(grepl("Cnig_chr_I.g907", gene))
df.AS %>% filter(grepl("Cnig_chr_I.g907", genes))
df.species_sex_class %>% filter(grepl("Cnig_chr_I.g907", rownames(df.species_sex_class)))

# fog-2 not found in C. briggase

# fog-3 WBGene00033342, I.g2603.t3__Cnig_chr_I.g2184.t1
df.inherit2 %>% filter(grepl("Cnig_chr_I.g2184", gene))
df.AS %>% filter(grepl("Cnig_chr_I.g2184", genes))
df.species_sex_class %>% filter(grepl("Cnig_chr_I.g2184", rownames(df.species_sex_class)))

# fem-1 WBGene00039054, IV.g11433.t1__Cnig_chr_IV.g13364.t1
df.inherit2 %>% filter(grepl("Cnig_chr_IV.g13364", gene))
df.AS %>% filter(grepl("Cnig_chr_IV.g13364", genes))
df.species_sex_class %>% filter(grepl("Cnig_chr_IV.g13364", rownames(df.species_sex_class)))

# fem-2 WBGene00000328, III.g7416.t1__Cnig_chr_III.g8341.t1
df.inherit2 %>% filter(grepl("Cnig_chr_III.g8341", gene))
df.AS %>% filter(grepl("Cnig_chr_III.g8341", genes))
df.species_sex_class %>% filter(grepl("Cnig_chr_III.g8341", rownames(df.species_sex_class)))

# fem-3 WBGene00000329, IV.g11153.t1__Cnig_chr_IV.g13029.t1
df.inherit2 %>% filter(grepl("Cnig_chr_IV.g13029", gene))
df.AS %>% filter(grepl("Cnig_chr_IV.g13029", genes))
df.species_sex_class %>% filter(grepl("Cnig_chr_IV.g13029", rownames(df.species_sex_class)))

df.sexdet.sex_species = df.species_sex_class %>% filter(grepl("Cnig_chr_V.g20275", rownames(df.species_sex_class))) %>% mutate(gene="her-1")
df.sexdet.sex_species = rbind(df.sexdet.sex_species, df.species_sex_class %>% filter(grepl("Cnig_chr_III.g9780", rownames(df.species_sex_class))) %>% mutate(gene="tra-1"))
df.sexdet.sex_species = rbind(df.sexdet.sex_species, df.species_sex_class %>% filter(grepl("Cnig_chr_II.g5233", rownames(df.species_sex_class))) %>% mutate(gene="tra-2"))
df.sexdet.sex_species = rbind(df.sexdet.sex_species, df.species_sex_class %>% filter(grepl("Cnig_chr_IV.g12476", rownames(df.species_sex_class))) %>% mutate(gene="tra-3"))
df.sexdet.sex_species = rbind(df.sexdet.sex_species, df.species_sex_class %>% filter(grepl("Cnig_chr_I.g907", rownames(df.species_sex_class))) %>% mutate(gene="fog-1"))
df.sexdet.sex_species = rbind(df.sexdet.sex_species, df.species_sex_class %>% filter(grepl("Cnig_chr_I.g2184", rownames(df.species_sex_class))) %>% mutate(gene="fog-3"))
df.sexdet.sex_species = rbind(df.sexdet.sex_species, df.species_sex_class %>% filter(grepl("Cnig_chr_IV.g13364", rownames(df.species_sex_class))) %>% mutate(gene="fem-1"))
df.sexdet.sex_species = rbind(df.sexdet.sex_species, df.species_sex_class %>% filter(grepl("Cnig_chr_III.g8341", rownames(df.species_sex_class))) %>% mutate(gene="fem-2"))
df.sexdet.sex_species = rbind(df.sexdet.sex_species, df.species_sex_class %>% filter(grepl("Cnig_chr_IV.g13029", rownames(df.species_sex_class))) %>% mutate(gene="fem-3"))


rlogdata = read.csv("tables/MDS_plot/rlog-transformed_exprdata.csv", row.names=1)
rlogdata.means = data.frame(CbrF=rowMeans(rlogdata[,1:3]),
                            CbrM=rowMeans(rlogdata[,4:6]),
                            HF1F=rowMeans(rlogdata[,7:9]),
                            HF1M=rowMeans(rlogdata[,10:12]),
                            CniF=rowMeans(rlogdata[,13:15]),
                            CniM=rowMeans(rlogdata[,16:18]))
scaledata = t(scale(t(rlogdata.means)))

tmp = data.frame(expression=scaledata["V.g17207__Cnig_chr_V.g20275",], samples=colnames(scaledata)) %>%
    mutate(sex=c("female","male")[ factor(substr(samples, 4, 4)) ], species=substr(samples, 1, 3), gene="her-1")
tmp = rbind(tmp, data.frame(expression=scaledata["III.g8535.t1__Cnig_chr_III.g9780.t2",], samples=colnames(scaledata)) %>%
    mutate(sex=c("female","male")[ factor(substr(samples, 4, 4)) ], species=substr(samples, 1, 3), gene="tra-1"))
tmp = rbind(tmp, data.frame(expression=scaledata["II.g4896.t1__Cnig_chr_II.g5233.t3",], samples=colnames(scaledata)) %>%
    mutate(sex=c("female","male")[ factor(substr(samples, 4, 4)) ], species=substr(samples, 1, 3), gene="tra-2"))
tmp = rbind(tmp, data.frame(expression=scaledata["IV.g10811.t1__Cnig_chr_IV.g12476.t2",], samples=colnames(scaledata)) %>%
    mutate(sex=c("female","male")[ factor(substr(samples, 4, 4)) ], species=substr(samples, 1, 3), gene="tra-3"))
tmp = rbind(tmp, data.frame(expression=scaledata["I.g1335__Cnig_chr_I.g907",], samples=colnames(scaledata)) %>%
    mutate(sex=c("female","male")[ factor(substr(samples, 4, 4)) ], species=substr(samples, 1, 3), gene="fog-1"))
tmp = rbind(tmp, data.frame(expression=scaledata["I.g2603.t3__Cnig_chr_I.g2184.t1",], samples=colnames(scaledata)) %>%
    mutate(sex=c("female","male")[ factor(substr(samples, 4, 4)) ], species=substr(samples, 1, 3), gene="fog-3"))
tmp = rbind(tmp, data.frame(expression=scaledata["IV.g11433.t1__Cnig_chr_IV.g13364.t1",], samples=colnames(scaledata)) %>%
    mutate(sex=c("female","male")[ factor(substr(samples, 4, 4)) ], species=substr(samples, 1, 3), gene="fem-1"))
tmp = rbind(tmp, data.frame(expression=scaledata["III.g7416.t1__Cnig_chr_III.g8341.t1",], samples=colnames(scaledata)) %>%
    mutate(sex=c("female","male")[ factor(substr(samples, 4, 4)) ], species=substr(samples, 1, 3), gene="fem-2"))
tmp = rbind(tmp, data.frame(expression=scaledata["IV.g11153.t1__Cnig_chr_IV.g13029.t1",], samples=colnames(scaledata)) %>%
    mutate(sex=c("female","male")[ factor(substr(samples, 4, 4)) ], species=substr(samples, 1, 3), gene="fem-3"))
df.sexdet.expr = tmp

tra_fog_fem = ggplot(df.sexdet.expr, aes(x=species, y=expression, group=sex, color=sex)) +
    geom_line(size=2) +
    scale_x_discrete(limits=c("Cbr","HF1","Cni")) +
    facet_rep_wrap(~gene) +
    scale_color_manual(values=cols.sex[c(1,3)], name=NULL) +
    xlab("") +
    background_grid()

tmp = df.AS %>% filter(grepl("Cnig_chr_V.g20275", genes)) %>% mutate(gene="her-1")
tmp = rbind(tmp, df.AS %>% filter(grepl("Cnig_chr_III.g9780", genes)) %>% mutate(gene="tra-1"))
tmp = rbind(tmp, df.AS %>% filter(grepl("Cnig_chr_II.g5233", genes)) %>% mutate(gene="tra-2"))
tmp = rbind(tmp, df.AS %>% filter(grepl("Cnig_chr_IV.g12476", genes)) %>% mutate(gene="tra-3"))
tmp = rbind(tmp, df.AS %>% filter(grepl("Cnig_chr_I.g907", genes)) %>% mutate(gene="fog-1"))
tmp = rbind(tmp, df.AS %>% filter(grepl("Cnig_chr_I.g2184", genes)) %>% mutate(gene="fog-3"))
tmp = rbind(tmp, df.AS %>% filter(grepl("Cnig_chr_IV.g13364", genes)) %>% mutate(gene="fem-1"))
tmp = rbind(tmp, df.AS %>% filter(grepl("Cnig_chr_III.g8341", genes)) %>% mutate(gene="fem-2"))
tmp = rbind(tmp, df.AS %>% filter(grepl("Cnig_chr_IV.g13029", genes)) %>% mutate(gene="fem-3"))
tmp = tmp %>% mutate(sex=sub("\\(M\\) hermaphrodite","male",sex))
df.sexdet.cis_trans = tmp

tra_fog_fem_cis_trans = ggplot(df.sexdet.cis_trans, aes(logFC.sp, logFC.ase, color=sex, label=gene)) +
    geom_point(data=df.AS, aes(logFC.sp, logFC.ase), color="grey80", alpha=0.1, inherit.aes=F) +
	geom_hline(yintercept=0, linetype=2) +
	geom_vline(xintercept=0, linetype=2) +
    scale_x_continuous(expand=c(0,0), limits=c(-4,4)) +
    scale_y_continuous(expand=c(0,0), limits=c(-4,2)) +
    geom_abline(slope=1, linetype=2) +
    geom_point(show.legend=F, size=2, color="black") +
	geom_label_repel(show.legend=T, size=4, min.segment.length=0, nudge_x=.2, key_glyph = "point") +
    #scale_shape_manual(values=c(10,11,12,13,15,17,18)) +
	scale_color_manual(values=cols.sex, name="") +
    annotate(geom="text", x=-3, y=.2, label="trans-only") +
    annotate(geom="text", x=1.6, y=-3.8, label="cis-trans (compensatory)") +
    annotate(geom="text", x=-2.5, y=-1.8, label="cis-only") +
	ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cnig"),"/",italic("Cbri"),")]", sep=''))) +
	xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cnig"),"/",italic("Cbri"),")]", sep=''))) +
	background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
	#facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	theme(strip.background = element_rect(fill="grey90", color=NA),
        legend.position=c(.7,.25)) #+
	#theme(strip.background = element_blank(), strip.text=element_blank(),
	#			axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

plot_grid(tra_fog_fem, tra_fog_fem_cis_trans, ncol=2, rel_widths=c(.5,.4))
ggsave("figures/Suppl_Fig_sexdet.png", device="png")
ggsave("figures/Suppl_Fig_sexdet.pdf", device="pdf", useDingbats=F)
#
#
#
# sp_mod_enrich = ggplot(sperm_models_enrich, aes(x=models, y=log2(enrichment), fill=sex)) +
# 	geom_bar(stat="identity", show.legend=F) +
# 	geom_text(aes(y=n_nudge_y, label=n)) +
# 	geom_text(aes(y=log2(enrichment)/2, label=sig), size=5, show.legend=F, color="white") +
# 	scale_fill_manual(values=cols.sex) +
# 	scale_x_discrete(labels=as.character(sperm_models_enrich$models2)) +
# 	scale_y_continuous(breaks=seq(-2, 4, 1)) +
# 	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
# 	xlab("") +
# 	ylab(expression(paste("gene enrichment "," [log"[2]," OR]",sep=""))) +
# 	theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
#
# bp_counts_chr = ggplot(sperm_class_male_chr_count, aes(x=chromosome, y=count, fill=class)) +
# 	geom_bar(stat="identity", show.legend=F, alpha=0.7) +
# 	scale_fill_manual(values=cols.class) +
# 	ylab("number of genes") +
# 	background_grid(major="y", size.major = 0.2, colour.major = "grey75")
#
# rects1 = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-1.7,1.5,-1.7), y2=c(1.5,-1.7,1.5))
#
# sp_enrich_class = ggplot(sperm_class_male_enrich_chr, aes(x=chr.pseudo, y=log2(enrichment), fill=class)) +
#  	geom_rect(data=rects1, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
# 	geom_bar(stat="identity", position=position_dodge(), alpha=0.7, show.legend=F) +
# 	ylim(-1.7, 1.7) +
# 	scale_x_continuous(breaks=1:6, labels=sperm_class_male_enrich_chr$chromosome[1:6]) +
# 	scale_fill_manual(values=cols.class) +
# 	background_grid(major="y", minor="y", size.major = 0.2, colour.major = "grey75", size.minor = 0.2, colour.minor = "grey75") +
# 	ylab(expression(paste("",sep=""))) +
# 	xlab("") +
# 	geom_text(aes(y=log2(enrichment)/2, label=sig), size=5, position=position_dodge(width=0.9), show.legend=F, color="white") +
# 	theme(axis.text.x=element_blank())
#
# rects2 = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-2,1.5,-2), y2=c(1.5,-2,1.5))
#
# sp_enrich_type = ggplot(sperm_type_male_enrich_chr, aes(x=chr.pseudo, y=log2(enrichment), fill=type)) +
#  	geom_rect(data=rects2, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
# 	geom_bar(stat="identity", position=position_dodge(), alpha=0.7, show.legend=F) +
# 	ylim(-2, 2) +
# 	scale_x_continuous(breaks=1:6, labels=sperm_type_male_enrich_chr$chromosome[1:6]) +
# 	scale_fill_manual(values=cols.type) +
# 	background_grid(major="y", minor="y", size.major = 0.2, colour.major = "grey75", size.minor = 0.2, colour.minor = "grey75") +
# 	ylab(expression(paste("gene enrichment "," [log"[2]," OR]",sep=""))) +
# 	xlab("chromosome") +
# 	geom_text(aes(y=log2(enrichment)/2, label=sig), size=5, position=position_dodge(width=0.9), show.legend=F, color="white")
#
# row3 = plot_grid(NULL,
#                     sp_mod_enrich, #+ theme(plot.margin = unit(c(t=0, r=0, b=0, l=0.5), "cm")),
# 					bp_counts_chr,
# 					plot_grid(sp_enrich_class + theme(plot.margin = unit(c(t=0, r=0, b=-0.3, l=0), "cm")),
# 										sp_enrich_type + theme(plot.margin = unit(c(t=-0.3, r=0, b=0, l=0), "cm")),
# 										ncol=1, align = "v"),
#                     NULL,
# 					nrow=1, rel_widths=c(0.2,0.4,0.3,0.4,0.2))
#
# plot_grid(row1, row2, ncol=1, rel_heights=c(0.55,0.45))
#
# ggsave("Fig4_species_sex_and_spermatogenic.pdf", device="pdf")

# Xsperm_molevol = read.csv("Xchr_sperm_cis_trans_molevol.csv", row.name=1)
# Xsperm_molevol$type_male = factor(Xsperm_molevol$type_male, type.levels)
#
# sp_dnds_pos = ggplot(Xsperm_molevol, aes(x=position/1e6, y=Ka/Ks_ENCc, color=type_male, group=type_male)) +
#     geom_rect(aes(xmin=8, xmax=10, ymin=min(Ka/Ks_ENCc), ymax=max(Ka/Ks_ENCc)), fill="grey90", color=NA) +
# 	geom_point(show.legend=F, alpha=0.7, shape=16) +
# 	geom_line(show.legend=F, alpha=0.7) +
# 	scale_x_continuous(breaks=seq(0,20,1), label=c(0,rep("",4),5,rep("",4),10,rep("",4),15,rep("",4),20)) +
# 	xlab("") +
#     ylab(expression(paste(italic("K")[a],"/",italic("K")[s],"\'"))) +
# 	scale_color_manual(values=cols.type, name="") +
# 	background_grid(major="y", size.major = 0.2, colour.major = "grey75")
#
# sp_ka_pos = ggplot(Xsperm_molevol, aes(x=position/1e6, y=Ka, color=type_male, group=type_male)) +
#     geom_rect(aes(xmin=8, xmax=10, ymin=min(Ka), ymax=max(Ka)), fill="grey90", color=NA) +
# 	geom_point(show.legend=F, alpha=0.7, shape=16) +
# 	geom_line(show.legend=F, alpha=0.7) +
# 	scale_x_continuous(breaks=seq(0,20,1), label=c(0,rep("",4),5,rep("",4),10,rep("",4),15,rep("",4),20)) +
# 	xlab("") +
#     ylab(expression(paste(italic("K")[a]))) +
# 	scale_color_manual(values=cols.type, name="") +
# 	background_grid(major="y", size.major = 0.2, colour.major = "grey75")
#
# sp_prop_pos = ggplot(Xsperm_molevol, aes(x=(position-500)/1e6, y=1-prop_cons_500bp_upstream_5bp_win, color=type_male, group=type_male)) +
#     geom_rect(aes(xmin=8, xmax=10, ymin=min(1-prop_cons_500bp_upstream_5bp_win, na.rm=T), ymax=max(1-prop_cons_500bp_upstream_5bp_win, na.rm=T)), fill="grey90", color=NA) +
#     geom_point(show.legend=F, alpha=0.7, shape=16) +
# 	geom_line(show.legend=F, alpha=0.7) +
# 	scale_x_continuous(breaks=seq(0,20,1), label=c(0,rep("",4),5,rep("",4),10,rep("",4),15,rep("",4),20)) +
# 	xlab("") +
# 	ylab(expression(paste("1-",italic("p")[cons]))) +
# 	scale_color_manual(values=cols.type, name="") +
# 	background_grid(major="y", size.major = 0.2, colour.major = "grey75")
#
# x_dens_pos = ggplot(filter(Xsperm_molevol, position < 20*1e6), aes(x=position/1e6, color=type_male)) +
#     geom_rect(aes(xmin=8, xmax=10, ymin=0, ymax=0.1), fill="grey90", color=NA) +
# 	geom_line(stat="density", show.legend=F, alpha=0.7) +
# 	scale_x_continuous(limits=c(0,20), breaks=seq(0,20,1), label=c(0,rep("",4),5,rep("",4),10,rep("",4),15,rep("",4),20)) +
# 	xlab("genome postion (Mbp)") +
# 	ylab("density") +
# 	scale_color_manual(values=cols.type, name="") +
# 	background_grid(major="y", size.major = 0.2, colour.major = "grey75")
#
# plot_grid(sp_dnds_pos + theme(plot.margin = unit(c(t=0, r=0, b=-0.3, l=0), "cm")),
#     sp_ka_pos + theme(plot.margin = unit(c(t=-0.15, r=0, b=-0.15, l=0), "cm")),
#     sp_prop_pos + theme(plot.margin = unit(c(t=-0.15, r=0, b=-0.15, l=0), "cm")),
#     x_dens_pos + theme(plot.margin = unit(c(t=-0.3, r=0, b=0, l=0), "cm")),
#     align="v", ncol=1, labels=c("A","B","C","D"), label_y=c(1,1.1,1.1,1.1))
#
# ggsave("sperm_cis_trans_molevol_Xchr.png")
# ggsave("sperm_cis_trans_molevol_Xchr.pdf", device="pdf")

# supplementary

cis_trans_count_sexbias = read.csv("cis_trans_count_sexbias.csv")
ggplot(cis_trans_count_sexbias, aes(x=effect, y=count, fill=effect)) +
    geom_bar(stat="identity", show.legend=F) +
    xlab("") +
    ylab("number of genes") +
    facet_wrap(~sex) +
    scale_fill_manual(values=cols.type[4:3])

ggsave("suppl_cis_trans_count_sexbias.png")
ggsave("suppl_cis_trans_count_sexbias.pdf", device="pdf")

cis_trans.sexbias = read.csv("cis_trans.sexbias.sex.csv")

cis_trans.sexbias = df.aes.inherit.merged %>% dplyr::select(logFC.ase, logFC.sp, sex, gene) %>% mutate(logFC.trans = logFC.sp-logFC.ase, effect = c("trans","cis")[ factor(abs(logFC.ase) > abs(logFC.trans)) ])
tmp = pivot_wider(cis_trans.sexbias, id_cols=gene, values_from=c(logFC.ase, logFC.sp, logFC.trans, effect), names_from=sex)

corr1 = ggplot(tmp, aes(logFC.ase_female, logFC.ase_male)) +
    geom_hline(yintercept=0, linetype=2) +
    geom_vline(xintercept=0, linetype=2) +
    geom_point(alpha=0.4) +
    geom_abline(slope=1, linetype=2) +
    geom_smooth(method="lm") +
    xlim(-15,15) +
    ylim(-15,15) +
    labs(x="logFC cis (female)", y="logFC cis (male)")
corr2 = ggplot(tmp, aes(logFC.trans_female, logFC.trans_male)) +
    geom_hline(yintercept=0, linetype=2) +
    geom_vline(xintercept=0, linetype=2) +
    geom_point(alpha=0.4) +
    geom_abline(slope=1, linetype=2) +
    xlim(-15,15) +
    ylim(-15,15) +
    geom_smooth(method="lm") +
    labs(x="logFC trans (female)", y="logFC trans (male)")
plot_grid(corr1, corr2)
ggsave("suppl_corr_sex_cis_trans_magnitude.png")


[sample(1:nrow(cis_trans.sexbias), replace=T),]

boot.sexbias = data.frame()
for (i in 1:500)
    boot.sexbias = rbind(boot.sexbias, (cis_trans.sexbias %>% filter((sex == "male" | sex == "female") & !is.na(effect)) %>% slice_sample(n = 22230, replace=T) %>% group_by(sex, effect) %>% count() )$n)
cis_trans.sexbias.count = as.data.frame(cis_trans.sexbias %>% filter((sex == "male" | sex == "female") & !is.na(effect)) %>% group_by(sex, effect) %>% count())
cis_trans.sexbias.count = cbind(cis_trans.sexbias.count, t(apply(boot.sexbias, 2, quantile, c(0.05,0.95))))

p1 = ggplot(cis_trans.sexbias %>% filter((sex == "male" | sex == "female") & !is.na(effect)), aes(x=logFC.ase, y=logFC.trans, color=effect)) +
    geom_point(alpha=0.3, show.legend=F) +
    facet_rep_wrap(~sex, ncol=1) +
    scale_color_manual(values=cols.type[4:3], name="") +
    geom_vline(xintercept=0, linetype=2) +
    geom_hline(yintercept=0, linetype=2) +
    xlim(-15,15) +
    ylim(-15,15) +
    xlab(expression(paste(italic(cis)," (log"[2]," fold-change)"))) +
    ylab(expression(paste(italic(trans)," (log"[2]," fold-change)"))) +
    background_grid() +
    theme(strip.background=element_blank(), strip.text=element_blank(),
            legend.position=c(0.1,0.1))
p2 = ggplot(cis_trans.sexbias.count, aes(x=effect, y=n, fill=effect)) +
    geom_bar(stat="identity", width=.75, show.legend=F) +
    geom_linerange(aes(ymin=`5%`, ymax=`95%`)) +
    scale_fill_manual(values=cols.type[4:3], name="") +
    facet_rep_wrap(~sex, ncol=1, strip.position="right") +
    xlab("") +
    ylab("number of genes") +
    background_grid(major="y", color.major="grey90") +
    theme(strip.background=element_rect(fill="grey75"))

plot_grid(p1, p2, ncol=2, rel_widths=c(1,0.7))
ggsave("figures/suppl_cis_trans.sexbias.png")

# supplementary

Cbr_allele = df.aes.inherit.merged %>%
    filter(sex == "male" & (class == "C. briggsae dominant" |
                    class == "C. nigoni dominant" | class == "overdominant" |
                    class == "underdominant" | class == "additive") & chromosome != "X") %>%
    mutate(Cbr = logFC.ase < 0) %>%
    group_by(class, type, Cbr) %>%
    count() %>%
    filter(type == "cis-only" | type == "cis + trans (enhancing)" | type == "cis x trans (compensatory)" | type == "cis-trans (compensatory)")

ggplot(Cbr_allele, aes(x=type, y=n, fill=Cbr)) +
    geom_bar(stat="identity", position="fill", width=0.7) +
    scale_fill_manual(values=cols.class[5:4], name="allele", labels=c("C. nigoni","C. briggsae")) +
    scale_x_discrete(labels=c(
        expression(paste(italic("cis"), "-only")),
        expression(paste(italic("cis + trans")," (enhancing)")),
        expression(paste(italic("cis x trans")," (compensatory)")),
        expression(paste(italic("cis-trans")," (compensatory)")))) +
    scale_y_continuous(expand=c(0,0,0,0), labels=c("0","","0.5","","1")) +
    xlab("") +
    ylab("proportion of genes") +
    background_grid(major="x", color.major="grey75") +
    theme(legend.text=element_text(face="italic")) +
    theme(legend.position=c(.8,.3)) +
    coord_flip() +
    facet_wrap(~class, scales="free")

ggsave("figures/suppl_dom_transgr_ASE.png")
ggsave("figures/suppl_dom_transgr_ASE.pdf", device="pdf")
