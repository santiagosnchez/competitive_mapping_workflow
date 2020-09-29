# libraries
library(ggplot2)
library(cowplot)
library(lemon)
library(ggridges)
library(ggstance)
library(ggsci)
library(tidyr)
library(dplyr)
library(gtable)
library(ggplotify)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(viridis)
theme_set(theme_cowplot())

# colors

brewer.pal(10, "Paired") -> cols
cols = cols[-c(7,8)]
cols.spp = cols[c(6,2,8)]
cols.sex = c("#80CDC1","#018571","#BEBEBE")
cols.heatmap = c(cols.spp[1],"grey",cols.sex[1],"black",cols.sex[2],cols.spp[2],"grey",cols.sex[3])
cols.class = c("black","grey",cols[ c(4,6,2,8,7) ])
cols.type = c("black","grey",pal_jama()(5))
cols.type[3] = "#09BC8A"
cols.type[7]= "#C49991"
cols.type = cols.type[c(1,2,4,3,5,6,7)]
cols.type[4] = "#476A6F"
cols.heatmap1 = cols.heatmap[c(1,2,2)]

# reaction norms

df.models.cent = read.csv("tables/df.centroides_and_CI.species_sex_reaction_norms.csv")

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
    scale_color_manual(values=cols.sex[1:2])

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
    scale_color_manual(values=cols.sex[1:2])

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
    scale_color_manual(values=cols.sex[1:2])

row1 = plot_grid(cons + theme(plot.margin = unit(c(t = 0, r = 9, b = 0, l = 0), "cm")),
        male_biased + theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "cm")),
        female_biased + theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "cm")),
        nrow=3)

row1
ggsave("figures/Fig5a_species_sex_reaction_norms.pdf", device="pdf", useDingbats=F)
ggsave("figures/Fig5a_species_sex_reaction_norms.png")

# read class count and expression data

df.species_sex_class = read.csv("tables/df.species_sex.csv", row.names=1)
df.species_sex_class.expr_type.counts = read.csv("tables/df.species_sex_class.expr_type.counts.csv")
df.species_sex_class.expr_class.counts = read.csv("tables/df.species_sex_class.expr_class.counts.csv")
df.models2 = read.csv("tables/species_sex_models_heatmap.csv")

class.levels = c("no change","ambiguous","additive","C. briggsae dominant", "C. nigoni dominant","overdominant","underdominant")
type.levels = c("conserved","ambiguous","trans-only","cis-only","cis + trans (enhancing)","cis x trans (compensatory)","cis-trans (compensatory)")
models.levels = sort(unique(as.character(df.species_sex_class$models)))[c(8,5,13,7,4,12,3,11,6,2,10,1,9)]
models.levels2 = sort(unique(as.character(df.species_sex_class$models2)))[c(8,5,13,7,4,12,3,11,6,2,10,1,9)]
models.counts = sort(unique(as.character(df.species_sex_class$counts)))[c(8,5,13,7,4,12,3,11,6,2,10,1,9)]

df.species_sex_class.expr_class.counts$class = factor(df.species_sex_class.expr_class.counts$class, class.levels)
df.species_sex_class.expr_type.counts$type = factor(df.species_sex_class.expr_type.counts$type, type.levels)
df.species_sex_class.expr_class.counts$models = factor(df.species_sex_class.expr_class.counts$models, models.levels)
df.species_sex_class.expr_type.counts$models = factor(df.species_sex_class.expr_type.counts$models, models.levels)
df.species_sex_class.expr_class.counts$models2 = factor(df.species_sex_class.expr_class.counts$models2, models.levels2)
df.species_sex_class.expr_type.counts$models2 = factor(df.species_sex_class.expr_type.counts$models2, models.levels2)
df.species_sex_class.expr_class.counts$sex2 = sapply(strsplit(as.character(df.species_sex_class.expr_class.counts$models), " "), function(x) x[2] )
df.species_sex_class.expr_type.counts$sex2 = sapply(strsplit(as.character(df.species_sex_class.expr_type.counts$models), " "), function(x) x[2] )
df.species_sex_class.expr_class.counts$sex2 = factor(df.species_sex_class.expr_class.counts$sex2, c("sex-neutral","male-biased","female-biased"))
df.species_sex_class.expr_type.counts$sex2 = factor(df.species_sex_class.expr_type.counts$sex2, c("sex-neutral","male-biased","female-biased"))
df.species_sex_class$models = factor(df.species_sex_class$models, models.levels)
df.species_sex_class$models2 = factor(df.species_sex_class$models2, models.levels2)
df.species_sex_class$sex = factor(df.species_sex_class$sex, c("sex-neutral","male-biased","female-biased"))
df.models2$models = factor(df.models2$models, models.levels)
df.models2$sex = factor(df.models2$sex, c("sex-neutral","male-biased","female-biased"))
models.counts = unique(as.character(df.models2$counts))[ order(unique(df.models2$models)) ]
df.models2$counts = factor(as.character(df.models2$counts), models.counts)

# plots

models_grid = ggplot(df.models2, aes(x=key, y=counts, fill=value)) +
	geom_tile(color="black", show.legend=F) +
	scale_fill_manual(values=cols.heatmap) +
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
	geom_density_ridges(aes(fill=sex), scale=2, panel_scaling=F, color=NA, alpha=0.5, show.legend=F) +
	geom_boxploth(alpha=0.5, outlier.shape=NA, width=0.5) +
	scale_fill_manual(values=rev(cols.sex)) +
	scale_x_continuous(position="top", limits=c(0,10), labels=c("0","","5","","10")) +
	#scale_y_discrete(labels=models.levels2) +
	xlab(expression(paste("log"[2]," expression divergence"))) +
	ylab("") +
    facet_rep_wrap(~sex, strip.position="right", drop=T, scales="free_y", ncol=1) +
	background_grid(major="x", size.major = 0.2, colour.major = "grey75", minor="x", size.minor = 0.2, colour.minor = "grey75") +
	theme(strip.background=element_blank(), strip.text=element_blank())

aa_div = ggplot(df.species_sex_class, aes(y=models2, x=Ka/Ks_ENCc)) +
	geom_density_ridges(aes(fill=sex), scale=1, panel_scaling=F, color=NA, alpha=0.5, show.legend=F) +
	geom_boxploth(alpha=0.5, outlier.shape=NA, width=0.5) +
	scale_fill_manual(values=rev(cols.sex)) +
	scale_x_continuous(position="top", limits=c(0,0.75), breaks=c(0,0.2,0.4,0.6), labels=c("0","0.2","0.4","0.6")) +
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
    mutate(chromosome2=c("auto","X")[ factor(chromosome == "X")]) %>%
    group_by(sex, models2, chromosome2) %>%
    summarise(n2=sum(n)))
x = enrichment(t(matrix(tmp1$n2, ncol=2, nrow=13, byrow=T)), odds.ratio=T)
y = enrichment(t(matrix(tmp1$n2, ncol=2, nrow=13, byrow=T)), odds.ratio=F)
tmp1 = as.data.frame(tmp1 %>%
    add_column(enrichment=as.vector(x)) %>%
    add_column(p.value=as.vector(y)) %>%
    group_by(sex, models2) %>%
    summarize(chromosome2, n2, enrichment, p.value, rel_enrich=enrichment/sum(enrichment)))
tmp2 =  as.data.frame(df.species_sex_class %>%
    filter(!(domain == "tipA" | domain == "tipB")) %>%
    mutate(domain2 = gsub("A$|B$","", as.character(domain))) %>%
    group_by(sex, models2, domain2) %>%
    summarize(n=n()))
x = enrichment(t(matrix(tmp2$n, ncol=2, nrow=13, byrow=T)), odds.ratio=T)
y = enrichment(t(matrix(tmp2$n, ncol=2, nrow=13, byrow=T)), odds.ratio=F)
tmp2 = as.data.frame(tmp2 %>%
    add_column(enrichment=as.vector(x)) %>%
    add_column(p.value=as.vector(y)) %>%
    group_by(sex, models2) %>%
    summarize(domain2, n, enrichment, p.value, rel_enrich=enrichment/sum(enrichment)))

df.enrich_arms_X = tmp2 %>%
    filter(domain2 == "arm") %>%
    rename(n_arm = n, rel_enrich.arm=rel_enrich, p.value_arm=p.value)
df.enrich_arms_X = df.enrich_arms_X %>%
    add_column(tmp1 %>%
        filter(chromosome2 == "X") %>%
        rename(n_X = n2, rel_enrich.X=rel_enrich, p.value_X=p.value) %>%
        dplyr::select(chromosome2,n_X, rel_enrich.X, p.value_X))
df.enrich_arms_X2 = df.enrich_arms_X %>%
    dplyr::select(rel_enrich.arm, rel_enrich.X, sex, models2) %>%
    pivot_longer(c(rel_enrich.arm, rel_enrich.X, -sex, -models2), names_to="arm_X", values_to="relative_enrichment")
df.enrich_arms_X2$arm_X = gsub("rel_enrich.arm","chromosome arms", df.enrich_arms_X2$arm_X)
df.enrich_arms_X2$arm_X = gsub("rel_enrich.X","X", df.enrich_arms_X2$arm_X)
df.enrich_arms_X2$models2 = factor(df.enrich_arms_X2$models2, models.levels2)
df.enrich_arms_X2$sex = factor(df.enrich_arms_X2$sex, c("sex-neutral","male-biased","female-biased"))

arms_X = ggplot(df.enrich_arms_X2, aes(x=relative_enrichment, y=models2, fill=sex)) +
    geom_bar(stat="identity", show.legend=F, alpha=0.7) +
    scale_x_continuous(limits=c(0,1), expand=c(0,0), position="top", labels=c("0","","0.5","","1")) +
    facet_rep_grid(sex~arm_X, scales="free_y") +
    scale_fill_manual(values=rev(cols.sex)) +
    labs(x="relative enrichment", y="") +
    theme(axis.text.y=element_blank()) +
    background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
    theme(strip.background = element_rect(fill="grey90", color=NA))

row3 = plot_grid(expr_div + theme(plot.margin = unit(c(t=0.2, r=0, b=0, l=-0.3), "cm")),
					aa_div + theme(plot.margin = unit(c(t=0.15, r=0, b=0, l=-0.1), "cm")),
                    arms_X + theme(plot.margin = unit(c(t=0.3, r=0, b=0, l=-0.2), "cm")),
					rel_widths=c(0.25,0.25,0.5), nrow=1)
row3

plot_grid(row2, row3, ncol=1)
ggsave("figures/Fig5b2_species_sex_expr_div_dnds_counts.pdf", device="pdf", useDingbats=F)
ggsave("figures/Fig5b2_species_sex_expr_div_dnds_counts.png")

# spermatogenic genes

sperm_models_enrich = read.csv("sperm_species_sex_models_enrich.csv")
sperm_class_male_chr_count = read.csv("sperm_class_male_chr_count.csv")
sperm_type_male_enrich_chr = read.csv("sperm_type_male_enrich_chr.csv")
sperm_class_male_enrich_chr = read.csv("sperm_class_male_enrich_chr.csv")

sperm_models_enrich$models = factor(sperm_models_enrich$models, models.levels)
sperm_models_enrich$n_nudge_y = log2(sperm_models_enrich$enrichment) + (log2(sperm_models_enrich$enrichment)/abs(log2(sperm_models_enrich$enrichment)) * 0.3)
sperm_class_male_chr_count$class = factor(sperm_class_male_chr_count$class, class.levels)
sperm_type_male_enrich_chr$type = factor(sperm_type_male_enrich_chr$type, type.levels)
sperm_type_male_enrich_chr$chr.pseudo = (1:6)[ factor(sperm_type_male_enrich_chr$chromosome) ]
sperm_class_male_enrich_chr$class = factor(sperm_class_male_enrich_chr$class, class.levels)
sperm_class_male_enrich_chr$chr.pseudo = (1:6)[ factor(sperm_class_male_enrich_chr$chromosome) ]

sperm_class_male_enrich_chr$sig = ""
sperm_type_male_enrich_chr$sig = ""
sperm_models_enrich$sig = ""
sperm_class_male_enrich_chr$sig[ sperm_class_male_enrich_chr$p.value < 0.05 & log2(sperm_class_male_enrich_chr$enrichment) > 0.5 ] = "*"
sperm_type_male_enrich_chr$sig[ sperm_type_male_enrich_chr$p.value < 0.05 & log2(sperm_type_male_enrich_chr$enrichment) > 0.5 ] = "*"
sperm_models_enrich$sig[ sperm_models_enrich$p.value < 0.05 & log2(sperm_models_enrich$enrichment) > 0.5 ] = "*"

brewer.pal(10, "Paired") -> cols
cols = cols[-c(7,8)]
cols.spp = cols[c(6,2,8)]
cols.sex = c("#80CDC1","#018571","#BEBEBE")
cols.heatmap = c(cols.spp[1],"grey",cols.sex[1],cols.spp[3],cols.sex[2],cols.spp[2],"grey",cols.sex[3])
cols.class = c("black","grey",cols[ c(4,6,2,8,7) ])
cols.type = c("black","grey",pal_jama()(4))
cols.type[3] = "#09BC8A"
cols.heatmap1 = cols.heatmap[c(1,2,2)]

## plots



# spermatogenic genes

sp_mod_enrich = ggplot(sperm_models_enrich, aes(x=models, y=log2(enrichment), fill=sex)) +
	geom_bar(stat="identity", show.legend=F) +
	geom_text(aes(y=n_nudge_y, label=n)) +
	geom_text(aes(y=log2(enrichment)/2, label=sig), size=5, show.legend=F, color="white") +
	scale_fill_manual(values=cols.sex) +
	scale_x_discrete(labels=as.character(sperm_models_enrich$models2)) +
	scale_y_continuous(breaks=seq(-2, 4, 1)) +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	xlab("") +
	ylab(expression(paste("gene enrichment "," [log"[2]," OR]",sep=""))) +
	theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))

bp_counts_chr = ggplot(sperm_class_male_chr_count, aes(x=chromosome, y=count, fill=class)) +
	geom_bar(stat="identity", show.legend=F, alpha=0.7) +
	scale_fill_manual(values=cols.class) +
	ylab("number of genes") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75")

rects1 = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-1.7,1.5,-1.7), y2=c(1.5,-1.7,1.5))

sp_enrich_class = ggplot(sperm_class_male_enrich_chr, aes(x=chr.pseudo, y=log2(enrichment), fill=class)) +
 	geom_rect(data=rects1, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
	geom_bar(stat="identity", position=position_dodge(), alpha=0.7, show.legend=F) +
	ylim(-1.7, 1.7) +
	scale_x_continuous(breaks=1:6, labels=sperm_class_male_enrich_chr$chromosome[1:6]) +
	scale_fill_manual(values=cols.class) +
	background_grid(major="y", minor="y", size.major = 0.2, colour.major = "grey75", size.minor = 0.2, colour.minor = "grey75") +
	ylab(expression(paste("",sep=""))) +
	xlab("") +
	geom_text(aes(y=log2(enrichment)/2, label=sig), size=5, position=position_dodge(width=0.9), show.legend=F, color="white") +
	theme(axis.text.x=element_blank())

rects2 = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-2,1.5,-2), y2=c(1.5,-2,1.5))

sp_enrich_type = ggplot(sperm_type_male_enrich_chr, aes(x=chr.pseudo, y=log2(enrichment), fill=type)) +
 	geom_rect(data=rects2, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
	geom_bar(stat="identity", position=position_dodge(), alpha=0.7, show.legend=F) +
	ylim(-2, 2) +
	scale_x_continuous(breaks=1:6, labels=sperm_type_male_enrich_chr$chromosome[1:6]) +
	scale_fill_manual(values=cols.type) +
	background_grid(major="y", minor="y", size.major = 0.2, colour.major = "grey75", size.minor = 0.2, colour.minor = "grey75") +
	ylab(expression(paste("gene enrichment "," [log"[2]," OR]",sep=""))) +
	xlab("chromosome") +
	geom_text(aes(y=log2(enrichment)/2, label=sig), size=5, position=position_dodge(width=0.9), show.legend=F, color="white")

row3 = plot_grid(NULL,
                    sp_mod_enrich, #+ theme(plot.margin = unit(c(t=0, r=0, b=0, l=0.5), "cm")),
					bp_counts_chr,
					plot_grid(sp_enrich_class + theme(plot.margin = unit(c(t=0, r=0, b=-0.3, l=0), "cm")),
										sp_enrich_type + theme(plot.margin = unit(c(t=-0.3, r=0, b=0, l=0), "cm")),
										ncol=1, align = "v"),
                    NULL,
					nrow=1, rel_widths=c(0.2,0.4,0.3,0.4,0.2))

plot_grid(row1, row2, ncol=1, rel_heights=c(0.55,0.45))

ggsave("Fig4_species_sex_and_spermatogenic.pdf", device="pdf")

Xsperm_molevol = read.csv("Xchr_sperm_cis_trans_molevol.csv", row.name=1)
Xsperm_molevol$type_male = factor(Xsperm_molevol$type_male, type.levels)

sp_dnds_pos = ggplot(Xsperm_molevol, aes(x=position/1e6, y=Ka/Ks_ENCc, color=type_male, group=type_male)) +
    geom_rect(aes(xmin=8, xmax=10, ymin=min(Ka/Ks_ENCc), ymax=max(Ka/Ks_ENCc)), fill="grey90", color=NA) +
	geom_point(show.legend=F, alpha=0.7, shape=16) +
	geom_line(show.legend=F, alpha=0.7) +
	scale_x_continuous(breaks=seq(0,20,1), label=c(0,rep("",4),5,rep("",4),10,rep("",4),15,rep("",4),20)) +
	xlab("") +
    ylab(expression(paste(italic("K")[a],"/",italic("K")[s],"\'"))) +
	scale_color_manual(values=cols.type, name="") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75")

sp_ka_pos = ggplot(Xsperm_molevol, aes(x=position/1e6, y=Ka, color=type_male, group=type_male)) +
    geom_rect(aes(xmin=8, xmax=10, ymin=min(Ka), ymax=max(Ka)), fill="grey90", color=NA) +
	geom_point(show.legend=F, alpha=0.7, shape=16) +
	geom_line(show.legend=F, alpha=0.7) +
	scale_x_continuous(breaks=seq(0,20,1), label=c(0,rep("",4),5,rep("",4),10,rep("",4),15,rep("",4),20)) +
	xlab("") +
    ylab(expression(paste(italic("K")[a]))) +
	scale_color_manual(values=cols.type, name="") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75")

sp_prop_pos = ggplot(Xsperm_molevol, aes(x=(position-500)/1e6, y=1-prop_cons_500bp_upstream_5bp_win, color=type_male, group=type_male)) +
    geom_rect(aes(xmin=8, xmax=10, ymin=min(1-prop_cons_500bp_upstream_5bp_win, na.rm=T), ymax=max(1-prop_cons_500bp_upstream_5bp_win, na.rm=T)), fill="grey90", color=NA) +
    geom_point(show.legend=F, alpha=0.7, shape=16) +
	geom_line(show.legend=F, alpha=0.7) +
	scale_x_continuous(breaks=seq(0,20,1), label=c(0,rep("",4),5,rep("",4),10,rep("",4),15,rep("",4),20)) +
	xlab("") +
	ylab(expression(paste("1-",italic("p")[cons]))) +
	scale_color_manual(values=cols.type, name="") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75")

x_dens_pos = ggplot(filter(Xsperm_molevol, position < 20*1e6), aes(x=position/1e6, color=type_male)) +
    geom_rect(aes(xmin=8, xmax=10, ymin=0, ymax=0.1), fill="grey90", color=NA) +
	geom_line(stat="density", show.legend=F, alpha=0.7) +
	scale_x_continuous(limits=c(0,20), breaks=seq(0,20,1), label=c(0,rep("",4),5,rep("",4),10,rep("",4),15,rep("",4),20)) +
	xlab("genome postion (Mbp)") +
	ylab("density") +
	scale_color_manual(values=cols.type, name="") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75")

plot_grid(sp_dnds_pos + theme(plot.margin = unit(c(t=0, r=0, b=-0.3, l=0), "cm")),
    sp_ka_pos + theme(plot.margin = unit(c(t=-0.15, r=0, b=-0.15, l=0), "cm")),
    sp_prop_pos + theme(plot.margin = unit(c(t=-0.15, r=0, b=-0.15, l=0), "cm")),
    x_dens_pos + theme(plot.margin = unit(c(t=-0.3, r=0, b=0, l=0), "cm")),
    align="v", ncol=1, labels=c("A","B","C","D"), label_y=c(1,1.1,1.1,1.1))

ggsave("sperm_cis_trans_molevol_Xchr.png")
ggsave("sperm_cis_trans_molevol_Xchr.pdf", device="pdf")

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

boot.sexbias = data.frame()
for (i in 1:1000)
    boot.sexbias = rbind(boot.sexbias, (cis_trans.sexbias[sample(1:nrow(cis_trans.sexbias), replace=T),] %>% group_by(bias, effect2) %>% count() )$n)
cis_trans.sexbias.count = as.data.frame(cis_trans.sexbias %>% group_by(bias, effect2) %>% count())
cis_trans.sexbias.count = cbind(cis_trans.sexbias.count, t(apply(boot.sexbias, 2, quantile, c(0.05,0.95))))

p1 = ggplot(cis_trans.sexbias, aes(x=logFC_F1_cis, y=logFC_trans, color=effect2)) +
    geom_point(show.legend=T) +
    facet_rep_wrap(~bias, ncol=1) +
    scale_color_manual(values=cols.type[4:3], name="") +
    geom_vline(xintercept=0, linetype=2) +
    geom_hline(yintercept=0, linetype=2) +
    xlab(expression(paste(italic(cis)," (log"[2]," fold-change)"))) +
    ylab(expression(paste(italic(trans)," (log"[2]," fold-change)"))) +
    theme(strip.background=element_blank(), strip.text=element_blank(),
            legend.position=c(0.1,0.1))
p2 = ggplot(cis_trans.sexbias.count, aes(x=effect2, y=n, fill=effect2)) +
    geom_bar(stat="identity", width=.75, show.legend=F) +
    geom_linerange(aes(ymin=`5%`, ymax=`95%`)) +
    scale_fill_manual(values=cols.type[4:3], name="") +
    facet_rep_wrap(~bias, ncol=1, strip.position="right") +
    xlab("") +
    ylab("number of genes") +
    background_grid(major="y", color.major="grey90") +
    theme(strip.background=element_rect(fill="grey75"))

plot_grid(p1, p2, ncol=2, rel_widths=c(1,0.7))
ggsave("suppl_cis_trans.sexbias.png")

# supplementary

Cbr_allele = df.species_sex_class %>%
    filter((class_male == "C. briggsae dominant" | class_male == "C. nigoni dominant" | class_male == "overdominant" | class_male == "underdominant") & chromosome != "X") %>%
    mutate(Cbr_allele = logFC_F1_cis_males < 0) %>%
    group_by(class_male, type_male, Cbr_allele) %>%
    count() %>%
    filter(type_male == "cis only" | type_male == "cis-trans (enhancing)" | type_male == "cis-trans (compensatory)")

ggplot(Cbr_allele, aes(x=type_male, y=n, fill=Cbr_allele)) +
    geom_bar(stat="identity", width=0.7) +
    scale_fill_manual(values=cols.class[5:4], name="allele", labels=c("C. nigoni","C. briggsae")) +
    scale_x_discrete(labels=c(
        expression(paste(italic("cis"), " only")),
        expression(paste(italic("cis-trans")," (enhancing)")),
        expression(paste(italic("cis-trans")," (compensatory)")))) +
    scale_y_continuous(expand=c(0,0,0,20)) +
    xlab("") +
    ylab("number of genes") +
    background_grid(major="x", color.major="grey75") +
    theme(legend.text=element_text(face="italic")) +
    coord_flip() +
    facet_wrap(~class_male, scales="free")

ggsave("suppl_dom_transgr_ASE.png")
ggsave("suppl_dom_transgr_ASE.pdf", device="pdf")
