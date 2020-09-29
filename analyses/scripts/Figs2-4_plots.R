
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

#### Plots

# load expression inheritance data
class.levels = c("no change","ambiguous","C. briggsae dominant", "C. nigoni dominant","additive","overdominant","underdominant")
class.legend = c("no change","ambiguous",expression(paste(italic("C. briggsae"), " dominant")),
                expression(paste(italic("C. nigoni"), " dominant")),"additive","overdominant","underdominant")
sex.levels = c("female","male","(F) hermaphrodite","(M) hermaphrodite")
df.inherit2 = read.csv("tables/expr_inheritance_logFC_F1_parents.csv")
df.inherit3 = read.csv("tables/expr_inheritance_chr_enrich.csv")
df.inherit2$class = factor(df.inherit2$class, levels=class.levels)
df.inherit3$class = factor(df.inherit3$class, levels=class.levels)
df.inherit2$sex = factor(df.inherit2$sex, levels=sex.levels)
df.inherit3$sex = factor(df.inherit3$sex, levels=sex.levels)

# load cis-trans regulatory data
type.levels = c("conserved","ambiguous","trans-only","cis-only", "cis + trans (enhancing)", "cis x trans (compensatory)", "cis-trans (compensatory)")
type.legend = c("conserved","ambiguous",expression(paste(italic("trans"), "-only")), expression(paste(italic("cis"),"-only")),
                expression(paste(italic("cis + trans"), " (enhancing)")), expression(paste(italic("cis x trans"), " (compensatory)")),
                expression(paste(italic("cis-trans"), " (compensatory)")))
df.AS = read.csv("tables/df.cis_trans.regulation_type.auto.csv")
df.AS2 = read.csv("tables/df.cis_trans.counts.chr_enrich.csv")
df.AS$type = factor(df.AS$type, levels=type.levels)
df.AS2$type = factor(df.AS2$type, levels=type.levels)
df.AS$sex = factor(df.AS$sex, levels=sex.levels)
df.AS2$sex = factor(df.AS2$sex, levels=sex.levels)

# get colors
# species

cols.inherit = rev(c("#A68481","#4F3B39","#6F7862","#B379AF","#4C5E91"))
cols.class = c("black","grey",cols.inherit)
# brewer.pal(10, "Paired") -> cols
# cols = cols[-c(7,8)]
# cols.class = c("black","grey",cols[ c(4,6,2,8,7) ])
cols.type1 = c("#D4C41C","#39B8E3","#3B9133","#A3627D","#DBA9BE")
cols.type = c("black","grey",cols.type1)
# cols.type[3] = "#09BC8A"
# cols.type[7]= "#C49991"
# cols.type = cols.type[c(1,2,4,3,5,6,7)]
# cols.type[4] = "#476A6F"

############################
## expression inheritance ##
############################

# biplot of expression difference between parents and F1
# colored by expression inheritance class
pA = ggplot(df.inherit2 %>% filter(sex == "female" | sex == "male"), aes(x=logFC.F1.vs.briggsae, y=logFC.F1.vs.nigoni, color=class)) +
	geom_point(show.legend=F, alpha=0.3) +
	scale_color_manual(values=cols.class) +
	geom_hline(yintercept=0, linetype=2) +
	geom_vline(xintercept=0, linetype=2) +
	xlab(expression(paste("log"[2],"(F1/",italic("Cbri"),")", sep=''))) + # bold("expr. diff."),
	ylab(expression(paste("log"[2],"(F1/",italic("Cnig"),")", sep=''))) + # bold("expr. diff."),
	background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
	facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	theme(strip.background = element_rect(fill="grey90", color=NA)) +
	theme(strip.background = element_blank(), strip.text=element_blank(),
				axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))
#p1

# F1 expression distance from origin (0,0)
pB = ggplot(df.inherit2 %>% filter(sex == "female" | sex == "male"), aes(y=class, x=F1_dist_from_zero, fill=class)) +
	geom_boxploth(alpha=0.7, outlier.shape=NA, width=0.5, show.legend=F) +
	geom_density_ridges(alpha=0.7, color=NA, show.legend=F) +
	scale_y_discrete(limits=rev(class.levels)) +
	scale_fill_manual(values=cols.class) +
	facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	xlab("F1 expression distance") +
	ylab("") +
	background_grid(major="x", minor="x", size.major = 0.2, colour.major = "grey75", size.minor = 0.2, colour.minor = "grey75") +
	xlim(0,7) +
	theme(strip.background = element_blank(), strip.text=element_blank(),
	 			axis.text.y=element_blank(),
	 			axis.ticks.y=element_blank())

# stacked barplot with gene counts for each class per chromosome
pC = ggplot(df.inherit3 %>% filter(sex == "male" | sex == "female"), aes(x=chromosome, y=n, fill=class)) +
	geom_bar(stat="identity", show.legend=F, alpha=0.7) +
    scale_y_continuous(breaks=seq(0,3000,500)) +
	scale_fill_manual(values=cols.class) +
	ylab("number of genes") +
	facet_rep_wrap(~sex, ncol=1, strip.position="right", scales="free") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	theme(strip.background = element_rect(fill="grey90", color=NA)) +
	theme(strip.background = element_blank(), strip.text=element_blank())
#p2

# for pD
# grey shade for each chromosome
rectsD = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-3.1,1.5,-3.1), y2=c(1.5,-3.1,1.5))

# gene enrichment of each expression inheritance class
pD = ggplot(df.inherit3 %>% filter(sex == "male" | sex == "female"), aes(x=chr.pseudo, y=log2(enrichment), fill=class)) +
 	geom_rect(data=rectsD, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
	geom_bar(stat="identity", position=position_dodge(), alpha=0.7, show.legend=T) +
	scale_y_continuous(breaks=seq(-3,2,1)) +
	scale_x_continuous(breaks=1:6, labels=df.inherit3$chromosome[1:6]) +
	scale_fill_manual(values=cols.class, labels=class.legend) +
	facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	background_grid(major="y", minor="y", size.major = 0.2, colour.major = "grey75", size.minor = 0.2, colour.minor = "grey75") +
	ylab(expression(paste("gene enrichment ","[log"[2]," odds ratio]",sep=""))) +
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
pE = ggplot(df.AS %>% filter(sex == "male" | sex == "female"), aes(logFC.sp, logFC.ase, color=type)) +
	geom_point(alpha=0.3, show.legend=F) +
	scale_color_manual(values=cols.type) +
	geom_hline(yintercept=0, linetype=2) +
	geom_vline(xintercept=0, linetype=2) +
	ylab(expression(paste(bold("ASE"), " [log"[2],"(",italic("Cnig"),"/",italic("Cbri"),")]", sep=''))) +
	xlab(expression(paste(bold("P"), " [log"[2],"(",italic("Cnig"),"/",italic("Cbri"),")]", sep=''))) +
	background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
	facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	theme(strip.background = element_rect(fill="grey90", color=NA)) +
	theme(strip.background = element_blank(), strip.text=element_blank(),
				axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

# boxplots of expression divergence between parent species
# for each regulation type
pF = ggplot(df.AS %>% filter(sex == "male" | sex == "female"), aes(y=type, x=abs(logFC.sp), fill=type)) +
	geom_boxploth(alpha=0.7, outlier.shape=NA, width=0.5, show.legend=F) +
	geom_density_ridges(alpha=0.7, color=NA, show.legend=F) +
	scale_y_discrete(limits=rev(type.levels)) +
	scale_fill_manual(values=cols.type) +
	facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	xlab(expression(paste("|log"[2],"-FC expression divergence|", sep=""))) +
	ylab("") +
	background_grid(major="x", minor="x", size.major = 0.2, colour.major = "grey75", size.minor = 0.2, colour.minor = "grey75") +
	xlim(0,7) +
	theme(strip.background = element_blank(), strip.text=element_blank(),
	 			axis.text.y=element_blank(),
	 			axis.ticks.y=element_blank())

# stacked barplot of counts of each regulation type per chromosome
pG = ggplot(df.AS %>% filter(sex == "male" | sex == "female"), aes(x=chromosome, fill=type)) +
	geom_bar(show.legend=F, alpha=0.7) +
	scale_fill_manual(values=cols.type) +
    scale_y_continuous(breaks=seq(0,3000,500)) +
	ylab("number of genes") +
	facet_rep_wrap(~sex, ncol=1, strip.position="right", scales="free_y") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	theme(strip.background = element_rect(fill="grey90", color=NA)) +
	theme(strip.background = element_blank(), strip.text=element_blank())

# for pH
# grey shade for each chromosome
rectsH = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-1,1,-1), y2=c(1,-1,1))

# gene enrichment of cis and trans regulation types in each chromosome
pH = ggplot(df.AS2 %>% filter(sex == "male" | sex == "female"), aes(x=chr.pseudo, y=log2(enrichment), fill=type)) +
 	geom_rect(data=rectsH, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
	geom_bar(stat="identity", position=position_dodge(), alpha=0.7, show.legend=T) +
	scale_y_continuous(breaks=seq(-1,1,0.5), labels=c("-1","","0","","1")) +
	scale_x_continuous(breaks=1:6, labels=df.AS2$chromosome[1:6]) +
	scale_fill_manual(values=cols.type, labels=type.legend) +
	facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	ylab(expression(paste("gene enrichment "," [log"[2]," odds ratio]",sep=""))) +
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
plot_grid(row1, row2, nrow=2)
# save to file
ggsave("figures/Fig2_F1_expression_inheritance.cis-trans.pdf", device="pdf", useDingbats=F)
ggsave("figures/Fig2_F1_expression_inheritance.png")

# save to file
# row1
# ggsave("figures/Fig2_F1_expression_inheritance.pdf", device="pdf", useDingbats=F)
# ggsave("figures/Fig2_F1_expression_inheritance.png")
# row2
# ggsave("figures/Fig3_cis-trans.pdf", device="pdf", useDingbats=F)
# ggsave("figures/Fig3_cis-trans.png")


# figure 3

class.vs.type = read.csv("tables/df.cis_trans.inherit.mat.csv")
class.vs.type$class = factor(class.vs.type$class, class.levels[ 3:7 ])
class.vs.type$type = factor(class.vs.type$type, type.levels[ 3:7 ])
class.vs.type$sex = factor(class.vs.type$sex, sex.levels)

f4_1 = ggplot(class.vs.type %>% filter(sex == "female" | sex == "male"), aes(x=class, y=type, fill=n)) +
    geom_tile() +
    xlab("") + ylab("") +
    scale_fill_viridis(direction=-1, name="number of \ngenes") +
    facet_rep_grid(sex~chromosome) +
    scale_y_discrete(labels=type.legend[3:7], limits=type.levels[3:7]) +
    scale_x_discrete(labels=class.legend[3:7][c(4:5,1:3)], limits=class.levels[3:7][c(4:5,1:3)]) +
    theme(strip.background = element_rect(fill="grey90", color=NA),
        strip.text.y = element_text(size=10),
        legend.title=element_text(size=10),
        legend.text =element_text(size=10),
        axis.text.y = element_text(size=10),
        #axis.text.x = element_blank())
		axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=10))
f4_2 = ggplot(class.vs.type %>% filter(sex == "(F) hermaphrodite" | sex == "(M) hermaphrodite"), aes(x=class, y=type, fill=n)) +
    geom_tile() +
    xlab("") + ylab("") +
    scale_fill_viridis(direction=-1, name="number of \ngenes") +
    facet_rep_grid(sex~chromosome) +
    scale_y_discrete(labels=type.legend[3:7], limits=type.levels[3:7]) +
    scale_x_discrete(labels=class.legend[3:7][c(4:5,1:3)], limits=class.levels[3:7][c(4:5,1:3)]) +
    theme(strip.background = element_rect(fill="grey90", color=NA),
        #strip.background.x = element_blank(),
        #strip.text.x = element_blank(),
        strip.text.y = element_text(size=10),
        #strip.background = element_rect(fill="grey90", color=NA),
        legend.title=element_text(size=10),
        legend.text =element_text(size=10),
        axis.text.y = element_text(size=10),
        #axis.text.x = element_blank())
		axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=10))

grob = as_grob(f4_1)
grob$grobs[[ which(grob$layout$name %in% "panel-2-6") ]] = nullGrob()
grob$grobs[[ which(grob$layout$name %in% "axis-l-2-6") ]] = nullGrob()
grob$grobs[[ which(grob$layout$name %in% "axis-b-6-1") ]] = nullGrob()
grob$layout[ which(grob$layout$name %in% "axis-b-6"), c("t","b")] = grob$layout[ which(grob$layout$name %in% "axis-b-6"), c("t","b")] - 3
f4_1 = as.ggplot(grob)

grob = as_grob(f4_2)
grob$grobs[[ which(grob$layout$name %in% "panel-2-6") ]] = nullGrob()
grob$grobs[[ which(grob$layout$name %in% "axis-l-2-6") ]] = nullGrob()
grob$grobs[[ which(grob$layout$name %in% "axis-b-6-1") ]] = nullGrob()
grob$layout[ which(grob$layout$name %in% "axis-b-6"), c("t","b")] = grob$layout[ which(grob$layout$name %in% "axis-b-6"), c("t","b")] - 3
f4_2 = as.ggplot(grob)

f4_1
ggsave("figures/Fig3_matrix_class_type_counts.pdf", device="pdf", useDingbats=F)
ggsave("figures/Fig3_matrix_class_type_counts.png")

f4_2
ggsave("figures/Suppl_Fig_matrix_class_type_counts_herm.pdf", device="pdf", useDingbats=F)
ggsave("figures/Suppl_Fig_matrix_class_type_counts_herm.png")

plot_grid(f4_1 + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
        f4_2 + theme(plot.margin = unit(c(-.7, 0, 0, 0), "cm")), rel_heights=c(0.5,0.6), nrow=2)

# ggsave("figures/Fig3_matrix_class_type_counts.pdf", device="pdf", useDingbats=F)
# ggsave("figures/Fig3_matrix_class_type_counts.png")

# Supplementary fig. S4

class_and_type.auto = read.csv("tables/df.cis_trans.inherit.merge.csv")
additive = na.omit(filter(class_and_type.auto, class == "additive"))
additive.count =  as.data.frame(additive %>% group_by(sex, type) %>% count())
boot = data.frame()
for (i in 1:1000){
    boot = rbind(boot, as.data.frame(additive %>% filter(sex == "female") %>% sample_frac(replace=TRUE) %>% group_by(sex,type) %>% count()))
    boot = rbind(boot, as.data.frame(additive %>% filter(sex == "male") %>% sample_frac(replace=TRUE) %>% group_by(sex,type) %>% count()))
}
additive.count = cbind(additive.count %>% filter(sex == "female" | sex == "male"), as.data.frame(boot %>% group_by(sex, type) %>% summarize(upper = quantile(n, 0.95), lower = quantile(n, 0.05)))[,c(3,4)])
ggplot(additive.count, aes(x=type, y=n, fill=type)) +
    geom_bar(stat="identity", show.legend=F, alpha=0.7) +
    geom_linerange(aes(ymin=lower, ymax=upper)) +
    scale_fill_manual(values=cols.type[c(2,5,7,4,3)], name="") +
    scale_y_continuous(expand=c(0,0)) +
    facet_rep_wrap(~sex) +
    xlab("") +
    ylab("number of genes") +
    coord_flip() +
    background_grid(major="x") +
    theme(axis.ticks.y=element_blank())

ggsave("figures/suppl_figS4_cis_in_additive.png")
ggsave("figures/suppl_figS4_cis_in_additive.pdf", device="pdf", useDingbats=F)
