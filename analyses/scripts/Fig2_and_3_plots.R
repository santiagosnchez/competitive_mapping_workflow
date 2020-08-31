
# libraries
library(ggplot2)
library(cowplot)
library(lemon)
library(ggridges)
library(ggstance)
library(ggsci)
library(tidyr)
library(dplyr)
library(gridExtra)
library(RColorBrewer)
library(viridis)
theme_set(theme_cowplot())

#### Plots

# load expression inheritance data
class.levels = c("no change","ambiguous","additive","C. briggsae dominant", "C. nigoni dominant","overdominant","underdominant")
class.legend = c("no change","ambiguous","additive",expression(paste(italic("C. briggsae"), " dominant")),
                expression(paste(italic("C. nigoni"), " dominant")),"overdominant","underdominant")
df.inherit2 = read.csv("expr_inheritance_logFC_F1_parents.csv")
df.inherit3 = read.csv("expr_inheritance_chr_enrich.csv")
df.inherit2$class = factor(df.inherit2$class, levels=class.levels)
df.inherit3$class = factor(df.inherit3$class, levels=class.levels)

# load cis-trans regulatory data
type.levels = c("conserved","ambiguous","trans only","cis only", "cis-trans (enhancing)", "cis-trans (compensatory)")
type.legend = c("conserved","ambiguous",expression(paste(italic("trans"), " only")), expression(paste(italic("cis")," only")),
                expression(paste(italic("cis-trans"), " (enhancing)")), expression(paste(italic("cis-trans"), " (compensatory)")))
df.AS = read.csv("df.cis_trans.regulation_type.csv")
df.AS2 = read.csv("df.cis_trans.counts.chr_enrich.csv")
df.AS$type = factor(df.AS$type, levels=type.levels)
df.AS2$type = factor(df.AS2$type, levels=type.levels)

# get colors
# species
brewer.pal(10, "Paired") -> cols
cols = cols[-c(7,8)]
cols.class = c("black","grey",cols[ c(4,6,2,8,7) ])
cols.type = c("black","grey",pal_jama()(4))
cols.type[3] = "#09BC8A"

############################
## expression inheritance ##
############################

# biplot of expression difference between parents and F1
# colored by expression inheritance class
pA = ggplot(df.inherit2, aes(x=logFC.F1.vs.briggsae, y=logFC.F1.vs.nigoni, color=class)) +
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
pB = ggplot(df.inherit2, aes(y=class, x=F1_dist_from_zero, fill=class)) +
	geom_boxploth(alpha=0.7, outlier.shape=NA, width=0.5, show.legend=F) +
	geom_density_ridges(alpha=0.7, color=NA, show.legend=F) +
	scale_y_discrete(limits=rev(class.levels)) +
	scale_fill_manual(values=cols.class) +
	facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	xlab("F1 expression distance") +
	ylab("") +
	background_grid(major="x", minor="x", size.major = 0.2, colour.major = "grey75", size.minor = 0.2, colour.minor = "grey75") +
	xlim(0,10) +
	theme(strip.background = element_blank(), strip.text=element_blank(),
	 			axis.text.y=element_blank(),
	 			axis.ticks.y=element_blank())

# stacked barplot with gene counts for each class per chromosome
pC = ggplot(df.inherit3, aes(x=chromosome, y=n, fill=class)) +
	geom_bar(stat="identity", show.legend=F, alpha=0.7) +
	scale_fill_manual(values=cols.class) +
	ylab("number of genes") +
	facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	theme(strip.background = element_rect(fill="grey90", color=NA)) +
	theme(strip.background = element_blank(), strip.text=element_blank())
#p2

# for pD
# grey shade for each chromosome
rectsD = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-3.7,2,-3.7), y2=c(2,-3.7,2))

# gene enrichment of each expression inheritance class
pD = ggplot(df.inherit3, aes(x=chr.pseudo, y=log2(enrichment), fill=class)) +
 	geom_rect(data=rectsD, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
	geom_bar(stat="identity", position=position_dodge(), alpha=0.7, show.legend=T) +
	scale_y_continuous(breaks=-3:2) +
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

##########################################
## cis trans expression regulation type ##
##########################################

# biplot of allele-specific expression in F1 (y-axis) vs
# expression difference between parents (cis + trans)
pE = ggplot(df.AS, aes(logFC_parents, logFC_F1_cis, color=type)) +
	geom_point(alpha=0.7, show.legend=F) +
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
pF = ggplot(df.AS, aes(y=type, x=abs(logFC_parents), fill=type)) +
	geom_boxploth(alpha=0.7, outlier.shape=NA, width=0.5, show.legend=F) +
	geom_density_ridges(alpha=0.7, color=NA, show.legend=F) +
	scale_y_discrete(limits=rev(type.levels)) +
	scale_fill_manual(values=cols.type) +
	facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	xlab("|expression divergence|") +
	ylab("") +
	background_grid(major="x", minor="x", size.major = 0.2, colour.major = "grey75", size.minor = 0.2, colour.minor = "grey75") +
	xlim(0,10) +
	theme(strip.background = element_blank(), strip.text=element_blank(),
	 			axis.text.y=element_blank(),
	 			axis.ticks.y=element_blank())

# stacked barplot of counts of each regulation type per chromosome
pG = ggplot(df.AS2, aes(x=chromosome, y=n, fill=type)) +
	geom_bar(stat="identity", show.legend=F, alpha=0.7) +
	scale_fill_manual(values=cols.type) +
	ylab("number of genes") +
	facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	theme(strip.background = element_rect(fill="grey90", color=NA)) +
	theme(strip.background = element_blank(), strip.text=element_blank())

# for pH
# grey shade for each chromosome
rectsH = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-1.5,1.5,-1.5), y2=c(1.5,-1.5,1.5))

# gene enrichment of cis and trans regulation types in each chromosome
pH = ggplot(df.AS2, aes(x=chr.pseudo, y=log2(enrichment), fill=type)) +
 	geom_rect(data=rectsH, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
	geom_bar(stat="identity", position=position_dodge(), alpha=0.7, show.legend=T) +
	scale_y_continuous(breaks=-1:1) +
	scale_x_continuous(breaks=1:6, labels=df.AS2$chromosome[1:6]) +
	scale_fill_manual(values=cols.type, labels=type.legend) +
	facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	background_grid(major="y", minor="y", size.major = 0.2, colour.major = "grey75", size.minor = 0.2, colour.minor = "grey75") +
	ylab(expression(paste("gene enrichment "," [log"[2]," odds ratio]",sep=""))) +
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
legend.type = as_grob(get_legend(pH))
pH = pH + theme(legend.position = "none")

# layout plots into a grid
row1 = plot_grid(pA + ggtitle("A") + theme(plot.margin = unit(c(0, 0, 0.2, 0), "cm")),
								pB + ggtitle("B") + theme(plot.margin = unit(c(0, 0.5, 0.2, 0), "cm")),
								pC + ggtitle("C") + theme(plot.margin = unit(c(0, 0.5, 0.2, 0), "cm")),
								pD + ggtitle("D") + theme(plot.margin = unit(c(0, 0.5, 0.2, 0), "cm")),
								legend.class, #+ theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
								nrow=1, align="both", rel_widths=c(0.4,0.4,0.35,0.7,0.3))
row2 = plot_grid(pE + ggtitle("E") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
								pF + ggtitle("F") + theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")),
								pG + ggtitle("G") + theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")),
								pH + ggtitle("H") + theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")),
								legend.type, # + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
								nrow=1, align="both", rel_widths=c(0.4,0.4,0.35,0.7,0.3))
plot_grid(row1, row2, nrow=2)

# save to file
ggsave("Fig2_F1_expression_inheritance_and_cis-trans.pdf", device="pdf")
ggsave("Fig2_F1_expression_inheritance_and_cis-trans.png")

# figure 3

class.vs.type = read.csv("class_vs_type.csv")
class.vs.type$class = factor(class.vs.type$class, class.levels)
class.vs.type$type = factor(class.vs.type$type, type.levels)

ggplot(class.vs.type, aes(x=type, y=class, fill=n)) +
    geom_tile() +
    xlab("") + ylab("") +
    scale_fill_viridis(direction=-1, name="number of \ngenes") +
    facet_rep_grid(sex~chromosome) +
    scale_x_discrete(labels=type.legend) +
    scale_y_discrete(labels=class.legend) +
    theme(strip.background = element_rect(fill="grey90", color=NA),
        legend.title=element_text(size=10),
        legend.text =element_text(size=10),
        axis.text.y = element_text(size=10),
		axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=10))

ggsave("Fig3_matrix_class_type_counts.pdf", device="pdf")
ggsave("Fig3_matrix_class_type_counts.png")

# Supplementary fig. S4

class_and_type.auto = read.csv("class_and_type.auto.csv", row.name=1)
additive = filter(class_and_type.auto, class == "additive")
additive.count =  as.data.frame(additive %>% group_by(sex, type) %>% count())
boot = data.frame()
for (i in 1:1000){
    boot = rbind(boot, as.data.frame(additive %>% filter(sex == "female") %>% sample_frac(replace=TRUE) %>% group_by(sex,type) %>% count()))
    boot = rbind(boot, as.data.frame(additive %>% filter(sex == "male") %>% sample_frac(replace=TRUE) %>% group_by(sex,type) %>% count()))
}
additive.count = cbind(additive.count, as.data.frame(boot %>% group_by(sex, type) %>% summarize(upper = quantile(n, 0.95), lower = quantile(n, 0.05)))[,c(3,4)])
ggplot(additive.count, aes(x=type, y=n, fill=type)) +
    geom_bar(stat="identity", show.legend=F) +
    geom_linerange(aes(ymin=lower, ymax=upper)) +
    scale_fill_manual(values=rev(cols.type[c(3,5,4,2)]), name="") +
    scale_y_continuous(expand=c(0,0)) +
    facet_rep_wrap(~sex) +
    xlab("") +
    ylab("number of genes") +
    coord_flip() +
    theme(axis.ticks.y=element_blank())

ggsave("suppl_figS4_cis_in_additive.png")
ggsave("suppl_figS4_cis_in_additive.pdf", device="pdf", useDingbats=F)
