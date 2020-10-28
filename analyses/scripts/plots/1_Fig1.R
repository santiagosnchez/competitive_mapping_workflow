# fig1
library(ggplot2)
library(cowplot)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(edgeR)
library(lemon)
library(ggstance)
library(ggrepel)
library(ggsci)
library(ggtext)
theme_set(theme_cowplot())

# get colors
# species
# brewer.pal(10, "Paired") -> cols
# cols = cols[-c(7,8)]
# cols.spp = cols[c(6,2,8)]
# cols.spp2 = cols[c(5,6,1,2)]

#cols.inherit_old = rev(c("#A68481","#FCA79F","#6F7862","#B379AF","#4C5E91"))
cols.inherit = c("#4C5E91","#B379AF","#3F784C","#C17817","#533E2D")
cols.spp = cols.inherit[1:3]

# sex
cols.sex = c("#80CDC1","#018571","#A68481","#BEBEBE")

# read data
df.mds = read.csv("tables/expression_dist.mds_coords.csv", row.names=1) # MDS distances
expr.div = read.csv("tables/df.expr_div.per_sex.csv") # expression divergence (logFC)
ase.div = read.csv("tables/df.ase.per_sex.csv") # allele-specific expression (logFC)
df.expr_div.chr.enrich = read.csv("tables/df.enrichment.exprdiv.DE.sex.csv")
df.ase.chr.enrich = read.csv("tables/df.enrichment.ase.DE.sex.csv")
df.DE.sex.enrich = read.csv("tables/DE_per_sex_chr_enrichment.csv")
df.DE.sex.counts = read.csv("tables/df.counts.perc.DE.sex.csv")
dds.rlog = read.csv("tables/rlog-transformed_exprdata.csv", row.names=1)

# process DESeq rlog data for females
# Dosage Compensation
dds.rlog.females = dds.rlog[, grep("_F", colnames(dds.rlog))]
dds.rlog.females.m = data.frame(Cbr=rowMeans(dds.rlog.females[,1:3]), F1=rowMeans(dds.rlog.females[,4:6]), Cni=rowMeans(dds.rlog.females[,7:9]))
dds.rlog.females.m$chromosome = sub("\\..*","", rownames(dds.rlog.females.m))
dds.rlog.females.m$autoX = c("autosomes","X")[ factor(dds.rlog.females.m$chromosome == "X") ]
dds.rlog.females.m = dds.rlog.females.m %>% group_by(chromosome) %>%
    mutate(pseudo.pos = seq_along(chromosome)) %>%
    mutate(rel.pos = scales:::rescale(pseudo.pos, c(0,1)))
df.dds.rlog.females.m = pivot_longer(dds.rlog.females.m, c(-chromosome, -rel.pos, -pseudo.pos, -autoX), names_to="species", values_to="log_norm_expr")
df.dds.rlog.females.m.sum = df.dds.rlog.females.m %>%
    group_by(species, autoX) %>%
    summarize(log_norm_expr.m=mean(log_norm_expr), upper=mean(log_norm_expr)+sd(log_norm_expr), lower=mean(log_norm_expr)-sd(log_norm_expr))

#df.DE.sex = read.csv("df.DE.sex.csv", row.names=1) # sex-bias differential expression
#df.DE.sex2 = gather(df.DE.sex, species, sex, -chromosome) # clean DF for plotting

# expression distance auto vs x

df.mds$sex = as.character(df.mds$sex)
df.mds$sex[1:3] = "H"
mds.auto.plt = ggplot(df.mds, aes(x=dim1, y=dim2, color=sex, shape=species)) +
	geom_point(size=5, alpha=0.7) +
	geom_vline(xintercept=0, linetype=2) +
	geom_hline(yintercept=0, linetype=2) +
	scale_shape_manual(values=c(16,17,15), labels=c(
		expression(italic("C. briggsae")),
		expression(italic("C. nigoni")),
		"F1")) +
	scale_y_continuous(limits=c(-3,3)) +
	scale_x_continuous(limits=c(-4,4)) +
	scale_color_manual(values=cols.sex, labels=c("female","hermaphrodite","male")) +
	ylab(expression(paste("log"[2]," expression (dim 2)", sep=""))) +
	xlab(expression(paste("log"[2]," expression (dim 1)",sep=""))) +
	background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
    #annotate(geom = "curve", x=2.5, y=-3.55, xend=0.6, yend=-2.7, curvature = -0.2, arrow = arrow(length=unit(2, "mm"))) +
    #annotate(geom = "text", x=2.6, y=-3.55, label="H", hjust=0) +
	theme(strip.background = element_rect(fill="grey90", color=NA)) +
	theme(axis.title=element_text(size=12), axis.text=element_text(size=10)) +
	theme(legend.text.align = 0, legend.title=element_blank(), legend.text=element_text(size=10))
	#theme(panel.border=element_blank(), axis.line=element_line()) +
	#coord_capped_cart(bottom='both', left='both') +
	#facet_rep_wrap(~chr)
#mds.auto.plt

dc.plt = ggplot(df.dds.rlog.females.m, aes(x=autoX, y=log_norm_expr, fill=species)) +
    geom_violin(draw_quantiles=0.5, color=NA, alpha=0.7, show.legend=T) +
    geom_point(data=df.dds.rlog.females.m.sum, aes(y=log_norm_expr.m), position=position_dodge(.9), color="black", show.legend=F) +
    geom_linerange(data=df.dds.rlog.females.m.sum, aes(y=log_norm_expr.m, ymax=upper, ymin=lower), position=position_dodge(.9)) +
    labs(x="", y=expression(paste("log"[2], "normalized counts (rlog)"))) +
    background_grid(major="y") +
    scale_fill_manual(values=cols.spp, name="", labels=c(
		expression(italic("C. briggsae")),
		expression(italic("C. nigoni")),
		"F1")) +
    theme(legend.text.align = 0, legend.title=element_blank(), legend.text=element_text(size=10)) +
    facet_wrap(~"female/hermaphrodite data only")


plot_grid(mds.auto.plt, dc.plt, ncol=1)
# supplementary
ggsave("figures/suppl_FigS1_expr_MDS_2.pdf", device="pdf", useDingbats=FALSE)
ggsave("figures/suppl_FigS1_expr_MDS_2.png", device="png")

## expression divergence histogram

exprdiff = rbind(data.frame(logFC = expr.div$logFC, species_DE = expr.div$species_DE, sp.vs.allele="expression divergence", sex = expr.div$sex),
    data.frame(logFC = ase.div$logFC.ase, species_DE = ase.div$species_DE, sp.vs.allele="ASE", sex = ase.div$sex))
exprdiff$species_DE = as.character(exprdiff$species_DE)
exprdiff$species_DE[ grepl("no sig", exprdiff$species_DE) ] = "no sig"

exprdiff.plt = ggplot(exprdiff %>% filter(sex != "hermaphrodite"), aes(logFC, fill=species_DE)) +
    geom_vline(xintercept=0) +
	geom_histogram(bins=50, show.legend=F, alpha=0.7) +
	scale_y_continuous(breaks=seq(0,4000,500)) +
	scale_fill_manual(values=c(cols.spp[c(1,2)],"grey"))  +
	#theme(legend.position=c(0.02,0.9), legend.text=element_text(size=11)) +
	#theme(legend.box.background=element_rect(fill="white", color="white")) +
	xlab(expression(paste(log[2]," fold-change",sep=""),"between species")) +
	ylab("number of genes") +
	theme(axis.title=element_text(size=12), axis.text=element_text(size=10), plot.margin=unit(c(1,1,2.5,1),"lines")) +
	#theme(legend.text.align = 0) +
	background_grid(major="y", minor="none", size.minor = 0.2, size.major = 0.2, colour.major = "grey70", colour.minor = "grey70") +
  facet_rep_grid(sex~sp.vs.allele) +
  theme(strip.background.y = element_blank(), strip.text.y=element_blank(), strip.background.x = element_rect(fill="grey90", color=NA)) +
  xlim(-10,10)

## barplot of between species DE enrichment across chromosomes

df.div.chr.enrich = rbind(cbind(df.expr_div.chr.enrich, sp.vs.allele="expression divergence"), cbind(df.ase.chr.enrich, sp.vs.allele="ASE"))

# add text column for significance asterisk
df.div.chr.enrich$sig = ""
df.div.chr.enrich$sig[ df.div.chr.enrich$p.value < 0.05 & abs(log2(df.div.chr.enrich$enrichment)) > 0.5 ] = "*"

# chromosome delimitation rectangles
rectsDE = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-2.1,1.5,-2.1), y2=c(1.5,-2.1,1.5))

enrich.div.chr.plt  = ggplot(df.div.chr.enrich %>% filter(sex != "hermaphrodite"), aes(x=chr.pseudo, y=log2(enrichment), fill=species_DE)) +
    geom_rect(data=rectsDE, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
    #geom_hline(yintercept=0) +
	geom_bar(stat="identity", position=position_dodge(width=0.9), show.legend=T, size=3, alpha=0.7) +
	geom_text(aes(y=log2(enrichment)/2, label=sig), size=5, position=position_dodge(width=0.9), show.legend=F) +
	scale_fill_manual(values=c(cols.spp[c(1,2)],"grey"), name="", labels=c(
		expression(italic("C. briggsae")),
		expression(italic("C. nigoni")),
		"no DE")) +
    scale_y_continuous(breaks=seq(-2,2,0.5), limits=c(-2,2), expand=c(0,0), labels=c("-2","","-1","","0","","1","","2")) +
	scale_x_continuous(breaks=1:6, labels=as.character(df.div.chr.enrich[1:6,"chromosome"])) +
	xlab("chromosome") +
	ylab(expression(atop("enrichment of genes with higher expression",paste("[log"[2]," odds ratio]",sep="")))) +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	facet_rep_grid(sex~sp.vs.allele) +
	theme(axis.title=element_text(size=12), axis.text=element_text(size=10)) +
    theme(legend.title=element_blank()) + #,
        # legend.box.background=element_rect(fill="white", color=NA)) +
	theme(strip.background = element_rect(fill="grey90", color=NA),
				axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)),
				legend.text.align = 0, legend.position=c(0.83,0.40),
				legend.text=element_text(size=10))
#enrich.div.chr.plt
## barplot of sex-bias expression between species and F1s

# add text offset
df.DE.sex.counts = df.DE.sex.counts %>% filter(sex != "no expression")
df.DE.sex.counts$text_offset = c(500,-400,-500,500,-500,500,500,-500,500,-500)

sex.de.plt = ggplot(df.DE.sex.counts, aes(x=species, y=count, color=sex, group=sex)) +
	#geom_bar(stat='identity', position=position_dodge(), alpha=0.7) +
	geom_point(shape=16, size=3, show.legend=F) +
	geom_line(show.legend=F) +
	geom_text(aes(y=count+text_offset, label=paste0(round(perc,0),"%")), color="black", size=3, hjust=0.5, show.legend=F) +
	scale_color_manual(values=cols.sex, name="", labels=c("female-biased","hermaphrodite","male-biased","sex-neutral")) +
	scale_y_continuous(limits=c(0,7000), expand=c(0,0),
											breaks=seq(0,7000,1000),
											labels=c("0","","2000","","4000","","6000","")) +
	scale_x_discrete(limits=c("Cbr","F1","Cni"),
		labels=c(
		expression(italic("C. briggsae")),
		"F1",
		expression(italic("C. nigoni")))) +
	ylab("number of genes") + xlab("") +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	theme(axis.title=element_text(size=12), axis.text=element_text(size=10)) +
  theme(legend.text.align = 0,
		#legend.title=element_text(face="bold", size=10),
		legend.text=element_text(size=10))
sex.de.plt

## barplot of sex-bias enrichment per chromosome

# make a strip labeller
species.labels = c("C. briggsae","C. nigoni","F1")
names(species.labels) = c("Cbr","Cni","F1")

# add text column for significance asterisk
df.DE.sex.enrich$sig = ""
df.DE.sex.enrich$sig[ df.DE.sex.enrich$p.value < 0.05 & abs(log2(df.DE.sex.enrich$enrichment)) > 0.5 ] = "*"

# chromosome delimitation rectangles
rectsSex = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-2,2,-2), y2=c(2,-2,2))

enrich.sex.plt  = ggplot(df.DE.sex.enrich, aes(x=chr.pseudo, y=log2(enrichment), fill=sex)) +
    geom_rect(data=rectsSex, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
    #geom_hline(yintercept=0) +
    geom_bar(stat="identity", position=position_dodge(width=0.9), show.legend=T, size=3, alpha=0.7) +
    geom_text(aes(y=log2(enrichment)/2, label=sig), size=5, fontface="bold", position=position_dodge(width=0.9), show.legend=F) +
    scale_fill_manual(values=cols.sex, name="", labels=c("female","hermaphrodite","male","sex-neutral")) +
    scale_y_continuous(breaks=seq(-2,2,0.5), limits=c(-2,2,1), expand=c(0,0), labels=c("-2","","-1","","0","","1","","2")) +
	scale_x_continuous(breaks=1:6, labels=as.character(df.DE.sex.enrich[1:6,"chromosome"])) +
	xlab("chromosome") +
	ylab(expression(atop("enrichment of sex-biased genes",paste("[log"[2]," odds ratio]",sep="")))) +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	facet_rep_wrap(~species, nrow=1, labeller=labeller(species = species.labels)) +
	theme(axis.title=element_text(size=12), axis.text=element_text(size=10)) +
	theme(strip.background = element_rect(fill="grey90", color=NA),
				strip.text=element_text(face="italic"),
				axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)),
                legend.text.align = 0, legend.position=c(0.35,0.25), legend.text=element_text(size=10))
enrich.sex.plt
## plot grid

plot_grid(
    plot_grid(exprdiff.plt + ggtitle("A") + theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "cm")),
              enrich.div.chr.plt + ggtitle("B") + theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = .5), "cm")),
              align="h", rel_widths=c(0.4,0.6)),
    plot_grid(sex.de.plt + ggtitle("C") + theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "cm")),
              enrich.sex.plt + ggtitle("D") + theme(plot.margin = unit(c(t = 0, r = 0, b = 0, l = .5), "cm")),
              align="h", rel_widths=c(0.3,0.7), ncol=2),
          nrow=2, rel_heights=c(1,0.7))

# plot_grid(plot_grid(
#   plot_grid(mds.auto.plt + ggtitle("A") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
#             sex.de.plt + ggtitle("D") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
#             ncol=1),
#   plot_grid(expression_div_hist + ggtitle("B") + theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")),
#             enrich.DE.chr.plt + ggtitle("C") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
#             align='h', rel_widths=c(0.4,0.6)),
#             ncol=2, align="h", rel_widths=c(0.3,0.6)),
#   enrich.sex.plt + ggtitle("E") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
#   nrow=2, align="v", rel_heights=c(0.6,0.4))

## save image

ggsave("figures/Fig1_expression_divergence_DE_sex.pdf", device="pdf", useDingbats=FALSE)
ggsave("figures/Fig1_expression_divergence_DE_sex.png", device="png")
