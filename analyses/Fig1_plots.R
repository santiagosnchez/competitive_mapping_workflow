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
brewer.pal(10, "Paired") -> cols
cols = cols[-c(7,8)]
cols.spp = cols[c(6,2,8)]
cols.spp2 = cols[c(5,6,1,2)]
# sex
cols.sex = c("#80CDC1","#018571","#BEBEBE")
cols.sex2 = colorRampPalette(cols.sex)(6)

# read data
df.mds = read.csv("expression_dist.mds_coords.csv", row.names=1) # MDS distances
tt.br.ni = read.csv("df.expr_div.per_sex.csv") # expression divergence (logFC)
df.DE.chr.enrich = read.csv("df.enrichment.DE.sex.csv")
df.DE.sex.enrich = read.csv("DE_per_sex_chr_enrichment.csv")
df.DE.sex.counts = read.csv("df.counts.perc.DE.sex.csv")
#df.DE.sex = read.csv("df.DE.sex.csv", row.names=1) # sex-bias differential expression
#df.DE.sex2 = gather(df.DE.sex, species, sex, -chromosome) # clean DF for plotting

# expression distance auto vs x

mds.auto.plt = ggplot(df.mds, aes(x=dim1, y=dim2, color=species, shape=sex)) +
	geom_point(size=3, alpha=0.7) +
	geom_vline(xintercept=0, linetype=2) +
	geom_hline(yintercept=0, linetype=2) +
	scale_color_manual(values=cols.spp, labels=c(
		expression(italic("C. briggsae")),
		expression(italic("C. nigoni")),
		"F1")) +
	scale_y_continuous(limits=c(-5.1,5.1)) +
	scale_x_continuous(limits=c(-5.1,5.1)) +
	scale_shape_manual(values=c(1,2,0), labels=c("F/H","M"), drop=FALSE) +
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

## expression divergence histogram

expression_div_hist = ggplot(tt.br.ni, aes(logFC, fill=species_DE)) +
	geom_histogram(bins=50, show.legend=F, alpha=0.7) +
	scale_y_continuous(breaks=seq(0,2500,500)) +
	scale_fill_manual(values=cols.spp2)  +
	#theme(legend.position=c(0.02,0.9), legend.text=element_text(size=11)) +
	#theme(legend.box.background=element_rect(fill="white", color="white")) +
	xlab(expression(atop(paste(log[2]," expression divergence",sep=""),"between species"))) +
	ylab("number of genes") +
	theme(axis.title=element_text(size=12), axis.text=element_text(size=10), plot.margin=unit(c(1,1,2.5,1),"lines")) +
	#theme(legend.text.align = 0) +
	background_grid(major="y", minor="none", size.minor = 0.2, size.major = 0.2, colour.major = "grey70", colour.minor = "grey70") +
  facet_rep_wrap(~sex, strip.position="right", nrow=2) +
  theme(strip.background = element_blank(), strip.text=element_blank()) +
  xlim(-10,10)

## barplot of between species DE enrichment across chromosomes

# add text column for significance asterisk
df.DE.chr.enrich$sig = ""
df.DE.chr.enrich$sig[ df.DE.chr.enrich$p.value < 0.05 & abs(log2(df.DE.chr.enrich$enrichment)) > 0.5 ] = "*"

# chromosome delimitation rectangles
rectsDE = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-2.1,1.5,-2.1), y2=c(1.5,-2.1,1.5))

enrich.DE.chr.plt  = ggplot(df.DE.chr.enrich, aes(x=chr.pseudo, y=log2(enrichment), fill=species_DE)) +
  geom_rect(data=rectsDE, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
	geom_bar(stat="identity", position=position_dodge(width=0.9), show.legend=T, size=3, alpha=0.7) +
	geom_text(aes(y=log2(enrichment)/2, label=sig), size=5, position=position_dodge(width=0.9), show.legend=F) +
	scale_fill_manual(values=c(cols.spp[c(1,2)],"grey"), name="", labels=c(
		expression(italic("C. briggsae")),
		expression(italic("C. nigoni")),
		"no DE")) +
  scale_y_continuous(breaks=seq(-2,1.5,1), limits=c(-2.1,1.5)) +
	scale_x_continuous(breaks=1:6, labels=as.character(df.DE.chr.enrich[1:6,"chromosome"])) +
	xlab("chromosome") +
	ylab(expression(atop("enrichment of genes with higher expression",paste("[log"[2]," odds ratio]",sep="")))) +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	facet_rep_wrap(~sex, ncol=1, strip.position="right") +
	theme(axis.title=element_text(size=12), axis.text=element_text(size=10)) +
	theme(strip.background = element_rect(fill="grey90", color=NA),
				axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)),
				legend.text.align = 0, legend.position=c(0.075,0.15),
				legend.text=element_text(size=10))

## barplot of sex-bias expression between species and F1s

# add text offset
df.DE.sex.counts$text_offset = c(-500,500,500,-500,500,500,-500,500,-500)

sex.de.plt = ggplot(df.DE.sex.counts, aes(x=species, y=n, color=sex, group=sex)) +
	#geom_bar(stat='identity', position=position_dodge(), alpha=0.7) +
	geom_point(shape=16, size=3, show.legend=T) +
	geom_line(show.legend=F) +
	geom_text(aes(y=n+text_offset, label=paste0(round(perc,0),"%")), color="black", size=2.5, hjust=0.5, show.legend=F) +
	scale_color_manual(values=cols.sex, name="", labels=c("female-biased","male-biased","sex-neutral")) +
	scale_y_continuous(limits=c(0,7300), expand=c(0,0),
											breaks=seq(0,7000,1000),
											labels=c("","1000","","3000","","5000","","7000")) +
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

## barplot of sex-bias enrichment per chromosome

# make a strip labeller
species.labels = c("C. briggsae","C. nigoni","F1")
names(species.labels) = c("Cbr","Cni","F1")

# add text column for significance asterisk
df.DE.sex.enrich$sig = ""
df.DE.sex.enrich$sig[ df.DE.sex.enrich$p.value < 0.05 & abs(log2(df.DE.sex.enrich$enrichment)) > 0.5 ] = "*"

# chromosome delimitation rectangles
rectsSex = data.frame(chr.pseudo=c(1,3,5), x1=c(0.5,2.5,4.5), x2=c(1.5,3.5,5.5), y1=c(-2.75,2,-2.75), y2=c(2,-2.75,2))

enrich.sex.plt  = ggplot(df.DE.sex.enrich, aes(x=chr.pseudo, y=log2(enrichment), fill=sex)) +
  geom_rect(data=rectsSex, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), fill="grey90", alpha=0.5, inherit.aes=F) +
	geom_bar(stat="identity", position=position_dodge(width=0.9), show.legend=F, size=3, alpha=0.7) +
	geom_text(aes(y=log2(enrichment)/2, label=sig), size=5, fontface="bold", position=position_dodge(width=0.9), show.legend=F) +
	scale_fill_manual(values=cols.sex) +
  scale_y_continuous(breaks=seq(-2,2,1), limits=c(-2.75,2)) +
	scale_x_continuous(breaks=1:6, labels=as.character(df.DE.sex.enrich[1:6,"chromosome"])) +
	xlab("chromosome") +
	ylab(expression(atop("enrichment of sex-biased genes",paste("[log"[2]," odds ratio]",sep="")))) +
	background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
	facet_rep_wrap(~species, nrow=1, labeller=labeller(species = species.labels)) +
	theme(axis.title=element_text(size=12), axis.text=element_text(size=10)) +
	theme(strip.background = element_rect(fill="grey90", color=NA),
				strip.text=element_text(face="italic"),
				axis.title.y=element_text(margin=margin(t = 0, r = -5, b = 0, l = 0)))

## plot grid

plot_grid(plot_grid(
  plot_grid(mds.auto.plt + ggtitle("A") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
            sex.de.plt + ggtitle("D") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
            ncol=1),
  plot_grid(expression_div_hist + ggtitle("B") + theme(plot.margin = unit(c(0, 0.5, 0, 0), "cm")),
            enrich.DE.chr.plt + ggtitle("C") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
            align='h', rel_widths=c(0.4,0.6)),
            ncol=2, align="h", rel_widths=c(0.3,0.6)),
  enrich.sex.plt + ggtitle("E") + theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
  nrow=2, align="v", rel_heights=c(0.6,0.4))

## save image

ggsave("Fig1_expression_divergence_DE_sex.pdf", device="pdf", useDingbats=FALSE)
ggsave("Fig1_expression_divergence_DE_sex.png", device="png")
