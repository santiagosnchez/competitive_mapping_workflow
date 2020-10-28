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
library(ggpointdensity)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(viridis)
theme_set(theme_cowplot())
source("enrichment_fun.R")


cis_trans.molevol.data = read.csv("tables/df.cis_divergence.molevol.csv")
cis_trans.molevol.data$domain2 = as.character(cis_trans.molevol.data$domain)
cis_trans.molevol.data$domain2[ grepl("arm",cis_trans.molevol.data$domain2) ] = "arms"
cis_trans.molevol.data$domain2[ grepl("tip",cis_trans.molevol.data$domain2) ] = "tips"
type.levels = c("conserved","ambiguous","trans-only","cis-only","cis + trans (enhancing)","cis x trans (compensatory)","cis-trans (compensatory)")

arms_centers.counts.male = as.data.frame(cis_trans.molevol.data %>% filter(sex == "male" & domain2 != "tips") %>% group_by(domain2, type) %>% count() )
arms_centers.counts.female = as.data.frame(cis_trans.molevol.data %>% filter(sex == "female" & domain2 != "tips") %>% group_by(domain2, type) %>% count() )
arms_centers.type.enrich = data.frame(sex=rep(c("female","male"), each=14),
                                      type=rep(levels(arms_centers.counts.male$type),2),
                                      domain=rep(levels(arms_centers.counts.male$domain2),2),
                                      enrichment=c(enrichment(matrix(arms_centers.counts.female$n, nrow=2, ncol=7), odds.ratio=T)[1,],
                                                  enrichment(matrix(arms_centers.counts.male$n, nrow=2, ncol=7), odds.ratio=T)[1,]),
                                      p.value=c(enrichment(matrix(arms_centers.counts.female$n, nrow=2, ncol=7), odds.ratio=F)[1,],
                                                  enrichment(matrix(arms_centers.counts.male$n, nrow=2, ncol=7), odds.ratio=F)[1,]))
arms_centers.type.enrich$type = factor(arms_centers.type.enrich$type, type.levels)

arms_centers.type.enrich$sig = ""
arms_centers.type.enrich$sig[ arms_centers.type.enrich$p.value < 0.05 ] = "*"

# colors
brewer.pal(10, "Paired") -> cols
cols = cols[-c(7,8)]
cols.class = c("black","grey",cols[ c(4,6,2,8,7) ])
cols.type = c("black","grey",pal_jama()(4))
cols.type[3] = "#09BC8A"

ggplot(arms_centers.type.enrich, aes(x=type, y=log2(enrichment), fill=type)) +
  geom_bar(position=position_dodge(), stat="identity", alpha=0.7) +
  geom_text(aes(y=log2(enrichment)/2, label=sig), size=5, show.legend=F, color="white", nudge_x=-0.2) +
  scale_fill_manual(values=cols.type) +
  ylab(expression(paste("enrichment in arms ","[log"[2]," OR]",sep=""))) +
  xlab("") +
  facet_rep_wrap(~sex, ncol=1, strip.position="right") +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
  theme(strip.background = element_rect(fill="grey90", color=NA),
				legend.title=element_blank(),
				legend.text=element_text(size=10),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip()

ggplot(cis_trans.molevol.data, aes(y=type, x=prop_cons_500bp_upstream_5bp_win, fill=type)) +
  geom_density_ridges(color=NA, alpha=0.7, show.legend=F) +
  geom_boxploth(show.legend=F, outlier.shape=NA, width=0.5, apha=0.7) +
  facet_rep_wrap(~sex, ncol=1) +
  xlim(0,1) +
  scale_fill_manual(values=cols.type) +
  #scale_y_discrete(limits=c("center","arms")) +
  xlab(expression(paste(italic("K")[a], " "))) +
  ylab("") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(axis.text.y=element_blank()) +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75")

## plots

brewer.pal(10, "Paired") -> cols
cols = cols[-c(7,8)]
cols.arms_centers = cols[c(4,8)]

# exclude genes that do not have significant ASE difference

pA = ggplot(filter(cis_trans.molevol.data %>% filter(domain2 != "tips"), sex == "female"),
      aes(x=position/1e6, y=dn)) +
  geom_point(aes(color=domain2), alpha=0.2, show.legend=F) +
  geom_smooth(show.legend=F, color="black") +
  facet_rep_wrap(~chromosome, nrow=1, scales="free_x") +
  scale_color_manual(values=cols.arms_centers) +
  #scale_y_discrete(limits=c("center","arms")) +
  ylab(expression(paste(italic("K")[a]))) +
  xlab("genomic position (Mbp)") +
  #theme(strip.background = element_blank(), strip.text=element_blank()) +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
  ylim(0,0.2)

pA2 = ggplot(filter(cis_trans.molevol.data %>% filter(domain2 != "tips"), sex == "female"),
      aes(x=position/1e6, y=ds_ENCc)) +
  geom_point(aes(color=domain2), alpha=0.2, show.legend=F) +
  geom_smooth(show.legend=F, color="black") +
  facet_rep_wrap(~chromosome, nrow=1, scales="free_x") +
  scale_color_manual(values=cols.arms_centers) +
  #scale_y_discrete(limits=c("center","arms")) +
  ylab(expression(paste(italic("K")[s],"\'"))) +
  xlab("genomic position (Mbp)") +
  theme(strip.background = element_blank(), strip.text=element_blank()) +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
  ylim(0,0.75)

pB = ggplot(filter(cis_trans.molevol.data %>% filter(domain2 != "tips"), sex == "female"),
      aes(x=position/1e6, y=1-prop_cons_500bp_upstream_5bp_win)) +
  geom_point(aes(color=domain2), alpha=0.2, show.legend=F) +
  geom_smooth(show.legend=F, color="black") +
  facet_rep_wrap(~chromosome, nrow=1, scales="free_x") +
  scale_color_manual(values=cols.arms_centers) +
  #scale_y_discrete(limits=c("center","arms")) +
  ylab(expression(paste("1-",italic("p")[cons]))) +
  xlab("") +
  #theme(strip.background = element_blank(), strip.text=element_blank()) +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
  ylim(0,1)

pC = ggplot(cis_trans.molevol.data %>% filter(!grepl("hermaphrodite", sex)) %>%
                    filter(domain2 != "tips") %>%
                    filter(type != "ambiguous" & type != "conserved" & type != "trans only" ),
      aes(x=position/1e6, y=abs(logFC.ase))) +
  geom_point(aes(color=domain2), alpha=0.2, show.legend=F) +
  geom_smooth(show.legend=F, color="black") +
  facet_rep_grid(sex~chromosome, drop=TRUE, scales="free_x") +
  scale_color_manual(values=cols.arms_centers) +
  #scale_y_discrete(limits=c("center","arms")) +
  ylab(expression(paste("log"[2], " ASE (",italic("cis"),")"))) +
  xlab("genomic position (Mbp)") +
  theme(strip.background.x = element_blank(), strip.text.x=element_blank()) +
  theme(strip.background.y = element_rect(fill="grey90", color=NA)) +
  #theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  background_grid(major="y", size.major = 0.2, colour.major = "grey75") +
  ylim(0,4)

# remove the empty X-chr facet for males
grob = as_grob(pC)
grob$grobs[[ which(grob$layout$name %in% "panel-2-6") ]] = nullGrob()
grob$grobs[[ which(grob$layout$name %in% "axis-l-2-6") ]] = nullGrob()
grob$grobs[[ which(grob$layout$name %in% "axis-b-6-1") ]] = nullGrob()
grob$layout[ which(grob$layout$name %in% "axis-b-6"), c("t","b")] = grob$layout[ which(grob$layout$name %in% "axis-b-6"), c("t","b")] - 3
pC = as.ggplot(grob)

pD1 = ggplot(cis_trans.molevol.data %>% filter(domain2 != "tips") %>% filter(!grepl("hermaphrodite", sex)),
      aes(y=domain2, x=sqrt(abs(logFC.ase)), fill=domain2)) +
  geom_density_ridges(color=NA, alpha=0.7, show.legend=F) +
  geom_boxploth(show.legend=F, outlier.shape=NA, width=0.5, alpha=0.7) +
  facet_rep_wrap(~sex) +
  scale_fill_manual(values=cols.arms_centers) +
  #scale_y_discrete(limits=c("center","arms")) +
  xlab(expression(paste("log"[2], " ", italic("cis")," [sqrt norm]"))) +
  ylab("") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
  xlim(0,4)

pD2 = ggplot(cis_trans.molevol.data %>% filter(domain2 != "tips") %>% filter(!grepl("hermaphrodite", sex)),
      aes(y=domain2, x=sqrt(abs(logFC.sp)-abs(logFC.ase)), fill=domain2)) +
  geom_density_ridges(color=NA, alpha=0.7, show.legend=F) +
  geom_boxploth(show.legend=F, outlier.shape=NA, width=0.5, alpha=0.7) +
  facet_rep_wrap(~sex) +
  scale_fill_manual(values=cols.arms_centers) +
  #scale_y_discrete(limits=c("center","arms")) +
  xlab(expression(paste("log"[2], " ", italic("trans"), " [sqrt norm]"))) +
  ylab("") +
  theme(strip.background = element_rect(fill="grey90", color=NA)) +
  theme(axis.text.y=element_blank()) +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75") +
  xlim(0,4)

part1 = plot_grid(
          pB + theme(plot.margin = unit(c(t=0.2, r=0.64, b=0, l=0), "cm")),
          pC + theme(plot.margin = unit(c(t=-0.1, r=-0.2, b=0, l=0.29), "cm")),
          plot_grid(pD1 + theme(plot.margin = unit(c(t=0, r=0.4, b=0, l=0), "cm")),
                    pD2 + theme(plot.margin = unit(c(t=0, r=1, b=0, l=-0.6), "cm"))),
          rel_heights = c(0.3,0.6,0.4), ncol=1, labels=c("A","B","C"))

# correlation
p1 = ggplot(filter(cis_trans.molevol.data, sex == "female"), aes(y=dn/ds_ENCc, x=abs(logFC.ase))) +
    geom_pointdensity(show.legend=F) +
    geom_smooth(method="lm", color="black") +
    ylab(expression(paste(italic("K")[a],"/",italic("K")[s],"\'",sep=""))) +
    xlab("") +
    background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
    facet_rep_wrap(~sex, ncol=1, strip.position="right") +
    scale_color_viridis(option="viridis") +
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    ylim(0,1)
p2 = ggplot(filter(cis_trans.molevol.data, sex == "female"), aes(y=dn, x=abs(logFC.ase))) +
    geom_pointdensity(show.legend=F) +
    geom_smooth(method="lm", color="black") +
    ylab(expression(paste(italic("K")[a]))) +
    xlab(expression(paste("log"[2], " ",italic("cis")))) +
    background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
    scale_color_viridis(option="viridis") +
    facet_rep_wrap(~sex, ncol=1, strip.position="right") +
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    ylim(0,0.4)
p3 = ggplot(filter(cis_trans.molevol.data, sex == "female"), aes(y=1-prop_cons_500bp_upstream_5bp_win, x=abs(logFC.ase))) +
    geom_pointdensity(show.legend=F) +
    geom_smooth(method="lm", color="black") +
    ylab(expression(paste("1-",italic("p")[cons]))) +
    xlab(expression(paste("log"[2], " ",italic("cis")))) +
    background_grid(major="xy", size.major = 0.2, colour.major = "grey75") +
    facet_rep_wrap(~sex, ncol=1, strip.position="right") +
    theme(strip.background = element_blank(), strip.text = element_blank()) +
    theme(axis.title.x=element_text(hjust=-0.35)) +
    #theme(strip.background = element_rect(fill="grey90", color=NA)) +
    scale_color_viridis(option="viridis")
part2 = plot_grid(p1, p3, nrow=1, labels = c("A","B"), rel_widths=c(0.45,0.5))
part2
#plot_grid(part1, part2, rel_heights=c(c(0.65, 0.35)), ncol=1)

ggsave("figures/Fig7_cis_molevol_corr.png")
ggsave("figures/Fig7_cis_molevol_corr.pdf", device="pdf", useDingbats=F)

# supplementary fig

tmp = as.data.frame( cis_trans.molevol.data %>% filter(sex == "male" & chromosome == "X"))
tmp2 = as.data.frame( tmp %>% filter(type == "cis only" | type == "trans only"))
tmp$type = factor(tmp$type, type.levels)
tmp2$type = factor(tmp2$type, type.levels)


p1 = ggplot(tmp, aes(x=dn/ds_ENCc, y = type, fill=type)) +
  geom_density_ridges(alpha=0.5, show.legend=F, color=NA) +
  geom_boxploth(alpha=0.7, show.legend=F) +
  scale_fill_manual(values=cols.type) +
  scale_x_log10(limits=c(0.001, 1)) +
  scale_y_discrete(expand=c(0.2,0)) +
  annotation_logticks(side="b") +
  ylab("") +
  xlab(expression(paste(italic("K")[a],"/",italic("K")[s],"\'"," [log"[10],"]", sep=""))) +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75")

p2 = ggplot(tmp, aes(x=1-prop_cons_500bp_upstream_5bp_win, y = type, fill=type)) +
  geom_density_ridges(alpha=0.5, show.legend=F, color=NA) +
  geom_boxploth(alpha=0.7, show.legend=F) +
  scale_fill_manual(values=cols.type) +
  ylab("") +
  xlab(expression(paste("1-",italic("p")[cons]))) +
  background_grid(major="x", size.major = 0.2, colour.major = "grey75")

plot_grid(p1, p2, ncol=1, labels=c("A","B"))

ggsave("suppl_cis_trans_X_male_dnds_prop_more.png")
ggsave("suppl_cis_trans_X_male_dnds_prop_more.pdf", device="pdf")

# supplementary

plot_grid(pA,# + theme(plot.margin = unit(c(t=0, r=0, b=1, l=0), "cm")),
          pA2,# + theme(plot.margin = unit(c(t=-1, r=0, b=0, l=0), "cm")),
          ncol=1, align="v")

ggsave("suppl_ka_ks_along_chr.png")
