library(qiime2R)
require(ggplot2)
library(btools)
library(vegan)
library(phyloseq)
library(dplyr)
library(gplots)
library(QsRutils)
library(data.table)
library(SRS)
library(reshape2)
library(matrixStats)
library(ggbeeswarm)
install.packages('plotmath')
librar
library(cmprsk)
library(ggstance)
library(ggpubr)
library(devtools)
library(BiodiversityR)
theme_set(theme_bw())
site_colors<-c('#5B4BB7', '#B79A4C', '#B75B4B')

#loading physeq object w/ full dataset #490 samples (bad samples (<3k) already discarded)
physeq<-qza_to_phyloseq(
  features="Files/ASV_table_final_goodsamples.qza",
  tree="Files/rooted-tree.qza",
  "Files/taxonomy_vs.qza",
  metadata = "Files/mapping_tmo_0.tsv")

#first lets work with tmo samples only (phase-wise)
physeq_tmo<-subset_samples(physeq,time!="MI"&time!="MF")

# primeiro, normalizaremos por SRS {3578}
physeq_SRS_tmo<-physeq_tmo
SRS_output <- SRS(data = as.data.frame(otu_table(physeq_SRS_tmo)), Cmin = min(sample_sums(physeq_tmo)), set_seed = 123) #running SRS
rownames(SRS_output)<-rownames(as.data.frame(otu_table(physeq_SRS_tmo))) #reassigning features names as rownames
otu_table(physeq_SRS_tmo)<-otu_table(SRS_output,taxa_are_rows = T)


#calculating important diversity metrics #Simpson only
div_SRS_tmo<-cbind(as.data.frame(sample_data(physeq_SRS_tmo)),
                   estimate_richness(physeq_SRS_tmo,measures=c("Simpson")))
                   
div_SRS_tmo$time<-factor(div_SRS_tmo$time, levels=c("P","A","E","E30","E75")) #fixar ordem de coletas
div_SRS_tmo[div_SRS_tmo=="BF"]<-"SB"
div_SRS_tmo[div_SRS_tmo=="FCG"]<-"GCF"
div_SRS_tmo[div_SRS_tmo=="MO"]<-"OM"


##################PLOTS#######################
#div x phase (site-separated)
F2A<-ggplot(data = div_SRS_tmo) + 
  aes(x = time, y = Simpson, color = site, group = time)+
  geom_boxplot(alpha=0.5, outlier.size = 1) + facet_grid(~ site,scales = "free_y") + 
  labs(y = "Diversity", x = 'Timepoint') + guides(color=FALSE)+
  scale_color_manual(values=site_colors)+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits = c(0,1.05), expand = c(0,0))+
  stat_compare_means(ref.group = "P",method = "wilcox.test", size = 3.2, label.y.npc = 0.95, 
                     aes(label = ifelse(..p.format..<0.05,..p.signif..,ifelse(..p.format..<0.1,
                                                                              round(as.numeric(..p.format..),2),''))))+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))

#regressão da diversidade por sítio e % dos pacientes enxertados por dia em relação à infusão
tte<-as.data.frame(read.csv("Files/DadosClinicosTMO_09.05.2021.csv"))[c("Patient.acronym","Time.to.engraftment")]
colnames(tte)<-c("patient","Time.to.engraftment")
cuminc_tte<-cuminc(tte$Time.to.engraftment,rep(1,each=31))
cuminc_tte<-data.frame("time"=as.vector(cuminc_tte[["1 1"]][["time"]]),"est"=as.vector(cuminc_tte[["1 1"]][["est"]]))
cuminc_tte<-rbind(cuminc_tte,c(100,1),c(-30,0))

div_SRS_tmo_tte<-merge(div_SRS_tmo,tte)
div_SRS_tmo_tte$days.to.I<-div_SRS_tmo_tte$Time.to.engraftment+div_SRS_tmo_tte$days.to.E
div_SRS_tmo_tte

FS6B<-ggplot(data = div_SRS_tmo_tte)+
  geom_smooth(aes(x=days.to.I,y=Simpson,color=site,group=site, fill = site))+
  scale_color_manual(values=site_colors)+
  scale_fill_manual(values=site_colors)+
  geom_line(data = cuminc_tte, aes(x=time,y=est*0.4+0.6), color = "black")+
  coord_cartesian(xlim=c(-15,95),ylim=c(0.6,1), expand = T)+
  scale_y_continuous(name = "Diversity",sec.axis = sec_axis(trans=~.*1/0.4-0.6/0.4,name="% of engrafted patients"))+
  labs(x = "Days from transplant", color = "Site")+guides(fill=F)+theme(panel.grid = element_blank())

############################################################################
#now we will compute all relevant ecological metrics
#calculations as in Orwin & Wardle (2004) + modifications
div_SRS.eco<-reshape(div_SRS_tmo, idvar = c("patient", "site"), timevar = "time", direction = "wide")

#lets compute resistance (only possible for patients with E & P samples)
div_SRS.eco$Simpson.rs<-1-2*abs(div_SRS.eco$Simpson.P-div_SRS.eco$Simpson.E)/
                            (div_SRS.eco$Simpson.P+abs(div_SRS.eco$Simpson.P-div_SRS.eco$Simpson.E))

#lets compute resilience considering E30
div_SRS.eco$Simpson.rl<-NA
for (i in 1:nrow(div_SRS.eco)){
    if (is.na(div_SRS.eco[["Simpson.E30"]][i])==F){
      div_SRS.eco[["Simpson.rl"]][i]<-(((2*abs(div_SRS.eco[["Simpson.P"]][i]-div_SRS.eco[["Simpson.E"]][i])/
                                            (abs(div_SRS.eco[["Simpson.P"]][i]-div_SRS.eco[["Simpson.E"]][i])+abs(div_SRS.eco[["Simpson.P"]][i]-div_SRS.eco[["Simpson.E30"]][i])))-1)/
                                          (div_SRS.eco[["days.to.E.E30"]][i]))
    }
}
#lets compute stability considering the last available sample (only possible for patients with E & P & (E30|E75) samples)
div_SRS.eco$Simpson.s<-NA
for (i in 1:nrow(div_SRS.eco)){
    if (is.na(div_SRS.eco[["Simpson.E30"]][i])==F){
      div_SRS.eco[["Simpson.s"]][i]<-1-(4*(abs(div_SRS.eco[["Simpson.P"]][i]-div_SRS.eco[["Simpson.E30"]][i])*abs(div_SRS.eco[["Simpson.P"]][i]-div_SRS.eco[["Simpson.E"]][i])))/
        ((div_SRS.eco[["Simpson.P"]][i]+abs(div_SRS.eco[["Simpson.P"]][i]-div_SRS.eco[["Simpson.E"]][i]))*(abs(div_SRS.eco[["Simpson.E"]][i]-div_SRS.eco[["Simpson.E30"]][i])+abs(div_SRS.eco[["Simpson.P"]][i]-div_SRS.eco[["Simpson.E"]][i])))
    }
}

#div at P x ecometrics (per site)
FS11A<-ggarrange(
ggscatter(data = div_SRS.eco, y = "Simpson.P", x = "Simpson.rs", facet.by = "site", ylab = "Diversity (P)",
          xlab = "Resistance", cor.coef = T, cor.method = "spearman", cor.coef.coord = c(0.05,0.51),
          add = 'reg.line', conf.int = F, alpha = 0.5, ggtheme = theme_bw(), cor.coef.size = 3.2)+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(color = 'black'),
        strip.background = element_rect(color="white", fill="white", linetype="solid")),
ggscatter(data = div_SRS.eco, y = "Simpson.P", x = "Simpson.rl", facet.by = "site", ylab = "Diversity (P)",
          xlab = "Resilience", cor.coef = T, cor.method = "spearman", cor.coef.coord = c(-0.04,0.51),
          add = 'reg.line', conf.int = F, alpha = 0.5, ggtheme = theme_bw(), cor.coef.size = 3.2)+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(color = 'black'),
        strip.background = element_rect(color="white", fill="white", linetype="solid")),
ggscatter(data = div_SRS.eco, y = "Simpson.P", x = "Simpson.s", facet.by = "site", ylab = "Diversity (P)",
          xlab = "Stability", cor.coef = T, cor.method = "spearman", cor.coef.coord = c(-0.9,0.51),
          add = 'reg.line', conf.int = F, alpha = 0.5, ggtheme = theme_bw(), cor.coef.size = 3.2)+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(color = 'black'),
        strip.background = element_rect(color="white", fill="white", linetype="solid")),
nrow = 3)


#now lets compare eco measures across sites
div_SRS.eco.melt<-melt(div_SRS.eco,id.vars=colnames(div_SRS.eco)[1:12],
                       variable.name="div_eco_metric",value.name="div_eco_measure")
div_SRS.eco.melt$div_eco_metric<-as.character(div_SRS.eco.melt$div_eco_metric)
div_SRS.eco.melt[div_SRS.eco.melt=="Simpson.rs"]<-"Resistance~(P%->%E)"
div_SRS.eco.melt[div_SRS.eco.melt=="Simpson.rl"]<-"Resilience~(E%->%E30)"
div_SRS.eco.melt[div_SRS.eco.melt=="Simpson.s"]<-"Stability~(P%->%E30)"
div_SRS.eco.melt$div_eco_metric<-factor(div_SRS.eco.melt$div_eco_metric, 
                                        levels = c("Resistance~(P%->%E)","Resilience~(E%->%E30)","Stability~(P%->%E30)"))

F2B<-ggplot(data = div_SRS.eco.melt) + 
  aes(x = site, y = div_eco_measure, color = site, group = site) + 
  geom_boxplot(alpha=0.5, outlier.size = 1)+
  #geom_quasirandom(size=1)+
  facet_wrap(. ~ div_eco_metric,scales = "free_y", ncol=3, labeller = label_parsed) + 
  labs(y = "Measure", x = 'Site') + guides(color=FALSE)+
  scale_y_continuous(expand = expand_scale(mult = c(0.04, 0.08)))+
  scale_color_manual(values=site_colors)+
  stat_compare_means(comparisons = list(c("OM","GCF"),c("OM","SB"),c("SB","GCF")),size = 3.2,
                     label="p.format",method = "wilcox.test",paired = F,step.increase = 0.09, tip.length = 0.005)+
theme(panel.grid = element_blank(),
      strip.background = element_rect(color="white", fill="white", linetype="solid"),
      strip.text=element_text(size=8),
      panel.border = element_blank(),axis.line = element_line(color = 'black'))
F2B

#RSxRL
F6A<-ggplot(div_SRS.eco, aes(y=Simpson.rs,x=Simpson.rl,fill=Simpson.s))+
  geom_point(size=3, color = "black",alpha=0.8, shape = 21)+facet_wrap(.~site)+
  scale_fill_gradient2(low="white",mid="yellow",high="red", midpoint = 0.07498193, breaks = c(-0.5,0,0.5))+
  scale_shape_manual(values=c(21,23,24))+
  scale_y_continuous(limits = c(0,1),breaks = c(0,0.2,0.4,0.6,0.8,1))+
  scale_x_continuous(limits=c(-0.045,0.045), breaks = c(-0.04,-0.02,0,0.02,0.04))+
labs(shape="Site",fill="  Stability\n(P -> E30)",x=expr(paste("Resilience (E"%->%"E30)")),y=expr(paste("Resistance (P"%->%"E)")))+
  theme(panel.grid = element_blank(),strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))

F6A

#saving alpha metrics to an external file
merge(dcast(div_SRS.eco,patient~site, value.var = c("Simpson.rs")) %>% `colnames<-`(c("patient","DB.rs", "GCF.rs", "OM.rs")),
dcast(div_SRS.eco,patient~site, value.var = c("Simpson.rl")) %>% `colnames<-`(c("patient","DB.rl", "GCF.rl", "OM.rl")),
dcast(div_SRS.eco,patient~site, value.var = c("Simpson.rl")) %>% `colnames<-`(c("patient","DB.rl", "GCF.rl", "OM.rl")))

write.csv(as.data.frame(data.table::dcast(div_SRS.eco,patient~site,value.var = c("Simpson.rs","Simpson.rl","Simpson.s","Simpson.P",
                                                         "Simpson.A","Simpson.E","Simpson.E30","Simpson.E75"))),
          file = "Files/300522_microbiota-variables.csv", row.names = F)

#math expressions
FS6A<-ggplot(data=data.frame(x=0,y=0))+geom_point(aes(x=x,y=y),alpha=0)+ylim(-3,3)+
  annotate('text', x = 0, y = 2,label = "italic(RS) == 1 - frac(2 * abs(alpha[p] - alpha[e]) , alpha[p] + abs(alpha[p] - alpha[e]))",parse = T, size=5)+
  annotate('text',x = 0, y = 3, label = "Resistance~(P%->%E)", size = 4, color = "darkred", parse = T)+
  annotate('text', x = 0, y = 0,label = "italic(RL) == (frac(2 * abs(alpha[p] - alpha[e]) ,abs(alpha[p] - alpha[e]) +  abs(alpha[p] - alpha[e30])) - 1)%/%t",parse = T, size=5)+
  annotate('text',x = 0, y = 1, label = "Resilience~(E%->%E30)", size = 4, color = "darkred", parse = T)+
  annotate('text', x = 0, y = -2,label = "italic(S) == 1 - frac(4 * abs(alpha[p] - alpha[e]) * abs(alpha[p] - alpha[e30]), (alpha[p] + abs(alpha[p] - alpha[e])) * (abs(alpha[p] - alpha[e]) + abs(alpha[p] - alpha[e30])))",parse = T, size=5)+
  annotate('text',x = 0, y = -1, label = "Stability~(P%->%E30)", size = 4, color = "darkred", parse = T)+
  theme(panel.grid = element_blank(), panel.border = element_blank(),axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank())

FS6A
