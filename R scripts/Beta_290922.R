library(qiime2R)
require(ggplot2)
library(btools)
library(vegan)
library(phyloseq)
library(geometry)
library(gplots)
library(cowplot)
library(QsRutils)
library(tidyverse)
library(dplyr)
library(ggpubr)
library(reshape2)
library(ggbeeswarm)
library(ggstance)
library(ggplotify)
library(ggpubr)
library(devtools)
library(pairwiseAdonis)
library(mctoolsr)
source("plotDistances.R")
library(SRS)
options(scipen=999)
theme_set(theme_bw())
site_colors<-c('#5B4BB7', '#B79A4C', '#B75B4B')
pacientes<-c("AAM"="#1","AFLF"="#2","AMSCP"="#3","BMF"="#4","CEFC"="#5","DBHB"="#6",
             "ECPL"="#7","ENS"="#8","FGTC"="#9","FM"="#10","IABC"="#11","ILSM"="#12",
             "JPM"="#13","JVS"="#14","LEMG"="#15","LFCS"="#16","MEFD"="#17","MFML"="#18",
             "MMMC"="#19","MSBS"="#20","OBR"="#21","OSLG"="#22","OSN"="#23","PHOS"="#24",
             "RJT"="#25","RLS"="#26","SMG"="#27","TPN"="#28","VFF"="#29","VGS"="#30","VRM"="#31")

#loading physeq object w/ full dataset #490 samples (bad samples (<3k) already discarded)
physeq<-qza_to_phyloseq(
  features="Files/ASV_table_final_goodsamples.qza",
  tree="Files/rooted-tree.qza",
  "Files/taxonomy_vs.qza",
  metadata ="Files/mapping_tmo_0.tsv")

metadata<-as.data.frame(sample_data(physeq))
metadata$time<-factor(metadata$time, levels = c('P','A','E','E30','E75'))
sample_data(physeq)<-metadata
#first lets work with tmo samples only (phase-wise)
physeq_tmo<-subset_samples(physeq,time!="MI"&time!="MF")

# primeiro, normalizaremos por SRS {3578}
physeq_SRS_tmo<-physeq_tmo
#set.seed(123)
SRS_output <- SRS(data = as.data.frame(otu_table(physeq_SRS_tmo)), Cmin = min(sample_sums(physeq_tmo)), set_seed=123) #running SRS
rownames(SRS_output)<-rownames(as.data.frame(otu_table(physeq_SRS_tmo))) #reassigning features names as rownames
otu_table(physeq_SRS_tmo)<-otu_table(SRS_output,taxa_are_rows = T)
#separate per site
physeq_SRS_tmo_mo<-subset_samples(physeq_SRS_tmo,site=="MO")
physeq_SRS_tmo_bf<-subset_samples(physeq_SRS_tmo,site=="BF")
physeq_SRS_tmo_fcg<-subset_samples(physeq_SRS_tmo,site=="FCG")

ord_mo_phase <- rbiom::unifrac(otu_table(physeq_SRS_tmo_mo),weighted=T, tree = phy_tree(physeq_SRS_tmo_mo))
ord_bf_phase <- rbiom::unifrac(otu_table(physeq_SRS_tmo_bf),weighted=T, tree = phy_tree(physeq_SRS_tmo_bf))
ord_fcg_phase <- rbiom::unifrac(otu_table(physeq_SRS_tmo_fcg),weighted=T, tree = phy_tree(physeq_SRS_tmo_fcg))

#plot distances to P centroid
dist_mo_pcentroid_wu<-data.frame("dist"=betadisper_P(ord_mo_phase, 
                                 as(sample_data(physeq_SRS_tmo_mo), "data.frame")$time, type="centroid")[["distances"]])
dist_mo_pcentroid_wu$sampleid<-rownames(dist_mo_pcentroid_wu)
metadata<-as(sample_data(physeq_SRS_tmo_mo),"data.frame")
metadata$sampleid<-rownames(metadata)
dist_mo_pcentroid_wu<-merge(dist_mo_pcentroid_wu,metadata,by="sampleid")

dist_bf_pcentroid_wu<-data.frame("dist"=betadisper_P(ord_bf_phase, 
                                                     as(sample_data(physeq_SRS_tmo_bf), "data.frame")$time, type="centroid")[["distances"]])
dist_bf_pcentroid_wu$sampleid<-rownames(dist_bf_pcentroid_wu)
metadata<-as(sample_data(physeq_SRS_tmo_bf),"data.frame")
metadata$sampleid<-rownames(metadata)
dist_bf_pcentroid_wu<-merge(dist_bf_pcentroid_wu,metadata,by="sampleid")

dist_fcg_pcentroid_wu<-data.frame("dist"=betadisper_P(ord_fcg_phase, 
                                                     as(sample_data(physeq_SRS_tmo_fcg), "data.frame")$time, type="centroid")[["distances"]])
dist_fcg_pcentroid_wu$sampleid<-rownames(dist_fcg_pcentroid_wu)
metadata<-as(sample_data(physeq_SRS_tmo_fcg),"data.frame")
metadata$sampleid<-rownames(metadata)
dist_fcg_pcentroid_wu<-merge(dist_fcg_pcentroid_wu,metadata,by="sampleid")

#bind all sites
dist_pcentroid_wu<-rbind(dist_mo_pcentroid_wu,dist_bf_pcentroid_wu,dist_fcg_pcentroid_wu)
#fix factor levels for phase
dist_pcentroid_wu$time<-factor(dist_pcentroid_wu$time, levels=c("P","A","E","E30","E75"))

dist_pcentroid_wu[dist_pcentroid_wu=="BF"]<-"SB"
dist_pcentroid_wu[dist_pcentroid_wu=="FCG"]<-"GCF"
dist_pcentroid_wu[dist_pcentroid_wu=="MO"]<-"OM"
#plot distances to P centroid
F2C<-ggplot(data=subset(dist_pcentroid_wu),aes(x = time, y = dist, color = site, group = time)) + 
  geom_boxplot(alpha=0.5, outlier.size = 1) + facet_wrap(.~site)+
  #geom_quasirandom(alpha=0.5, size=1)+
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1), expand = c(0,0))+
  labs(y = "Distance to P centroid", x = 'Timepoint') + guides(color=FALSE)+
  scale_color_manual(values=site_colors)+
  stat_compare_means(ref.group = "P",method = "wilcox.test", size = 3.2, label.y.npc = 0.95, 
                     aes(label = ifelse(..p.format..<0.05,..p.signif..,ifelse(..p.format..<0.1,
                                                                             round(as.numeric(..p.format..),2),''))))+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))

F2C
#plot distances from centroid
dist_mo_centroid_wu<-data.frame("dist"=betadisper(ord_mo_phase, 
                                                     as(sample_data(physeq_SRS_tmo_mo), "data.frame")$time, type="centroid")[["distances"]])
dist_mo_centroid_wu$sampleid<-rownames(dist_mo_centroid_wu)
metadata<-as(sample_data(physeq_SRS_tmo_mo),"data.frame")
metadata$sampleid<-rownames(metadata)
dist_mo_centroid_wu<-merge(dist_mo_centroid_wu,metadata,by="sampleid")

dist_bf_centroid_wu<-data.frame("dist"=betadisper(ord_bf_phase, 
                                                     as(sample_data(physeq_SRS_tmo_bf), "data.frame")$time, type="centroid")[["distances"]])
dist_bf_centroid_wu$sampleid<-rownames(dist_bf_centroid_wu)
metadata<-as(sample_data(physeq_SRS_tmo_bf),"data.frame")
metadata$sampleid<-rownames(metadata)
dist_bf_centroid_wu<-merge(dist_bf_centroid_wu,metadata,by="sampleid")

dist_fcg_centroid_wu<-data.frame("dist"=betadisper(ord_fcg_phase, 
                                                      as(sample_data(physeq_SRS_tmo_fcg), "data.frame")$time, type="centroid")[["distances"]])
dist_fcg_centroid_wu$sampleid<-rownames(dist_fcg_centroid_wu)
metadata<-as(sample_data(physeq_SRS_tmo_fcg),"data.frame")
metadata$sampleid<-rownames(metadata)
dist_fcg_centroid_wu<-merge(dist_fcg_centroid_wu,metadata,by="sampleid")

#bind all sites
dist_centroid_wu<-rbind(dist_mo_centroid_wu,dist_bf_centroid_wu,dist_fcg_centroid_wu)
#fix factor levels for phase
dist_centroid_wu$time<-factor(dist_centroid_wu$time, levels=c("P","A","E","E30","E75"))

dist_centroid_wu[dist_centroid_wu=="BF"]<-"SB"
dist_centroid_wu[dist_centroid_wu=="FCG"]<-"GCF"
dist_centroid_wu[dist_centroid_wu=="MO"]<-"OM"

#plot distances from centroid
F6C<-ggplot(data=subset(dist_centroid_wu),aes(x = time, y = dist, color = site, group = time)) + 
  geom_boxplot(alpha=0.5, outlier.size = 1) + facet_wrap(.~site)+
  #geom_quasirandom(alpha=0.5, size=1)+
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.2,0.4,0.6,0.8,1), expand = c(0,0))+
  labs(y = "Distance to centroid", x = 'Timepoint') + guides(color=FALSE)+
  scale_color_manual(values=site_colors)+
  stat_compare_means(ref.group = "P",method = "wilcox.test", size = 3.2, label.y.npc = 0.95, 
                    aes(label = ifelse(..p.format..<0.05,..p.signif..,ifelse(..p.format..<0.1,
                                                          round(as.numeric(..p.format..),2),''))))+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))
F6C
#plot distances withing group per phase (ie, literal beta-diversity)
metadata<-as.data.frame(sample_data(physeq_SRS_tmo_mo))
metadata$sampleid<-rownames(as.data.frame(sample_data(physeq_SRS_tmo_mo)))
sample_data(physeq_SRS_tmo_mo)<-metadata
dist_mo_phase_wu<-plotDistances(physeq_SRS_tmo_mo,m="wunifrac",plot=F,s="sampleid",d="time")
dist_mo_phase_wu_within_phase<-subset(dist_mo_phase_wu,Type1==Type2)
dist_mo_phase_wu_within_phase<-subset(dist_mo_phase_wu_within_phase, !duplicated(subset(dist_mo_phase_wu_within_phase, 
                                                                                        select=c(value))))
dist_mo_phase_wu_within_phase$site<-"OM"

metadata<-as.data.frame(sample_data(physeq_SRS_tmo_bf))
metadata$sampleid<-rownames(as.data.frame(sample_data(physeq_SRS_tmo_bf)))
sample_data(physeq_SRS_tmo_bf)<-metadata
dist_bf_phase_wu<-plotDistances(physeq_SRS_tmo_bf,m="wunifrac",plot=F,s="sampleid",d="time")
dist_bf_phase_wu_within_phase<-subset(dist_bf_phase_wu,Type1==Type2)
dist_bf_phase_wu_within_phase<-subset(dist_bf_phase_wu_within_phase, !duplicated(subset(dist_bf_phase_wu_within_phase, 
                                                                                        select=c(value))))
dist_bf_phase_wu_within_phase$site<-"SB"

metadata<-as.data.frame(sample_data(physeq_SRS_tmo_fcg))
metadata$sampleid<-rownames(as.data.frame(sample_data(physeq_SRS_tmo_fcg)))
sample_data(physeq_SRS_tmo_fcg)<-metadata
dist_fcg_phase_wu<-plotDistances(physeq_SRS_tmo_fcg,m="wunifrac",plot=F,s="sampleid",d="time")
dist_fcg_phase_wu_within_phase<-subset(dist_fcg_phase_wu,Type1==Type2)
dist_fcg_phase_wu_within_phase<-subset(dist_fcg_phase_wu_within_phase, !duplicated(subset(dist_fcg_phase_wu_within_phase, 
                                                                                        select=c(value))))
dist_fcg_phase_wu_within_phase$site<-"GCF"

#bind all sites
dist_wu_within_phase<-rbind(dist_mo_phase_wu_within_phase,dist_bf_phase_wu_within_phase,dist_fcg_phase_wu_within_phase)
#fix factor levels for phase
dist_wu_within_phase$Type1<-factor(dist_wu_within_phase$Type1, levels=c("P","A","E","E30","E75"))

F6D<-ggplot(data=dist_wu_within_phase,aes(x = Type1, y = value, group = Type1)) +
  geom_violin(aes(fill=site),size=0.25, color = 'white')+
  geom_boxplot(color="black",width = 0.2, outlier.shape = NA) + 
  facet_wrap(.~site)+
  scale_y_continuous(limits = c(0,1.3), breaks = c(0,0.2,0.4,0.6,0.8,1,1.2), expand = c(0,0))+
  labs(y = "Distance within group", x = 'Timepoint') + guides(fill=FALSE)+
  scale_fill_manual(values=site_colors)+
  stat_compare_means(ref.group = "P",method = "wilcox.test", size = 3.2, label.y.npc = 0.95, 
                     aes(label = ifelse(..p.format..<0.05,..p.signif..,ifelse(..p.format..<0.1,
                                                                              round(as.numeric(..p.format..),2),''))))+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))

F6D
#boxplot com distância média intrapaciente entre sítios por fase
metadata<-as(sample_data(physeq_SRS_tmo),"data.frame")
metadata$sampleid<-rownames(metadata)
sample_data(physeq_SRS_tmo)<-metadata
dist_wu<-plotDistances(physeq_SRS_tmo,m="wunifrac",plot=F,s="sampleid",d="time")
colnames(dist_wu)[1]<-"sampleid"
dist_wu<-merge(metadata,dist_wu)
colnames(dist_wu)<-c("sampleid1","patient1","site1","time1","days.to.E1","sampleid","value","Type1","Type2")
dist_wu<-merge(metadata,dist_wu)
colnames(dist_wu)<-c("sampleid2","patient2","site2","time2","days.to.E2","sampleid1","patient1","site1","time1","days.to.E1","value","Type1","Type2")
dist_wu<-dist_wu[c("patient2","site2","time2","patient1","site1","time1","value")]
dist_wu_subset<-subset(dist_wu,patient1==patient2&time1==time2)
#subset(data.frame(table(dist_wu_subset$value)), Freq == 2)
dist_wu_subset<-subset(dist_wu_subset, !duplicated(subset(dist_wu_subset, select=c(value))))
#dist_wu_subset<-rbind(dist_wu_subset,c('AFLF','BF','E30','AFLF','FCG','E30',5.347207e-05))
dist_wu_subset$value<-as.numeric(dist_wu_subset$value)

dist_wu_subset_min<-data.frame(patient = character(0), time = character(0), min_dist = numeric(0))
for (pt in unique(dist_wu_subset$patient1)){
  for (tm in unique(dist_wu_subset$time1)){
    if(nrow(subset(dist_wu_subset,patient1==pt&time1==tm))==3){
    dist_wu_subset_min<-rbind(dist_wu_subset_min,c(pt,tm,min(subset(dist_wu_subset,patient1==pt&time1==tm)$value)))
    }
  }
}

colnames(dist_wu_subset_min)<-c("patient","time","min_dist")
dist_wu_subset_min<-subset(dist_wu_subset_min,min_dist!='Inf')
dist_wu_subset_min$min_dist<-as.numeric(dist_wu_subset_min$min_dist)
dist_wu_subset_min$time<-factor(dist_wu_subset_min$time,levels = c('P','A','E','E30','E75'))
dist_wu_subset$time1<-factor(dist_wu_subset$time1,levels = c('P','A','E','E30','E75'))

F1C<-ggplot(data = dist_wu_subset_min) + 
  aes(x = time, y = min_dist, group = time) + 
  geom_boxplot(alpha = 0.5, color = 'black', outlier.size = 1) +
  #geom_quasirandom(alpha=0.5, size=1)+
  labs(y = "Min. dist. between sites", x = 'Timepoint') + guides(color=FALSE)+
  theme(panel.grid.minor = element_blank())+
  scale_y_continuous(expand = c(0,0),limits=c(0,0.6),breaks=c(0,0.2,0.4,0.6))+
  stat_compare_means(ref.group = "P",method = "wilcox.test", size = 3.2, label.y.npc = 0.95, 
                     aes(label = ifelse(..p.format..<0.05,..p.signif..,ifelse(..p.format..<0.1,
                                                                              round(as.numeric(..p.format..),2),''))))+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))

F1C

###13/12/21 - boxplot da distância intrapaciente entre coletas subsequentes (facetada por sítio)
metadata<-as(sample_data(physeq_SRS_tmo_mo),"data.frame")
metadata$sampleid<-rownames(metadata)
sample_data(physeq_SRS_tmo_mo)<-metadata
dist_mo_phase_wu<-plotDistances(physeq_SRS_tmo_mo,m="wunifrac",plot=F,s="sampleid",d="time")

#pick patients names and remove non-intrapatient comparisons
dist_mo_phase_ind<-merge(metadata[c("patient","sampleid")],dist_mo_phase_wu,by.x = "sampleid",by.y = "Var1")
colnames(dist_mo_phase_ind)<-c("sampleid","patient1","Var2","value","Type1","Type2")
dist_mo_phase_ind<-merge(metadata[c("patient","sampleid")],dist_mo_phase_ind,by.x = "sampleid",by.y = "Var2")
colnames(dist_mo_phase_ind)<-c("sampleid0","patient0","sampleid1","patient1","value","Type1","Type2")

dist_mo_phase_ind_intra<-subset(dist_mo_phase_ind,patient0==patient1)
dist_mo_phase_ind_intra_sq<-subset(dist_mo_phase_ind_intra,(Type2=="P"&Type1=="A")|(Type2=="A"&Type1=="E")|(Type2=="E"&Type1=="E30")|(Type2=="E30"&Type1=="E75"))
dist_mo_phase_ind_intra_sq$site<-"MO"

metadata<-as(sample_data(physeq_SRS_tmo_bf),"data.frame")
metadata$sampleid<-rownames(metadata)
sample_data(physeq_SRS_tmo_bf)<-metadata
dist_bf_phase_wu<-plotDistances(physeq_SRS_tmo_bf,m="wunifrac",plot=F,s="sampleid",d="time")

#pick patients names and remove non-intrapatient comparisons
dist_bf_phase_ind<-merge(metadata[c("patient","sampleid")],dist_bf_phase_wu,by.x = "sampleid",by.y = "Var1")
colnames(dist_bf_phase_ind)<-c("sampleid","patient1","Var2","value","Type1","Type2")
dist_bf_phase_ind<-merge(metadata[c("patient","sampleid")],dist_bf_phase_ind,by.x = "sampleid",by.y = "Var2")
colnames(dist_bf_phase_ind)<-c("sampleid0","patient0","sampleid1","patient1","value","Type1","Type2")

dist_bf_phase_ind_intra<-subset(dist_bf_phase_ind,patient0==patient1)
dist_bf_phase_ind_intra_sq<-subset(dist_bf_phase_ind_intra,(Type2=="P"&Type1=="A")|(Type2=="A"&Type1=="E")|(Type2=="E"&Type1=="E30")|(Type2=="E30"&Type1=="E75"))
dist_bf_phase_ind_intra_sq$site<-"BF"

metadata<-as(sample_data(physeq_SRS_tmo_fcg),"data.frame")
metadata$sampleid<-rownames(metadata)
sample_data(physeq_SRS_tmo_fcg)<-metadata
dist_fcg_phase_wu<-plotDistances(physeq_SRS_tmo_fcg,m="wunifrac",plot=F,s="sampleid",d="time")

#pick patients names and remove non-intrapatient comparisons
dist_fcg_phase_ind<-merge(metadata[c("patient","sampleid")],dist_fcg_phase_wu,by.x = "sampleid",by.y = "Var1")
colnames(dist_fcg_phase_ind)<-c("sampleid","patient1","Var2","value","Type1","Type2")
dist_fcg_phase_ind<-merge(metadata[c("patient","sampleid")],dist_fcg_phase_ind,by.x = "sampleid",by.y = "Var2")
colnames(dist_fcg_phase_ind)<-c("sampleid0","patient0","sampleid1","patient1","value","Type1","Type2")

dist_fcg_phase_ind_intra<-subset(dist_fcg_phase_ind,patient0==patient1)
dist_fcg_phase_ind_intra_sq<-subset(dist_fcg_phase_ind_intra,(Type2=="P"&Type1=="A")|(Type2=="A"&Type1=="E")|(Type2=="E"&Type1=="E30")|(Type2=="E30"&Type1=="E75"))
dist_fcg_phase_ind_intra_sq$site<-"FCG"

dist_all_phase_ind_intra_sq<-rbind(dist_mo_phase_ind_intra_sq, dist_bf_phase_ind_intra_sq, dist_fcg_phase_ind_intra_sq)
dist_all_phase_ind_intra_sq$subseq<-interaction(dist_all_phase_ind_intra_sq$Type2,dist_all_phase_ind_intra_sq$Type1)
dist_all_phase_ind_intra_sq$subseq<-factor(dist_all_phase_ind_intra_sq$subseq, levels = c("P.A","A.E","E.E30","E30.E75"))

#fazer o ajuste pelo tempo entre coletas
metadata<-as(sample_data(physeq_SRS_tmo),"data.frame")
dist_all_phase_ind_intra_sq$interval<-0
for (i in 1:nrow(dist_all_phase_ind_intra_sq)){
  days.between<-subset(metadata,site==dist_all_phase_ind_intra_sq[["site"]][i]&patient==dist_all_phase_ind_intra_sq[["patient0"]][i]&
                         time==dist_all_phase_ind_intra_sq[["Type1"]][i])$days.to.E-
    subset(metadata,site==dist_all_phase_ind_intra_sq[["site"]][i]&patient==dist_all_phase_ind_intra_sq[["patient0"]][i]&
             time==dist_all_phase_ind_intra_sq[["Type2"]][i])$days.to.E
  print(days.between)
  dist_all_phase_ind_intra_sq[["interval"]][i]<-days.between
}

dist_all_phase_ind_intra_sq[dist_all_phase_ind_intra_sq=="BF"]<-"SB"
dist_all_phase_ind_intra_sq[dist_all_phase_ind_intra_sq=="FCG"]<-"GCF"
dist_all_phase_ind_intra_sq[dist_all_phase_ind_intra_sq=="MO"]<-"OM"

F6B<-ggplot(data=dist_all_phase_ind_intra_sq,aes(x = subseq, y = value/log(interval), color = site)) + 
  geom_boxplot(alpha = 0.5, outlier.size = 1) + facet_wrap(.~site)+
  #geom_point(alpha=0.5, size=1)+
  guides(color = 'none')+
  scale_color_manual(values = site_colors)+
  scale_y_continuous(expand = c(0,0), breaks = c(0,0.2,0.4,0.6), limits = c(0,0.6))+
  scale_x_discrete(labels = c("P.A"="P-A","A.E"="A-E","E.E30"="E-E30","E30.E75"="E30-E75"))+
  #geom_line(aes(group=patient0),size=0.1,alpha=0.5)+
  labs(y = "Compositional drift", x = 'Interval')+
  stat_summary(aes(group=1),fun = "median",geom="line",size=1,alpha=1)+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
F6B

#plot effect size difference between Pvsotherphase per site (calculated elsewhere)
F_phase_P<-data.frame("site"=c("SB","SB","SB","SB","GCF","GCF","GCF","GCF","OM","OM","OM","OM"),
                    "F."=c(3.61,6.21,3.23,5.90,2.55,3.49,3.80,5.23,3.79,4.73,1.81,1.52),
                    "P"=c(0.011,0.001,0.004,0.002,0.042,0.005,0.004,0.004,0.01,0.003,0.121,0.175),
                    "time"=c("P vs A","P vs E","P vs E30","P vs E75"))
F_phase_P$P.signif<-ifelse(F_phase_P$P<.05,F_phase_P$P.signif<-"Yes",F_phase_P$P.signif<-"No")

F2D<-ggplot(F_phase_P, aes(x=time,y=F.,color=site,group=site,shape=P.signif))+geom_line()+
  geom_point(size=3.2)+labs(x="Comparison",y="PERMANOVA F",shape="P < .05",color="Site")+
  scale_color_manual(values=site_colors)+
  scale_shape_manual(values=c(1,8))+
  scale_y_continuous(limits = c(0,6.5), expand = c(0,0))+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))
F2D

#LETS COMPARE SITES#####
ord_tmo<-rbiom::unifrac(otu_table(physeq_SRS_tmo),weighted=T, tree = phy_tree(physeq_SRS_tmo))
adonis(ord_tmo ~ site, as(sample_data(physeq_SRS_tmo), "data.frame"))
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# site        2     3.040 1.52014   11.59 0.05037  0.001 ***
#   Residuals 437    57.315 0.13115         0.94963           
# Total     439    60.355                 1.00000 

# pairwise.adonis(ord_tmo,as(sample_data(physeq_SRS_tmo), "data.frame")$site)
# pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1  MO vs BF  1 0.6143586  4.769275 0.01617969   0.003      0.009   *
#   2 MO vs FCG  1 1.8636365 13.979507 0.04583753   0.001      0.003   *
#   3 BF vs FCG  1 2.0757047 15.805728 0.05118340   0.001      0.003   *
#now lets compare sites



#o mesmo, mas em cada fase separadamente
ord_tmo_phase<-rbiom::unifrac(otu_table(subset_samples(physeq_SRS_tmo,time=="P")),weighted=T, tree = phy_tree(subset_samples(physeq_SRS_tmo,time=="P")))
adonis(ord_tmo_phase ~ site, as(sample_data(subset_samples(physeq_SRS_tmo,time=="P")), "data.frame"))
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# site       2    1.1945 0.59724  6.3893 0.12555  0.001 ***
#   Residuals 89    8.3193 0.09348         0.87445           
# Total     91    9.5138                 1.00000
pairwise.adonis(ord_tmo_phase,as(sample_data(subset_samples(physeq_SRS_tmo,time=="P")), "data.frame")$site)
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1  MO vs BF  1 0.2706044 2.982168 0.04811332   0.016      0.048   .
# 2 MO vs FCG  1 0.8728765 9.690889 0.14107968   0.001      0.003   *
# 3 BF vs FCG  1 0.6474280 6.506061 0.09782659   0.001      0.003   *

ord_tmo_phase<-rbiom::unifrac(otu_table(subset_samples(physeq_SRS_tmo,time=="A")),weighted=T, tree = phy_tree(subset_samples(physeq_SRS_tmo,time=="A")))
adonis(ord_tmo_phase ~ site, as(sample_data(subset_samples(physeq_SRS_tmo,time=="A")), "data.frame"))
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# site       2    0.8261 0.41303  3.0119 0.06273  0.005 **
# Residuals 90   12.3419 0.13713         0.93727          
# Total     92   13.1679                 1.00000
pairwise.adonis(ord_tmo_phase,as(sample_data(subset_samples(physeq_SRS_tmo,time=="A")), "data.frame")$site)
# pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1  MO vs BF  1 0.1069915 0.7891616 0.01298195   0.509      1.000    
# 2 MO vs FCG  1 0.6373036 4.4975769 0.06973249   0.004      0.012   .
# 3 BF vs FCG  1 0.4947859 3.6891334 0.05792406   0.005      0.015   
ord_tmo_phase<-rbiom::unifrac(otu_table(subset_samples(physeq_SRS_tmo,time=="E")),weighted=T, tree = phy_tree(subset_samples(physeq_SRS_tmo,time=="E")))
adonis(ord_tmo_phase ~ site, as(sample_data(subset_samples(physeq_SRS_tmo,time=="E")), "data.frame"))
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# site       2    0.8269 0.41346  2.4899 0.05414   0.02 *
#   Residuals 87   14.4465 0.16605         0.94586         
# Total     89   15.2734                 1.00000 
pairwise.adonis(ord_tmo_phase,as(sample_data(subset_samples(physeq_SRS_tmo,time=="E")), "data.frame")$site)
# pairs Df  SumsOfSqs   F.Model          R2 p.value p.adjusted sig
# 1  MO vs BF  1 0.05681239 0.3361893 0.005863475   0.902      1.000    
# 2 MO vs FCG  1 0.58832370 3.7063406 0.060064178   0.009      0.027   .
# 3 BF vs FCG  1 0.58641242 3.4412259 0.055111440   0.013      0.039   .
ord_tmo_phase<-rbiom::unifrac(otu_table(subset_samples(physeq_SRS_tmo,time=="E30")),weighted=T, tree = phy_tree(subset_samples(physeq_SRS_tmo,time=="E30")))
adonis(ord_tmo_phase ~ site, as(sample_data(subset_samples(physeq_SRS_tmo,time=="E30")), "data.frame"))
# Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# site       2    0.4848 0.24239  1.8199 0.04016  0.076 .
# Residuals 87   11.5878 0.13319         0.95984         
# Total     89   12.0726                 1.00000 
pairwise.adonis(ord_tmo_phase,as(sample_data(subset_samples(physeq_SRS_tmo,time=="E30")), "data.frame")$site)
# pairs Df SumsOfSqs  F.Model         R2 p.value p.adjusted sig
# 1  MO vs BF  1 0.2944604 2.216273 0.03680522   0.060      0.180    
# 2 MO vs FCG  1 0.2031679 1.396796 0.02351635   0.201      0.603    
# 3 BF vs FCG  1 0.2295511 1.892997 0.03160631   0.076      0.228  
ord_tmo_phase<-rbiom::unifrac(otu_table(subset_samples(physeq_SRS_tmo,time=="E75")),weighted=T, tree = phy_tree(subset_samples(physeq_SRS_tmo,time=="E75")))
adonis(ord_tmo_phase ~ site, as(sample_data(subset_samples(physeq_SRS_tmo,time=="E75")), "data.frame"))
# Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
# site       2    0.6000 0.300022  3.7538 0.09443  0.003 **
#   Residuals 72    5.7546 0.079925         0.90557          
# Total     74    6.3546                  1.00000 
pairwise.adonis(ord_tmo_phase,as(sample_data(subset_samples(physeq_SRS_tmo,time=="E75")), "data.frame")$site)
# 1  MO vs BF  1 0.3134952 4.448210 0.08481147   0.007      0.021   .
# 2 MO vs FCG  1 0.2080866 2.272351 0.04520081   0.069      0.207    
# 3 BF vs FCG  1 0.3784845 4.869584 0.09210558   0.004      0.012   .

metadata<-as.data.frame(sample_data(physeq_SRS_tmo))
metadata$site<-factor(metadata$site, levels = c('FCG','MO','BF'))
sample_data(physeq_SRS_tmo)<-metadata

F1A<-ggarrange(phyloseq::plot_ordination(physeq_SRS_tmo, ordinate(physeq_SRS_tmo, "PCoA", 
                                                             distance = rbiom::unifrac(otu_table(physeq_SRS_tmo), 
                                                                                       weighted=T, tree = phy_tree(physeq_SRS_tmo))),
                                    type="samples", color = "site",shape = "time") + scale_shape_manual(values=c(19,NA,NA,NA,NA))+
            geom_point(size=1.5, alpha = 0.8) + facet_grid(. ~ "P") + labs(color="Site",x=" ", y=" ")+
            scale_color_manual(values=site_colors,labels=c("BF"="SB","MO"="OM","FCG"="GCF"))+scale_fill_manual(values=site_colors,labels=c("BF"="SB","MO"="OM","FCG"="GCF"))+
            scale_x_continuous(limits = c(-0.75,0.75),breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75), expand = c(0,0))+scale_y_continuous(limits = c(-0.5,0.5),breaks = c(-0.5,-0.25,0,0.25,0.5), expand = c(0,0))+
            stat_ellipse(data=subset(phyloseq::plot_ordination(physeq_SRS_tmo, ordinate(physeq_SRS_tmo, "PCoA", 
                                                                                        distance = rbiom::unifrac(otu_table(physeq_SRS_tmo), 
                                                                                                                  weighted=T, tree = phy_tree(physeq_SRS_tmo))),
                                                               type="samples", color = "site",shape = "time")$data,time=="P"),geom = "polygon", level=0.95, alpha=0.1, aes(fill=site))+guides(fill=F,shape=F)+theme(panel.grid = element_blank(),strip.background = element_rect(color="white", fill="white", linetype="solid"),
                                                                                                                                                                                                                    panel.border = element_blank(),axis.line = element_line(color = 'black'),plot.margin = margin(t=-5,b=0)),#+
          #annotate("text", x = -Inf, y = -Inf, color = "black",label = "Global: P = 0.001; F = 6.39\nBD vs. FCG: P = 0.002; F = 6.51\nBD vs. MO: P = 0.018; F = 2.98\nFCG vs. MO: P = 0.001; F = 9.69", vjust = 0, hjust = 0, size=3),
          
          phyloseq::plot_ordination(physeq_SRS_tmo, ordinate(physeq_SRS_tmo, "PCoA", 
                                                             distance = rbiom::unifrac(otu_table(physeq_SRS_tmo), 
                                                                                       weighted=T, tree = phy_tree(physeq_SRS_tmo))),
                                    type="samples", color = "site",shape = "time") + scale_shape_manual(values=c(NA,19,NA,NA,NA))+
            geom_point(size=1.5, alpha = 0.8) + facet_grid(. ~ "A") + labs(color="Site",x=" ", y=" ")+
            scale_color_manual(values=site_colors,labels=c("BF"="SB","MO"="OM","FCG"="GCF"))+scale_fill_manual(values=site_colors,labels=c("BF"="SB","MO"="OM","FCG"="GCF"))+
            scale_x_continuous(limits = c(-0.75,0.75),breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75), expand = c(0,0))+scale_y_continuous(limits = c(-0.5,0.5),breaks = c(-0.5,-0.25,0,0.25,0.5), expand = c(0,0))+
            stat_ellipse(data=subset(phyloseq::plot_ordination(physeq_SRS_tmo, ordinate(physeq_SRS_tmo, "PCoA", 
                                                                                        distance = rbiom::unifrac(otu_table(physeq_SRS_tmo), 
                                                                                                                  weighted=T, tree = phy_tree(physeq_SRS_tmo))),
                                                               type="samples", color = "site",shape = "time")$data,time=="A"),geom = "polygon", level=0.95, alpha=0.1, aes(fill=site))+guides(fill=F,shape=F)+theme(panel.grid = element_blank(),strip.background = element_rect(color="white", fill="white", linetype="solid"),
                                                                                                                                                                                                                    panel.border = element_blank(),axis.line = element_line(color = 'black'),plot.margin = margin(t=-5,b=0)),#+
          #annotate("text", x = -Inf, y = -Inf, color = "black",label = "Global: P = 0.005; F = 3.01\nBD vs. FCG: P = 0.009; F = 3.69\nBD vs. MO: P = 0.52; F = 0.79\nFCG vs. MO: P = 0.004; F = 4.50", vjust = 0, hjust = 0, size=3),
          phyloseq::plot_ordination(physeq_SRS_tmo, ordinate(physeq_SRS_tmo, "PCoA", 
                                                             distance = rbiom::unifrac(otu_table(physeq_SRS_tmo), 
                                                                                       weighted=T, tree = phy_tree(physeq_SRS_tmo))),
                                    type="samples", color = "site",shape = "time") + scale_shape_manual(values=c(NA,NA,19,NA,NA))+
            geom_point(size=1.5, alpha = 0.8) + facet_grid(. ~ "E") + labs(color="Site",x=" ", y="Axis #2 [14%]")+
            scale_x_continuous(limits = c(-0.75,0.75),breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75), expand = c(0,0))+scale_y_continuous(limits = c(-0.5,0.5),breaks = c(-0.5,-0.25,0,0.25,0.5), expand = c(0,0))+
            scale_color_manual(values=site_colors,labels=c("BF"="SB","MO"="OM","FCG"="GCF"))+scale_fill_manual(values=site_colors,labels=c("BF"="SB","MO"="OM","FCG"="GCF"))+
            stat_ellipse(data=subset(phyloseq::plot_ordination(physeq_SRS_tmo, ordinate(physeq_SRS_tmo, "PCoA", 
                                                                                        distance = rbiom::unifrac(otu_table(physeq_SRS_tmo), 
                                                                                                                  weighted=T, tree = phy_tree(physeq_SRS_tmo))),
                                                               type="samples", color = "site",shape = "time")$data,time=="E"),geom = "polygon", level=0.95, alpha=0.1, aes(fill=site))+guides(fill=F,shape=F)+theme(panel.grid = element_blank(),strip.background = element_rect(color="white", fill="white", linetype="solid"),
                                                                                                                                                                                                                    panel.border = element_blank(),axis.line = element_line(color = 'black'),plot.margin = margin(t=-5,b=0)),#+
          #annotate("text", x = -Inf, y = -Inf, color = "black",label = "Global: P = 0.02; F = 2.49\nBD vs. FCG: P = 0.021; F = 3.44\nBD vs. MO: P = 0.89; F = 0.34\nFCG vs. MO: P = 0.01; F = 3.71", vjust = 0, hjust = 0, size=3),
          phyloseq::plot_ordination(physeq_SRS_tmo, ordinate(physeq_SRS_tmo, "PCoA", 
                                                             distance = rbiom::unifrac(otu_table(physeq_SRS_tmo), 
                                                                                       weighted=T, tree = phy_tree(physeq_SRS_tmo))),
                                    type="samples", color = "site",shape = "time") + scale_shape_manual(values=c(NA,NA,NA,19,NA))+
            geom_point(size=1.5, alpha = 0.8) + facet_grid(. ~ "E30") + labs(color="Site",x=" ", y=" ")+
            scale_color_manual(values=site_colors,labels=c("BF"="SB","MO"="OM","FCG"="GCF"))+scale_fill_manual(values=site_colors,labels=c("BF"="SB","MO"="OM","FCG"="GCF"))+
            scale_x_continuous(limits = c(-0.75,0.75),breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75), expand = c(0,0))+scale_y_continuous(limits = c(-0.5,0.5),breaks = c(-0.5,-0.25,0,0.25,0.5), expand = c(0,0))+
            stat_ellipse(data=subset(phyloseq::plot_ordination(physeq_SRS_tmo, ordinate(physeq_SRS_tmo, "PCoA", 
                                                                                        distance = rbiom::unifrac(otu_table(physeq_SRS_tmo), 
                                                                                                                  weighted=T, tree = phy_tree(physeq_SRS_tmo))),
                                                               type="samples", color = "site",shape = "time")$data,time=="E30"),geom = "polygon", level=0.95, alpha=0.1, aes(fill=site))+guides(fill=F,shape=F)+theme(panel.grid = element_blank(),strip.background = element_rect(color="white", fill="white", linetype="solid"),
                                                                                                                                                                                                                      panel.border = element_blank(),axis.line = element_line(color = 'black'),plot.margin = margin(t=-5,b=0)),#+
          #annotate("text", x = -Inf, y = -Inf, color = "black",label = "Global: P = 0.076; F = 1.82\nBD vs. FCG: P = 0.078; F = 1.89\nBD vs. MO: P = 0.064; F = 2.22\nFCG vs. MO: P = 0.2; F = 1.40", vjust = 0, hjust = 0, size=3),
          phyloseq::plot_ordination(physeq_SRS_tmo, ordinate(physeq_SRS_tmo, "PCoA", 
                                                             distance = rbiom::unifrac(otu_table(physeq_SRS_tmo), 
                                                                                       weighted=T, tree = phy_tree(physeq_SRS_tmo))),
                                    type="samples", color = "site",shape = "time") + scale_shape_manual(values=c(NA,NA,NA,NA,19))+
            geom_point(size=1.5, alpha = 0.8) + facet_grid(. ~ "E75") + labs(color="Site",x="Axis #1 [41.5%]", y=" ")+
            scale_color_manual(values=site_colors,labels=c("BF"="SB","MO"="OM","FCG"="GCF"))+scale_fill_manual(values=site_colors,labels=c("BF"="SB","MO"="OM","FCG"="GCF"))+
            scale_x_continuous(limits = c(-0.75,0.75),breaks = c(-0.75,-0.5,-0.25,0,0.25,0.5,0.75), expand = c(0,0))+scale_y_continuous(limits = c(-0.5,0.5),breaks = c(-0.5,-0.25,0,0.25,0.5), expand = c(0,0))+
            stat_ellipse(data=subset(phyloseq::plot_ordination(physeq_SRS_tmo, ordinate(physeq_SRS_tmo, "PCoA", 
                                                                                        distance = rbiom::unifrac(otu_table(physeq_SRS_tmo), 
                                                                                                                  weighted=T, tree = phy_tree(physeq_SRS_tmo))),
                                                               type="samples", color = "site",shape = "time")$data,time=="E75"),geom = "polygon", level=0.95, alpha=0.1, aes(fill=site))+guides(fill=F,shape=F)+theme(panel.grid = element_blank(),strip.background = element_rect(color="white", fill="white", linetype="solid"),
                                                                                                                                                                                                                      panel.border = element_blank(),axis.line = element_line(color = 'black'),plot.margin = margin(t=-5,b=0)),#+
          #annotate("text", x = -Inf, y = -Inf, color = "black",label = "Global: P = 0.003; F = 3.75\nBD vs. FCG: P = 0.004; F = 4.87\nBD vs. MO: P = 0.002; F = 4.45\nFCG vs. MO: P = 0.069; F = 2.27", vjust = 0, hjust = 0, size=3),
          nrow = 5, common.legend = T)
F1A

#plot effect size difference between sites per phase (F vs phase)
F_phase<-data.frame("comparison"=c("GCF vs SB","OM vs SB","GCF vs OM","GCF vs SB","OM vs SB","GCF vs OM","GCF vs SB","OM vs SB","GCF vs OM","GCF vs SB","OM vs SB","GCF vs OM","GCF vs SB","OM vs SB","GCF vs OM"),
                    "F."=c(6.51,2.98,9.69,3.69,0.79,4.50,3.44,0.34,3.71,1.89,2.22,1.40,4.87,4.45,2.27),
                    "P"=c(0.002,0.018,0.001,0.009,0.52,0.004,0.021,0.89,0.01,0.078,0.064,0.199,0.004,0.002,0.069),
                    "time"=c("P","P","P","A","A","A","E","E","E","E30","E30","E30","E75","E75","E75"))
F_phase$P.signif<-ifelse(F_phase$P<.05,F_phase$P.signif<-"Yes",F_phase$P.signif<-"No")
F_phase$time<-factor(F_phase$time, levels=c("P","A","E","E30","E75"))
F_phase$comparison<-factor(F_phase$comparison, levels=c("GCF vs OM","GCF vs SB","OM vs SB"))

F1B<-ggplot(F_phase, aes(x=time,y=F.,group=comparison,shape=P.signif, color = comparison))+
  geom_line()+scale_shape_manual(values=c(1,8))+
  scale_y_continuous(expand=c(0,0), breaks = c(0,2,6,8,10),limits = c(0,10))+
  scale_color_manual(values = site_colors)+
  geom_point(size=3.2)+labs(x="Timepoint",y="PERMANOVA F",shape="P < .05",color="Comparison")+
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),axis.line = element_line(color = 'black'))
F1B


#plot beta-volatility (distance to P)#######################################
#plot beta-volatility (distance to P)
#MO
metadata<-as(sample_data(physeq_SRS_tmo_mo),"data.frame")
metadata$sampleid<-rownames(metadata)
sample_data(physeq_SRS_tmo_mo)<-metadata
dist_mo_phase_wu<-plotDistances(physeq_SRS_tmo_mo,m="wunifrac",plot=F,s="sampleid",d="time")

#pick patients names and remove non-intrapatient comparisons
dist_mo_phase_ind<-merge(metadata[c("patient","sampleid")],dist_mo_phase_wu,by.x = "sampleid",by.y = "Var1")
colnames(dist_mo_phase_ind)<-c("sampleid","patient1","Var2","value","Type1","Type2")
dist_mo_phase_ind<-merge(metadata[c("patient","sampleid")],dist_mo_phase_ind,by.x = "sampleid",by.y = "Var2")
colnames(dist_mo_phase_ind)<-c("sampleid0","patient0","sampleid1","patient1","value","Type1","Type2")
dist_mo_phase_ind<-subset(dist_mo_phase_ind,patient0==patient1 & Type1=="P")

# #AAM dont have P and AMSCP don have post-engr samples
# #lets remove AMSCP only, since AAM has already been removed
dist_mo_phase_ind<-subset(dist_mo_phase_ind,patient0!="AMSCP")

#now we can use these values to classify recoverers and non-recoverers, 
#but first we need to include columns with the final distance
dist_mo_phase_ind_class<-dist_mo_phase_ind
dist_mo_phase_ind_class$fdis.e30<-NA
dist_mo_phase_ind_class$fdis.e30.class<-NA

for (i in 1:nrow(dist_mo_phase_ind)){
  dist_mo_phase_ind_class[["fdis.e30"]][i]<-subset(dist_mo_phase_ind,patient0==dist_mo_phase_ind[["patient0"]][i]&Type2=="E30")$value
  if(dist_mo_phase_ind_class[["fdis.e30"]][i]<0.5){dist_mo_phase_ind_class[["fdis.e30.class"]][i]<-"R"}else{
    dist_mo_phase_ind_class[["fdis.e30.class"]][i]<-"NR"}
  }

dist_mo_phase_ind_class <- mutate_if(dist_mo_phase_ind_class, is.character, as.factor) #char to factors

#length(unique(subset(dist_bf_phase_ind_class,fdis.e30.class=="R"&value>=0.5)$patient0)) #n de R que utltrapassam 0.5 nas fases anteriores

#BF
metadata<-as(sample_data(physeq_SRS_tmo_bf),"data.frame")
metadata$sampleid<-rownames(metadata)
sample_data(physeq_SRS_tmo_bf)<-metadata
dist_bf_phase_wu<-plotDistances(physeq_SRS_tmo_bf,m="wunifrac",plot=F,s="sampleid",d="time")

#pick patients names and remove non-intrapatient comparisons
dist_bf_phase_ind<-merge(metadata[c("patient","sampleid")],dist_bf_phase_wu,by.x = "sampleid",by.y = "Var1")
colnames(dist_bf_phase_ind)<-c("sampleid","patient1","Var2","value","Type1","Type2")
dist_bf_phase_ind<-merge(metadata[c("patient","sampleid")],dist_bf_phase_ind,by.x = "sampleid",by.y = "Var2")
colnames(dist_bf_phase_ind)<-c("sampleid0","patient0","sampleid1","patient1","value","Type1","Type2")
dist_bf_phase_ind<-subset(dist_bf_phase_ind,patient0==patient1 & Type1=="P")

#AMSCP don have post-engr samples
#lets remove AMSCP only
dist_bf_phase_ind<-subset(dist_bf_phase_ind,patient0!="AMSCP")

#now we can use these values to classify recoverers and non-recoverers, 
#but first we need to include columns with the final distance
dist_bf_phase_ind_class<-dist_bf_phase_ind
dist_bf_phase_ind_class$fdis.e30<-NA
dist_bf_phase_ind_class$fdis.e30.class<-NA

for (i in 1:nrow(dist_bf_phase_ind)){
  dist_bf_phase_ind_class[["fdis.e30"]][i]<-subset(dist_bf_phase_ind,patient0==dist_bf_phase_ind[["patient0"]][i]&Type2=="E30")$value
  if(dist_bf_phase_ind_class[["fdis.e30"]][i]<0.5){dist_bf_phase_ind_class[["fdis.e30.class"]][i]<-"R"}else{
    dist_bf_phase_ind_class[["fdis.e30.class"]][i]<-"NR"}
}

dist_bf_phase_ind_class <- mutate_if(dist_bf_phase_ind_class, is.character, as.factor) #char to factors

#FCG
metadata<-as(sample_data(physeq_SRS_tmo_fcg),"data.frame")
metadata$sampleid<-rownames(metadata)
sample_data(physeq_SRS_tmo_fcg)<-metadata
dist_fcg_phase_wu<-plotDistances(physeq_SRS_tmo_fcg,m="wunifrac",plot=F,s="sampleid",d="time")

#pick patients names and remove non-intrapatient comparisons
dist_fcg_phase_ind<-merge(metadata[c("patient","sampleid")],dist_fcg_phase_wu,by.x = "sampleid",by.y = "Var1")
colnames(dist_fcg_phase_ind)<-c("sampleid","patient1","Var2","value","Type1","Type2")
dist_fcg_phase_ind<-merge(metadata[c("patient","sampleid")],dist_fcg_phase_ind,by.x = "sampleid",by.y = "Var2")
colnames(dist_fcg_phase_ind)<-c("sampleid0","patient0","sampleid1","patient1","value","Type1","Type2")
dist_fcg_phase_ind<-subset(dist_fcg_phase_ind,patient0==patient1 & Type1=="P")

#AMSCP don have post-engr samples
#lets remove AMSCP only
dist_fcg_phase_ind<-subset(dist_fcg_phase_ind,patient0!="AMSCP")

#now we can use these values to classify recoverers and non-recoverers, 
#but first we need to include columns with the final distance
dist_fcg_phase_ind_class<-dist_fcg_phase_ind
dist_fcg_phase_ind_class$fdis.e30<-NA
dist_fcg_phase_ind_class$fdis.e30.class<-NA
for (i in 1:nrow(dist_fcg_phase_ind)){
  dist_fcg_phase_ind_class[["fdis.e30"]][i]<-subset(dist_fcg_phase_ind,patient0==dist_fcg_phase_ind[["patient0"]][i]&Type2=="E30")$value
  if(dist_fcg_phase_ind_class[["fdis.e30"]][i]<0.5){dist_fcg_phase_ind_class[["fdis.e30.class"]][i]<-"R"}else{
    dist_fcg_phase_ind_class[["fdis.e30.class"]][i]<-"NR"}
}

dist_fcg_phase_ind_class <- mutate_if(dist_fcg_phase_ind_class, is.character, as.factor) #char to factors

dist_mo_phase_ind_class$site<-"OM"
dist_bf_phase_ind_class$site<-"SB"
dist_fcg_phase_ind_class$site<-"GCF"

dist_from_P<-rbind(dist_mo_phase_ind_class,dist_bf_phase_ind_class,dist_fcg_phase_ind_class)

#plot distances to P per patient (Beta-vol)
F6E<-ggplot(data=subset(dist_from_P, Type2%in%c('E','E30')),aes(color=site, x = Type2, y = value, group = patient0))+
  geom_line(alpha=0.5,size=0.3) + facet_wrap(.~site)+guides(color = F)+
  stat_summary(aes(group=site),fun.y = "median",geom="line",size=1.8,alpha=0.6)+
  labs(y = "Distance to P", x = 'Timepoint') + 
  scale_x_discrete(expand=c(0.1,0.1))+
  #geom_hline(yintercept=0.5,linetype="dashed", color = "black")+
  scale_color_manual(values=site_colors)+
  scale_y_continuous(limits = c(0,1.11), breaks = c(0,0.25,0.5,0.75,1), expand = c(0,0))+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),axis.line = element_line(color = 'black'),
        strip.background = element_rect(color="white", fill="white", linetype="solid"))
F6E

#plot heatmap with all classifications
dist_from_P_class<-dist_from_P[!duplicated(dist_from_P[,c("patient0","site")]),]
dist_from_P_class$patient0 <- with(dist_from_P_class,factor(patient0,levels = rev(sort(unique(patient0)))))
dist_from_P_class<-dist_from_P_class[c("patient0","site","fdis.e30.class")]
dist_from_P_class.melt<-melt(dist_from_P_class,id.vars=c("patient0","site"),
                                variable.name="div_beta_eco_metric",value.name="div_beta_eco_measure_class")
dist_from_P_class.melt$patient0<-as.character(dist_from_P_class.melt$patient0)
dist_from_P_class.melt$patient0<-factor(dist_from_P_class.melt$patient0,
                                        levels = sort(unique(dist_from_P_class.melt$patient0)))

F7B<-ggplot(dist_from_P_class.melt)+aes(x=patient0,y=site,fill=div_beta_eco_measure_class)+
  geom_tile(color="white")+coord_equal()+
  scale_fill_manual(values=c('#616161', '#008d78'))+
  scale_x_discrete(labels=pacientes)+labs(x="Patient",y="Site",fill="Class")+
  theme(legend.justification = "top", panel.grid = element_blank(), 
        axis.text.x=element_text(angle=90,hjust=1,vjust =0.5))

#manipular um pouco os dfs para salvar r/nr e fdis.e30
#dist_from_P_class_tosave<-dist_from_P[!duplicated(dist_from_P[,c("patient0","site")]),]
#dist_from_P_class_tosave$patient0 <- with(dist_from_P_class_tosave,factor(patient0,levels = rev(sort(unique(patient0)))))
#dist_from_P_class_tosave<-dist_from_P_class_tosave[c("patient0","site","fdis.e30.class","fdis.e30")]

#dist_from_P_class_tosave<-merge(dcast(dist_from_P_class_tosave,patient0~site, value.var = c("fdis.e30.class")) %>% `colnames<-`(c("patient","DB.fdis.e30.class", "GCF.fdis.e30.class", "OM.fdis.e30.class")),
#dcast(dist_from_P_class_tosave,patient0~site, value.var = c("fdis.e30")) %>% `colnames<-`(c("patient","DB.fdis.e30", "GCF.fdis.e30", "OM.fdis.e30")))
#dist_from_P_class_tosave$patient<-as.character(dist_from_P_class_tosave$patient)
#dist_from_P_class_tosave<-rbind(dist_from_P_class_tosave,c("AMSCP",NA,NA,NA,NA,NA,NA))

#write.csv(dist_from_P_class_tosave[order(dist_from_P_class_tosave$patient),],
#          file = "Files/tempfile.csv", row.names = F) #chamei de tempfile pois vou excluir apos juntar manualmente com arquivo microbiota-variables que foi gerado no script da alpha

##########16/04 - pcoas individuais por paciente, ordenados na grade pela distância final. obviamento separando por sítio, e tanto considerando fdis do e30 e do e75

#função que gera um vetor de tamanho n repleto da entrada all, exceto na posição k, em que entra exc
values_vector<-function(k,all,exc,list){
  vec<-rep(all,length(list)) #inicializa o vetor
  names(vec)<-list
  vec[k]<-exc
  return(vec)
}

#BF - E30
#recalculate fdis
metadata<-as(sample_data(physeq_SRS_tmo_bf),"data.frame")
metadata$sampleid<-rownames(metadata)
sample_data(physeq_SRS_tmo_bf)<-metadata
#calc for wunifrac-curtis beta-div metric
dist_bf_phase_wu<-plotDistances(p = physeq_SRS_tmo_bf,m = "wunifrac",s = "sampleid",d = "time",plot = F)
dist_bf_phase_wu$beta_metric<-"wu"

#bind
dist_bf_phase<-dist_bf_phase_wu
#fix factor levels for phase
dist_bf_phase$Type2<-factor(dist_bf_phase$Type2, levels=c("P","A","E","E30","E75"))

#pick patients names and remove non-intrapatient comparisons
dist_bf_phase_ind<-merge(metadata[c("patient","sampleid")],dist_bf_phase,by.x = "sampleid",by.y = "Var1")
colnames(dist_bf_phase_ind)<-c("sampleid","patient1","Var2","value","Type1","Type2","beta_metric")
dist_bf_phase_ind<-merge(metadata[c("patient","sampleid")],dist_bf_phase_ind,by.x = "sampleid",by.y = "Var2")
colnames(dist_bf_phase_ind)<-c("sampleid0","patient0","sampleid1","patient1","value","Type1","Type2","beta_metric")
dist_bf_phase_ind<-subset(dist_bf_phase_ind,patient0==patient1 & Type1=="P")

ordered_patient_list_bf.e30<-subset(dist_bf_phase_ind,Type2=="E30")[order(subset(dist_bf_phase_ind,Type2=="E30")$value),]$patient0
physeq_SRS_tmo_bf_subset<-subset_samples(physeq_SRS_tmo_bf,time!="E75"&patient%in%ordered_patient_list_bf.e30)

p_list<-list()
i=1
for (pt in ordered_patient_list_bf.e30){
  print(pt)
p_list[[i]]<-phyloseq::plot_ordination(physeq_SRS_tmo_bf_subset, ordinate(physeq_SRS_tmo_bf_subset, "PCoA", 
               distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_bf_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_bf_subset))),
                                       type="samples", color = "time", shape = "patient") +
  #geom_point(shape = 1, size = 3,colour = "black") + 
  geom_point(size = 3, alpha = 1)+ 
  geom_path(data=subset(phyloseq::plot_ordination(physeq_SRS_tmo_bf_subset, ordinate(physeq_SRS_tmo_bf_subset, "PCoA", 
            distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_bf_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_bf_subset))),
            type="samples", color = "time")$data,patient==pt)[order(match(subset(phyloseq::plot_ordination(physeq_SRS_tmo_bf_subset, 
            ordinate(physeq_SRS_tmo_bf_subset, "PCoA",distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_bf_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_bf_subset))),
            type="samples", color = "time")$data,patient==pt)$time,c("P","A","E","E30","E75"))),],
            aes(x=Axis.1,y=Axis.2,group=patient),alpha=1,color=ifelse(i<=23, '#008d78', '#616161'),
            arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")),lwd=0.8)+
  labs(color="Timepoint",title = as.character(pacientes[pt]))+
  scale_color_manual(values=c('#CD7700',"gray70","gray62","gray54"))+
  scale_shape_manual(values=values_vector(pt,all=NA,exc=19,ordered_patient_list_bf.e30), guide = F)+
  theme(panel.grid = element_blank(), plot.title = element_text(size=9), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank())
i=i+1
}

FS12C<-annotate_figure(ggarrange(plotlist = p_list, nrow=5,ncol = 6, common.legend = T), 
                fig.lab.pos = "top.left", fig.lab.size = 14, fig.lab = "SB")

#MO - E30
#recalculate fdis
metadata<-as(sample_data(physeq_SRS_tmo_mo),"data.frame")
metadata$sampleid<-rownames(metadata)
sample_data(physeq_SRS_tmo_mo)<-metadata
#calc for wunifrac-curtis beta-div metric
dist_mo_phase_wu<-plotDistances(p = physeq_SRS_tmo_mo,m = "wunifrac",s = "sampleid",d = "time",plot = F)
dist_mo_phase_wu$beta_metric<-"wu"

#bind
dist_mo_phase<-dist_mo_phase_wu
#fix factor levels for phase
dist_mo_phase$Type2<-factor(dist_mo_phase$Type2, levels=c("P","A","E","E30","E75"))

#pick patients names and remove non-intrapatient comparisons
dist_mo_phase_ind<-merge(metadata[c("patient","sampleid")],dist_mo_phase,by.x = "sampleid",by.y = "Var1")
colnames(dist_mo_phase_ind)<-c("sampleid","patient1","Var2","value","Type1","Type2","beta_metric")
dist_mo_phase_ind<-merge(metadata[c("patient","sampleid")],dist_mo_phase_ind,by.x = "sampleid",by.y = "Var2")
colnames(dist_mo_phase_ind)<-c("sampleid0","patient0","sampleid1","patient1","value","Type1","Type2","beta_metric")
dist_mo_phase_ind<-subset(dist_mo_phase_ind,patient0==patient1 & Type1=="P")

ordered_patient_list_mo.e30<-subset(dist_mo_phase_ind,Type2=="E30")[order(subset(dist_mo_phase_ind,Type2=="E30")$value),]$patient0
physeq_SRS_tmo_mo_subset<-subset_samples(physeq_SRS_tmo_mo,time!="E75"&patient%in%ordered_patient_list_mo.e30)

p_list<-list()
i=1
for (pt in ordered_patient_list_mo.e30){
  print(pt)
  p_list[[i]]<-phyloseq::plot_ordination(physeq_SRS_tmo_mo_subset, ordinate(physeq_SRS_tmo_mo_subset, "PCoA", 
                                                                            distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_mo_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_mo_subset))),
                                         type="samples", color = "time", shape = "patient") +
    #geom_point(shape = 1, size = 3,colour = "black") + 
    geom_point(size = 3, alpha = 1)+ 
    geom_path(data=subset(phyloseq::plot_ordination(physeq_SRS_tmo_mo_subset, ordinate(physeq_SRS_tmo_mo_subset, "PCoA", 
                                                                                       distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_mo_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_mo_subset))),
                                                    type="samples", color = "time")$data,patient==pt)[order(match(subset(phyloseq::plot_ordination(physeq_SRS_tmo_mo_subset, 
                                                                                                                                                   ordinate(physeq_SRS_tmo_mo_subset, "PCoA",distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_mo_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_mo_subset))),
                                                                                                                                                   type="samples", color = "time")$data,patient==pt)$time,c("P","A","E","E30","E75"))),],
              aes(x=Axis.1,y=Axis.2,group=patient),alpha=1,color=ifelse(i<=20,'#008d78', '#616161'),
              arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")),lwd=0.8)+
    labs(color="Timepoint",title = as.character(pacientes[pt]))+
    scale_color_manual(values=c('#CD7700',"gray70","gray62","gray54"))+
    scale_shape_manual(values=values_vector(pt,all=NA,exc=19,ordered_patient_list_mo.e30), guide = F)+
    theme(panel.grid = element_blank(), plot.title = element_text(size=9), axis.title = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank())
  i=i+1
}

FS12B<-annotate_figure(ggarrange(plotlist = p_list, nrow=5,ncol = 6, common.legend = T), 
                fig.lab.pos = "top.left", fig.lab.size = 14, fig.lab = "OM")
i=1
pt='JVS'
F7Aa<-phyloseq::plot_ordination(physeq_SRS_tmo_mo_subset, ordinate(physeq_SRS_tmo_mo_subset, "PCoA", 
                                                                   distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_mo_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_mo_subset))),
                                type="samples", color = "time", shape = "patient") +
  #geom_point(shape = 1, size = 3,colour = "black") + 
  geom_point(size = 3, alpha = 1)+ 
  geom_path(data=subset(phyloseq::plot_ordination(physeq_SRS_tmo_mo_subset, ordinate(physeq_SRS_tmo_mo_subset, "PCoA", 
                                                                                     distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_mo_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_mo_subset))),
                                                  type="samples", color = "time")$data,patient==pt)[order(match(subset(phyloseq::plot_ordination(physeq_SRS_tmo_mo_subset, 
                                                                                                                                                 ordinate(physeq_SRS_tmo_mo_subset, "PCoA",distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_mo_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_mo_subset))),
                                                                                                                                                 type="samples", color = "time")$data,patient==pt)$time,c("P","A","E","E30","E75"))),],
            aes(x=Axis.1,y=Axis.2,group=patient),alpha=1,color=ifelse(i<=20,'#008d78', '#616161'),
            arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")),lwd=0.8)+
  labs(color="Timepoint",title = as.character(pacientes[pt]))+
  scale_color_manual(values=c('#CD7700',"gray70","gray62","gray54"))+
  scale_shape_manual(values=values_vector(pt,all=NA,exc=19,ordered_patient_list_mo.e30), guide = F)+
  theme(panel.grid = element_blank(), plot.title = element_text(size=9), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank())+
  annotate("text",x=0.5,y=-0.35,hjust=1,vjust=0,label="Recoverer")
F7Aa

i=26
pt='OSN'
F7Ab<-phyloseq::plot_ordination(physeq_SRS_tmo_mo_subset, ordinate(physeq_SRS_tmo_mo_subset, "PCoA", 
                                                                   distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_mo_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_mo_subset))),
                                type="samples", color = "time", shape = "patient") +
  #geom_point(shape = 1, size = 3,colour = "black") + 
  geom_point(size = 3, alpha = 1)+ 
  geom_path(data=subset(phyloseq::plot_ordination(physeq_SRS_tmo_mo_subset, ordinate(physeq_SRS_tmo_mo_subset, "PCoA", 
                                                                                     distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_mo_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_mo_subset))),
                                                  type="samples", color = "time")$data,patient==pt)[order(match(subset(phyloseq::plot_ordination(physeq_SRS_tmo_mo_subset, 
                                                                                                                                                 ordinate(physeq_SRS_tmo_mo_subset, "PCoA",distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_mo_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_mo_subset))),
                                                                                                                                                 type="samples", color = "time")$data,patient==pt)$time,c("P","A","E","E30","E75"))),],
            aes(x=Axis.1,y=Axis.2,group=patient),alpha=1,color=ifelse(i<=20,'#008d78', '#616161'),
            arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")),lwd=0.8)+
  labs(color="Timepoint",title = as.character(pacientes[pt]))+
  scale_color_manual(values=c('#CD7700',"gray70","gray62","gray54"))+
  scale_shape_manual(values=values_vector(pt,all=NA,exc=19,ordered_patient_list_mo.e30), guide = F)+
  theme(panel.grid = element_blank(), plot.title = element_text(size=9), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.text = element_blank())+
  annotate("text",x=0.5,y=-0.35,hjust=1,vjust=0,label="Non-recoverer")
F7Ab
  
F7A<-ggarrange(F7Aa, F7Ab, nrow = 1, ncol = 2, common.legend = T)
F7A


#FCG - E30
#recalculate fdis
metadata<-as(sample_data(physeq_SRS_tmo_fcg),"data.frame")
metadata$sampleid<-rownames(metadata)
sample_data(physeq_SRS_tmo_fcg)<-metadata
#calc for wunifrac-curtis beta-div metric
dist_fcg_phase_wu<-plotDistances(p = physeq_SRS_tmo_fcg,m = "wunifrac",s = "sampleid",d = "time",plot = F)
dist_fcg_phase_wu$beta_metric<-"wu"

#bind
dist_fcg_phase<-dist_fcg_phase_wu
#fix factor levels for phase
dist_fcg_phase$Type2<-factor(dist_fcg_phase$Type2, levels=c("P","A","E","E30","E75"))

#pick patients names and remove non-intrapatient comparisons
dist_fcg_phase_ind<-merge(metadata[c("patient","sampleid")],dist_fcg_phase,by.x = "sampleid",by.y = "Var1")
colnames(dist_fcg_phase_ind)<-c("sampleid","patient1","Var2","value","Type1","Type2","beta_metric")
dist_fcg_phase_ind<-merge(metadata[c("patient","sampleid")],dist_fcg_phase_ind,by.x = "sampleid",by.y = "Var2")
colnames(dist_fcg_phase_ind)<-c("sampleid0","patient0","sampleid1","patient1","value","Type1","Type2","beta_metric")
dist_fcg_phase_ind<-subset(dist_fcg_phase_ind,patient0==patient1 & Type1=="P")

ordered_patient_list_fcg.e30<-subset(dist_fcg_phase_ind,Type2=="E30")[order(subset(dist_fcg_phase_ind,Type2=="E30")$value),]$patient0
physeq_SRS_tmo_fcg_subset<-subset_samples(physeq_SRS_tmo_fcg,time!="E75"&patient%in%ordered_patient_list_fcg.e30)

p_list<-list()
i=1
for (pt in ordered_patient_list_fcg.e30){
  print(pt)
  p_list[[i]]<-phyloseq::plot_ordination(physeq_SRS_tmo_fcg_subset, ordinate(physeq_SRS_tmo_fcg_subset, "PCoA", 
                                                                            distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_fcg_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_fcg_subset))),
                                         type="samples", color = "time", shape = "patient") +
    #geom_point(shape = 1, size = 3,colour = "black") + 
    geom_point(size = 3, alpha = 1)+ 
    geom_path(data=subset(phyloseq::plot_ordination(physeq_SRS_tmo_fcg_subset, ordinate(physeq_SRS_tmo_fcg_subset, "PCoA", 
                                                                                       distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_fcg_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_fcg_subset))),
                                                    type="samples", color = "time")$data,patient==pt)[order(match(subset(phyloseq::plot_ordination(physeq_SRS_tmo_fcg_subset, 
                                                                                                                                                   ordinate(physeq_SRS_tmo_fcg_subset, "PCoA",distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_fcg_subset), weighted=T, tree=phy_tree(physeq_SRS_tmo_fcg_subset))),
                                                                                                                                                   type="samples", color = "time")$data,patient==pt)$time,c("P","A","E","E30","E75"))),],
              aes(x=Axis.1,y=Axis.2,group=patient),alpha=1,color=ifelse(i<=23,'#008d78', '#616161'),
              arrow = arrow(type = "open", angle = 30, length = unit(0.1, "inches")),lwd=0.8)+
    labs(color="Timepoint",title = as.character(pacientes[pt]))+
    scale_color_manual(values=c('#CD7700',"gray70","gray62","gray54"))+
    scale_shape_manual(values=values_vector(pt,all=NA,exc=19,ordered_patient_list_fcg.e30), guide = F)+
    theme(panel.grid = element_blank(), plot.title = element_text(size=9), axis.title = element_blank(),
          axis.ticks = element_blank(), axis.text = element_blank())
  i=i+1
}

FS12A<-annotate_figure(ggarrange(plotlist = p_list, nrow=5,ncol = 6, common.legend = T), 
                fig.lab.pos = "top.left", fig.lab.size = 14, fig.lab = "GCF")

#####04/05/21 correlação entre div at P e fdis (precisa ter gerado os dfs usados abaixo antes)
div_dist<-read.csv("Files/300522_microbiota-variables.csv")

# cor.test(bf_fdis.e30_div.p$fdis.e30,bf_fdis.e30_div.p$Simpson.P,method = "spearman") #p-value = 0.413
# cor.test(mo_fdis.e30_div.p$fdis.e30,mo_fdis.e30_div.p$Simpson.P,method = "spearman") #p-value = 0.016
# cor.test(fcg_fdis.e30_div.p$fdis.e30,fcg_fdis.e30_div.p$Simpson.P,method = "spearman") #p-value = 0.51
div_dist_fcg<-div_dist[c("Simpson.P_GCF","GCF.fdis.e30")]
colnames(div_dist_fcg)<-c('Simpson.P','fdis.e30')
div_dist_fcg$site<-'GCF'
div_dist_om<-div_dist[c("Simpson.P_OM","OM.fdis.e30")]
colnames(div_dist_om)<-c('Simpson.P','fdis.e30')
div_dist_om$site<-'OM'
div_dist_db<-div_dist[c("Simpson.P_DB","DB.fdis.e30")]
colnames(div_dist_db)<-c('Simpson.P','fdis.e30')
div_dist_db$site<-'SB'

div_dist_all<-rbind(rbind(div_dist_db,div_dist_fcg),div_dist_om)

F8B<-ggscatter(div_dist_all,y="Simpson.P",x="fdis.e30", add = "reg.line", conf.int = F, facet.by = 'site', alpha=0.5,
          cor.coef = TRUE, cor.method = "spearman",xlab = "Distance to P (E30)", ylab = "Diversity (P)", ylim=c(0.5,1), 
          xlim=c(0,1), cor.coef.coord = c(0.03,0.52),ggtheme = theme_bw(), cor.coef.size = 3.2)+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(color = 'black'),
         strip.background = element_rect(color="white", fill="white", linetype="solid"))
F8B

#####16/08/22 - composition stability
###02/10/21 - volatilidade da beta (convex hull volumes/areas)
find_comp_stab<-function(physeq){
 result<-data.frame(patient = as.character(), site = as.character(), area = as.numeric())
 ord <- rbiom::unifrac(otu_table(physeq),weighted=T, tree = phy_tree(physeq))
 pcoa_result <- data.frame(phyloseq::plot_ordination(physeq, ordinate(physeq, "PCoA",distance = ord),type="samples")["data"])
 for (site in unique(pcoa_result$data.site)){
   pcoa_result_st<-subset(pcoa_result, data.site == site)
   for (patient in unique(pcoa_result_st$data.patient)){
     pcoa_result_pt<-subset(pcoa_result_st, data.patient == patient)
     area<-as.numeric(convhulln(as.matrix(pcoa_result_pt[,1:2]), output.options = 'FA')$area) 
     result<-rbind(result,c(patient,site,area))
   }
 }
 colnames(result)<-c('patient','site','area')
 result$area<-1-as.numeric(result$area)
 return(result)
}


comp_stab_bf<-find_comp_stab(subset_samples(physeq_SRS_tmo,time%in%c('P','E','E30')&
                          site=="BF"&patient!="AMSCP"&patient!="RJT"))
colnames(comp_stab_bf)[3]<-"comp_S.DB"

comp_stab_fcg<-find_comp_stab(subset_samples(physeq_SRS_tmo,time%in%c('P','E','E30')&
                                         site=="FCG"&patient!="AMSCP"))
colnames(comp_stab_fcg)[3]<-"comp_S.GCF"

comp_stab_mo<-find_comp_stab(subset_samples(physeq_SRS_tmo,time%in%c('P','E','E30')&
                                          site=="MO"&patient!="AAM"&patient!="AMSCP"&patient!="CEFC"&patient!="DBHB"))
colnames(comp_stab_mo)[3]<-"comp_S.OM"
 
comp_stab_eachsite<-merge(merge(comp_stab_bf[c(1,3)],comp_stab_fcg[c(1,3)], by = 'patient',all = T), 
                       comp_stab_mo[c(1,3)], by = 'patient',all = T)
comp_stab_eachsite<-rbind(comp_stab_eachsite,c("AMSCP",NA,NA,NA))

comp_stab_eachsite.melt<-reshape2::melt(comp_stab_eachsite,c("patient"))
comp_stab_eachsite.melt$site<-NA
comp_stab_eachsite.melt$site[comp_stab_eachsite.melt$variable=="comp_S.DB"]<-"SB"
comp_stab_eachsite.melt$site[comp_stab_eachsite.melt$variable=="comp_S.GCF"]<-"GCF"
comp_stab_eachsite.melt$site[comp_stab_eachsite.melt$variable=="comp_S.OM"]<-"OM"
comp_stab_eachsite.melt$value<-as.numeric(comp_stab_eachsite.melt$value)

FS6C<-ggplot(data = comp_stab_eachsite.melt) + 
  aes(x = site, y = value, color = site, group = site) + 
  geom_boxplot(alpha=0.5, outlier.size = 1)+
  #geom_quasirandom(size=1)+
  labs(y = 'Compositional stability', x = 'Site') + guides(color='none')+
  scale_y_continuous(expand = c(0,0),limits = c(-1,1.3))+
  scale_color_manual(values = site_colors)+
  stat_compare_means(comparisons = list(c("OM","GCF"),c("OM","SB"),c("SB","GCF")),size = 3.2,
                     label="p.format",method = "wilcox.test",paired = F,step.increase = 0.09, tip.length = 0.005)+ #Mann-whitney
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))

microbiota_vars<-read.csv("Files/300522_microbiota-variables.csv")
#microbiota_vars<-merge(microbiota_vars,comp_stab_eachsite)
#write.csv(microbiota_vars, file = "Files/300522_microbiota-variables.csv")


#comp stability x diversity at P
comp.stab.p.db<-microbiota_vars[c('patient','Simpson.P_DB',"comp_S.DB")]
colnames(comp.stab.p.db)<-c('patient',"div_P","comp_S")
comp.stab.p.db$site<-"SB"
comp.stab.p.gcf<-microbiota_vars[c('patient','Simpson.P_GCF',"comp_S.GCF")]
colnames(comp.stab.p.gcf)<-c('patient',"div_P","comp_S")
comp.stab.p.gcf$site<-"GCF"
comp.stab.p.om<-microbiota_vars[c('patient','Simpson.P_OM',"comp_S.OM")]
colnames(comp.stab.p.om)<-c('patient',"div_P","comp_S")
comp.stab.p.om$site<-"OM"
comp.stab.p<-rbind(rbind(comp.stab.p.db,comp.stab.p.gcf),comp.stab.p.om)
comp.stab.p$comp_S<-as.numeric(comp.stab.p$comp_S)

FS11B<-ggscatter(data = comp.stab.p, y = "div_P", x = "comp_S", facet.by = "site", ylab = "Diversity (P)",alpha=0.5,
          xlab = "Compositional stability", cor.coef = T, cor.method = "spearman", cor.coef.coord = c(-1,0.51),
          add = 'reg.line', conf.int = F,ggtheme = theme_bw(),cor.coef.size = 3.2)+
  theme(panel.grid = element_blank(),panel.border = element_blank(),axis.line = element_line(color = 'black'),
        strip.background = element_rect(color="white", fill="white", linetype="solid"))

#########################################################
FS11C<-ggarrange(phyloseq::plot_ordination(physeq_SRS_tmo_fcg, ordinate(physeq_SRS_tmo_fcg, "PCoA", 
                                                                 distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_fcg), weighted=T, tree=phy_tree(physeq_SRS_tmo_fcg))),
                                    type="samples", color = "time") +
            geom_point(size=1.5, alpha = 0.8) + facet_grid(. ~ "GCF") + labs(color="Timepoint")+
            scale_color_manual(values=c('#f3b300','#b36600','#000000','#376387','#509dc2'))+
            scale_fill_manual(values=c('#f3b300','#b36600','#000000','#376387','#509dc2'))+
            stat_ellipse(geom = "polygon", level=0.95, alpha=0, aes(fill=time))+guides(fill=F)+
            theme(panel.grid = element_blank(),strip.background = element_rect(color="white", fill="white", linetype="solid"),
                   panel.border = element_blank(),axis.line = element_line(color = 'black')),
          phyloseq::plot_ordination(physeq_SRS_tmo_mo, ordinate(physeq_SRS_tmo_mo, "PCoA", 
                                                                distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_mo), weighted=T, tree=phy_tree(physeq_SRS_tmo_mo))),
                                    type="samples", color = "time") +
            geom_point(size=1.5, alpha = 0.8) + facet_grid(. ~ "OM") + labs(color="Timepoint")+
            scale_color_manual(values=c('#f3b300','#b36600','#000000','#376387','#509dc2'))+
            scale_fill_manual(values=c('#f3b300','#b36600','#000000','#376387','#509dc2'))+
            stat_ellipse(geom = "polygon", level=0.95, alpha=0, aes(fill=time))+guides(fill=F)+
            theme(panel.grid = element_blank(),strip.background = element_rect(color="white", fill="white", linetype="solid"),
                   panel.border = element_blank(),axis.line = element_line(color = 'black')),
          phyloseq::plot_ordination(physeq_SRS_tmo_bf, ordinate(physeq_SRS_tmo_bf, "PCoA", 
                                                                distance = rbiom::unifrac(otu_table(physeq_SRS_tmo_bf), weighted=T, tree=phy_tree(physeq_SRS_tmo_bf))),
                                    type="samples", color = "time") +
            geom_point(size=1.5, alpha = 0.8) + facet_grid(. ~ "SB") + labs(color="Timepoint")+
            scale_color_manual(values=c('#f3b300','#b36600','#000000','#376387','#509dc2'))+
            scale_fill_manual(values=c('#f3b300','#b36600','#000000','#376387','#509dc2'))+
            stat_ellipse(geom = "polygon", level=0.95, alpha=0, aes(fill=time))+guides(fill=F)+
            theme(panel.grid = element_blank(),strip.background = element_rect(color="white", fill="white", linetype="solid"),
                   panel.border = element_blank(),axis.line = element_line(color = 'black')),ncol=3,common.legend = T)
