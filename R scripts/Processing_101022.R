require(ggplot2)
theme_set(theme_bw())
library(forcats)
library(scales)
require(gridExtra)
require(grid)
require(egg)
require(lemon)
library(data.table)
require(reshape)
require(purrr)
require(ggpubr)
require(cowplot)
require(dplyr)
library(ggstance)
library(Rmisc)
library(dplyr)
library(magick)
library(png)
library(tidyr)
require(permutations)
require(survminer)
require(survival)
require(cmprsk)
library(forcats)
library(openxlsx)
library(ggnewscale)
library(ggrepel)

metadata <- as.data.frame(read.csv("Files/mapping_tmo_0.tsv", header = TRUE, sep = "\t"))
demux <- as.data.frame(read.csv("Files/demux-paired-end.csv", header = TRUE))
colnames(demux)<-c("sample.id","nseqs")
demux<-merge(metadata,demux)
demux_tmo <- subset(demux,time!="MI"&time!="MF") #tmo indica que mi e mf foram removidas
sum(demux_tmo$nseqs) #53253725
median(demux_tmo$nseqs) #104230.5
min(demux_tmo$nseqs) #2059
max(demux_tmo$nseqs) #502409
quantile(demux_tmo$nseqs, probs = seq(0,1,by=0.01))

FS2B<-ggplot(demux_tmo,aes(x=nseqs))+
  geom_vline(xintercept = median(demux_tmo$nseqs),color="grey40",linetype="dashed")+
  geom_histogram(color="#B17B4E",fill="grey97",breaks = seq(0,510000,by=10000),alpha=0.8)+
  geom_col(aes(x=nseqs,y=1),color="#B17B4E",fill="grey97",width = 1)+
  labs(x="N of raw reads",y="N of samples")+
  scale_x_continuous(limits=c(0,510000),breaks = seq(0,510000,by=50000),expand=c(0.007,0.007))+
  scale_y_continuous(limits=c(0,45),breaks=seq(0,45,by=5),expand=c(0.01,0.01))+
  theme(panel.grid = element_blank())

metadata <- as.data.frame(read.csv("Files/mapping_tmo_0.tsv", header = TRUE, sep = "\t"))
final <- as.data.frame(read.csv("Files/nseqs_final.csv", header = FALSE))
colnames(final)<-c("sample.id","nseqs")
final<-merge(metadata,final)
final_tmo <- subset(final,time!="MI"&time!="MF") #tmo indica que mi e mf foram removidas
sum(final_tmo$nseqs) #31343619
median(final_tmo$nseqs) #63075.5
min(final_tmo$nseqs) #87
max(final_tmo$nseqs) #310082
quantile(final_tmo$nseqs, probs = seq(0,1,by=0.01))

FS2E<-ggplot(final_tmo,aes(x=nseqs))+
  geom_vline(xintercept = median(final_tmo$nseqs),color="grey40",linetype="dashed")+
  geom_histogram(color="#4E84B1",fill="grey97",breaks = seq(0,510000,by=10000),alpha=0.8)+
  geom_col(aes(x=nseqs,y=1),color="#4E84B1",fill="grey97",width = 1)+
  labs(x="N of retained reads",y="N of samples")+
  scale_x_continuous(limits=c(0,510000),breaks = seq(0,510000,by=50000),expand=c(0.007,0.007))+
  scale_y_continuous(limits=c(0,45),breaks=seq(0,45,by=5),expand=c(0.01,0.01))+
  theme(panel.grid = element_blank())

#calculate n of asvs generated for whole dataset
physeq<-qza_to_phyloseq(
  features="Files/ASV_table_final_goodsamples.qza",
  tree="Files/rooted-tree.qza",
  "Files/taxonomy_vs.qza",
  metadata = "Files/mapping_tmo_0.tsv")

physeq #4212 ASVs considerando também as MI e MF
physeq_tmo<-subset_samples(physeq,time!="MI"&time!="MF")
prune_taxa(taxa_sums(physeq_tmo) > 0, physeq_tmo) #4046 ASVs sem considerar MI e MF

###PIZZA
dadaed <- as.data.frame(read.csv("Files/demux-dadaed.tsv", header = TRUE,sep="\t"))
#colnames(dadaed)<-c("sample.id","nseqs")
dadaed<-merge(metadata,dadaed)
dadaed_tmo <- subset(dadaed,time!="MI"&time!="MF")

chim <- as.data.frame(read.csv("Files/demux-dadaed-chimeraremoved.csv",header = FALSE))
colnames(chim)<-c("sample.id","nseqs")
chim<-merge(metadata,chim)
chim_tmo <- subset(chim,time!="MI"&time!="MF")

dadaed_chim_tmo<-merge(dadaed_tmo,chim_tmo)

final_pc=sum(final_tmo$nseqs)/sum(demux_tmo$nseqs) #58.85714% #AO FIM
q_filtered_pc=1-(sum(dadaed_chim_tmo$filtered)/sum(demux_tmo$nseqs)) #0.3551806 q-filtered
non_denoised_pc=1-(sum(dadaed_chim_tmo$denoised)/sum(demux_tmo$nseqs))-q_filtered_pc
non_merged_pc=1-(sum(dadaed_chim_tmo$merged)/sum(demux_tmo$nseqs))-q_filtered_pc
non_bimeric_pc=1-(sum(dadaed_chim_tmo$non.chimeric)/sum(demux_tmo$nseqs))-q_filtered_pc-non_merged_pc
non_chimeric_pc=1-(sum(dadaed_chim_tmo$nseqs)/sum(demux_tmo$nseqs))-q_filtered_pc-non_merged_pc-non_bimeric_pc
non_bacterial_pc=1-final_pc-q_filtered_pc-non_merged_pc-non_bimeric_pc-non_chimeric_pc

pie_df <- data.frame(
  group = c("fq_filtered_pc", "enon_merged_pc","dnon_bimeric_pc","cnon_chimeric_pc","bnon_bacterial_pc","afinal_pc"),
  value = c(q_filtered_pc, non_merged_pc,non_bimeric_pc,non_chimeric_pc,non_bacterial_pc,final_pc)
)
pie_df$pos = c(0.1775903,0.3650224,0.3792219,0.3915219,0.4054464,0.7057143)

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    axis.text.x=element_blank()
  )
FS2C<-ggplot(pie_df, aes(x="", y=value, fill=group))+
  geom_bar(width = 1, stat = "identity")+geom_label_repel(aes(x = 1.4, y = pos, label = paste0(as.character(round(value*100,1)),"%")), 
                                                        nudge_x = .3, 
                                                        segment.size = .7, 
                                                        show.legend = FALSE, color = "black")+
  coord_polar("y", start=0)+blank_theme+labs(title="Raw reads proportion",fill="")+
  scale_fill_manual(values=c("#709cc0","grey90","grey82","grey74","grey66","grey58"),
                    labels=c("afinal_pc"="Retained","bnon_bacterial_pc"="Non-bacterial","dnon_bimeric_pc"= "Bimeras",
                             "cnon_chimeric_pc" = "Chimeras","enon_merged_pc"="Not merged","fq_filtered_pc"="Low-quality"),
                    guide = guide_legend(reverse = TRUE))

median_df<-data.frame(
  group = c("step0", "step1","step2","step3","step4","step5","step6"),
  value = c(median(demux_tmo$nseqs), median(dadaed_chim_tmo$filtered),
            median(dadaed_chim_tmo$denoised),median(dadaed_chim_tmo$merged),
            median(dadaed_chim_tmo$non.chimeric),median(dadaed_chim_tmo$nseqs),
            median(final_tmo$nseqs))
)
options(scipen=10000)
FS2D<-ggplot(median_df,aes(x=group,y=value))+geom_col(aes(fill=group),width=0.7)+labs(x="Pipeline step",y="Median N of reads per sample")+
  scale_x_discrete(labels=c("step0"= "Raw reads","step1"="After step 1","step2"="After step 2","step3"="After step 3",
                            "step4"="After step 4","step5"="After step 5","step6" = "After step 6 \n (retained) "))+
  geom_text(aes(y=value+3000,label=value),size=3.1)+
  scale_y_continuous(breaks=seq(0,100000,by=20000),expand=c(0.01,0.01),limits = c(0,110000))+
  scale_fill_manual(values=c("#B17B4E","grey58","grey66","grey74","grey82","grey90","#4E84B1"))+
  guides(fill=FALSE)+theme(panel.grid = element_blank(),axis.text.x=element_text(angle=20,hjust=0.5,vjust =0.8))