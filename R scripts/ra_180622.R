library(qiime2R)
require(ggplot2)
library(btools)
library(vegan)
library(scales)
library(phyloseq)
library(gplots)
library(dplyr)
library(QsRutils)
library(reshape2)
library(ggrepel)
library(ggn)
library(ggbeeswarm)
library(ggnewscale)
library(ggstance)
library(ggpubr)
library(cowplot)
library(pals)
library(RColorBrewer)
library(pairwiseAdonis)
library(mctoolsr)
library(metagMisc)
library(SRS)
source("plotDistances.R")
theme_set(theme_bw())
pacientes<-c("AAM"="#1","AFLF"="#2","AMSCP"="#3","BMF"="#4","CEFC"="#5","DBHB"="#6",
             "ECPL"="#7","ENS"="#8","FGTC"="#9","FM"="#10","IABC"="#11","ILSM"="#12",
             "JPM"="#13","JVS"="#14","LEMG"="#15","LFCS"="#16","MEFD"="#17","MFML"="#18",
             "MMMC"="#19","MSBS"="#20","OBR"="#21","OSLG"="#22","OSN"="#23","PHOS"="#24",
             "RJT"="#25","RLS"="#26","SMG"="#27","TPN"="#28","VFF"="#29","VGS"="#30","VRM"="#31")
labeller <- function(variable,value){
  return(pacientes[value])
}

coherent_taxa_names<-function(physeq){
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "d__", replacement = "")
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "uncultured", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "metagenome", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "human_gut", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = "_sp.", replacement = NA)
  tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = " sp.", replacement = NA)
  
  tax.clean <- data.frame(tax_table(physeq))
  tax.clean[is.na(tax.clean)] <- ""
  for (i in 1:nrow(tax.clean)){
    if (tax.clean[i,2] == ""){
      kingdom <- paste("Kingdom_", tax.clean[i,1], sep = "")
      tax.clean[i, 2:7] <- kingdom
    } else if (tax.clean[i,3] == ""){
      phylum <- paste("Phylum_", tax.clean[i,2], sep = "")
      tax.clean[i, 3:7] <- phylum
    } else if (tax.clean[i,4] == ""){
      class <- paste0("Class_", tax.clean[i,3], sep = "")
      tax.clean[i, 4:7] <- class
    } else if (tax.clean[i,5] == ""){
      order <- paste("Order_", tax.clean[i,4], sep = "")
      tax.clean[i, 5:7] <- order
    } else if (tax.clean[i,6] == ""){
      family <- paste("Family_", tax.clean[i,5], sep = "")
      tax.clean[i, 6:7] <- family
    } else if (tax.clean[i,7] == ""){
      tax.clean$Species[i] <- paste("Genus",tax.clean$Genus[i], sep = "_")
    }
  }
  tax_table(physeq) <- as.matrix(tax.clean)
  #tax_table(physeq)[, colnames(tax_table(physeq))] <- gsub(tax_table(physeq)[, colnames(tax_table(physeq))],     pattern = " ", replacement = "")
  return(physeq)
}

coherent_asv_names<-function(physeq){
  asvs_list<-paste0("ASV",1:length(as.data.frame(physeq@tax_table)[,1]),"_",
                    as.data.frame(physeq@tax_table)['Species'][1:length(as.data.frame(physeq@tax_table)[,1]),])
  
  new_tax_table<-as.data.frame(physeq@tax_table)
  new_tax_table$ASV<-asvs_list
  tax_table(physeq)<-as.matrix(new_tax_table)
  return(physeq)
  }

#loading physeq object w/ full dataset #490 samples (bad samples (<3k) already discarded)
physeq<-coherent_taxa_names(qza_to_phyloseq(
  features="Files/ASV_table_final_goodsamples.qza",
  tree="Files/rooted-tree.qza",
  "Files/taxonomy_vs.qza",
  metadata = "Files/mapping_tmo_0.tsv"))

#first lets work with tmo samples only (phase-wise)
physeq_tmo<-subset_samples(physeq,time!="MI"&time!="MF")

physeq_ra_tmo<-transform_sample_counts(physeq_tmo, function(x) x / sum(x) ) #make relative abundances

# #first lets work with tmo samples only (phase-wise)
# physeq_SRS_tmo_mo<-subset_samples(physeq_SRS_tmo,site=="MO")
# physeq_SRS_tmo_bf<-subset_samples(physeq_SRS_tmo,site=="BF")
# physeq_SRS_tmo_fcg<-subset_samples(physeq_SRS_tmo,site=="FCG")

ra_bump<-function(physeq_ra1, physeq_ra2, physeq_ra3, taxrank = "Genus", min_ra=0.01, min_n_samples = 1, xlab = "Fase"){
  physeq_ra1<-phyloseq_average(physeq_ra1, group = "time",avg_type = "arithmetic")
  physeq_ra1<-tax_glom(physeq_ra1, taxrank = taxrank)
  physeq_ra1<-prune_taxa(genefilter_sample(physeq_ra1, filterfun_sample(function(x) x >= min_ra), 
                                           A = min_n_samples),physeq_ra1)
  ranking1<-as.data.frame(otu_table(physeq_ra1))
  ranking1<-as.data.frame(sapply(1-ranking1, rank))
  ranking1<-cbind(as.data.frame(physeq_ra1@tax_table@.Data)[taxrank],ranking1)
  ranking1.melt<-reshape2::melt(ranking1)
  colnames(ranking1.melt)<-c("taxa","time","ranking")
  ranking1.melt$site<-'SB'
  
  physeq_ra2<-phyloseq_average(physeq_ra2, group = "time",avg_type = "arithmetic")
  physeq_ra2<-tax_glom(physeq_ra2, taxrank = taxrank)
  physeq_ra2<-prune_taxa(genefilter_sample(physeq_ra2, filterfun_sample(function(x) x >= min_ra), 
                                           A = min_n_samples),physeq_ra2)
  ranking2<-as.data.frame(otu_table(physeq_ra2))
  ranking2<-as.data.frame(sapply(1-ranking2, rank))
  ranking2<-cbind(as.data.frame(physeq_ra2@tax_table@.Data)[taxrank],ranking2)
  ranking2.melt<-reshape2::melt(ranking2)
  colnames(ranking2.melt)<-c("taxa","time","ranking")
  ranking2.melt$site<-'GCF'
  
  physeq_ra3<-phyloseq_average(physeq_ra3, group = "time",avg_type = "arithmetic")
  physeq_ra3<-tax_glom(physeq_ra3, taxrank = taxrank)
  physeq_ra3<-prune_taxa(genefilter_sample(physeq_ra3, filterfun_sample(function(x) x >= min_ra), 
                                           A = min_n_samples),physeq_ra3)
  ranking3<-as.data.frame(otu_table(physeq_ra3))
  ranking3<-as.data.frame(sapply(1-ranking3, rank))
  ranking3<-cbind(as.data.frame(physeq_ra3@tax_table@.Data)[taxrank],ranking3)
  ranking3.melt<-reshape2::melt(ranking3)
  colnames(ranking3.melt)<-c("taxa","time","ranking")
  ranking3.melt$site<-'OM'
  
  ranking.melt<-rbind(ranking1.melt,ranking2.melt,ranking3.melt)
  ranking.melt$time<-factor(ranking.melt$time, levels=c("P","A","E","E30","E75"))
  
  ggplot(ranking.melt,aes(x=time,y=ranking,group=taxa,color=taxa))+
    geom_line(size=2)+geom_point(size=3,pch=21,fill="white")+
    #geom_label_repel(data=subset(ranking.melt,ranking<=10)[order(subset(ranking.melt,ranking<=10)[,'taxa'],subset(ranking.melt,ranking<=10)[,'time']),][!duplicated(subset(ranking.melt,ranking<=10)[order(subset(ranking.melt,ranking<=10)[,'taxa'],subset(ranking.melt,ranking<=10)[,'time']),]$taxa),],
    #                   aes(label=taxa,y=ranking,x=time),color="black",size=3,nudge_x = -0.5,label.padding = 0,label.size = 0,alpha=0.6)+
    scale_y_reverse(breaks=seq(1,10,1), labels = c("1"=' 1st',"2"=' 2nd',"3"=' 3rd',"4"=' 4th',"5"=' 5th',"6"=' 6th',"7"=' 7th',"8"=' 8th',"9"=' 9th',"10"=' 10th'))+coord_cartesian(ylim=c(10, 1))+
    scale_color_manual(values=c(glasbey(n=32)))+
    labs(x="Timepoint",y="Mean RA ranking")+
    theme(panel.grid = element_blank(), axis.title.y = element_text(size = 9),
          legend.title = element_text(size=14), legend.position = 'none',
          strip.background = element_rect(color="white", fill="white", linetype="solid"))+guides(color=F)+
    scale_x_discrete(expand = c(0.05,0.05))+facet_wrap(.~site, nrow = 1)
}

ra_persample_avg<-function(physeq_ra1,physeq_ra2,physeq_ra3, taxrank = "Genus", min_ra=0.01, min_n_samples = 1, xlab = "Phase"){
  physeq_ra1<-phyloseq_average(physeq_ra1, group = "time",avg_type = "arithmetic")
  physeq_ra1<-tax_glom(physeq_ra1, taxrank = taxrank)
  physeq_ra1<-prune_taxa(genefilter_sample(physeq_ra1, filterfun_sample(function(x) x >= min_ra), 
                                           A = min_n_samples),physeq_ra1)
  data1<-as.data.frame(otu_table(physeq_ra1))
  data1$taxa<-as.data.frame(physeq_ra1@tax_table@.Data)[[taxrank]]
  data1<-rbind(data1,c(1-colSums(data1[,1:(length(data1)-1)]),'Other'))
  data1.melt<-reshape2::melt(data1, c("taxa"))
  data1.melt$value<-as.numeric(data1.melt$value)
  data1.melt$site<-"SB"
  
  physeq_ra2<-phyloseq_average(physeq_ra2, group = "time",avg_type = "arithmetic")
  physeq_ra2<-tax_glom(physeq_ra2, taxrank = taxrank)
  physeq_ra2<-prune_taxa(genefilter_sample(physeq_ra2, filterfun_sample(function(x) x >= min_ra), 
                                           A = min_n_samples),physeq_ra2)
  data2<-as.data.frame(otu_table(physeq_ra2))
  data2$taxa<-as.data.frame(physeq_ra2@tax_table@.Data)[[taxrank]]
  data2<-rbind(data2,c(1-colSums(data2[,1:(length(data2)-1)]),'Other'))
  data2.melt<-reshape2::melt(data2, c("taxa"))
  data2.melt$value<-as.numeric(data2.melt$value)
  data2.melt$site<-"GCF"
  
  physeq_ra3<-phyloseq_average(physeq_ra3, group = "time",avg_type = "arithmetic")
  physeq_ra3<-tax_glom(physeq_ra3, taxrank = taxrank)
  physeq_ra3<-prune_taxa(genefilter_sample(physeq_ra3, filterfun_sample(function(x) x >= min_ra), 
                                           A = min_n_samples),physeq_ra3)
  data3<-as.data.frame(otu_table(physeq_ra3))
  data3$taxa<-as.data.frame(physeq_ra3@tax_table@.Data)[[taxrank]]
  data3<-rbind(data3,c(1-colSums(data3[,1:(length(data3)-1)]),'Other'))
  data3.melt<-reshape2::melt(data3, c("taxa"))
  data3.melt$value<-as.numeric(data3.melt$value)
  data3.melt$site<-"OM"
  
  data.melt<-rbind(data1.melt,data2.melt,data3.melt)
  ng<-length(unique(data.melt$taxa))-1
  
  print(ng)
  
  data.melt$taxa<-factor(data.melt$taxa, levels = c(sort(unique(subset(data.melt,taxa!='Other')$taxa)),'Other'))
  data.melt$variable<-factor(data.melt$variable, levels = c('P','A','E','E30','E75'))
  ggplot(data.melt, aes(x=variable, y=value)) + 
    geom_bar(data=data.melt,
             aes(x=variable, y=value, fill=taxa), stat="identity", color = "black", size = 0.1)+
    scale_fill_manual(values=c(glasbey(n=ng),'gray50'), guide = guide_legend(ncol=1,bycol=TRUE))+
    theme(panel.grid = element_blank(), legend.position = "right",
          legend.text = element_text(size = 10), legend.key.size = unit(0.85, 'lines'),
          legend.title = element_text(size=14),
          strip.background = element_rect(color="white", fill="white", linetype="solid"))+
    labs(x=xlab, y = 'Mean RA (%)', fill = taxrank)+
    scale_y_continuous(expand = c(0,0), limits = c(0,1.0001))+
    scale_x_discrete(expand = c(0,0))+
    facet_wrap(.~site, nrow = 1)
  
}

F3A<-ra_persample_avg(physeq_ra1=subset_samples(physeq_ra_tmo,site=="BF"),
                 physeq_ra2=subset_samples(physeq_ra_tmo,site=="FCG"),
                 physeq_ra3=subset_samples(physeq_ra_tmo,site=="MO"),min_ra = 0.02)+
  guides(fill=guide_legend(nrow=6,byrow = T))+
  theme(legend.position = "top",legend.text = element_text(size=8.5), 
        legend.key.size = unit(0.4, 'cm'),axis.title.x = element_blank())

F3B<-ra_bump(physeq_ra1=subset_samples(physeq_ra_tmo,site=="BF"),
        physeq_ra2=subset_samples(physeq_ra_tmo,site=="FCG"),
        physeq_ra3=subset_samples(physeq_ra_tmo,site=="MO"), min_ra = 0.02)

FS7A<-ra_persample_avg(physeq_ra1=subset_samples(physeq_ra_tmo,site=="BF"),
                   physeq_ra2=subset_samples(physeq_ra_tmo,site=="FCG"),
                   physeq_ra3=subset_samples(physeq_ra_tmo,site=="MO"),min_ra = 0.02, taxrank = "Phylum")+
    guides(fill=guide_legend(nrow=2,byrow = T))+
    theme(legend.position = "top",legend.text = element_text(size=8.5), 
          legend.key.size = unit(0.4, 'cm'),axis.title.x = element_blank())

FS7B<-ra_bump(physeq_ra1=subset_samples(physeq_ra_tmo,site=="BF"),
          physeq_ra2=subset_samples(physeq_ra_tmo,site=="FCG"),
          physeq_ra3=subset_samples(physeq_ra_tmo,site=="MO"), min_ra = 0.02, taxrank = "Phylum")

FS7C<-ra_persample_avg(physeq_ra1=subset_samples(physeq_ra_tmo,site=="BF"),
                   physeq_ra2=subset_samples(physeq_ra_tmo,site=="FCG"),
                   physeq_ra3=subset_samples(physeq_ra_tmo,site=="MO"),min_ra = 0.02, taxrank = "Class")+
    guides(fill=guide_legend(nrow=3,byrow = T))+
    theme(legend.position = "top",legend.text = element_text(size=8.5), 
          legend.key.size = unit(0.4, 'cm'),axis.title.x = element_blank())

FS7D<-ra_bump(physeq_ra1=subset_samples(physeq_ra_tmo,site=="BF"),
          physeq_ra2=subset_samples(physeq_ra_tmo,site=="FCG"),
          physeq_ra3=subset_samples(physeq_ra_tmo,site=="MO"), min_ra = 0.02, taxrank = "Class")

FS7E<-ra_persample_avg(physeq_ra1=subset_samples(physeq_ra_tmo,site=="BF"),
                   physeq_ra2=subset_samples(physeq_ra_tmo,site=="FCG"),
                   physeq_ra3=subset_samples(physeq_ra_tmo,site=="MO"),min_ra = 0.02, taxrank = "Order")+
    guides(fill=guide_legend(nrow=4,byrow = T))+
    theme(legend.position = "top",legend.text = element_text(size=8.5), 
          legend.key.size = unit(0.4, 'cm'),axis.title.x = element_blank())

FS7F<-ra_bump(physeq_ra1=subset_samples(physeq_ra_tmo,site=="BF"),
          physeq_ra2=subset_samples(physeq_ra_tmo,site=="FCG"),
          physeq_ra3=subset_samples(physeq_ra_tmo,site=="MO"), min_ra = 0.02, taxrank = "Order")

FS7G<-ra_persample_avg(physeq_ra1=subset_samples(physeq_ra_tmo,site=="BF"),
                   physeq_ra2=subset_samples(physeq_ra_tmo,site=="FCG"),
                   physeq_ra3=subset_samples(physeq_ra_tmo,site=="MO"),min_ra = 0.02, taxrank = "Family")+
    guides(fill=guide_legend(nrow=7,byrow = T))+
    theme(legend.position = "top",legend.text = element_text(size=8.5), 
          legend.key.size = unit(0.4, 'cm'),axis.title.x = element_blank())

FS7H<-ra_bump(physeq_ra1=subset_samples(physeq_ra_tmo,site=="BF"),
          physeq_ra2=subset_samples(physeq_ra_tmo,site=="FCG"),
          physeq_ra3=subset_samples(physeq_ra_tmo,site=="MO"), min_ra = 0.02, taxrank = "Family")
