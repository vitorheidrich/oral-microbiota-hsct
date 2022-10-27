library(nloptr)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggbeeswarm)
theme_set(theme_bw())
library(stringr)
library(phyloseq)
library(qiime2R)
library(ggpubr)
library(tidytext)
library(pals)
site_colors<-c('#5B4BB7', '#B79A4C', '#B75B4B')

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
#################################################
physeq<-coherent_taxa_names(qza_to_phyloseq(
  features="Files/ASV_table_final_goodsamples.qza",
  tree="Files/rooted-tree.qza",
  "Files/taxonomy_vs.qza",
  metadata = "Files/mapping_tmo_0.tsv"))

metadata_filtered<-subset(as.data.frame(sample_data(physeq)), time!= "MI" & time != "MF")

physeq_tmo<-subset_samples(physeq,time!="MI"&time!="MF")
physeq_ra_tmo<-transform_sample_counts(physeq_tmo, function(x) x / sum(x) )

physeq_ra_tmo_g<-tax_glom(physeq_ra_tmo,taxrank = "Genus")

df_blooms<-data.frame(patient = "a", site = "a", taxa = "a", time= "a")

for (sitio in c("BF","MO","FCG")){
  for(paciente in unique(as(sample_data(physeq_ra_tmo_g),"data.frame")$patient)){
    for(tempo in c("A","E","E30","E75")){
      test_presence_df_tm<-subset(metadata_filtered,site==sitio&time==tempo)$patient
      test_presence_df_P<-subset(metadata_filtered,site==sitio&time=="P")$patient
      if (paciente %in% test_presence_df_tm & paciente %in% test_presence_df_P){
        #print(paste(sitio,tempo,paciente))
        sampleid_tm<-rownames(subset(metadata_filtered,site==sitio&time==tempo&patient==paciente))
        sampleid_P<-rownames(subset(metadata_filtered,site==sitio&time=="P"&patient==paciente))
        sampleid_P<-str_replace(sampleid_P,pattern="sa",replacement="")
        #print(paste(sampleid_tm,sampleid_P))
        blooming_taxa<-genefilter_sample(prune_samples(sampleid_tm,physeq_ra_tmo_g), 
                                     filterfun_sample(function(x) x >= 0.3), A=1)&
        genefilter_sample(prune_samples(sampleid_P,physeq_ra_tmo_g), 
                        filterfun_sample(function(x) x < 0.01), A=1)
        blooming_taxa_df<-as.data.frame(blooming_taxa)
        if(TRUE %in% blooming_taxa_df$blooming_taxa){
          blooming_taxa_code<-rownames(subset(blooming_taxa_df,blooming_taxa==T))
          blooming_taxa_name<-as.data.frame(tax_table(physeq_ra_tmo_g))[blooming_taxa_code,]$Genus
          for (i in 1:length(blooming_taxa_name)){
            if(any(df_blooms$patient == paciente & df_blooms$site == sitio & df_blooms$taxa == blooming_taxa_name[i])==F){
              print(blooming_taxa_name)
              df_blooms<-rbind(df_blooms,c(paciente,sitio,blooming_taxa_name[i],tempo))
            }
          }
        }
      }
    }
  }
} 
df_blooms<-df_blooms[-1,]

length(unique(df_blooms$patient)) #27 patients w bloom overall
length(unique(df_blooms$taxa)) #22 genera w bloom overall
length(unique(subset(df_blooms, site=="BF")$patient)) #23 patients w bloom BF
length(unique(subset(df_blooms, site=="FCG")$patient)) #14 patients w bloom FCG
length(unique(subset(df_blooms, site=="MO")$patient)) #16 patients w bloom MO
length(unique(subset(df_blooms, site=="BF")$taxa)) #16 genera w bloom BF
length(unique(subset(df_blooms, site=="FCG")$taxa)) #14 genera w bloom FCG
length(unique(subset(df_blooms, site=="MO")$taxa)) #16 genera w bloom MO
sort(table(df_blooms$taxa)) #most frequent overall: Lactobacillus (12), Enterococcus (10), Staphylococcus (8)
sort(table(subset(df_blooms, site=="BF")$taxa)) #most frequent BF: Enterococcus (7), Lactobacillus (6)
sort(table(subset(df_blooms, site=="FCG")$taxa)) #most frequent FCG: Lactobacillus (4), Staphylococcus (4)
sort(table(subset(df_blooms, site=="MO")$taxa)) #most frequent MO: Prevotella (3)
sort(table(subset(df_blooms, site=="BF")$patient)) #n per patient BF: 3x3, 6x2, 14x1 = 35
sort(table(subset(df_blooms, site=="FCG")$patient)) #n per patient FCG: 2x3, 6x2, 6x1 = 24
sort(table(subset(df_blooms, site=="MO")$patient)) #n per patient MO: 1x3, 4x2, 11x1 = 22 (TOTAL = 81)
sort(table(df_blooms$time)) #most frequent times overall: E (43, 53%), E30 (21, 26%), A (12, 15%), E75 (5, 6%)
sort(table(subset(df_blooms, site=="BF")$time)) #most frequent times BF: E (22, 63%), E30 (8, 23%), A (4, 11%), E75 (1, 3%)
sort(table(subset(df_blooms, site=="FCG")$time)) #most frequent times FCG: E (12, 50%), E30 (5, 21%), A (3, 12%), E75 (4, 17%)
sort(table(subset(df_blooms, site=="MO")$time)) #most frequent times MO: E (9, 41%), E30 (8, 36%), A (5, 23%), E75 (0, 0%)
df_blooms %>% #concomitant bloom of the same genus in all sites
  group_by(patient, taxa, time) %>% #JPM Abiotrophia E30
  filter(n() >= 3)                  #AAM Stenotrophomonas E
                                    #SMG Lactobacillus E

dacol<-read.csv("Files/Coletas_TMO_dcast.csv")
dacol<-dacol[c("patient","dacolA","dacolE","dacolE30","dacolE75")]
colnames(dacol)<-c("patient","A","E","E30","E75")
dacol.melt<-reshape2::melt(dacol,c("patient"))
colnames(dacol.melt)<-c("patient","time","day")

write.csv(merge(df_blooms,dacol.melt),file = "Files/260522_blooms_list.csv")

# write.csv(data.frame(bloom_BF = sort(unique(as(sample_data(physeq_ra_tmo_g),"data.frame")$patient)) %in% 
#   sort(unique(subset(df_blooms, site=="BF")$patient)),
#   bloom_FCG = sort(unique(as(sample_data(physeq_ra_tmo_g),"data.frame")$patient)) %in% 
#     sort(unique(subset(df_blooms, site=="FCG")$patient)),
#   bloom_MO = sort(unique(as(sample_data(physeq_ra_tmo_g),"data.frame")$patient)) %in% 
#     sort(unique(subset(df_blooms, site=="MO")$patient))),file = "moch.csv")

bacteremia_bloom<-read.csv("Files/Bacteremia-vs-Bloom_270522.csv")
fisher.test(table(bacteremia_bloom[c("bacteremia_durantecoletas","bloom_BF")]))
fisher.test(table(bacteremia_bloom[c("bacteremia_durantecoletas","bloom_FCG")]))
fisher.test(table(bacteremia_bloom[c("bacteremia_durantecoletas","bloom_MO")]))
fisher.test(table(bacteremia_bloom[c("bacteremia_durantecoletas","bloom")]))
fisher.test(table(bacteremia_bloom[c("bacteremia_durantecoletas","bloom_all")]))

fisher.test(table(bacteremia_bloom[c("infresp_durantecoletas","bloom_BF")]))
fisher.test(table(bacteremia_bloom[c("infresp_durantecoletas","bloom_FCG")]))
fisher.test(table(bacteremia_bloom[c("infresp_durantecoletas","bloom_MO")]))
fisher.test(table(bacteremia_bloom[c("infresp_durantecoletas","bloom")]))
fisher.test(table(bacteremia_bloom[c("infresp_durantecoletas","bloom_all")]))

#############################################################################
df_blooms<-read.csv("Files/260522_blooms_list.csv")
df_blooms %>% group_by(site, taxa) %>% count() 
df_blooms[df_blooms=="BF"]<-"SB"
df_blooms[df_blooms=="FCG"]<-"GCF"
df_blooms[df_blooms=="MO"]<-"OM"

F4D<-ggplot(data=df_blooms %>% group_by(site, taxa) %>% count(),
       aes(x=n,y=reorder_within(taxa, n, site),fill=site))+
  geom_bar(stat = "identity", position = "stack")+
  labs(x = "N of events")+guides(fill="none")+
  scale_x_continuous(breaks=c(0,2,4,6,8),limits = c(0,8))+scale_y_discrete(labels = function(x) sub('__.*', '', x))+
  facet_wrap(~site, scales = "free_y")+
  theme(panel.grid = element_blank(), panel.border = element_blank(),axis.line = element_line(color = 'black'),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        axis.title.y = element_blank())+
  scale_fill_manual(values=site_colors)

F4D

F4A<-ggplot(data=df_blooms %>% group_by(site) %>% count() %>% as.data.frame() %>% mutate(prop = prop.table(n)),
       aes(y=prop,x=reorder(site, -n)))+
  geom_bar(stat = "identity", position = "stack", fill = 'black')+
  labs(y = "% of events", x = 'Site')+guides(fill="none")+
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5),limits = c(0,0.5), expand = c(0,0))+
  scale_x_discrete(labels = function(x) sub('__.*', '', x))+
  facet_wrap(~'Site')+
  theme(panel.grid = element_blank(), panel.border = element_blank(),axis.line = element_line(color = 'black'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title.x = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"))
F4A

F4B<-ggplot(data=df_blooms %>% group_by(time) %>% count() %>% as.data.frame() %>% mutate(prop = prop.table(n)),
            aes(y=prop,x=reorder(time, -n)))+
  geom_bar(stat = "identity", position = "stack", fill = 'black')+
  labs(y = "% of events", x = 'Timepoint')+guides(fill="none")+
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5),limits = c(0,0.55), expand = c(0,0))+
  scale_x_discrete(labels = function(x) sub('__.*', '', x))+
  facet_wrap(~'Timepoint')+
  theme(panel.grid = element_blank(), panel.border = element_blank(),axis.line = element_line(color = 'black'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"))
F4B

F4C<-ggplot(data=df_blooms %>% group_by(taxa) %>% count() %>% as.data.frame() %>% mutate(prop = prop.table(n)),
            aes(y=prop,x=reorder(taxa, -n)))+
  geom_bar(stat = "identity", position = "stack", fill = 'black')+
  labs(y = "% of events", x = 'Genus')+guides(fill="none")+
  scale_y_continuous(breaks=c(0,0.05, 0.1, 0.15, 0.2),limits = c(0,0.155), expand = c(0,0))+
  scale_x_discrete(labels = function(x) sub('__.*', '', x))+
  facet_wrap(~'Genus')+
  theme(panel.grid = element_blank(), panel.border = element_blank(),axis.line = element_line(color = 'black'),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), axis.title = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"))
F4C

#############################################################################
df_blooms<-read.csv("Files/260522_blooms_list.csv")
df_blooms %>% group_by(site, taxa) %>% count() 
df_blooms[df_blooms=="BF"]<-"DB"
df_blooms[df_blooms=="FCG"]<-"GCF"
df_blooms[df_blooms=="MO"]<-"OM"

df_blooms_woe75<-subset(df_blooms,time!="E75") #removendo blooms mais tardios

bloom_class_db_any<-as.data.frame(table(unique(subset(df_blooms_woe75,site=="DB")$patient))) %>% rename('patient'='Var1',"DB_bloom"='Freq') %>% mutate(across(everything(), as.character))
bloom_class_gcf_any<-as.data.frame(table(unique(subset(df_blooms_woe75,site=="GCF")$patient))) %>% rename('patient'='Var1',"GCF_bloom"='Freq') %>% mutate(across(everything(), as.character))
bloom_class_om_any<-as.data.frame(table(unique(subset(df_blooms_woe75,site=="OM")$patient))) %>% rename('patient'='Var1',"OM_bloom"='Freq') %>% mutate(across(everything(), as.character))
bloom_class_db_lactob<-as.data.frame(table(unique(subset(df_blooms_woe75,site=="DB"&taxa=="Lactobacillus")$patient))) %>% rename('patient'='Var1',"DB_bloom_lactob"='Freq') %>% mutate(across(everything(), as.character))
bloom_class_db_entero<-as.data.frame(table(unique(subset(df_blooms_woe75,site=="DB"&taxa=="Enterococcus")$patient))) %>% rename('patient'='Var1',"DB_bloom_entero"='Freq') %>% mutate(across(everything(), as.character))

bloom_classification<-full_join(full_join(full_join(full_join(bloom_class_db_any,bloom_class_gcf_any,by = "patient"), 
                                     bloom_class_om_any,by = "patient"), bloom_class_db_lactob,by = "patient"), bloom_class_db_entero,by = "patient")%>% replace(is.na(.), 0)
bloom_classification[1,4]<-NA #pt AAM não tem amostra OM P, entao nao da pra saber se tem bloom
#colar pacientes que não tem bloom
setdiff(unique(read.csv("Files/mapping_tmo_0.tsv",sep="\t")$patient), bloom_classification$patient)
#"AMSCP" "BMF"   "IABC"  "OBR"   "PHOS"  "RLS"
bloom_classification<-rbind(bloom_classification,c("AMSCP",0,0,0,0,0))
bloom_classification<-rbind(bloom_classification,c("BMF",0,0,0,0,0))
bloom_classification<-rbind(bloom_classification,c("IABC",0,0,0,0,0))
bloom_classification<-rbind(bloom_classification,c("OBR",0,0,0,0,0))
bloom_classification<-rbind(bloom_classification,c("PHOS",0,0,0,0,0))
bloom_classification<-rbind(bloom_classification,c("RLS",0,0,0,0,0))

#load atb classification
atb_classification<-read.csv('Files/300622_ATB_classification.csv')

bloom_atb<-merge(bloom_classification,atb_classification)

#MANN-WHITNEY test against lot and dot
wilcox.test(subset(bloom_atb,DB_bloom==1)$lot,subset(bloom_atb,DB_bloom==0)$lot)#p = 0.03649
median(subset(bloom_atb,DB_bloom==1)$lot, na.rm = T) #17.5
median(subset(bloom_atb,DB_bloom==0)$lot, na.rm = T) #11

wilcox.test(subset(bloom_atb,DB_bloom==1)$dot,subset(bloom_atb,DB_bloom==0)$dot)#p = 0.02718
median(subset(bloom_atb,DB_bloom==1)$dot, na.rm = T) #28
median(subset(bloom_atb,DB_bloom==0)$dot, na.rm = T) #15

wilcox.test(subset(bloom_atb,DB_bloom_entero==1)$lot,subset(bloom_atb,DB_bloom_entero==0)$lot)#p = 0.1472
median(subset(bloom_atb,DB_bloom_entero==1)$lot, na.rm = T) #23
median(subset(bloom_atb,DB_bloom_entero==0)$lot, na.rm = T) #14

wilcox.test(subset(bloom_atb,DB_bloom_entero==1)$dot,subset(bloom_atb,DB_bloom_entero==0)$dot)#p = 0.2014
median(subset(bloom_atb,DB_bloom_entero==1)$dot, na.rm = T) #28
median(subset(bloom_atb,DB_bloom_entero==0)$dot, na.rm = T) #21

wilcox.test(subset(bloom_atb,DB_bloom_lactob==1)$lot,subset(bloom_atb,DB_bloom_lactob==0)$lot)#p = 0.1017
median(subset(bloom_atb,DB_bloom_lactob==1)$lot, na.rm = T) #23
median(subset(bloom_atb,DB_bloom_lactob==0)$lot, na.rm = T) #14

wilcox.test(subset(bloom_atb,DB_bloom_lactob==1)$dot,subset(bloom_atb,DB_bloom_lactob==0)$dot)#p = 0.07305
median(subset(bloom_atb,DB_bloom_lactob==1)$dot, na.rm = T)#31
median(subset(bloom_atb,DB_bloom_lactob==0)$dot, na.rm = T)#20

wilcox.test(subset(bloom_atb,GCF_bloom==1)$lot,subset(bloom_atb,GCF_bloom==0)$lot)#p = 0.0156
median(subset(bloom_atb,GCF_bloom==1)$lot, na.rm = T) #24.5
median(subset(bloom_atb,GCF_bloom==0)$lot, na.rm = T) #12.5

wilcox.test(subset(bloom_atb,GCF_bloom==1)$dot,subset(bloom_atb,GCF_bloom==0)$dot)#p = 0.004478
median(subset(bloom_atb,GCF_bloom==1)$dot, na.rm = T) #33
median(subset(bloom_atb,GCF_bloom==0)$dot, na.rm = T) #12.5

wilcox.test(subset(bloom_atb,OM_bloom==1)$lot,subset(bloom_atb,OM_bloom==0)$lot) #p = 0.02227
median(subset(bloom_atb,OM_bloom==1)$lot, na.rm = T) #17.5
median(subset(bloom_atb,OM_bloom==0)$lot, na.rm = T) #10

wilcox.test(subset(bloom_atb,OM_bloom==1)$dot,subset(bloom_atb,OM_bloom==0)$dot) #p = 0.01869
median(subset(bloom_atb,OM_bloom==1)$dot, na.rm = T) #28
median(subset(bloom_atb,OM_bloom==0)$dot, na.rm = T) #14

bloom_atb.melt<-melt(bloom_atb,c('patient','dot','lot','carbapenems','glycopeptides','penicillins','cephalosporins',
                                 'carbapenems_cat','glycopeptides_cat','penicillins_cat','cephalosporins_cat'))
colnames(bloom_atb.melt)<-c('patient','dot','lot','carbapenems','glycopeptides','penicillins','cephalosporins',
                            'carbapenems_cat','glycopeptides_cat','penicillins_cat','cephalosporins_cat','bloom_site','bloom_value')
bloom_atb.melt.melt<-melt(bloom_atb.melt,c('patient','carbapenems','glycopeptides','penicillins','cephalosporins',
                                           'carbapenems_cat','glycopeptides_cat','penicillins_cat','cephalosporins_cat','bloom_site','bloom_value'))
bloom_atb.melt.melt$bloom_site<-as.character(bloom_atb.melt.melt$bloom_site)
bloom_atb.melt.melt$variable<-as.character(bloom_atb.melt.melt$variable)
bloom_atb.melt.melt$value<-as.numeric(bloom_atb.melt.melt$value)

bloom_atb.melt.melt[bloom_atb.melt.melt=="DB_bloom"]<-"SB"
bloom_atb.melt.melt[bloom_atb.melt.melt=="GCF_bloom"]<-"GCF"
bloom_atb.melt.melt[bloom_atb.melt.melt=="OM_bloom"]<-"OM"
bloom_atb.melt.melt[bloom_atb.melt.melt=="dot"]<-"DOT"
bloom_atb.melt.melt[bloom_atb.melt.melt=="lot"]<-"LOT"
bloom_atb.melt.melt$bloom_value[bloom_atb.melt.melt$bloom_value=="0"]<-"No"
bloom_atb.melt.melt$bloom_value[bloom_atb.melt.melt$bloom_value=="1"]<-"Yes"

bloom_atb.melt.melt$variable<-factor(bloom_atb.melt.melt$variable,levels=c('LOT','DOT'))

F5B<-ggplot(data = subset(bloom_atb.melt.melt,bloom_site!='DB_bloom_entero'&bloom_site!='DB_bloom_lactob'&is.na(bloom_value)==F)) + 
  aes(x = bloom_value, y = value, color = bloom_site, group = bloom_value) + 
  geom_boxplot(alpha=0.5, outlier.size = 1)+
  facet_grid(variable ~ bloom_site,scales = "free_y") + 
  labs(y = "Antibiotic metric", x = 'Bloom') + guides(color=FALSE)+
  scale_color_manual(values=site_colors)+
  scale_y_continuous(expand = expand_scale(mult = c(0.01, 0.12)))+
  stat_compare_means(comparisons = list(c("No","Yes")),label="p.signif",method = "wilcox.test",paired = F, tip.length = 0.005) +#Mann-whitney with P as reference
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))

#fisher test against atb classes
fisher.test(table(bloom_atb[c("DB_bloom","carbapenems")])) #p = 0.4172
fisher.test(table(bloom_atb[c("DB_bloom","glycopeptides")])) #p = 0.2098
fisher.test(table(bloom_atb[c("DB_bloom","penicillins")])) #p = 0.6378
fisher.test(table(bloom_atb[c("DB_bloom","cephalosporins")])) #p = 0.6431

fisher.test(table(bloom_atb[c("DB_bloom_entero","carbapenems")])) #p = 1
fisher.test(table(bloom_atb[c("DB_bloom_entero","glycopeptides")])) #p = 1
fisher.test(table(bloom_atb[c("DB_bloom_entero","penicillins")])) #p = 1
fisher.test(table(bloom_atb[c("DB_bloom_entero","cephalosporins")])) #p = 0.1434

fisher.test(table(bloom_atb[c("DB_bloom_lactob","carbapenems")])) #p = 1
fisher.test(table(bloom_atb[c("DB_bloom_lactob","glycopeptides")])) #p = 0.05683
fisher.test(table(bloom_atb[c("DB_bloom_lactob","penicillins")])) #p = 0.2901
fisher.test(table(bloom_atb[c("DB_bloom_lactob","cephalosporins")])) #p = 1

fisher.test(table(bloom_atb[c("GCF_bloom","carbapenems")])) #p = 0.4425
fisher.test(table(bloom_atb[c("GCF_bloom","glycopeptides")])) #p = 0.006825
fisher.test(table(bloom_atb[c("GCF_bloom","penicillins")])) #p = 0.08371
fisher.test(table(bloom_atb[c("GCF_bloom","cephalosporins")])) #p = 1

fisher.test(table(bloom_atb[c("OM_bloom","carbapenems")])) #p = 0.1426
fisher.test(table(bloom_atb[c("OM_bloom","glycopeptides")])) #p = 0.0667
fisher.test(table(bloom_atb[c("OM_bloom","penicillins")])) #p = 1
fisher.test(table(bloom_atb[c("OM_bloom","cephalosporins")])) #p = 1

#DB_bloom pval adjustment for number of atb classes
p.adjust(c(0.4172,0.2098,0.6378,0.6431), method = "bonferroni") #1.0000 0.8392 1.0000 1.0000
#GCF_bloom pval adjustment for number of atb classes
p.adjust(c(0.4425,0.006285,0.08371,1), method = "bonferroni") #1.00000 0.02514 0.33484 1.00000
#OM_bloom pval adjustment for number of atb classes
p.adjust(c(0.1426,0.0667,1,1), method = "bonferroni") #0.5704 0.2668 1.0000 1.0000

############################################################################
pacientes<-c("AAM"="#1","AAM"="#2","AMSCP"="#3","BMF"="#4","CEFC"="#5","DBHB"="#6",
             "ECPL"="#7","ENS"="#8","FGTC"="#9","FM"="#10","IABC"="#11","ILSM"="#12",
             "JPM"="#13","JVS"="#14","LEMG"="#15","LFCS"="#16","MEFD"="#17","MFML"="#18",
             "MMMC"="#19","MSBS"="#20","OBR"="#21","OSLG"="#22","OSN"="#23","PHOS"="#24",
             "RJT"="#25","RLS"="#26","SMG"="#27","TPN"="#28","VFF"="#29","VGS"="#30","VRM"="#31")
labeller <- function(variable,value){
  return(pacientes[value])
}

#ler as datas de coleta
#ler as datas de coleta
dacol<-read.csv("Files/Coletas_TMO.tsv", sep = "\t")
dacol<-subset(dacol, sample%in%c(1,2,3,4,5)) #remover MI e MF

#formatar as datas para o mesmo formato utilizado nos dados de ATB
dacol$data <- format(strptime(dacol$data, format = "%d.%m.%Y"), "%Y-%m-%d")
dacol<-dcast(dacol, ID ~ sample, value.var = "data")#unmelt
colnames(dacol)<-c("patient","dacolP","dacolA","dacolE","dacolE30","dacolE75")

daatb<-merge(dacol,read.csv("Files/ATB_data_220622.csv"))

#reshape to long data format
daatb<-reshape(daatb, direction='long', 
               varying=colnames(daatb)[9:ncol(daatb)], 
               timevar='atb',
               times = paste0("atb",seq(1,15,1)),
               v.names=c('beg','code','name',"end"),
               idvar='patient')
daatb$atb<-NULL
colnames(daatb)[9:12]<-c("code","name","end","beg")
daatb[daatb==""]<-NA
#delete rows with NA atb names
daatb <- drop_na(data = daatb, name)

#corrigindo essas inconsistências
daatb[daatb=="cepime"]<-"cefepime"
daatb[daatb=="targo"]<-"teico"
daatb[daatb=="Poli B"]<-"polimixina"
daatb[daatb=="tazocin"]<-"tazo"
daatb[daatb=="claro"]<-"claritro"
daatb[daatb=="amoxicilina"]<-"amoxi"
daatb[daatb=="meropenem"]<-"mero"
#padronizando nomes para terem 4 letras
daatb$name <- str_replace_all(daatb$name, c("clavulin"="clav","amoxi"="amox","cefepime"="cefe","metro"="metr","ceftriaxona"="ceft",
                                            "vanco"="vanc","teico"="teic","cipro"="cipr","claritro"="clar","bactrim"="bact",
                                            "polimixina"="poli","dapto"="dapt", "amica"="amic"))

#calculando dias em relação à infusão
daatb$dacolP.inf<-as.numeric(as.Date(daatb$dacolP)-as.Date(daatb$dainf))
daatb$dacolA.inf<-as.numeric(as.Date(daatb$dacolA)-as.Date(daatb$dainf))
daatb$dacolE.inf<-as.numeric(as.Date(daatb$dacolE)-as.Date(daatb$dainf))
daatb$dacolE30.inf<-as.numeric(as.Date(daatb$dacolE30)-as.Date(daatb$dainf))
daatb$dacolE75.inf<-as.numeric(as.Date(daatb$dacolE75)-as.Date(daatb$dainf))
daatb$daengr.inf<-as.numeric(as.Date(daatb$daengr)-as.Date(daatb$dainf))
daatb$beg.inf<-as.numeric(as.Date(daatb$beg)-as.Date(daatb$dainf))
daatb$end.inf<-as.numeric(as.Date(daatb$end)-as.Date(daatb$dainf))
daatb$duration<-daatb$end.inf-daatb$beg.inf
daatb$last_col<-apply(daatb[c("dacolE.inf","dacolE30.inf","dacolE75.inf")] %>% replace(is.na(.), 0), MARGIN = 1, max)

#remover usos de ATB fora do período de coletas
daatb<-subset(daatb,beg.inf<=100)
daatb<-subset(daatb,name!='levo'&name!='bact') #remover atbs da profilaxia

#gerar vetor com nomes completos dos atbs e df com as classes
atb_names<-c("clav"="amoxicilin clavunalate", "tazo"="piperacilin tazobactam", "amox"="amoxicilin", "cefe"="cefepime", "mero"="meropenem", 
             "metr"="metronidazole", "ceft"="ceftriaxone", "vanc"="vancomycin", "teic"="teicoplanin", "cipr"="ciprofloxacin", 
             "levo"="levofloxacin", "doxi"="doxycycline", "ampi"="ampicilin", "clar"="clarithromycin", "bact"="sulfamethoxazole trimethoprim", 
             "erta"="ertapenem", "poli"="polymyxin b", "dapt"="daptomycin","line"="linezolid", "tige" ="tigecycline", "amic"="amikacin")
atb_classes<-data.frame(name = c("tazo","cefe","mero","bact","clav","ceft","cipr","vanc","levo","line","poli","amox","teic","tige",
                                 "clar","doxi","metr","ampi","erta","dapt",'amic'), 
                        class = c("Penicillins","Cephalosporins","Carbapenems","Sulfonamides","Penicillins","Cephalosporins","Quinolones",
                                  "Glycopeptides","Quinolones","Oxazolidinones","Polymixins","Glycylcyclines","Glycopeptides","Glycylcyclines",
                                  "Macrolides","Tetracyclines","Nitroimidazoles","Penicillins","Carbapenems","Cyclic lipopeptides","Aminoglycosides"))
daatb <- merge(daatb,atb_classes)


#plot timelines focused on patients
atb_timeline_AAM<-ggplot(subset(daatb,patient=="AAM"))+
  geom_vline(aes(xintercept = 0), size = 0.5, color = "red", linetype = 'longdash')+
  annotate(geom="text", x=-6.2, y = Inf, label="SC infusion", color="red", vjust = 2.9, size = 3.2)+
  geom_vline(aes(xintercept = daengr.inf), size = 0.5, color = "blue", linetype = 'longdash')+
  annotate(geom="text", x=19, y = Inf, label="SC engraftment", color="blue", vjust = 2.9, size = 3.2)+
  geom_vline(aes(xintercept = 76), size = 0.5, color = "black", linetype = 'longdash')+
  annotate(geom="text", x=72.5, y = Inf, label="Death", color="black", vjust = 2.9, size = 3.2)+
  geom_linerange(aes(y=name,xmin=beg.inf-0.05,xmax=end.inf+0.05, color=class), size = 4)+
  geom_text(aes(y=name,x=(beg.inf+end.inf)/2,label=name), color = "black", size = 4)+
  scale_x_continuous(breaks = seq(-40,100,20),limits=c(-47,87),expand = c(0,0.05))+
  facet_wrap(.~'Antibiotic usage timeline')+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  #scale_color_manual(values=c(glasbey()))+
  labs(x="Days from transplant", color = "Antibiotic class")+
  #coord_cartesian(xlim = c(-20,100))+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(), strip.background = element_rect(color="white", fill="white", linetype="solid"),
        legend.position = 'bottom')

atb_timeline_AAM

atb_timeline_AFLF<-ggplot(subset(daatb,patient=="AFLF"))+
  geom_vline(aes(xintercept = 0), size = 0.5, color = "red", linetype = 'longdash')+
  annotate(geom="text", x=-4.1, y = Inf, label="SC infusion", color="red", vjust = 6.5, size = 3.2)+
  geom_vline(aes(xintercept = daengr.inf), size = 0.5, color = "blue", linetype = 'longdash')+
  annotate(geom="text", x=7.4, y = Inf, label="SC engraftment", color="blue", vjust = 6.5, size = 3.2)+
  geom_vline(aes(xintercept = 55), size = 0.5, color = "black", linetype = 'longdash')+
  annotate(geom="text", x=52.6, y = Inf, label="Death", color="black", vjust = 6.5, size = 3.2)+
  geom_linerange(aes(y=name,xmin=beg.inf-0.05,xmax=end.inf+0.05, color=class), size = 4)+
  geom_text(aes(y=name,x=(beg.inf+end.inf)/2,label=name), color = "black", size = 4)+
  scale_x_continuous(breaks = seq(-40,100,20),limits=c(-24,64),expand = c(0,0))+
  guides(color=guide_legend(nrow=2,byrow=TRUE))+
  facet_wrap(.~'Antibiotic usage timeline')+
  labs(x="Days from transplant", color = "Antibiotic class")+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(), strip.background = element_rect(color="white", fill="white", linetype="solid"),
        legend.position = 'bottom')

atb_timeline_AFLF

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

physeq<-qza_to_phyloseq(
  features="Files/ASV_table_final_goodsamples.qza",
  tree="Files/rooted-tree.qza",
  "Files/taxonomy_vs.qza",
  metadata = "Files/mapping_tmo_0.tsv")

physeq_tmo<-subset_samples(physeq, time!="MI"& time!="MF")
physeq_ra_tmo<-transform_sample_counts(coherent_taxa_names(physeq_tmo), function(x) x / sum(x) )

physeq_ra_tmo_subset<-subset_samples(physeq_ra_tmo,patient%in%c('AAM','AFLF'))
metadata<-as.data.frame(sample_data(physeq_ra_tmo_subset))
metadata$days.to.I<-metadata$days.to.E+c(rep(15,12),rep(28,11))
sample_data(physeq_ra_tmo_subset)<-metadata

ra_perday<-function(physeq_ra, physeq_ra_all_sites, pat, taxrank = "Genus", time = "days.to.E", min_ra=0.01, min_n_samples = 1, xlab = "Days", xlim_custom = c(-25,25)){
  physeq_ra<-tax_glom(physeq_ra, taxrank = taxrank)
  physeq_ra_all_sites<-tax_glom(physeq_ra_all_sites, taxrank = taxrank)
  taxa_list_id<-genefilter_sample(physeq_ra_all_sites, filterfun_sample(function(x) x >= min_ra), 
                                  A = min_n_samples)
  taxa_list_id<-rownames(subset(as.data.frame(taxa_list_id),taxa_list_id=="TRUE"))
  taxa_list<-subset(as.data.frame(tax_table(physeq_ra_all_sites)),
                    rownames(as.data.frame(tax_table(physeq_ra_all_sites)))%in%taxa_list_id)[[taxrank]]
  print(taxa_list)
  #physeq_ra<-prune_taxa(genefilter_sample(physeq_ra_all_sites, filterfun_sample(function(x) x >= min_ra), 
  #                                        A = min_n_samples),physeq_ra)
  #return(physeq_ra_all_sites)
  
  data<-as.data.frame(otu_table(physeq_ra))
  rownames(data)<-as.data.frame(physeq_ra@tax_table@.Data)[[taxrank]]
  data$taxa<-as.data.frame(physeq_ra@tax_table@.Data)[[taxrank]]
  data<-subset(data,taxa%in%taxa_list)
  data.melt<-reshape2::melt(data, c("taxa"))
  metadata<-as.data.frame(sample_data(physeq_ra))
  data.melt$time<-0
  data.melt$phase<-"X"
  
  print('ok1')
  for (sample in unique(data.melt$variable)){
    data.melt[data.melt$variable==sample,4]<-subset(metadata, rownames(metadata)==sample)[[time]][1] #collecting time in days
    data.melt[data.melt$variable==sample,5]<-subset(metadata, rownames(metadata)==sample)$time[1] #collecting phase
  }
  
  data.melt<-arrange(data.melt, taxa)
  data.melt$taxa<-factor(data.melt$taxa, levels = sort(unique(data.melt$taxa)))
  
  data.melt.bar<-arrange(data.melt, taxa)
  data.melt$taxa<-factor(data.melt$taxa, levels = sort(unique(data.melt$taxa)))
  
  ggplot(data.melt, aes(x=time, y=value)) + 
    geom_area(aes(x=time, y=value, fill=taxa), alpha=0.4)+
    geom_bar(data=data.melt,
             aes(x=time, y=value, fill=taxa), stat="identity", color = "black",width = 1, size = 0.1)+
    scale_fill_manual(values=c(glasbey(n=32),as.character(alphabet2(n=26))), guide = guide_legend(ncol=7,bycol=TRUE, title = taxrank))+
    geom_text(aes(y=0.8,x=time-2, label=phase), size = 3, alpha = 0.3)+
    theme(panel.grid = element_blank(), legend.position = "bottom",
          legend.text = element_text(size =8), legend.key.size = unit(0.7, 'lines'))+
    labs(x=xlab, y = 'Relative abundance (%)',fill = "Genus")+
    scale_y_continuous(expand = c(0,0), limits = c(0,1))+
    scale_x_continuous(limits = xlim_custom, breaks = seq(-40,100,20),expand = c(0,0))
  
}

ra_AFLF_FCG<-ra_perday(subset_samples(physeq_ra_tmo_subset,site=="FCG"&patient=='AFLF'), 
                    subset_samples(physeq_ra_tmo_subset,patient=='AFLF'),
                    min_ra = 0.1, time = "days.to.I", xlab = "Days from transplant", taxrank = "Genus", xlim_custom = c(-24,64))+
  facet_wrap(.~"GCF")+theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), legend.position = 'top', strip.background = element_rect(color="white", fill="white", linetype="solid"))

ra_AFLF_OM<-ra_perday(subset_samples(physeq_ra_tmo_subset,site=="MO"&patient=='AFLF'), 
          subset_samples(physeq_ra_tmo_subset,patient=='AFLF'),
          min_ra = 0.1, time = "days.to.I", xlab = "Days from transplant", taxrank = "Genus", xlim_custom = c(-24,64))+
  facet_wrap(.~"OM")+theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = 'none', strip.background = element_rect(color="white", fill="white", linetype="solid"))

ra_AFLF_SB<-ra_perday(subset_samples(physeq_ra_tmo_subset,site=="BF"&patient=='AFLF'), 
          subset_samples(physeq_ra_tmo_subset,patient=='AFLF'),
          min_ra = 0.1, time = "days.to.I", xlab = "Days from transplant", taxrank = "Genus", xlim_custom = c(-24,64))+
  facet_wrap(.~"SB")+theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), legend.position = 'none', strip.background = element_rect(color="white", fill="white", linetype="solid"))

            
F5C<-annotate_figure(
  plot_grid(
  ra_AFLF_FCG,
  ra_AFLF_OM+theme(axis.title.y = element_text(size = 9)),
  ra_AFLF_SB,
  atb_timeline_AFLF,
  align = 'v', axis = 'lr', nrow = 4, rel_heights=c(0.67,0.485,0.485,0.8)),
  fig.lab = "Patient #2", fig.lab.size = 12, fig.lab.pos = 'top')

F5C

ra_AAM_FCG<-ra_perday(subset_samples(physeq_ra_tmo_subset,site=="FCG"&patient=='AAM'), 
                       subset_samples(physeq_ra_tmo_subset,patient=='AAM'),
                       min_ra = 0.1, time = "days.to.I", xlab = "Days from transplant", taxrank = "Genus", xlim_custom = c(-47,87))+
  facet_wrap(.~"GCF")+theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), legend.position = 'top', strip.background = element_rect(color="white", fill="white", linetype="solid"))

ra_AAM_OM<-ra_perday(subset_samples(physeq_ra_tmo_subset,site=="MO"&patient=='AAM'), 
                      subset_samples(physeq_ra_tmo_subset,patient=='AAM'),
                      min_ra = 0.1, time = "days.to.I", xlab = "Days from transplant", taxrank = "Genus", xlim_custom = c(-47,87))+
  facet_wrap(.~"OM")+theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), legend.position = 'none', strip.background = element_rect(color="white", fill="white", linetype="solid"))

ra_AAM_SB<-ra_perday(subset_samples(physeq_ra_tmo_subset,site=="BF"&patient=='AAM'), 
                      subset_samples(physeq_ra_tmo_subset,patient=='AAM'),
                      min_ra = 0.1, time = "days.to.I", xlab = "Days from transplant", taxrank = "Genus", xlim_custom = c(-47,87))+
  facet_wrap(.~"SB")+theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), legend.position = 'none', strip.background = element_rect(color="white", fill="white", linetype="solid"))


FS10<-annotate_figure(
  plot_grid(
    ra_AAM_FCG,
    ra_AAM_OM+theme(axis.title.y = element_text(size = 9)),
    ra_AAM_SB,
    atb_timeline_AAM,
    align = 'v', axis = 'lr', nrow = 4, rel_heights=c(0.67,0.485,0.485,0.8)),
  fig.lab = "Patient #1", fig.lab.size = 12, fig.lab.pos = 'top')

FS10

annotate_figure(ggarrange(ra_AAM, atb_timeline_AAM, heights = c(0.65,0.35), nrow = 2), fig.lab = "#1", fig.lab.size = 14)
annotate_figure(ggarrange(ra_AFLF, atb_timeline_AFLF, heights = c(0.65,0.35), nrow = 2), fig.lab = "#2", fig.lab.size = 14)
