library(nloptr)
library(reshape2)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggbeeswarm)
theme_set(theme_bw())
library(stringr)
library(phyloseq)
#install.packages('outliers')
library(outliers)
library(qiime2R)
library(ggpubr)
library(tidytext)
library(pals)

microbiota<-read.csv('Files/300522_microbiota-variables.csv')

#############BLOOM
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

bloom_rnr<-merge(bloom_classification,microbiota)

fisher.test(table(bloom_rnr[c('OM_bloom','OM.fdis.e30.class_bin')])) #p = 0.1296; OR = 0.2457774 (0.02010178-1.74315342)

#############ATB
atb<-read.csv('Files/300622_ATB_classification.csv')
atb_rnr<-merge(atb,microbiota)

fisher.test(table(atb_rnr[c('cephalosporins','OM.fdis.e30.class_bin')])) #p = 0.2089; OR = 3.058668  (0.40911-24.15968)
fisher.test(table(atb_rnr[c('carbapenems','OM.fdis.e30.class_bin')])) #p = 0.4118; OR = 0.3615939 (0.02956483-2.59608492)
fisher.test(table(atb_rnr[c('glycopeptides','OM.fdis.e30.class_bin')])) #p = 0.6942; OR = 0.6214524 (0.07792966-3.99682826)
fisher.test(table(atb_rnr[c('penicillins','OM.fdis.e30.class_bin')])) #p = 1; OR = 0.8791066 (0.09649062-11.87283340)

wilcox.test(subset(atb_rnr,OM.fdis.e30.class_bin==1)$lot,subset(atb_rnr,OM.fdis.e30.class_bin==0)$lot)#p = 0.2986
median(subset(atb_rnr,OM.fdis.e30.class_bin==1)$lot, na.rm = T) #14.5
IQR(subset(atb_rnr,OM.fdis.e30.class_bin==1)$lot, na.rm = T) #8.75
median(subset(atb_rnr,OM.fdis.e30.class_bin==0)$lot, na.rm = T) #19
IQR(subset(atb_rnr,OM.fdis.e30.class_bin==0)$lot, na.rm = T) #12

wilcox.test(subset(atb_rnr,OM.fdis.e30.class_bin==1)$dot,subset(atb_rnr,OM.fdis.e30.class_bin==0)$dot)#p = 0.4639
median(subset(atb_rnr,OM.fdis.e30.class_bin==1)$dot, na.rm = T) #21
IQR(subset(atb_rnr,OM.fdis.e30.class_bin==1)$dot, na.rm = T) #19.25
median(subset(atb_rnr,OM.fdis.e30.class_bin==0)$dot, na.rm = T) #22
IQR(subset(atb_rnr,OM.fdis.e30.class_bin==0)$dot, na.rm = T) #19.25

atb_rnr.melt<-melt(atb_rnr[c('dot','lot','OM.fdis.e30.class_bin')],c('OM.fdis.e30.class_bin'))
atb_rnr.melt$variable<-as.character(atb_rnr.melt$variable)
atb_rnr.melt$value<-as.numeric(atb_rnr.melt$value)

atb_rnr.melt[atb_rnr.melt=="dot"]<-"DOT"
atb_rnr.melt[atb_rnr.melt=="lot"]<-"LOT"
atb_rnr.melt$OM.fdis.e30.class_bin[atb_rnr.melt$OM.fdis.e30.class_bin=="0"]<-"NR"
atb_rnr.melt$OM.fdis.e30.class_bin[atb_rnr.melt$OM.fdis.e30.class_bin=="1"]<-"R"

atb_rnr.melt$variable<-factor(atb_rnr.melt$variable,levels=c('LOT','DOT'))

FS14A<-ggplot(data = subset(atb_rnr.melt,is.na(OM.fdis.e30.class_bin)==F)) + 
  aes(x = OM.fdis.e30.class_bin, y = value, color = OM.fdis.e30.class_bin, group = OM.fdis.e30.class_bin) + 
  geom_boxplot(alpha=0.5, outlier.size = 1) +
  #geom_quasirandom(size=1)+ 
  facet_wrap(. ~ variable,scales = "free_y") + 
  labs(y = "Days", x = 'OM recovery') + guides(color=FALSE)+
  scale_color_manual(values=c('#616161', '#008d78'))+
  scale_y_continuous(expand = expand_scale(mult = c(0.0, 0.06)))+
  stat_compare_means(comparisons = list(c("NR","R")),method = "wilcox.test",paired = F,
                     tip.length = 0.005, size = 3.2, label.y.npc = 0.95, aes(label=round(as.numeric(..p.format..),2))) +#Mann-whitney with P as reference
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))

#############CLINICAL PARAMETERS
clindata<-read.csv('Files/DadosClinicosTMO_09.05.2021.csv')
clindata$Patient<-NULL
clindata$patient<-clindata$Patient.acronym
clindata$Patient.acronym<-NULL
clindata_rnr<-merge(clindata,microbiota)

wilcox.test(subset(clindata_rnr,OM.fdis.e30.class_bin==1)$Age,subset(clindata_rnr,OM.fdis.e30.class_bin==0)$Age)#p = 0.6268
median(subset(clindata_rnr,OM.fdis.e30.class_bin==1)$Age, na.rm = T) #51.09651
IQR(subset(clindata_rnr,OM.fdis.e30.class_bin==1)$Age, na.rm = T) #15.12526
median(subset(clindata_rnr,OM.fdis.e30.class_bin==0)$Age, na.rm = T) #51.07461
IQR(subset(clindata_rnr,OM.fdis.e30.class_bin==0)$Age, na.rm = T) #21.4538

fisher.test(table(clindata_rnr[c('Sex','OM.fdis.e30.class_bin')])) #p = 0.6942
fisher.test(table(clindata_rnr[c('Underlying.disease','OM.fdis.e30.class_bin')])) #p = 0.4223
fisher.test(table(clindata_rnr[c('Pretransplant.comorbidity..HCT.CI.','OM.fdis.e30.class_bin')])) #p = 0.4758
fisher.test(table(clindata_rnr[c('Disease.risk.index','OM.fdis.e30.class_bin')])) #p = 1
fisher.test(table(clindata_rnr[c('Conditioning.intensity','OM.fdis.e30.class_bin')])) #p = 0.6942
fisher.test(table(clindata_rnr[c('Total.body.irradiation','OM.fdis.e30.class_bin')])) #p = 0.3962
fisher.test(table(clindata_rnr[c('T.cell.depletion','OM.fdis.e30.class_bin')])) #p = 0.427
fisher.test(table(clindata_rnr[c('Graft.source','OM.fdis.e30.class_bin')])) #p = 0.4311
fisher.test(table(clindata_rnr[c('Donor','OM.fdis.e30.class_bin')])) #p = 0.535

#############ALPHA DIVERSITY
microbiota.melt<-melt(microbiota[c('Simpson.P_OM','Simpson.A_OM','Simpson.E_OM','Simpson.E30_OM','OM.fdis.e30.class_bin')],c('OM.fdis.e30.class_bin'))
microbiota.melt$variable<-as.character(microbiota.melt$variable)
microbiota.melt$value<-as.numeric(microbiota.melt$value)

microbiota.melt[microbiota.melt=="Simpson.P_OM"]<-"P"
microbiota.melt[microbiota.melt=="Simpson.A_OM"]<-"A"
microbiota.melt[microbiota.melt=="Simpson.E_OM"]<-"E"
microbiota.melt[microbiota.melt=="Simpson.E30_OM"]<-"E30"
microbiota.melt$OM.fdis.e30.class_bin[microbiota.melt$OM.fdis.e30.class_bin=="0"]<-"NR"
microbiota.melt$OM.fdis.e30.class_bin[microbiota.melt$OM.fdis.e30.class_bin=="1"]<-"R"

microbiota.melt$variable<-factor(microbiota.melt$variable,levels=c('P','A','E','E30'))

F8A<-ggplot(data = subset(microbiota.melt,is.na(OM.fdis.e30.class_bin)==F)) + 
  aes(x = OM.fdis.e30.class_bin, y = value, color = OM.fdis.e30.class_bin, group = OM.fdis.e30.class_bin) + 
  geom_boxplot(alpha=0.5, outlier.size = 1) +
  #geom_quasirandom(size=1)+
  facet_wrap(. ~ variable, nrow = 1) + 
  labs(y = "Diversity", x = 'OM recovery') + guides(color=FALSE)+
  scale_color_manual(values=c('#616161', '#008d78'))+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), limits=c(0,1.05), expand = c(0,0))+
  stat_compare_means(ref.group = 'NR', method = "wilcox.test",paired = F, tip.length = 0.005,size = 3.2, label.y.npc = 0.95, bracket.size = 0.3,
                     aes(label = ifelse(..p.format..<0.05,..p.signif..,round(as.numeric(..p.format..),2)))) +#Mann-whitney with P as reference
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))

#############TIME TO ENGRAFTMENT
clindata<-read.csv('Files/DadosClinicosTMO_09.05.2021.csv')
clindata$Patient<-NULL
clindata$patient<-clindata$Patient.acronym
clindata$Patient.acronym<-NULL
clindata_rnr<-merge(clindata,microbiota)

FS14B<-ggplot(data = subset(clindata_rnr,is.na(OM.fdis.e30.class)==F)) + 
  aes(x = OM.fdis.e30.class, y = Time.to.engraftment, color = OM.fdis.e30.class, group = OM.fdis.e30.class) + 
  geom_boxplot(alpha=0.5, outlier.size = 1) +
  #geom_quasirandom(size=1)+
  labs(y = "TTE (days)", x = 'OM recovery') + guides(color=FALSE)+
  scale_color_manual(values=c('#616161', '#008d78'))+
  scale_y_continuous(limits = c(0,28.5),breaks=c(0,7,14,21,28), expand = c(0,0))+
  stat_compare_means(comparisons = list(c("NR","R")),method = "wilcox.test",paired = F,
                     tip.length = 0.005, size = 3.2, label.y.npc = 0.95, aes(label=round(as.numeric(..p.format..),2))) +#Mann-whitney with P as reference
  theme(panel.grid = element_blank(),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))


#############HEMOGRAMA
hemograma<-read.csv("Files/Dados_hemog_resiliencia_250522.csv")

# hemograma$diff_exp_act<-as.Date(hemograma$date_exp)-as.Date(hemograma$date_act)
# 
# median(subset(hemograma,time=='P')$diff_exp_act, na.rm = T) #0
# range(subset(hemograma,time=='P')$diff_exp_act, na.rm = T) #-3, 2
# median(subset(hemograma,time=='A')$diff_exp_act, na.rm = T) #0
# range(subset(hemograma,time=='A')$diff_exp_act, na.rm = T) #0, 0
# median(subset(hemograma,time=='E')$diff_exp_act, na.rm = T) #0
# range(subset(hemograma,time=='E')$diff_exp_act, na.rm = T) #0, 0
# median(subset(hemograma,time=='E30')$diff_exp_act, na.rm = T) #0
# range(subset(hemograma,time=='E30')$diff_exp_act, na.rm = T) #-3, 4
# median(subset(hemograma,time=='I3M')$diff_exp_act, na.rm = T) #0
# range(subset(hemograma,time=='I3M')$diff_exp_act, na.rm = T) #-6, 6
# median(subset(hemograma,time=='I6M')$diff_exp_act, na.rm = T) #0
# range(subset(hemograma,time=='I6M')$diff_exp_act, na.rm = T) #-40, 28
# median(subset(hemograma,time=='I9M')$diff_exp_act, na.rm = T) #0
# range(subset(hemograma,time=='I9M')$diff_exp_act, na.rm = T) #-38, 48
# median(subset(hemograma,time=='I12M')$diff_exp_act, na.rm = T) #0
# range(subset(hemograma,time=='I12M')$diff_exp_act, na.rm = T) #-74, 12

hemograma<-subset(hemograma,time!="A")
#detecting outliers and removing
grubbs.test(hemograma$rbc)
grubbs.test(hemograma$rbc,opposite = T)

grubbs.test(hemograma$pla)
hemograma$pla[hemograma$pla==1150000]<-NA
grubbs.test(hemograma$pla)
grubbs.test(hemograma$pla,opposite = T)

grubbs.test(hemograma$wbc)
hemograma$wbc[hemograma$wbc==15530]<-NA
grubbs.test(hemograma$wbc)
hemograma$wbc[hemograma$wbc==14670]<-NA
grubbs.test(hemograma$wbc)
hemograma$wbc[hemograma$wbc==14020]<-NA
grubbs.test(hemograma$wbc)
hemograma$wbc[hemograma$wbc==13610]<-NA
grubbs.test(hemograma$wbc)
grubbs.test(hemograma$wbc,opposite = T)

grubbs.test(hemograma$neu)
hemograma$neu[hemograma$neu==14010]<-NA
grubbs.test(hemograma$neu)
hemograma$neu[hemograma$neu==13100]<-NA
grubbs.test(hemograma$neu)
hemograma$neu[hemograma$neu==10860]<-NA
grubbs.test(hemograma$neu)
hemograma$neu[hemograma$neu==10200]<-NA
grubbs.test(hemograma$neu)
grubbs.test(hemograma$neu,opposite = T)

grubbs.test(hemograma$lym)
hemograma$lym[hemograma$lym==7570]<-NA
grubbs.test(hemograma$lym)
grubbs.test(hemograma$lym,opposite = T)

grubbs.test(hemograma$mon)
hemograma$mon[hemograma$mon==7200]<-NA
grubbs.test(hemograma$mon)
hemograma$mon[hemograma$mon==2300]<-NA
grubbs.test(hemograma$mon)
hemograma$mon[hemograma$mon==1810]<-NA
grubbs.test(hemograma$mon)
hemograma$mon[hemograma$mon==1480]<-NA
grubbs.test(hemograma$mon)
grubbs.test(hemograma$mon,opposite = T)

hemograma_rnr<-merge(hemograma,microbiota[c('patient','OM.fdis.e30.class')])

hemograma_rnr.melt<-reshape2::melt(hemograma_rnr[c('patient','time','rbc','wbc','neu','lym','mon','pla', 'OM.fdis.e30.class')],
                              c('patient','time','OM.fdis.e30.class'))
hemograma_rnr.melt$time<-factor(hemograma_rnr$time,levels=c('P','E','E30','I3M','I6M','I9M','I12M'))
hemograma_rnr.melt$variable<-as.character(hemograma_rnr.melt$variable)
hemograma_rnr.melt[hemograma_rnr.melt=="rbc"]<-"Erythrocytes"
hemograma_rnr.melt[hemograma_rnr.melt=="wbc"]<-"Leukocytes"
hemograma_rnr.melt[hemograma_rnr.melt=="neu"]<-"Neutrophils"
hemograma_rnr.melt[hemograma_rnr.melt=="lym"]<-"Lymphocytes"
hemograma_rnr.melt[hemograma_rnr.melt=="mon"]<-"Monocytes"
hemograma_rnr.melt[hemograma_rnr.melt=="pla"]<-"Platelets"
hemograma_rnr.melt$variable<-factor(hemograma_rnr.melt$variable,levels = c("Erythrocytes","Neutrophils",
                                                                           "Platelets","Lymphocytes",
                                                                           "Leukocytes","Monocytes"))
hemograma_rnr.melt$ref_low<-0
hemograma_rnr.melt$ref_up<-0
hemograma_rnr.melt$ref_low[hemograma_rnr.melt$variable=="Erythrocytes"]<-3900000
hemograma_rnr.melt$ref_up[hemograma_rnr.melt$variable=="Erythrocytes"]<-5700000
hemograma_rnr.melt$ref_low[hemograma_rnr.melt$variable=="Leukocytes"]<-3500
hemograma_rnr.melt$ref_up[hemograma_rnr.melt$variable=="Leukocytes"]<-10500
hemograma_rnr.melt$ref_low[hemograma_rnr.melt$variable=="Neutrophils"]<-1700
hemograma_rnr.melt$ref_up[hemograma_rnr.melt$variable=="Neutrophils"]<-7000
hemograma_rnr.melt$ref_low[hemograma_rnr.melt$variable=="Lymphocytes"]<-900
hemograma_rnr.melt$ref_up[hemograma_rnr.melt$variable=="Lymphocytes"]<-2900
hemograma_rnr.melt$ref_low[hemograma_rnr.melt$variable=="Monocytes"]<-300
hemograma_rnr.melt$ref_up[hemograma_rnr.melt$variable=="Monocytes"]<-900
hemograma_rnr.melt$ref_low[hemograma_rnr.melt$variable=="Platelets"]<-150000
hemograma_rnr.melt$ref_up[hemograma_rnr.melt$variable=="Platelets"]<-450000

F8C<-ggplot(data = subset(hemograma_rnr.melt, is.na(OM.fdis.e30.class)==F&time%in%c('P','E','E30'))) + 
  aes(x = time, y = value, color = OM.fdis.e30.class, group = interaction(time,OM.fdis.e30.class)) + 
  geom_hline(aes(yintercept=ref_low), color = "red", alpha = 0.5, linetype = "dotted")+
  geom_hline(aes(yintercept=ref_up), color = "red", alpha = 0.5, linetype = "dotted")+
  #geom_point(position = position_dodge(width=0.75),size=0.5,alpha=0.5) +
  geom_quasirandom(size=1, dodge.width = 0.5, alpha = 0.5)+
  facet_wrap(. ~ variable, nrow = 3, scales = "free") + 
  labs(y = "Blood count", x = 'Timepoint', color = 'OM recovery') +
  scale_color_manual(values=c('#616161', '#008d78'))+
  scale_y_continuous(expand = expand_scale(mult = c(0.04, 0.2)))+
  #scale_y_log10(expand = expand_scale(mult = c(0.04, 0.18)))+
  stat_compare_means(label="p.signif",method = "wilcox.test",paired = F, tip.length = 0.005, size = 5.5, color='black',
                     hide.ns = TRUE) +#Mann-whitney with P as reference
  theme(panel.grid = element_blank(), legend.position = 'top',
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))

FS14C<-ggplot(data = subset(hemograma_rnr.melt, is.na(OM.fdis.e30.class)==F&time%in%c('I3M','I6M','I9M','I12M'))) + 
  aes(x = time, y = value, color = OM.fdis.e30.class, group = interaction(time,OM.fdis.e30.class)) + 
  geom_hline(aes(yintercept=ref_low), color = "red", alpha = 0.5, linetype = "dotted")+
  geom_hline(aes(yintercept=ref_up), color = "red", alpha = 0.5, linetype = "dotted")+
  geom_hline(yintercept = 0, color = "red", alpha = 0, linetype = "dotted")+
  #geom_point(position = position_dodge(width=0.75),size=0.5,alpha=0.5) +
  geom_quasirandom(size=1, dodge.width = 0.5, alpha = 0.5)+
  facet_wrap(. ~ variable, nrow = 3, scales = "free") + 
  labs(y = "Blood count", x = 'Timepoint', color = 'OM recovery') +
  scale_color_manual(values=c('#616161', '#008d78'))+
  scale_y_continuous(expand = expand_scale(mult = c(0.04, 0.2)))+
  stat_summary(fun = median, geom = "line", size=0.75, alpha = 0.5, aes(group = OM.fdis.e30.class))+
  #scale_y_log10(expand = expand_scale(mult = c(0.04, 0.18)))+
  stat_compare_means(label="p.signif",method = "wilcox.test",paired = F, tip.length = 0.005, color = 'black', size = 5.5,
                     hide.ns = TRUE) +#Mann-whitney with P as reference
  theme(panel.grid = element_blank(), legend.position = 'top',
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.border = element_blank(),axis.line = element_line(color = 'black'))