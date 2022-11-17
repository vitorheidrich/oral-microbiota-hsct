library(ggplot2)
library(reshape2)
library(tidyr)
library(stringr)
library(plyr)
library(dplyr)
library(cowplot)
library(ggpubr)
library(vistime)
library(ggforce)
library(ggridges)
theme_set(theme_bw())

pacientes<-c("AAM"="#1","AFLF"="#2","AMSCP"="#3","BMF"="#4","CEFC"="#5","DBHB"="#6",
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
atb_names<-c("clav"="amoxicillin clavulanate", "tazo"="piperacillin tazobactam", "amox"="amoxicillin", "cefe"="cefepime", "mero"="meropenem", 
             "metr"="metronidazole", "ceft"="ceftriaxone", "vanc"="vancomycin", "teic"="teicoplanin", "cipr"="ciprofloxacin", 
             "levo"="levofloxacin", "doxi"="doxycycline", "ampi"="ampicillin", "clar"="clarithromycin", "bact"="sulfamethoxazole trimethoprim", 
             "erta"="ertapenem", "poli"="polymixin b", "dapt"="daptomycin","line"="linezolid", "tige" ="tigecycline", "amic"="amikacin")
atb_classes<-data.frame(name = c("tazo","cefe","mero","bact","clav","ceft","cipr","vanc","levo","line","poli","amox","teic","tige",
                                 "clar","doxi","metr","ampi","erta","dapt",'amic'), 
                        class = c("penicillins","cephalosporins","carbapenems","sulfonamides","penicillins","cephalosporins","quinolones",
                                  "glycopeptides","quinolones","oxazolidinones","polymixins","glycylcyclines","glycopeptides","glycylcyclines",
                                  "macrolides","tetracyclines","nitroimidazoles","penicillins","carbapenems","cyclic lipopeptides","aminoglycosides"))
daatb <- merge(daatb,atb_classes)
daatb_restrito<-subset(daatb,end.inf>=dacolP.inf&beg.inf<=dacolE30.inf)

#CONTINUAR DAQUI
#agora vamos calcular alguns números simples
length(unique(daatb_restrito$name)) #17 atbs utilizados no total
length(unique(daatb_restrito$class)) #12 atbs utilizados no total
#é preciso sempre lembrar que o pcte IABC não está representando nos dataframes
#mín de atbs usados: 0 (IABC)
max(as.numeric(table(daatb_restrito$patient))) #10
median(c(0,as.numeric(table(daatb_restrito$patient)))) #3
#median (range): 3 (0-10)
barplot(sort(c(0,as.numeric(table(daatb_restrito$patient))))) #view

length(interaction(daatb_restrito$patient,daatb_restrito$name)) #105 usos de atbs
length(unique(interaction(daatb_restrito$patient,daatb_restrito$name))) #86 usos inéditos de atbs (não conta quando o paciente utiliza pela segunda vez um atb)
#min de atbs diferentes usados: 0 (IABC)
max(table(gsub("\\..*","",unique(interaction(daatb_restrito$patient,daatb_restrito$name))))) #9
median(c(IABC=0,table(gsub("\\..*","",unique(interaction(daatb_restrito$patient,daatb_restrito$name)))))) #3
#median (range): 3 (0-9)
barplot(sort(c(IABC=0,table(gsub("\\..*","",unique(interaction(daatb$patient,daatb$name))))))) #view

#numero de pacientes que usam ao menos uma vez cada um dos atbs
table(gsub(".*\\.","",unique(interaction(daatb_restrito$patient,daatb_restrito$name))))
# ampi cefe ceft cipr clar clav dapt doxi erta line mero metr poli tazo teic tige vanc 
# 1   22    2    2    2    1    2    1    1    2   19    3    2    6    5    1   14

length(unique(interaction(daatb_restrito$patient,daatb_restrito$class))) #81 usos inéditos de classe (não conta quando o paciente utiliza pela segunda vez a classe)
max(table(gsub("\\..*","",unique(interaction(daatb_restrito$patient,daatb_restrito$class))))) #7
median(c(IABC=0,table(gsub("\\..*","",unique(interaction(daatb_restrito$patient,daatb_restrito$class)))))) #3
#median (range): 3 (0-7)

#numero de pacientes que usam ao menos uma vez um atb de dada classe
table(gsub(".*\\.","",unique(interaction(daatb_restrito$patient,daatb_restrito$class))))
#carbapenems      cephalosporins cyclic lipopeptides       glycopeptides      glycylcyclines 
#19                  22                   2                  18                   1 
#macrolides     nitroimidazoles      oxazolidinones          penicillins          polymixins 
#2                   3                   2                   7                   2 
#quinolones       tetracyclines 
#2                   1

#classes de interesse para classificações (usados por >20%): cephalosporins, carbapenems, glycopeptides, penicillins

#numero de clases que cada paciente usa
#min 0 (IABC)
max(table(gsub("\\..*","",unique(interaction(daatb_restrito$patient,daatb_restrito$class))))) #7
median(c(IABC=0,table(gsub("\\..*","",unique(interaction(daatb_restrito$patient,daatb_restrito$class)))))) #3
#median (range): 3 (0-7)
sum(table(gsub("\\..*","",unique(interaction(daatb_restrito$patient,daatb_restrito$class)))))
barplot(sort(c(IABC=0,table(gsub("\\..*","",unique(interaction(daatb_restrito$patient,daatb_restrito$class))))))) #view

#plot frequencia de uso de cada atb, com cores indicando as classes
freq.df<-as.data.frame(table(gsub(".*\\.","",unique(interaction(daatb_restrito$patient,daatb_restrito$name))))/31)
colnames(freq.df)<-c("name","freq")
freq.df<-merge(freq.df,atb_classes)
freq.df<-arrange(freq.df,class,name)
freq.df$name<-factor(freq.df$name, levels = freq.df$name)

rect.df<-cbind(as.data.frame(cumsum(table(freq.df$class))+0.5),c(0,head(as.numeric(cumsum(table(freq.df$class))),length(unique(freq.df$class))-1))+0.5)
colnames(rect.df)<-c("fim","ini")
rect.df
FS9A<-ggplot(freq.df)+
        geom_rect(data=rect.df, aes(ymin=ini,ymax=fim,xmin=0,xmax=1,fill=rownames(rect.df)), alpha=0.4)+
        geom_bar(aes(y=name,x=freq,fill=class),stat="identity", color="black")+
        geom_text(data=rect.df,aes(label=rownames(rect.df),y=ini+0.7,x=0.75),vjust=1.05)+
        labs(x="Proportion of patients", y = "Antibiotic agent")+geom_hline(yintercept=rect.df$fim, lty="longdash",size=0.3)+
        scale_y_discrete(labels=atb_names, expand = c(0,0.5))+scale_x_continuous(limits=c(0,1), expand = c(0,0), breaks = c(0,0.2,0.4,0.6,0.8,1))+
        theme(panel.grid = element_blank())+
  guides(fill=F)

#plot frequencia de uso de cada classe
freq.df<-as.data.frame(table(gsub(".*\\.","",unique(interaction(daatb_restrito$patient,daatb_restrito$class))))/31)
colnames(freq.df)<-c("class","freq")
#freq.df<-merge(freq.df,atb_classes)
freq.df<-arrange(freq.df,class)
freq.df$class<-factor(freq.df$class, levels = freq.df$class)

#rect.df<-cbind(as.data.frame(cumsum(table(freq.df$class))+0.5),c(0,head(as.numeric(cumsum(table(freq.df$class))),length(unique(freq.df$class))-1))+0.5)
#colnames(rect.df)<-c("fim","ini")
#rect.df
FS9B<-ggplot(freq.df)+
        #geom_rect(data=rect.df, aes(xmin=ini,xmax=fim,ymin=0,ymax=1,fill=rownames(rect.df)), alpha=0.4)+
        geom_bar(aes(y=class,x=freq,fill=class),stat="identity", color="black")+
        #geom_text(data=rect.df,aes(label=rownames(rect.df),x=ini+0.5,y=1),angle=90,hjust=1.05)+
        labs(x="Proportion of patients", y = "Antibiotic class")+
        #geom_vline(xintercept=rect.df$fim, lty="longdash",size=0.3)+
        scale_y_discrete(expand = c(0,0.5))+
        scale_x_continuous(limits=c(0,1), expand = c(0,0), breaks = c(0,0.2,0.4,0.6,0.8,1))+
        geom_vline(xintercept = 0.2, linetype = "longdash")+
        theme(panel.grid = element_blank())+guides(fill=F)

#tilemap com uso de cada ATB por paciente
tilemap.df<-distinct(daatb_restrito,name,patient, .keep_all= TRUE) #remover usos do mesmo atb por um paciente
tilemap.df<-tilemap.df[c("name","patient","class")] #keep just relevant columns
tilemap.df<-rbind(tilemap.df,c(NA,"IABC",NA)) #include patient IABC (no ATBs used)
tilemap.df<-melt(dcast(tilemap.df, patient ~ name, value.var = "name"),c("patient")) #reshape df
colnames(tilemap.df)<-c("patient","name","value")
tilemap.df<-distinct(merge(tilemap.df,daatb[c("name","class")],by="name"),name,patient, .keep_all= TRUE) #remerge to get atb classes
tilemap.df$value<-ifelse(is.na(tilemap.df$value)==T,"No","Yes")
tilemap.df<-arrange(tilemap.df,class,name)
tilemap.df$name<-factor(tilemap.df$name, levels = unique(tilemap.df$name))

F5A<-ggplot(tilemap.df, aes(x=patient,y=forcats::fct_rev(name),fill=value, color = class))+geom_tile(size=1)+
        scale_fill_manual(values=c("grey60","lemonchiffon"))+coord_fixed()+
        scale_x_discrete(labels = pacientes)+
        scale_y_discrete(labels = atb_names)+labs(x = "Patient", y = "Antibiotic agent", color = "Antibiotic class", fill = "Antibiotic use")+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              legend.text = element_text(size = 8.5), legend.key.size = unit(0.7, 'lines'))

#####03/11/21 - atb timelines per patient
median(daatb_restrito$duration) #7 dias
max(daatb_restrito$duration) #36 dias
min(daatb_restrito$duration) #1 dias


#plot timelines focused on patients
ggplot(daatb)+
        #geom_rect(aes(xmin=dacolP.inf,xmax=last_col,ymin=-Inf,ymax=+Inf), fill="gray90", alpha=0.4)+
        geom_vline(aes(xintercept = dacolP.inf), size = 0.5, linetype = "dashed", color = "red")+
        #geom_label(aes(y = 0.2, x = dacolP.inf), label = "P")+
        #geom_vline(aes(xintercept = dacolA.inf), size = 0.5, linetype = "dashed")+
        #geom_label(aes(y = 0.2, x = dacolA.inf), label = "A")+
        #geom_vline(aes(xintercept = dacolE.inf), size = 0.5, linetype = "dashed")+
        geom_vline(aes(xintercept = dacolE30.inf), size = 0.5, linetype = "dashed", color = "blue")+
        #geom_label(aes(y = 0.2, x = dacolE30.inf), label = "E30")+
        #geom_vline(aes(xintercept = dacolE75.inf), size = 0.5, linetype = "dashed")+
        #geom_label(aes(y = 0.2, x = dacolE75.inf), label = "E75")+
        geom_vline(aes(xintercept = 0), size = 1, color = "red")+
        #geom_label(aes(y = 0.2, x = 0), label = "I")+
        geom_vline(aes(xintercept = daengr.inf), size = 1, color = "blue")+
        #geom_label(aes(y = 0.2, x = dacolE.inf), label = "E")+
        geom_linerange(aes(y=name,xmin=beg.inf-0.05,xmax=end.inf+0.05, color=class), size = 4)+
        geom_text(aes(y=name,x=(beg.inf+end.inf)/2,label=name), color = "black", size = 4)+
        #scale_y_discrete(expand = c(0.1,0.1))+
        scale_x_continuous(breaks = seq(-75,151,25))+
        facet_col(.~patient, scales = "free_y", space = "free", labeller = labeller)+
        guides(color=guide_legend(nrow=2,byrow=TRUE))+
        labs(x="Days from transplant", color = "Antibiotic class")+
        coord_cartesian(xlim = c(-50,100))+
        theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
              panel.grid = element_blank(), strip.background =element_rect(fill="lightgray"), legend.position = "top") #1000x2000

#plot timelines focused on ATBs with ridgeline plot
FS1<-annotate_figure(plot_grid(
ggplot(as.data.frame(lapply(daatb, rep, daatb$duration))%>%
               group_by(name,patient) %>%
               mutate(taking = beg.inf-1+row_number()) %>%
               ungroup(), aes(x = taking, y = name, fill = class, height = ..count..)) +
        geom_density_ridges(rel_min_height = 0.01, alpha = 0.5, binwidth = 1, scale = 3, stat = "binline") +
        #scale_fill_manual(name = "Temp. [F]", option = "C") +
        labs(x="Days from transplant", y = "Antibiotic", fill = "Antibiotic class")+
        #theme_ipsum() +
        coord_cartesian(xlim = c(-50,100))+
        scale_y_discrete(labels = atb_names)+
        theme(legend.position="top",panel.spacing = unit(0.1, "lines"),strip.text.x = element_text(size = 8),
              panel.grid.minor = element_blank()),
ggplot(as.data.frame(lapply(daatb, rep, daatb$duration))%>%
               group_by(name,patient) %>%
               mutate(taking = beg.inf-1+row_number()) %>%
               ungroup(), aes(x = taking, y = class, fill = class, height = ..count..)) +
        geom_density_ridges(rel_min_height = 0.01, alpha = 0.5, binwidth = 1, scale = 3, stat = "binline") +
        #scale_fill_manual(name = "Temp. [F]", option = "C") +
        labs(x="Days from transplant", y = "Antibiotic class")+
        coord_cartesian(xlim = c(-50,100))+
        #theme_ipsum() +
        theme(legend.position="none",panel.spacing = unit(0.1, "lines"),strip.text.x = element_text(size = 8),
              panel.grid.minor = element_blank()),
ncol = 1, nrow = 2, align = "v", axis = "lr", rel_heights = c(2,1.1)),
fig.lab = 'ab', fig.lab.size = 18)

#Calcular DOT e LOT
#é preciso limitar a duração para os usos que ultrapassam os limites (P e E30)
daatb_restrito_cap<-daatb_restrito
for (i in 1:nrow(daatb_restrito_cap)){
  diff_E30<-daatb_restrito_cap[["end.inf"]][[i]]-daatb_restrito_cap[["dacolE30.inf"]][[i]]
  diff_P<-daatb_restrito_cap[["dacolP.inf"]][[i]]-daatb_restrito_cap[["beg.inf"]][[i]]
  if(diff_E30>0){
    daatb_restrito_cap[["end.inf"]][[i]]<-daatb_restrito_cap[["dacolE30.inf"]][[i]]
    daatb_restrito_cap[["duration"]][[i]]<-daatb_restrito_cap[["duration"]][[i]]-diff_E30
  }
  if(diff_P>0){
    daatb_restrito_cap[["beg.inf"]][[i]]<-daatb_restrito_cap[["dacolP.inf"]][[i]]
    daatb_restrito_cap[["duration"]][[i]]<-daatb_restrito_cap[["duration"]][[i]]-diff_P
  }
}

dot_lot<-data.frame(patient=character(),dot=numeric(),lot=numeric())
for (pt in unique(daatb_restrito_cap$patient)){
  daatb_restrito_cap_pt<-subset(daatb_restrito_cap,patient==pt)
  dot<-sum(daatb_restrito_cap_pt$duration)
  
dut<-c()
for (i in 1:nrow(daatb_restrito_cap_pt)){
  dut<-append(dut,as.Date(daatb_restrito_cap_pt[["beg"]][[i]])+
    seq(0,daatb_restrito_cap_pt[["duration"]][[i]]-1,1))
}
counter=0
for (day in as.Date(daatb_restrito_cap_pt[["dacolP"]][[1]])+
            seq(0,daatb_restrito_cap_pt[["dacolE30.inf"]][[1]]-daatb_restrito_cap_pt[["dacolP.inf"]][[1]])){
  if (day %in% dut){counter=counter+1}
  
}
dot_lot<-rbind(dot_lot,c(pt,dot,counter))
}
colnames(dot_lot)<-c('patient','dot','lot')

#Classificações de uso de classes de ATBs específicas
cepha_class<-as.data.frame(table(unique(subset(daatb_restrito,class=="cephalosporins")$patient))) %>% rename('patient'='Var1',"cephalosporins"='Freq') %>% mutate(across(everything(), as.character))
carba_class<-as.data.frame(table(unique(subset(daatb_restrito,class=="carbapenems")$patient))) %>% rename('patient'='Var1',"carbapenems"='Freq') %>% mutate(across(everything(), as.character))
glyco_class<-as.data.frame(table(unique(subset(daatb_restrito,class=="glycopeptides")$patient))) %>% rename('patient'='Var1',"glycopeptides"='Freq') %>% mutate(across(everything(), as.character))
penic_class<-as.data.frame(table(unique(subset(daatb_restrito,class=="penicillins")$patient))) %>% rename('patient'='Var1',"penicillins"='Freq') %>% mutate(across(everything(), as.character))

atb_classification<-join(join(join(join(dot_lot,carba_class,by = "patient"), 
               glyco_class,by = "patient"), penic_class,by = "patient"), cepha_class,by = "patient")%>% replace(is.na(.), 0) #joining dfs
atb_classification<-rbind(atb_classification,c("IABC",0,0,0,0,0,0))
atb_classification<-rbind(atb_classification,c("AMSCP",NA,NA,NA,NA,NA,NA,NA))
#min dot = 0 (IABC)
#min lot = 0
median(as.numeric(atb_classification$dot)) #median dot = 22
median(as.numeric(atb_classification$lot)) #median lot = 15.5
max(as.numeric(atb_classification$dot)) #max dot = 112
max(as.numeric(atb_classification$lot)) #max lot = 58
write.csv(atb_classification, 'Files/300622_ATB_classification.csv', row.names = F) #saving as csv

min(daatb_restrito_cap$duration) #1 dia
max(daatb_restrito_cap$duration) #33 dias
median(daatb_restrito_cap$duration) #7 dias
