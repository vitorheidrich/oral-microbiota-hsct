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
#install.packages("rJava")
#install.packages("xlsx")
library(openxlsx)
library(ggnewscale)
library(ggrepel)
pacientes<-c("AAM"="#1","AFLF"="#2","AMSCP"="#3","BMF"="#4","CEFC"="#5","DBHB"="#6",
             "ECPL"="#7","ENS"="#8","FGTC"="#9","FM"="#10","IABC"="#11","ILSM"="#12",
             "JPM"="#13","JVS"="#14","LEMG"="#15","LFCS"="#16","MEFD"="#17","MFML"="#18",
             "MMMC"="#19","MSBS"="#20","OBR"="#21","OSLG"="#22","OSN"="#23","PHOS"="#24",
             "RJT"="#25","RLS"="#26","SMG"="#27","TPN"="#28","VFF"="#29","VGS"="#30","VRM"="#31")

coletas <- as.data.frame(read.csv("Files/Coletas_TMO_dcast.csv", header = TRUE))
coletas <- coletas[c("patient","dacolP","dacolA","dacolE","dacolE30","dacolE75","daengr")]
coletas$dacolP <- as.Date(coletas$dacolP, format= "%Y-%m-%d")
coletas$dacolA <- as.Date(coletas$dacolA, format= "%Y-%m-%d")
coletas$dacolE <- as.Date(coletas$dacolE, format= "%Y-%m-%d")
coletas$dacolE30 <- as.Date(coletas$dacolE30, format= "%Y-%m-%d")
coletas$dacolE75 <- as.Date(coletas$dacolE75, format= "%Y-%m-%d")
coletas$daengr <- as.Date(coletas$daengr, format= "%Y-%m-%d")

coletas$P<-as.numeric(coletas$dacolP-coletas$daengr)
coletas$A<-as.numeric(coletas$dacolA-coletas$daengr)
coletas$E<-as.numeric(coletas$dacolE-coletas$daengr)
coletas$E30<-as.numeric(coletas$dacolE30-coletas$daengr)
coletas$E75<-as.numeric(coletas$dacolE75-coletas$daengr)

median(as.numeric(coletas$P),na.rm=TRUE) #-24
median(as.numeric(coletas$A),na.rm=TRUE) #-10
median(as.numeric(coletas$E),na.rm=TRUE) #0
median(as.numeric(coletas$E30),na.rm=TRUE) #31
median(as.numeric(coletas$E75),na.rm=TRUE) #74.5

min(as.numeric(coletas$P),na.rm=TRUE) #-45
min(as.numeric(coletas$A),na.rm=TRUE) #-25
min(as.numeric(coletas$E),na.rm=TRUE) #0
min(as.numeric(coletas$E30),na.rm=TRUE) #20
min(as.numeric(coletas$E75),na.rm=TRUE) #60

max(as.numeric(coletas$P),na.rm=TRUE) #-17
max(as.numeric(coletas$A),na.rm=TRUE) #-1
max(as.numeric(coletas$E),na.rm=TRUE) #8
max(as.numeric(coletas$E30),na.rm=TRUE) #45
max(as.numeric(coletas$E75),na.rm=TRUE) #131

coletas.melt<-reshape2::melt(coletas[c("patient","P","A","E","E30","E75")], c("patient"))
coletas.melt<-na.omit(coletas.melt)
coletas.melt$variable<-as.character(coletas.melt$variable)
coletas.melt<-rbind(coletas.melt,c("AAM","Death",48))
coletas.melt<-rbind(coletas.melt,c("AFLF","Death",40))
coletas.melt<-rbind(coletas.melt,c("AMSCP","Death",20))
coletas.melt<-rbind(coletas.melt,c("OBR","Death",141))
coletas.melt<-rbind(coletas.melt,c("VRM","Death",41))
coletas.melt$variable<-factor(coletas.melt$variable, levels = c("P","A","E","E30","E75","Death"))
coletas.melt$value<-as.numeric(coletas.melt$value)

FS3<-ggplot()+
  geom_point(data=subset(coletas.melt,variable!="Death"), aes(x=patient,y=value,color=variable),size=2)+
  geom_point(data=subset(coletas.melt,variable=="Death"), aes(x=patient,y=value,shape=variable),size=2)+
  scale_color_manual(values=c("red3","darkorange","yellow3","green3","dodgerblue3"))+
  scale_shape_manual(values = c(4))+
  coord_flip()+
  geom_hline(yintercept = c(-24),color="red3", linetype = 'longdash')+
  geom_hline(yintercept = c(-10),color="darkorange", linetype = 'longdash')+
  geom_hline(yintercept = c(-0),color="yellow3", linetype = 'longdash')+
  geom_hline(yintercept = c(31),color="green3", linetype = 'longdash')+
  geom_hline(yintercept = c(74.5),color="dodgerblue3", linetype = 'longdash')+
  labs(x="Patient",y="Days relative to engraftment",color="Timepoint",shape="")+
  scale_y_continuous(breaks=c(-45,-30,-15,0,15,30,45,60,75,90,105,120,135,150), limits=c(-45,143))+
  scale_x_discrete(labels = pacientes)+
  theme(panel.grid.minor = element_blank(), legend.position = 'top')
FS3
