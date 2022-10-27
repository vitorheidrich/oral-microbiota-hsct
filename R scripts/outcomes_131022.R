library(qiime2R)
library(ggplot2)
library(btools)
library(survminer)
library(scales)
library(survival)
library(cmprsk)
#install.packages(c("survival","cmprsk","survminer"))
library(vegan)
library(cowplot)
library(phyloseq)
library(dplyr)
library(plyr)
library(gplots)
source("plotDistances.R")
library(rbiom)
library(QsRutils)
library(SRS)
library(reshape2)
library(purrr)
library(matrixStats)
library(ggbeeswarm)
library(ggstance)
library(ggpubr)
library(magrittr)
theme_set(theme_bw())

clindata<-read.csv("Files/DadosClinicosTMO_09.05.2021.csv")
clindata$Patient<-NULL
clindata$patient<-clindata$Patient.acronym
clindata$Patient.acronym<-NULL
microdata<-read.csv("Files/300522_microbiota-variables.csv")
clinmicrodata<-merge(clindata,microdata)


#importante, pacientes com recaída antes da coleta final (seja qual for a que estivermos considerando), precisam ser excluídos da
#associação entre este desfecho ou PFS com as variáveis que incorporam a coleta final
#são eles: VRM, OBR e AMSCP.
#o paciente AMSCP não fez coletas após a enxertia, então ele já está excluído
######SUMMARIZE F&G#######
summarize_fg <-function(input_df,time_col,stat_col,beg,end,cencode_char,failcode_char){
  fg_results <- data.frame(matrix(nrow=end-beg+1,ncol=9))
  #print(rownames(input_df[beg:end]))
  rownames(fg_results)<-colnames(input_df[beg:end])
  colnames(fg_results)<-c("P","HR","HR 2.5","HR 97.5","N","NGroup1","NGroup0","%Group1","%Group0")
  input_df_copy<-input_df
  for (i in colnames(input_df[beg:end])){
    input_df<-input_df_copy
    if (i%in%c("DB.fdis.e30.class_bin","GCF.fdis.e30.class_bin","OM.fdis.e30.class_bin")&stat_col=="Relapse"){
      input_df<-subset(input_df,patient!="VRM"&patient!="OBR")
    }
    skip_to_next <- FALSE
    tryCatch(summary(crr(ftime = input_df[[time_col]], fstatus = input_df[[stat_col]], 
                         cov1 = input_df[[i]], cencode=cencode_char, failcode=failcode_char)), 
             error = function(e) { skip_to_next <<- TRUE})
    if (skip_to_next){
      fg_results[i,]=c("error","error","error","error","error","error","error","error","error")
      { next } 
    }
    crr_result=summary(crr(ftime = input_df[[time_col]], fstatus = input_df[[stat_col]], 
                           cov1 = input_df[[i]], cencode=cencode_char, failcode=failcode_char))
    fg_results[i,]=c(crr_result[["coef"]][5],round(crr_result[["conf.int"]][1],2),round(crr_result[["conf.int"]][3],2),
                     round(crr_result[["conf.int"]][4],2),crr_result[["n"]],nrow(subset(input_df,input_df[[i]]==1)),
                     nrow(subset(input_df,input_df[[i]]==0)),
                     round((nrow(subset(input_df,input_df[[i]]==1 & input_df[[stat_col]] == failcode_char))/nrow(subset(input_df,input_df[[i]]==1)))*100,0),
                     round((nrow(subset(input_df,input_df[[i]]==0 & input_df[[stat_col]] == failcode_char))/nrow(subset(input_df,input_df[[i]]==0)))*100,0))
    
  }
  return(fg_results)
}

TRD_FG_summary <- summarize_fg(clinmicrodata,time_col = 'Time.to.TRD', stat_col = "TRD",
                               beg=88, end=90, cencode_char = "Alive", failcode_char = "TRD")

Relapse_FG_summary <- summarize_fg(clinmicrodata,time_col = 'Time.to.relapse', stat_col = "Relapse",
                                   beg=88, end=90, cencode_char = "Alive", failcode_char = "Relapse")

###SUMMARIXE COX###
summarize_cox <-function(input_df,time_col,stat_col,beg,end){
  input_df_original<-input_df
  cox_results <- data.frame(matrix(nrow=end-beg+1,ncol=9))
  rownames(cox_results)<-colnames(input_df[beg:end])
  colnames(cox_results)<-c("P","HR","HR 2.5","HR 97.5","N","NGroup1","NGroup0","%Group1","%Group0")
  surv_object_survival <- Surv(time = input_df[[time_col]], 
                               event = as.numeric(input_df[[stat_col]]),type="right")
  if (stat_col=="Relapsezeroone"){
  surv_object_survival_ptremoved <- Surv(time = subset(input_df,patient!="VRM"&patient!="OBR")[[time_col]], 
                               event = as.numeric(subset(input_df,patient!="VRM"&patient!="OBR")[[stat_col]]),type="right")}
  for (i in colnames(input_df[beg:end])){
    if (i%in%c("DB.fdis.e30.class_bin","GCF.fdis.e30.class_bin","OM.fdis.e30.class_bin")&stat_col=="Relapsezeroone"){
      univ_formulas_survival <- sapply(i,
                                       function(x) as.formula(paste('surv_object_survival_ptremoved~', x)))
      input_df<-subset(input_df_original,patient!="VRM"&patient!="OBR")
    }else{
    input_df<-input_df_original
    univ_formulas_survival <- sapply(i,
                                     function(x) as.formula(paste('surv_object_survival~', x)))}
    univ_models_survival <- lapply(univ_formulas_survival, function(x){coxph(x, data = input_df)})
    cox_result <- lapply(univ_models_survival,
                         function(x){ 
                           return(summary(x))})
    
    cox_results[i,]=c(cox_result[[i]]$coefficients[5],round(cox_result[[i]]$conf.int[1],2),round(cox_result[[i]]$conf.int[3],2),
                      round(cox_result[[i]]$conf.int[4],2),cox_result[[i]][["n"]],nrow(subset(input_df,input_df[[i]]==1)),
                      nrow(subset(input_df,input_df[[i]]==0)),
                      round((nrow(subset(input_df,input_df[[i]]==1 & input_df[[stat_col]] == 1))/nrow(subset(input_df,input_df[[i]]==1)))*100,0),
                      round((nrow(subset(input_df,input_df[[i]]==0 & input_df[[stat_col]] == 1))/nrow(subset(input_df,input_df[[i]]==0)))*100,0))
    
  }
  return(cox_results)
}

OS_cox_summary <- summarize_cox(clinmicrodata,time_col = 'Time.to.death',
                                stat_col = 'Deathzeroone',beg=88, end=90)

PFS_cox_summary <- summarize_cox(clinmicrodata,time_col = 'Time.to.relapse',
                                 stat_col = 'Relapsezeroone',beg=88, end=90)

summarize_cox_tte <-function(input_df,time_col,beg,end){
  cox_results <- data.frame(matrix(nrow=end-beg+1,ncol=9))
  rownames(cox_results)<-colnames(input_df[beg:end])
  colnames(cox_results)<-c("P","HR","HR 2.5","HR 97.5","N","NGroup1","NGroup0","%Group1","%Group0")
  surv_object_survival <- Surv(time = input_df[[time_col]], 
                               event = rep(1,31),type="right")
  for (i in colnames(input_df[beg:end])){
      univ_formulas_survival <- sapply(i,
                                       function(x) as.formula(paste('surv_object_survival~', x)))
    univ_models_survival <- lapply(univ_formulas_survival, function(x){coxph(x, data = input_df)})
    cox_result <- lapply(univ_models_survival,
                         function(x){ 
                           return(summary(x))})
    
    cox_results[i,]=c(cox_result[[i]]$coefficients[5],round(cox_result[[i]]$conf.int[1],2),round(cox_result[[i]]$conf.int[3],2),
                      round(cox_result[[i]]$conf.int[4],2),cox_result[[i]][["n"]],nrow(subset(input_df,input_df[[i]]==1)),
                      nrow(subset(input_df,input_df[[i]]==0)),
                      100,
                      100)
    
  }
  return(cox_results)
}

TTE_cox_summary <- summarize_cox_tte(clinmicrodata,time_col = 'Time.to.engraftment',beg=85, end=87)

######HEATMAP WITH RESULTS###
Relapse_FG_summary$outcome<-"Relapse"
Relapse_FG_summary$site<-sub("_[^_]+$", "",(sub('.*\\. ?(\\w+)', '\\1', rownames(Relapse_FG_summary))))
Relapse_FG_summary$metric<-sub('\\.[^\\.]*$', '', rownames(Relapse_FG_summary))
TRD_FG_summary$outcome<-"TRD"
TRD_FG_summary$site<-sub("_[^_]+$", "",(sub('.*\\. ?(\\w+)', '\\1', rownames(TRD_FG_summary))))
TRD_FG_summary$metric<-sub('\\.[^\\.]*$', '', rownames(TRD_FG_summary))
OS_cox_summary$outcome<-"OS"
OS_cox_summary$site<-sub("_[^_]+$", "",(sub('.*\\. ?(\\w+)', '\\1', rownames(OS_cox_summary))))
OS_cox_summary$metric<-sub('\\.[^\\.]*$', '', rownames(OS_cox_summary))
PFS_cox_summary$outcome<-"PFS"
PFS_cox_summary$site<-sub("_[^_]+$", "",(sub('.*\\. ?(\\w+)', '\\1', rownames(PFS_cox_summary))))
PFS_cox_summary$metric<-sub('\\.[^\\.]*$', '', rownames(PFS_cox_summary))

outcomes_summary<-rbind(Relapse_FG_summary,TRD_FG_summary,PFS_cox_summary,OS_cox_summary)
rownames(outcomes_summary)<-seq(1,nrow(outcomes_summary),1)

outcomes_summary$signif.corrected<-ifelse(outcomes_summary$P*3<0.001,"***",
                                ifelse(outcomes_summary$P*3<0.01,"**",
                                       ifelse(outcomes_summary$P*3<0.05,"*","")))
outcomes_summary$Q<-ifelse(outcomes_summary$P*3<1,outcomes_summary$P*3,1)

ggplot(outcomes_summary,aes(x=outcome,y=metric))+geom_tile(aes(fill=HR),color="white")+
  geom_text(aes(label=round(Q,3)), size=5,color="#ff1573")+
  scale_fill_gradient2(low="#15d0ff", mid="white", high="#ffd015",trans="log", 
                       midpoint=0,limits=c(min(outcomes_summary$HR),1/min(outcomes_summary$HR)),
                       breaks=c(0.1,0.5,1,2,10)) +
  scale_y_discrete(labels=c("DB","GCF","OM"))+
  facet_wrap(.~"recovery")+coord_equal()+labs(x="Outcome",y="Site")+
  theme(panel.grid = element_blank())

write.csv(outcomes_summary, file = "outcomes_res_summary_070622.csv", row.names = F)

#cuminc and survival plots
plot_cuminc <- function(clin_data, site_phase, plottitle, pval, hr, hrl, hru, time_col, status_col, cen_code, categories = c("Low","High"), colores = c('red','dodgerblue3')){
  #fit cuminc
  cuminc_fit_relapse <- cuminc(ftime = clin_data[[time_col]], fstatus = clin_data[[status_col]],
                               group = clin_data[[site_phase]],cencode = cen_code)
  print(cuminc_fit_relapse)
  #organize data to plot with basic r
  cuminc_fit_relapse_plot_data<- 
    cuminc_fit_relapse %>%
    list_modify("Tests" = NULL) %>%
    map_df(`[`, c("time", "est"), .id = "id") %>%
    filter(id %in% c(paste(categories[1],status_col), paste(categories[2],status_col))) 
  
  cuminc_fit_relapse_plot_data[cuminc_fit_relapse_plot_data==paste(categories[1],status_col)]<-categories[1]
  cuminc_fit_relapse_plot_data[cuminc_fit_relapse_plot_data==paste(categories[2],status_col)]<-categories[2]
  
  print(cuminc_fit_relapse_plot_data,n=40)
  
  censors<-subset(clin_data,clin_data[[status_col]]==cen_code)
  censors$height<-NA
  for (i in 1:nrow(censors)){
    censors[i,c("height")]<-tail(subset(cuminc_fit_relapse_plot_data,id==censors[i,c(site_phase)]&time<=censors[i,c(time_col)])$est,1)
  }
  print(censors)
  #plot cuminc curves with basic r
  scramble_plot <- ggplot(palette = colores)+
    geom_step(data=cuminc_fit_relapse_plot_data, aes(x=time, y=est, color=id), lwd = 1.2) + 
    geom_point(data=censors, shape = 3,
               aes(x=.data[[time_col]],y=.data[["height"]],color=.data[[site_phase]]))+
    ylim(c(0,1))+
    theme_classic()+
    theme(plot.title = element_text(size = 14),
          legend.title = element_blank(),
          legend.position = c(0.1,0.65))+
    labs(x = "Months",
         y = paste("Cumulative incidence of",ifelse(status_col=="Relapse",'relapse',status_col)),
         title = plottitle)+
    annotate("text", x=0, y = 1, hjust = 0,
             label = paste0("HR = ",hr," (",hrl,"–",hru,"); P = ", pval), size = 3.2)+
    scale_x_continuous(breaks=seq(0,60,6),limits=c(0,62))+
    scale_color_manual(values = colores, labels = c(category1=categories[1],category2=categories[2]))
  
  #fit surv curves to get number at risk table
  scramble_fit <- survfit(as.formula(paste0("Surv(",time_col,", ifelse(",status_col," != 2,1,0)) ~ ", site_phase)),subset(clinmicrodata,patient!="OBR"&patient!="VRM"))
  #plot number at risk table
  numatrisk_plot <- ggsurvplot(
    scramble_fit, palette = colores,
    risk.table = TRUE,
    risk.table.y.text = FALSE,
    ylab = "Months",
    tables.theme = theme_cleantable(),
    risk.table.fontsize = 3,
    break.time.by = 6,
    xlim=c(0,62)
  )
  #plot cuminc curves and number at risk table arranged in a single plot
  final_plot <- cowplot::plot_grid(scramble_plot+theme(plot.title = element_text(size=11), axis.title = element_text(size = 10)),
                                   numatrisk_plot$table+theme(plot.title = element_text(size=10)),nrow=2,rel_heights = c(4,1),align='v',axis='b')
  return (final_plot)
}

FS13L<-plot_cuminc(clinmicrodata,
            site_phase <- "DB.fdis.e30.class", plottitle = "SB recovery", colores = c('#616161', '#008d78'),
            pval = 0.60,hr=0.64,hru = 3.50,hrl = 0.12,categories = c("NR","R"),
            time_col <- "Time.to.TRD",status_col <- "TRD", cen_code = "Alive")
FS13J<-plot_cuminc(clinmicrodata,
            site_phase <- "GCF.fdis.e30.class", plottitle = "GCF recovery", colores = c('#616161', '#008d78'),
            pval = 0.69,hr=0.71,hru = 3.66,hrl = 0.14,categories = c("NR","R"),
            time_col <- "Time.to.TRD",status_col <- "TRD", cen_code = "Alive")
FS13K<-plot_cuminc(clinmicrodata,
            site_phase <- "OM.fdis.e30.class", plottitle = "OM recovery", colores = c('#616161', '#008d78'),
            pval = 0.05,hr=0.19,hru = 1,hrl = 0.04,categories = c("NR","R"),
            time_col <- "Time.to.TRD",status_col <- "TRD", cen_code = "Alive")

#ALERTA: CASO PRECISE EXCLUIR PACIENTE PRECISA EXCLUIR TAMBÉM NO PLOT_CUMINC NA PARTE DE GERAR O GRÁFICO DE NUMBER AT RISK
FS13I<-plot_cuminc(subset(clinmicrodata,patient!="OBR"&patient!="VRM"),
            site_phase <- "DB.fdis.e30.class", plottitle = "SB recovery", colores = c('#616161', '#008d78'),
            pval = 0.56,hr=1.54,hru = 6.66,hrl = 0.36,categories = c("NR","R"),
            time_col <- "Time.to.relapse",status_col <- "Relapse", cen_code = "Alive")
FS13G<-plot_cuminc(subset(clinmicrodata,patient!="OBR"&patient!="VRM"),
            site_phase <- "GCF.fdis.e30.class", plottitle = "GCF recovery", colores = c('#616161', '#008d78'),
            pval = 0.52,hr=0.66,hru = 2.36,hrl = 0.18,categories = c("NR","R"),
            time_col <- "Time.to.relapse",status_col <- "Relapse", cen_code = "Alive")
FS13H<-plot_cuminc(subset(clinmicrodata,patient!="OBR"&patient!="VRM"),
            site_phase <- "OM.fdis.e30.class", plottitle = "OM recovery", colores = c('#616161', '#008d78'),
            pval = 0.011,hr=0.2,hru = 0.69,hrl = 0.06,categories = c("NR","R"),
            time_col <- "Time.to.relapse",status_col <- "Relapse", cen_code = "Alive")


plot_survival <- function(clindata, time_col, status_col, site_phase, plottitle, colores = c("red","dodgerblue3"), ylabel,
                          pval, hr, hrl, hru,categories=c("NR","R")){
  #drop na rows at diversity measure
  clindata <- clindata[!is.na(clindata[[site_phase]]),]
  #building survival object
  surv_object <- Surv(time = clindata[[time_col]], event = as.numeric(clindata[[status_col]]),type="right")
  #fit km curves
  surv_fit <- surv_fit(formula = as.formula(paste0("surv_object ~", site_phase)), data = clindata)
  
  #plotting survival curves
  p_surv <- ggsurvplot(surv_fit, data = clindata, size = 1, palette = colores,
                       conf.int = FALSE, pval = F, risk.table = TRUE, xlab = "Months", ylab = ylabel,
                       risk.table.heights = 0.25,
                       ggtheme = theme_bw(), 
                       break.time.by = 6, xlim=c(0,62), 
                       risk.table.y.text = FALSE,risk.table.title="Number at risk",
                       fontsize = 3, tables.theme = clean_theme()) #+
  #annotate("text", x=0, y = 1, hjust = 0,
  #         label = paste0("HR = ",hr," (",hrl,"-",hru,"); P = ", pval))
  
  final_plot <- cowplot::plot_grid(p_surv$plot+annotate("text", x=1, y = 0.05, hjust = 0,label = paste0("HR = ",hr," (",hrl,"–",hru,"); P = ", pval), size = 3.2)+
                                     labs(title=plottitle)+theme(axis.line = element_line(colour = "black"),plot.title = element_text(size=11),
                                                                 panel.grid.major = element_blank(), axis.title = element_text(size = 10),
                                                                 panel.grid.minor = element_blank(),
                                                                 panel.border = element_blank(),
                                                                 panel.background = element_blank(), legend.position = c(0.1,0.35),legend.title = element_blank())+
                                     scale_color_manual(values=colores,labels=categories),
                                   p_surv$table+theme(plot.title = element_text(size=10),legend.position = "none", panel.grid = element_blank(), axis.title = element_blank(),
                                                      panel.border = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank()),
                                   nrow=2,rel_heights = c(4,1),align='v',axis='b')
  return (final_plot)
}

FS13C<-plot_survival(clinmicrodata,time_col<-"Time.to.death",status_col<-"Deathzeroone",
              ylabel="Overall survival",pval=0.83,hr=0.87,hrl=0.24,hru=3.14,colores=c('#616161', '#008d78'),
              site_phase <- "DB.fdis.e30.class", plottitle="SB recovery")
FS13A<-plot_survival(clinmicrodata,time_col<-"Time.to.death",status_col<-"Deathzeroone",
              ylabel="Overall survival",pval=0.09,hr=0.39,hrl=0.13,hru=1.16,colores=c('#616161', '#008d78'),
              site_phase <- "GCF.fdis.e30.class", plottitle="GCF recovery")
FS13B<-plot_survival(clinmicrodata,time_col<-"Time.to.death",status_col<-"Deathzeroone",
              ylabel="Overall survival",pval=0.002,hr=0.17,hrl=0.05,hru=0.52,colores=c('#616161', '#008d78'),
              site_phase <- "OM.fdis.e30.class", plottitle="OM recovery")

FS13F<-plot_survival(subset(clinmicrodata,patient!="OBR"&patient!="VRM"),time_col<-"Time.to.relapse",status_col<-"Relapsezeroone",
              ylabel="Progression-free survival",pval=0.86,hr=1.15,hrl=0.23,hru=5.83,colores=c('#616161', '#008d78'),
              site_phase <- "DB.fdis.e30.class", plottitle="SB recovery")
FS13D<-plot_survival(subset(clinmicrodata,patient!="OBR"&patient!="VRM"),time_col<-"Time.to.relapse",status_col<-"Relapsezeroone",
              ylabel="Progression-free survival",pval=0.32,hr=0.51,hrl=0.13,hru=1.93,colores=c('#616161', '#008d78'),
              site_phase <- "GCF.fdis.e30.class", plottitle="GCF recovery")
FS13E<-plot_survival(subset(clinmicrodata,patient!="OBR"&patient!="VRM"),time_col<-"Time.to.relapse",status_col<-"Relapsezeroone",
              ylabel="Progression-free survival",pval=0.001,hr=0.06,hrl=0.01,hru=0.34,colores=c('#616161', '#008d78'),
              site_phase <- "OM.fdis.e30.class", plottitle="OM recovery")


####05/07/22 - associação entre parametros clínicos e desfechos
clindata_bin<-read.csv("Files/DadosClinicosTMO_09.05.2021.csv")
clindata_bin$Patient<-NULL
clindata_bin$patient<-clindata_bin$Patient.acronym
clindata_bin$Patient.acronym<-NULL
atbdata<-read.csv("Files/300622_ATB_classification.csv")
clindata_bin<-merge(clindata_bin,atbdata)

TRD_FG_summary_clin <- summarize_fg(clindata_bin,time_col = 'Time.to.TRD', stat_col = "TRD",
                               beg=38, end=55, cencode_char = "Alive", failcode_char = "TRD")

Relapse_FG_summary_clin <- summarize_fg(clindata_bin,time_col = 'Time.to.relapse', stat_col = "Relapse",
                                   beg=38, end=55, cencode_char = "Alive", failcode_char = "Relapse")

OS_cox_summary_clin <- summarize_cox(clindata_bin,time_col = 'Time.to.death',
                                stat_col = 'Deathzeroone',beg=38, end=55)

PFS_cox_summary_clin <- summarize_cox(clindata_bin,time_col = 'Time.to.relapse',
                                 stat_col = 'Relapsezeroone',beg=38, end=55)

####################análise multivariada####################
#juntar dados clinicos binarios com dados de microbiota binarios
microbiota<-read.csv("Files/300522_microbiota-variables.csv")
clindata_microbiota<-merge(clindata_bin,microbiota)

#OS
summary(coxph(formula = Surv(time = clindata_microbiota[['Time.to.death']], 
                                    event = as.numeric(clindata_microbiota[['Deathzeroone']]),type="right")~
           OM.fdis.e30.class_bin+Disease.risk.index_High_bin+Conditioning.intensity_Myeloablative_bin+dot, 
           data = clindata_microbiota))[c("n","coefficients","conf.int")]
# $n
# [1] 29
# 
# $coefficients
# coef  exp(coef)   se(coef)         z     Pr(>|z|)
# OM.fdis.e30.class_bin                    -2.45501100 0.08586225 0.71518160 -3.432710 0.0005975809
# Disease.risk.index_High_bin               1.84747462 6.34377882 0.70942057  2.604202 0.0092088403
# Conditioning.intensity_Myeloablative_bin -1.83355714 0.15984397 0.87103230 -2.105039 0.0352878959
# dot                                       0.04162637 1.04250489 0.02083646  1.997766 0.0457420840
# 
# $conf.int
# exp(coef) exp(-coef)  lower .95  upper .95
# OM.fdis.e30.class_bin                    0.08586225 11.6465617 0.02113677  0.3487915
# Disease.risk.index_High_bin              6.34377882  0.1576348 1.57938552 25.4804981
# Conditioning.intensity_Myeloablative_bin 0.15984397  6.2561009 0.02899160  0.8812929
# dot                                      1.04250489  0.9592281 1.00078796  1.0859608


#PFS
summary(coxph(formula = Surv(time = subset(clindata_microbiota,patient!="VRM"&patient!="OBR")[['Time.to.relapse']], 
                             event = as.numeric(subset(clindata_microbiota,patient!="VRM"&patient!="OBR")[['Relapsezeroone']]),type="right")~
                OM.fdis.e30.class_bin+Disease.risk.index_High_bin, 
              data = subset(clindata_microbiota,patient!="VRM"&patient!="OBR")))[c("n","coefficients","conf.int")]
# $n
# [1] 27
# 
# $coefficients
# coef  exp(coef)  se(coef)         z    Pr(>|z|)
# OM.fdis.e30.class_bin       -2.366171 0.09383935 0.8465167 -2.795185 0.005186997
# Disease.risk.index_High_bin  1.284852 3.61413370 0.7355965  1.746681 0.080692753
# 
# $conf.int
# exp(coef) exp(-coef)  lower .95  upper .95
# OM.fdis.e30.class_bin       0.09383935 10.6565108 0.01785783  0.4931071
# Disease.risk.index_High_bin 3.61413370  0.2766915 0.85479766 15.2807653

#Relapse
summary(crr(ftime = subset(clindata_microbiota,patient!="VRM"&patient!="OBR")$Time.to.relapse, 
                  fstatus = subset(clindata_microbiota,patient!="VRM"&patient!="OBR")$Relapse, 
                  cov1 = subset(clindata_microbiota,patient!="VRM"&patient!="OBR")[c("OM.fdis.e30.class_bin","Disease.risk.index_High_bin")], 
              cencode='Alive', failcode='Relapse'))[c("n","coef","conf.int")]
# $n
# [1] 27
# 
# $coef
# coef exp(coef)  se(coef)         z p-value
# OM.fdis.e30.class_bin       -1.683949 0.1856394 0.5561067 -3.028105  0.0025
# Disease.risk.index_High_bin  1.546302 4.6940816 0.5923147  2.610610  0.0090
# 
# $conf.int
# exp(coef) exp(-coef)       2.5%      97.5%
#   OM.fdis.e30.class_bin       0.1856394  5.3867888 0.06241819  0.5521142
# Disease.risk.index_High_bin 4.6940816  0.2130342 1.47018319 14.9875217

#TRD - just in case
summary(crr(ftime = clindata_microbiota$Time.to.TRD, 
            fstatus = clindata_microbiota$TRD, 
            cov1 = clindata_microbiota[c("OM.fdis.e30.class_bin","dot")], 
            cencode='Alive', failcode='TRD'))[c("n","coef","conf.int")]

# $n
# [1] 29
# 
# $coef
# coef exp(coef)   se(coef)         z p-value
# OM.fdis.e30.class_bin -1.34476343 0.2606014 0.93989182 -1.430764   0.150
# dot                    0.04671518 1.0478235 0.01994218  2.342531   0.019
# 
# $conf.int
# exp(coef) exp(-coef)       2.5%    97.5%
#   OM.fdis.e30.class_bin 0.2606014  3.8372786 0.04129898 1.644425
# dot                   1.0478235  0.9543592 1.00765840 1.089590

#forest plot for these associations
F7F<-ggplot(data=data.frame(var=c("OM recovery (R)","DRI (High)","Cond Int (Myeloablative)", "DOT"),
                       hr=c(0.08586225, 6.34377882, 0.15984397, 1.04250489),
                       hrl=c(0.02113677, 1.57938552, 0.02899160, 1.00078796),
                       hru=c(0.3487915, 25.4804981, 0.8812929, 1.0859608),
                       p=c(0.0005975809, 0.0092088403, 0.0352878959, 0.0457420840)), 
       aes(y=var, x=hr, xmin=hrl, xmax=hru, label=ifelse(p<0.001,"p < 0.001",paste("p =",round(p,3))),color = var))+
  geom_point(shape = 18, size = 2.5)+geom_errorbarh(height=.1)+scale_x_log10(breaks=c(0.04,0.2,1,5,25),limits=c(0.017,59), expand = c(0,0))+
  labs(x = "HR",title = "Overall survival")+
  scale_color_manual(values = c('#616161','#616161','#616161', '#008d78'))+
  geom_vline(xintercept = 1, linetype="longdash")+
  geom_label(nudge_y = 0.3, size=3, color = 'black',label.padding = unit(0, "lines"), label.size = 0)+
  theme(axis.title.y = element_blank(), panel.grid = element_blank(), legend.position = 'none',
        axis.title = element_text(size = 10), plot.title = element_text(size = 11))

F7G<-ggplot(data=data.frame(var=c("OM recovery (R)","DRI (High)"),
                       hr=c(0.09383935,3.61413370),
                       hrl=c(0.01785783,0.85479766),
                       hru=c(0.4931071,15.2807653),
                       p=c(0.005186997,0.080692753)), 
       aes(y=var, x=hr, xmin=hrl, xmax=hru, label=ifelse(p<0.001,"p < 0.001",paste("p =",round(p,3))),color = var))+
  geom_point(shape = 18, size = 2.5)+geom_errorbarh(height=.1)+scale_x_log10(breaks=c(0.04,0.2,1,5,25),limits=c(0.017,59), expand = c(0,0))+
  labs(x = "HR",title = "Progression-free survival")+
    scale_color_manual(values = c('#616161', '#008d78'))+
    geom_vline(xintercept = 1, linetype="longdash")+
    geom_label(nudge_y = 0.3, size=3, color = 'black',label.padding = unit(0, "lines"), label.size = 0)+
    theme(axis.title.y = element_blank(), panel.grid = element_blank(), legend.position = 'none',
          axis.title = element_text(size = 10), plot.title = element_text(size = 11))

F7H<-ggplot(data=data.frame(var=c("OM recovery (R)","DRI (High)"),
                       hr=c(0.1856394, 4.6940816),
                       hrl=c(0.06241819,1.47018319),
                       hru=c(0.5521142,14.9875217),
                       p=c(0.0025,0.0090)), 
            aes(y=var, x=hr, xmin=hrl, xmax=hru, label=ifelse(p<0.001,"p < 0.001",paste("p =",round(p,3))),color = var))+
  geom_point(shape = 18, size = 2.5)+geom_errorbarh(height=.1)+scale_x_log10(breaks=c(0.04,0.2,1,5,25),limits=c(0.017,59), expand = c(0,0))+
  labs(x = "HR",title = "Risk of relapse")+
  scale_color_manual(values = c('#616161', '#008d78'))+
  geom_vline(xintercept = 1, linetype="longdash")+
  geom_label(nudge_y = 0.3, size=3, color = 'black',label.padding = unit(0, "lines"), label.size = 0)+
  theme(axis.title.y = element_blank(), panel.grid = element_blank(), legend.position = 'none',
        axis.title = element_text(size = 10), plot.title = element_text(size = 11))
