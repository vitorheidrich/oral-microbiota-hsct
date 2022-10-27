library(dplyr)
library(nloptr)
library(ggplot2)
theme_set(theme_bw())
library(phyloseq)
library(qiime2R)
library(ggpubr)
library(tidyverse)
library(cowplot)
library(ANCOMBC)
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
coherent_asv_names<-function(physeq){
  asvs_list<-paste0("ASV",1:length(as.data.frame(physeq@tax_table)[,1]),"_",
                    as.data.frame(physeq@tax_table)['Species'][1:length(as.data.frame(physeq@tax_table)[,1]),])
  
  new_tax_table<-as.data.frame(physeq@tax_table)
  new_tax_table$ASV<-asvs_list
  tax_table(physeq)<-as.matrix(new_tax_table)
  return(physeq)}

#################################################
physeq<-coherent_asv_names(coherent_taxa_names(qza_to_phyloseq(
  features="Files/ASV_table_final_goodsamples.qza",
  tree="Files/rooted-tree.qza",
  "Files/taxonomy_vs.qza",
  metadata = "Files/mapping_tmo_0.tsv")))

physeq_tmo<-subset_samples(physeq,time!="MI"&time!="MF")

##################################################
plot_ancombc<-function(physeq,rank,group,revert_order=F, q_cutoff = 0.05, beta_cutoff = 0, color="black",x_lim = 5, facet_title = "ANCOM-BC"){
  physeq<-coherent_taxa_names(physeq)
  if (rank!="ASV"){
    physeq<-tax_glom(physeq, taxrank = rank)
  }
  #print(physeq)
  out = ancombc(phyloseq = physeq, 
                formula = group, p_adj_method = "holm", zero_cut = 0.75, lib_cut = 0, 
                group = group, struc_zero = TRUE, neg_lb = TRUE,
                max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
  res<-as.data.frame(out$res)
  colnames(res)<-names(out$res)
  res$taxa<-subset(as.data.frame(physeq@tax_table@.Data),
                   rownames(as.data.frame(physeq@tax_table@.Data))%in%rownames(res))[rownames(res), ][[rank]]
  res$star=ifelse(res$q_val==0, "z",
                  ifelse(res$q_val<.001, "***", 
                         ifelse(res$q_val<.01, "**",
                                ifelse(res$q_val<.05, "*", ""))))
  groups<-sort(unique(as.data.frame(sample_data(physeq))[[group]]))
  if(groups[1]=="BF"){groups[1]<-"SB"}
  if(groups[2]=="BF"){groups[2]<-"SB"}
  if(groups[1]=="FCG"){groups[1]<-"GCF"}
  if(groups[2]=="FCG"){groups[2]<-"GCF"}
  if(groups[1]=="MO"){groups[1]<-"OM"}
  if(groups[2]=="MO"){groups[2]<-"OM"}
  
  if(revert_order==T){
    groups<-sort(groups,decreasing = T)
    res$beta<--1*res$beta
  }
  
  res<-subset(res,q_val<=q_cutoff&abs(beta)>=beta)
  if (nrow(res)!=0){
    p <- ggplot(res,aes(x=beta,y=reorder(taxa,beta),color = factor(sign(beta))))+
      scale_colour_manual(values = c("-1"=color[1],"1"=color[2]))+guides(color="none")+
      geom_vline(xintercept = 0, linetype="dashed")+
      geom_errorbar(aes(xmin=beta-se, xmax=beta+se), width=.15)+
      geom_point(shape=18,size=4.5)+
      geom_text(data=subset(res,q_val<.25),aes(x=ifelse(beta>0,beta+se+0.1*(beta+se),beta-se-0.1*(abs(beta)+se)),
                                               y=reorder(taxa,beta), label = star),color='black',size=4)+
      scale_x_continuous(limits = c(-x_lim,x_lim))+
      labs(y=rank, x=expr(paste(!!groups[1] %<-%  'ANCOM-BC (',beta,')'  %->% !!groups[2])))+
      facet_wrap(.~eval(facet_title))+
      theme(axis.title.y = element_blank(), panel.grid.minor = element_blank(),
            strip.background = element_rect(color="white", fill="white", linetype="solid"))
    
    return(p)
  }else{return("no significant taxa found")}
}

#ANCOM-BC - compare phases with P separately per site
abc_bd_a<-plot_ancombc(subset_samples(physeq_tmo,time%in%c("P","A")&site=="BF"), 
                       rank='Genus', group = 'time', color = c(site_colors[1],site_colors[1]), x_lim = 7.5, revert_order = T, facet_title = "SB")
abc_bd_e<-plot_ancombc(subset_samples(physeq_tmo,time%in%c("P","E")&site=="BF"), 
                       rank='Genus', group = 'time', color = c(site_colors[1],site_colors[1]), x_lim = 7.5, revert_order = T, facet_title = "SB")
abc_bd_e30<-plot_ancombc(subset_samples(physeq_tmo,time%in%c("P","E30")&site=="BF"), 
                         rank='Genus', group = 'time', color = c(site_colors[1],site_colors[1]), x_lim = 7.5, revert_order = T, facet_title = "SB")
abc_bd_e75<-plot_ancombc(subset_samples(physeq_tmo,time%in%c("P","E75")&site=="BF"), 
                         rank='Genus', group = 'time', color = c(site_colors[1],site_colors[1]), x_lim = 7.5, revert_order = T, facet_title = "SB")
F3Cbf<-plot_grid(abc_bd_a, abc_bd_e, abc_bd_e30, abc_bd_e75, ncol = 1, nrow = 4, align = "v", axis = "lr", rel_heights = c(5,18,3,3)+4)

abc_fcg_a<-plot_ancombc(subset_samples(physeq_tmo,time%in%c("P","A")&site=="FCG"), 
                        rank='Genus', group = 'time', color = c(site_colors[2],site_colors[2]), x_lim = 6, revert_order = T, facet_title = "GCF")
abc_fcg_e<-plot_ancombc(subset_samples(physeq_tmo,time%in%c("P","E")&site=="FCG"), 
                        rank='Genus', group = 'time', color = c(site_colors[2],site_colors[2]), x_lim = 6, revert_order = T, facet_title = "GCF")
abc_fcg_e30<-plot_ancombc(subset_samples(physeq_tmo,time%in%c("P","E30")&site=="FCG"), 
                          rank='Genus', group = 'time', color = c(site_colors[2],site_colors[2]), x_lim = 6, revert_order = T, facet_title = "GCF")
abc_fcg_e75<-plot_ancombc(subset_samples(physeq_tmo,time%in%c("P","E75")&site=="FCG"), 
                          rank='Genus', group = 'time', color = c(site_colors[2],site_colors[2]), x_lim = 6, revert_order = T, facet_title = "GCF")

abc_mo_a<-plot_ancombc(subset_samples(physeq_tmo,time%in%c("P","A")&site=="MO"), 
                       rank='Genus', group = 'time', color = c(site_colors[3],site_colors[3]), x_lim = 8, revert_order = T, facet_title = "OM")
abc_mo_e<-plot_ancombc(subset_samples(physeq_tmo,time%in%c("P","E")&site=="MO"), 
                       rank='Genus', group = 'time', color = c(site_colors[3],site_colors[3]), x_lim = 8, revert_order = T, facet_title = "OM")
abc_mo_e30<-plot_ancombc(subset_samples(physeq_tmo,time%in%c("P","E30")&site=="MO"), 
                         rank='Genus', group = 'time', color = c(site_colors[3],site_colors[3]), x_lim = 8, revert_order = T, facet_title = "OM")
abc_mo_e75<-plot_ancombc(subset_samples(physeq_tmo,time%in%c("P","E75")&site=="MO"), 
                         rank='Genus', group = 'time', color = c(site_colors[3],site_colors[3]), x_lim = 8, revert_order = T, facet_title = "OM")


F3C<-plot_grid(abc_fcg_a+theme(axis.title.x = element_text(size = 10)), abc_mo_a+theme(axis.title.x = element_text(size = 10)), abc_bd_a+theme(axis.title.x = element_text(size = 10)),
               abc_fcg_e+theme(axis.title.x = element_text(size = 10)), abc_mo_e+theme(axis.title.x = element_text(size = 10)), abc_bd_e+theme(axis.title.x = element_text(size = 10)),
               abc_fcg_e30+theme(axis.title.x = element_text(size = 10)), abc_mo_e30+theme(axis.title.x = element_text(size = 10)), abc_bd_e30+theme(axis.title.x = element_text(size = 10)),
               abc_fcg_e75+theme(axis.title.x = element_text(size = 10)), abc_mo_e75+theme(axis.title.x = element_text(size = 10)), abc_bd_e75+theme(axis.title.x = element_text(size = 10)),
               ncol = 3, nrow = 4, align = "vh", axis = "lrbt", rel_heights = c(5,18,5,3)+4)

#ANCOM-BC - compare sites pairwise in each phase
abc_p_bf_fcg<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("BF","FCG")&time=="P"), 
                       rank='Genus', group = 'site', color = c(site_colors[3],site_colors[1]), x_lim = 6, revert_order = F, facet_title = "P")
abc_p_bf_mo<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("BF","MO")&time=="P"), 
                           rank='Genus', group = 'site', color = c(site_colors[3],site_colors[2]), x_lim = 6, revert_order = F, facet_title = "P")
abc_p_fcg_mo<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("FCG","MO")&time=="P"), 
                          rank='Genus', group = 'site', color = c(site_colors[1],site_colors[2]), x_lim = 6, revert_order = F, facet_title = "P")
abc_a_bf_fcg<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("BF","FCG")&time=="A"), 
                           rank='Genus', group = 'site', color = c(site_colors[3],site_colors[1]), x_lim = 6, revert_order = F, facet_title = "A")
abc_a_bf_mo<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("BF","MO")&time=="A"), 
                          rank='Genus', group = 'site', color = c(site_colors[3],site_colors[2]), x_lim = 6, revert_order = F, facet_title = "A")
abc_a_fcg_mo<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("FCG","MO")&time=="A"), 
                           rank='Genus', group = 'site', color = c(site_colors[1],site_colors[2]), x_lim = 6, revert_order = F, facet_title = "A")
abc_e_bf_fcg<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("BF","FCG")&time=="E"), 
                           rank='Genus', group = 'site', color = c(site_colors[3],site_colors[1]), x_lim = 6, revert_order = F, facet_title = "E")
abc_e_bf_mo<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("BF","MO")&time=="E"), 
                          rank='Genus', group = 'site', color = c(site_colors[3],site_colors[2]), x_lim = 6, revert_order = F, facet_title = "E")
abc_e_fcg_mo<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("FCG","MO")&time=="E"), 
                           rank='Genus', group = 'site', color = c(site_colors[1],site_colors[2]), x_lim = 6, revert_order = F, facet_title = "E")
abc_e30_bf_fcg<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("BF","FCG")&time=="E30"), 
                           rank='Genus', group = 'site', color = c(site_colors[3],site_colors[1]), x_lim = 6, revert_order = F, facet_title = "E30")
abc_e30_bf_mo<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("BF","MO")&time=="E30"), 
                          rank='Genus', group = 'site', color = c(site_colors[3],site_colors[2]), x_lim = 6, revert_order = F, facet_title = "E30")
abc_e30_fcg_mo<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("FCG","MO")&time=="E30"), 
                           rank='Genus', group = 'site', color = c(site_colors[1],site_colors[2]), x_lim = 6, revert_order = F, facet_title = "E30")
abc_e75_bf_fcg<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("BF","FCG")&time=="E75"), 
                           rank='Genus', group = 'site', color = c(site_colors[3],site_colors[1]), x_lim = 6, revert_order = F, facet_title = "E75")
abc_e75_bf_mo<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("BF","MO")&time=="E75"), 
                          rank='Genus', group = 'site', color = c(site_colors[3],site_colors[2]), x_lim = 6, revert_order = F, facet_title = "E75")
abc_e75_fcg_mo<-plot_ancombc(subset_samples(physeq_tmo,site%in%c("FCG","MO")&time=="E75"), 
                           rank='Genus', group = 'site', color = c(site_colors[1],site_colors[2]), x_lim = 6, revert_order = F, facet_title = "E75")

FS5<-plot_grid(abc_p_fcg_mo+theme(axis.title.x = element_text(size = 10)),abc_p_bf_fcg+theme(axis.title.x = element_text(size = 10)),abc_p_bf_mo+theme(axis.title.x = element_text(size = 10)),
               abc_a_fcg_mo+theme(axis.title.x = element_text(size = 10)),abc_a_bf_fcg+theme(axis.title.x = element_text(size = 10)),abc_a_bf_mo+theme(axis.title.x = element_text(size = 10)),
               abc_e_fcg_mo+theme(axis.title.x = element_text(size = 10)),abc_e_bf_fcg+theme(axis.title.x = element_text(size = 10)),abc_e_bf_mo+theme(axis.title.x = element_text(size = 10)),
               abc_e30_fcg_mo+theme(axis.title.x = element_text(size = 10)),abc_e30_bf_fcg+theme(axis.title.x = element_text(size = 10)),abc_e30_bf_mo+theme(axis.title.x = element_text(size = 10)),
               abc_e75_fcg_mo+theme(axis.title.x = element_text(size = 10)),abc_e75_bf_fcg+theme(axis.title.x = element_text(size = 10)),abc_e75_bf_mo+theme(axis.title.x = element_text(size = 10)),
          align = "vh", axis = "lrbt", ncol = 3, nrow = 5, rel_heights = c(11,2,1,2,4)+6)

#N of differential taxa for a given comparison
n_diff_ancom_bc<-plot_ancombc<-function(physeq,rank,group, q_cutoff = 0.05, beta_cutoff = 0){
  physeq<-coherent_taxa_names(physeq)
  if (rank!="ASV"){
    physeq<-tax_glom(physeq, taxrank = rank)
  }
  #print(physeq)
  out = ancombc(phyloseq = physeq, 
                formula = group, p_adj_method = "holm", zero_cut = 0.75, lib_cut = 0, 
                group = group, struc_zero = TRUE, neg_lb = TRUE,
                max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE)
  res<-as.data.frame(out$res)
  colnames(res)<-names(out$res)
  res$taxa<-subset(as.data.frame(physeq@tax_table@.Data),
                   rownames(as.data.frame(physeq@tax_table@.Data))%in%rownames(res))[rownames(res), ][[rank]]
  res$star=ifelse(res$q_val==0, "z",
                  ifelse(res$q_val<.001, "***", 
                         ifelse(res$q_val<.01, "**",
                                ifelse(res$q_val<.05, "*", ""))))
  groups<-sort(unique(as.data.frame(sample_data(physeq))[[group]]))
  if(groups[1]=="BF"){groups[1]<-"DB"}
  if(groups[2]=="BF"){groups[2]<-"DB"}
  if(groups[1]=="FCG"){groups[1]<-"GCF"}
  if(groups[2]=="FCG"){groups[2]<-"GCF"}
  if(groups[1]=="MO"){groups[1]<-"OM"}
  if(groups[2]=="MO"){groups[2]<-"OM"}
  
  res<-subset(res,q_val<=q_cutoff&abs(beta)>=beta)
  #print(res)
  return(nrow(res))
}

abc_bd_a_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,time%in%c("P","A")&site=="BF"), rank='Genus', group = 'time')
abc_bd_e_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,time%in%c("P","E")&site=="BF"), rank='Genus', group = 'time')
abc_bd_e30_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,time%in%c("P","E30")&site=="BF"), rank='Genus', group = 'time')
abc_bd_e75_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,time%in%c("P","E75")&site=="BF"), rank='Genus', group = 'time')

abc_gcf_a_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,time%in%c("P","A")&site=="FCG"), rank='Genus', group = 'time')
abc_gcf_e_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,time%in%c("P","E")&site=="FCG"), rank='Genus', group = 'time')
abc_gcf_e30_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,time%in%c("P","E30")&site=="FCG"), rank='Genus', group = 'time')
abc_gcf_e75_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,time%in%c("P","E75")&site=="FCG"), rank='Genus', group = 'time')

abc_om_a_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,time%in%c("P","A")&site=="MO"), rank='Genus', group = 'time')
abc_om_e_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,time%in%c("P","E")&site=="MO"), rank='Genus', group = 'time')
abc_om_e30_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,time%in%c("P","E30")&site=="MO"), rank='Genus', group = 'time')
abc_om_e75_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,time%in%c("P","E75")&site=="MO"), rank='Genus', group = 'time')

#ANCOM-BC - n of differential taxa per phase

abc_df<-data.frame(site=c(rep("SB",4),rep("GCF",4),rep("OM",4)),
                   time=c(rep(c("P vs A","P vs E","P vs E30","P vs E75"),3)),
                   n=c(abc_bd_a_n,abc_bd_e_n,abc_bd_e30_n,abc_bd_e75_n,
                       abc_gcf_a_n,abc_gcf_e_n,abc_gcf_e30_n,abc_gcf_e75_n,
                       abc_om_a_n,abc_om_e_n,abc_om_e30_n,abc_om_e75_n))

FS8<-ggplot(abc_df,aes(y=n,x=time,color=site,group=site))+
  geom_point(size=3.2, alpha = 0.75)+geom_line(alpha=0.75)+
  scale_color_manual(values=site_colors)+
  labs(x="Comparison",y="N differentially abundant genera", color = "Site")+
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.line = element_line(color = 'black'))+
  scale_y_continuous(limits=c(0,20),expand = c(0,0))
FS8

#ANCOM-BC - n of differential taxa in pairwise comparisons between sites
abc_p_bf_fcg_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("BF","FCG")&time=="P"),rank='Genus', group = "site")
abc_p_bf_mo_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("BF","MO")&time=="P"),rank='Genus', group = "site")
abc_p_fcg_mo_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("FCG","MO")&time=="P"),rank='Genus', group = "site")
abc_a_bf_fcg_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("BF","FCG")&time=="A"),rank='Genus', group = "site")
abc_a_bf_mo_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("BF","MO")&time=="A"),rank='Genus', group = "site")
abc_a_fcg_mo_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("FCG","MO")&time=="A"),rank='Genus', group = "site")
abc_e_bf_fcg_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("BF","FCG")&time=="E"),rank='Genus', group = "site")
abc_e_bf_mo_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("BF","MO")&time=="E"),rank='Genus', group = "site")
abc_e_fcg_mo_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("FCG","MO")&time=="E"),rank='Genus', group = "site")
abc_e30_bf_fcg_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("BF","FCG")&time=="E30"),rank='Genus', group = "site")
abc_e30_bf_mo_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("BF","MO")&time=="E30"),rank='Genus', group = "site")
abc_e30_fcg_mo_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("FCG","MO")&time=="E30"),rank='Genus', group = "site")
abc_e75_bf_fcg_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("BF","FCG")&time=="E75"),rank='Genus', group = "site")
abc_e75_bf_mo_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("BF","MO")&time=="E75"),rank='Genus', group = "site")
abc_e75_fcg_mo_n<-n_diff_ancom_bc(subset_samples(physeq_tmo,site%in%c("FCG","MO")&time=="E75"),rank='Genus', group = "site")

abc_df_site<-data.frame(comparison=c(rep("GCF vs SB",5),rep("OM vs SB",5),rep("GCF vs OM",5)),
                   time=c(rep(c('P','A','E','E30','E75'),3)),
                   n=c(abc_p_bf_fcg_n,abc_a_bf_fcg_n,abc_e_bf_fcg_n,abc_e30_bf_fcg_n,abc_e75_bf_fcg_n,
                       abc_p_bf_mo_n,abc_a_bf_mo_n,abc_e_bf_mo_n,abc_e30_bf_mo_n,abc_e75_bf_mo_n,
                       abc_p_fcg_mo_n,abc_a_fcg_mo_n,abc_e_fcg_mo_n,abc_e30_fcg_mo_n,abc_e75_fcg_mo_n))
abc_df_site$time<-factor(abc_df_site$time, levels = c('P','A','E','E30','E75'))

F1D<-ggplot(abc_df_site,aes(y=n,x=time,color=comparison,group=comparison))+
  geom_point(size=3.2, alpha = 0.75)+geom_line(alpha=0.75)+
  scale_color_manual(values=site_colors)+
  labs(x="Timepoint",y="N differentially abundant genera", color = "Comparison")+
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.line = element_line(color = 'black'))+
  scale_y_continuous(limits=c(0,12),expand = c(0,0))