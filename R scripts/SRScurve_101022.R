library(qiime2R)
require(ggplot2)
theme_set(theme_bw())
library(vegan)
library(phyloseq)
library(QsRutils)
library(SRS)


#loading physeq object w/ full dataset #495 samples
physeq<-qza_to_phyloseq(
  features="Files/ASV_table_final.qza",
  tree="Files/rooted-tree.qza",
  "Files/taxonomy_vs.qza",
  metadata = "Files/mapping_tmo_0.tsv")

physeq_tmo<-subset_samples(physeq, time!="MI"& time!="MF")

#primeiro, good's coverage (n de singletons na amostra/n de reads da amostra)
gc<-goods(t(otu_table(physeq_tmo)))
min(gc$goods)
#todas as amostras tem gc>=99.97584%

#e após remover amostras defectivas?
physeq_good<-qza_to_phyloseq(
  features="Files/ASV_table_final_goodsamples.qza",
  tree="Files/rooted-tree.qza",
  "Files/taxonomy_vs.qza",
  metadata = "Files/mapping_tmo_0.tsv")

physeq_tmo_good<-subset_samples(physeq_good, time!="MI"& time!="MF")

gc_good<-goods(t(otu_table(physeq_tmo_good)))
min(gc_good$goods)
#mesma coisa, todas as amostras tem gc>=99.97584%

#FS4
SRScurve(as.data.frame(otu_table(physeq_tmo)),max.sample.size = 15000, step =  50, metric = "simpson", label = F, 
         col=alpha(rgb(0,0,0), 0.25), ylab = "Diversity", xlab = "Sequencing depth",  xaxt = "n")
par(las = 1)
abline(v=3000, lty = "dotted")
axis(side = 1, at = seq(0, 15000, 500), labels = FALSE)
text(x = seq(0, 15000, 500), y = -0.1, labels = seq(0, 15000, 500), xpd = NA, srt = 45, adj=1, cex = 1)
