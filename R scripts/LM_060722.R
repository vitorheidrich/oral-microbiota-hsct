library(shiny)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(lattice)
library(DHARMa)

clindata_bin<-read.csv("Files/DadosClinicosTMO_09.05.2021.csv")
clindata_bin$Patient<-NULL
clindata_bin$patient<-clindata_bin$Patient.acronym
clindata_bin$Patient.acronym<-NULL
atbdata<-read.csv("Files/300622_ATB_classification.csv")
clindata_bin<-merge(clindata_bin,atbdata)
microbiota<-read.csv("Files/300522_microbiota-variables.csv")
clindata_microbiota<-merge(clindata_bin,microbiota)

# model_s_DB<-lm(Simpson.s_DB~Age+Sex+Underlying.disease+Pretransplant.comorbidity..HCT.CI.+Disease.risk.index+Conditioning.intensity+Total.body.irradiation+T.cell.depletion+Graft.source+Donor+carbapenems_cat+glycopeptides_cat+penicillins_cat+cephalosporins_cat+dot, data = clindata_microbiota)
# model_s_GCF<-lm(Simpson.s_GCF~Age+Sex+Underlying.disease+Pretransplant.comorbidity..HCT.CI.+Disease.risk.index+Conditioning.intensity+Total.body.irradiation+T.cell.depletion+Graft.source+Donor+carbapenems_cat+glycopeptides_cat+penicillins_cat+cephalosporins_cat+dot, data = clindata_microbiota)
# model_s_OM<-lm(Simpson.s_OM~Age+Sex+Underlying.disease+Pretransplant.comorbidity..HCT.CI.+Disease.risk.index+Conditioning.intensity+Total.body.irradiation+T.cell.depletion+Graft.source+Donor+carbapenems_cat+glycopeptides_cat+penicillins_cat+cephalosporins_cat+dot, data = clindata_microbiota)

model_s_DB<-lm(Simpson.s_DB~cephalosporins_cat+carbapenems_cat+glycopeptides_cat+penicillins_cat+dot, data = clindata_microbiota)
model_s_GCF<-lm(Simpson.s_GCF~cephalosporins_cat+carbapenems_cat+glycopeptides_cat+penicillins_cat+dot, data = clindata_microbiota)
model_s_OM<-lm(Simpson.s_OM~cephalosporins_cat+carbapenems_cat+glycopeptides_cat+penicillins_cat+dot, data = clindata_microbiota)

summary(model_s_DB)
# Call:
#   lm(formula = Simpson.s_DB ~ cephalosporins_cat + carbapenems_cat + 
#        glycopeptides_cat + penicillins_cat + dot, data = clindata_microbiota)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -1.14259 -0.10336  0.03761  0.18918  0.71718 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            0.713128   0.182892   3.899 0.000722 ***
#   cephalosporins_catYes  0.294667   0.176404   1.670 0.108393    
# carbapenems_catYes     0.102401   0.173917   0.589 0.561741    
# glycopeptides_catYes  -0.136454   0.184620  -0.739 0.467320    
# penicillins_catYes     0.188979   0.221867   0.852 0.403125    
# dot                   -0.008520   0.004054  -2.102 0.046735 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3736 on 23 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.3132,	Adjusted R-squared:  0.1639 
# F-statistic: 2.098 on 5 and 23 DF,  p-value: 0.1023

summary(model_s_GCF)
# Call:
#   lm(formula = Simpson.s_GCF ~ cephalosporins_cat + carbapenems_cat + 
#        glycopeptides_cat + penicillins_cat + dot, data = clindata_microbiota)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.48051 -0.05378  0.03459  0.09972  0.40539 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)            0.966019   0.100701   9.593  1.1e-09 ***
#   cephalosporins_catYes  0.042670   0.096906   0.440   0.6636    
# carbapenems_catYes     0.008822   0.095117   0.093   0.9269    
# glycopeptides_catYes  -0.089505   0.101106  -0.885   0.3848    
# penicillins_catYes     0.087212   0.122408   0.712   0.4830    
# dot                   -0.005674   0.002218  -2.558   0.0172 *  
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.2061 on 24 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.3975,	Adjusted R-squared:  0.2719 
# F-statistic: 3.166 on 5 and 24 DF,  p-value: 0.02463

summary(model_s_OM)
# Call:
#   lm(formula = Simpson.s_OM ~ cephalosporins_cat + carbapenems_cat + 
#        glycopeptides_cat + penicillins_cat + dot, data = clindata_microbiota)
# 
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.95956 -0.16965  0.09829  0.17650  0.44797 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)   
# (Intercept)            0.679418   0.185163   3.669  0.00143 **
#   cephalosporins_catYes  0.201127   0.185808   1.082  0.29133   
# carbapenems_catYes     0.168657   0.170562   0.989  0.33400   
# glycopeptides_catYes   0.253394   0.187687   1.350  0.19136   
# penicillins_catYes     0.165546   0.215811   0.767  0.45157   
# dot                   -0.016729   0.004576  -3.656  0.00147 **
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3443 on 21 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.448,	Adjusted R-squared:  0.3165 
# F-statistic: 3.408 on 5 and 21 DF,  p-value: 0.02071

model_comp.s_DB<-lm(comp_S.DB~cephalosporins_cat+carbapenems_cat+glycopeptides_cat+penicillins_cat+dot, data = clindata_microbiota)
model_comp.s_GCF<-lm(comp_S.GCF~cephalosporins_cat+carbapenems_cat+glycopeptides_cat+penicillins_cat+dot, data = clindata_microbiota)
model_comp.s_OM<-lm(comp_S.OM~cephalosporins_cat+carbapenems_cat+glycopeptides_cat+penicillins_cat+dot, data = clindata_microbiota)

summary(model_comp.s_DB)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.61686 -0.24623  0.00588  0.20900  0.71284 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)            0.392789   0.188773   2.081   0.0488 *
#   cephalosporins_catYes -0.367667   0.182076  -2.019   0.0553 .
# carbapenems_catYes     0.171253   0.179510   0.954   0.3500  
# glycopeptides_catYes  -0.462381   0.190557  -2.426   0.0235 *
#   penicillins_catYes    -0.179279   0.229002  -0.783   0.4417  
# dot                    0.004986   0.004184   1.192   0.2456  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3856 on 23 degrees of freedom
# (2 observations deleted due to missingness)
# Multiple R-squared:  0.3452,	Adjusted R-squared:  0.2029 
# F-statistic: 2.425 on 5 and 23 DF,  p-value: 0.06621
summary(model_comp.s_GCF)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.59830 -0.20612 -0.02995  0.20110  0.56194 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)  
# (Intercept)            0.446800   0.169742   2.632   0.0146 *
# cephalosporins_catYes -0.293657   0.163345  -1.798   0.0848 .
# carbapenems_catYes     0.020564   0.160329   0.128   0.8990  
# glycopeptides_catYes  -0.195220   0.170425  -1.145   0.2633  
# penicillins_catYes    -0.375559   0.206331  -1.820   0.0812 .
# dot                    0.003110   0.003738   0.832   0.4136  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3475 on 24 degrees of freedom
# (1 observation deleted due to missingness)
# Multiple R-squared:  0.2181,	Adjusted R-squared:  0.0552 
# F-statistic: 1.339 on 5 and 24 DF,  p-value: 0.2819
summary(model_comp.s_OM)
# Residuals:
#   Min       1Q   Median       3Q      Max 
# -0.71132 -0.16035 -0.04329  0.23098  0.59550 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)
# (Intercept)            0.062275   0.197815   0.315    0.756
# cephalosporins_catYes -0.065766   0.198505  -0.331    0.744
# carbapenems_catYes    -0.194800   0.182217  -1.069    0.297
# glycopeptides_catYes   0.017063   0.200512   0.085    0.933
# penicillins_catYes     0.156495   0.230558   0.679    0.505
# dot                    0.001931   0.004888   0.395    0.697
# 
# Residual standard error: 0.3678 on 21 degrees of freedom
# (4 observations deleted due to missingness)
# Multiple R-squared:  0.1064,	Adjusted R-squared:  -0.1063 
# F-statistic: 0.5002 on 5 and 21 DF,  p-value: 0.7726

