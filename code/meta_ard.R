#         #        #         #        #         #        #         #          # 
#                       Meta-analysis for Absolute Risk Differences
# Inputs  : ardse.csv 
# Outputs : Meta-analysis forest plots for ARD for each outcome by drugclass
#         #        #         #        #         #        #         #          # 
library(data.table)
library(tidyr)
library(dplyr)
library(scales)
library(stringr)
library(metafor)
library(meta)
# Set working directory
wd <- "~/Documents/Sys-Rev-Meta"


setwd(wd)

# Load data 
df.ard <- read.csv("data/ardse.csv")
# Labels for titles

outcomes <- unique(df.ard$outcome)
df.ard$outcome <- ifelse(df.ard$outcome == "allcauseMort", "All-cause Mortality",
                          df.ard$outcome)
df.ard$outcome <- ifelse(df.ard$outcome == "CVMort", "Cardiovascular Mortality",
                          df.ard$outcome)
df.ard$outcome <- ifelse(df.ard$outcome == "HospHF", "Hospitalization for Heart Failure",
                          df.ard$outcome)
df.ard$outcome <- ifelse(df.ard$outcome == "stroke", "Stroke",
                          df.ard$outcome)
df.ard$outcome <- ifelse(df.ard$outcome == "MI", "Myocardial Infarction",
                          df.ard$outcome)
df.ard$outcome <- ifelse(df.ard$outcome == "sustGFRdecl", "Composite Renal Outcome",
                          df.ard$outcome)
df.ard$Class <- ifelse(df.ard$type == "GLP1", "GLP-1RAs", "SGLT2i")

# Vector of outcomes

outcomes <- unique(df.ard$outcome)

df.ard <- subset(df.ard, subset = !(df.ard$trialname == "LEADER/SUSTAIN-6"))
# Meta analysis by outcome/drugclass (ard)
for(i in 1:length(outcomes)){
  
  
  # I: Subset data per outcome/drugclass
  d.o <- subset(df.ard, subset = (df.ard$outcome == outcomes[i]))
  d.g <- subset(d.o, subset = (Class == "GLP-1RAs"))
  d.s <- subset(d.o, subset = (Class == "SGLT2i"))
  # II. Metagen for each outcome/class
  # II.1  model.g for GLP-1RA  
  model.g <- metagen(TE           = 100*d.g$ard,
                     seTE         = 100*d.g$ardse2,
                     data         = d.g,
                     sm           = "RD",
                     studlab      = d.g$trialname,
                     fixed        = FALSE,
                     random       = TRUE,
                     method.tau   = "REML",
                     overall      = FALSE,
                     hakn         = TRUE,
                     #subgroup     = d.g$Class,
                     title        = "study weights"
  )
  # II.2 model.s for SGLT2i  
  model.s <- metagen(TE           = 100*d.s$ard,
                     seTE         = 100*d.s$ardse2,
                     data         = d.s,
                     sm           = "RD",
                     studlab      = d.s$trialname,
                     fixed        = FALSE,
                     random       = TRUE,
                     method.tau   = "REML",
                     overall      = TRUE,
                     hakn         = TRUE,
                     #subgroup     = d.s$Class,
                     title        = "study weights"
  )
  # Forest plots
  # III.a GLP-1Ra
  #png(paste0("_forest/meta_ard_glp1_", outcomes[i],".png"), 
  #     width = 6, height = 4, units = 'in', res = 300)    
  forest.meta(model.g, 
         #atransf  = log, 
         overall = T,
         studlab  = TRUE,
         fontsize = 9,
         xlim     = c(-0.3,0.2)*100,
         at       = c(-0.1,-0.1,0,0.1)*100,
         fs.axis  = 7,
         ref      = 0,
         level    = 0.95,
         weight.subgroup = "weight",
         smlab = "",
         col.diamond = "plum3",
         col.square  = "plum4",
         leftcols = c("studlab"),
         leftlabs = paste0(unique(d.o$outcome), "\nGLP-1RAs"),
         label = T,
         print.Q = T, print.pval.Q = T,
         print.I2 = T,
         print.tau2 = F,
         fs.hetstat = 7,
   #      text.addline1 =  unique(d.o$outcomen),
         subgroup=F,
         subgroup.name = "AAA",
         lab.NA = "",
         lab.NA.weight = "",
         lab.NA.effect = ""
  )
  #dev.off()
  # III.b SGLT2i  
  #png(paste0("_forest/meta_ard_sglt2_", outcomes[i],".png"), 
  #    width = 6, height = 4, units = 'in', res = 300) 
  forest(model.s, 
         #atransf  = log, 
         main ="a", top = 3,
         overall = T,
       #  studlab  = TRUE,
         fontsize = 9,
         xlim     = c(-0.3,0.2)*100,
         at       = c(-0.1,-0.1,0,0.1)*100,
         fs.axis  = 7,
         ref      = 0,
         level    = 0.95,
         weight.subgroup = "weight",
         col.diamond = "plum4",
         col.square  = "plum3",
         leftcols = c("studlab"),
         leftlabs ="\nSGLT2i",
         smlab = "",  
         label = T,
         print.Q = T, print.pval.Q = T,
         print.I2 = T,
         print.tau2 = F,
         fs.hetstat = 7,
  #       text.addline1 =  unique(d.o$outcomen),
         subgroup=F
  )
  #dev.off()
}

