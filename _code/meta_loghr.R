#         #        #         #        #         #        #         #          # 
#                       Meta-analysis and forest plots
# Inputs  : df.csv 
# Outputs : forest plots for logHR, outcome, by drug class
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

# A: Meta analysis and forest plots for log hazard ratios

# Load data ------
  df <- read.csv("data/df.csv")



# Soloist and Scored with divergent measure for HHF  
  df <- subset(df, !( (trialname == "SOLOIST-WHF"& df$outcome == "HospHeartFailure") 
                     |( df$trialname == "SCORED" & df$outcome == "HospHeartFailure")
                     )
               )
  # Pooled result for leader/sustain
  df <- subset(df, subset = !(df$trialname == "LEADER/SUSTAIN-6"))
  
# Compute some variables of interest ------

# Log variance for HR, based on CI
df$loghr  <- log(df$hr)

# Parmar et al to extract logvi from CI
df$logvi  <- ((log(df$uci) - log(df$lci)) / (2 * qnorm(.975) ))^2 
df$logsei <- sqrt(df$logvi)

# Outcome labels for titles 
df$outcomen <- df$outcome
df$outcomen <- ifelse(df$outcome == "allcauseMort", "All-cause Mortality",
                      df$outcomen)
df$outcomen <- ifelse(df$outcome == "CVMort", "Cardiovascular Mortality",
                      df$outcomen)
df$outcomen <- ifelse(df$outcome == "HospHeartFailure", "Hospitalization for Heart Failure",
                      df$outcomen)
df$outcomen <- ifelse(df$outcome == "stroke", "Stroke",
                      df$outcomen)
df$outcomen <- ifelse(df$outcome == "MI", "Myocardial Infarction",
                      df$outcomen)
df$outcomen <- ifelse(df$outcome == "sustGFRdecl", "Composite Renal Outcome",
                      df$outcomen)
outcomes <- unique(df$outcome) #vector of outcomes

df$Class <- ifelse(df$type == "GLP1", "GLP-1RAs", "SGLT2i")

# Meta analysis by outcome/drugclass (loghr)
l.mg <- list()
l.ms <- list()
for(i in 1:length(outcomes)){
  # I: Subset data per outcome/drugclass
  d.o <- subset(df, subset = (df$outcome == outcomes[i]))
  d.g <- subset(d.o, subset = (Class == "GLP-1RAs"))
  d.s <- subset(d.o, subset = (Class == "SGLT2i"))
  # II. Metagen for each outcome/class
  # II.1  model.g for GLP-1RA
  l.mg[[i]] <- metagen(TE         = d.g$loghr,
                     seTE         = d.g$logsei,
                     data         = d.g,
                     method.ci = "z",
                     sm           = "HR",
                     studlab      = d.g$trialname,
                     fixed        = FALSE,
                     random       = TRUE,
                     method.tau   = "REML",
                     overall      = FALSE,
                     hakn         = TRUE,
                     title        = "study weights"
  )
  # II.2 model.s for SGLT2i
  l.ms[[i]] <- metagen(TE           = d.s$loghr,
                     seTE         = d.s$logsei,
                     data         = d.s,
                     method.ci = "t",
                     sm           = "HR",
                     studlab      = d.s$trialname,
                     fixed        = FALSE,
                     random       = TRUE,
                     method.tau   = "REML",
                     overall      = TRUE,
                     hakn         = TRUE,
                     title        = "study weights"
  )
  # III. Forest plot
  # III.a GLP-1Ra
   png(paste0("~/Documents/_meta/_forest/meta_loghr_glp1_", outcomes[i],".png"), 
      width = 6, height = 4, units = 'in', res = 300)  
  par(mfrow = c(2,2))
  forest(l.mg[[i]], 
         atransf  = log, 
         overall = T,
         studlab  = TRUE,
         effect = TRUE,
         predict = FALSE,
         fontsize = 9,
         xlim     = c(0.2,2.6),
         at       = c(0.5,0.75,1,1.5),
         fs.axis  = 7,
         ref      = 1,
         level    = 0.95,
         weight.subgroup = "weight",
         smlab = "",
         col.diamond = "slategray2",
         col.square  = "slategray4",
         leftcols = c("studlab"),
         leftlabs = paste0(unique(d.o$outcomen), "\nGLP-1RAs"),
         label = T,
         print.Q = T, print.pval.Q = T,
         print.I2 = T,
         print.tau2 = F,
         fs.hetstat = 7,
         subgroup=F
         )
   dev.off()
  # III.b GLP-1Ra    
   png(paste0("~/Documents/_meta/_forest/meta_loghr_sglt2i_", outcomes[i],".png"), 
       width = 6, height = 4, units = 'in', res = 300)  

   

  forest(l.ms[[i]], 
         atransf  = log, 
         overall = T,
         studlab  = TRUE,
         fontsize = 9,
         xlim     = c(0.2,2.6),
         at       = c(0.5,0.75,1,1.5),
         fs.axis  = 7,
         ref      = 1,
         level    = 0.95,
         weight.subgroup = "weight",
         smlab = "",
         col.diamond = "slategray4",
         col.square  = "slategray2",
         leftcols = c("studlab"),
         leftlabs ="\nSGLT2i",
         label = T,
         print.Q = T, print.pval.Q = T,
         print.I2 = T,
         print.tau2 = F,
         fs.hetstat = 7,
  #      text.addline1 =  "SGLT2i",
         subgroup=F
  )
  dev.off()
  
}





       






