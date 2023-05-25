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
library(png)
library(jpeg)
library(ggimage)
library(patchwork)
library(magick)
library(ggpubr)
# Set working directory
 
  wd <- "~/Documents/Sys-Rev-Meta"
  setwd(wd)

# A: Meta analysis and forest plots for log hazard ratios

# Load data ------
  df <- read.csv("_data/df.csv")



# Soloist and Scored with divergent measure for HHF  
  df <- subset(df, !( (df$trialname == "SOLOIST-WHF"& df$outcome == "HospHeartFailure") 
                     )
               )
  df$hr <- ifelse((df$trialname == "SOLOIST-WHF"& df$outcome == "HospHeartFailure") ,
                   NA, 
                  df$hr)
  df$lci <- ifelse((df$trialname == "SOLOIST-WHF"& df$outcome == "HospHeartFailure") ,
                   NA, 
                  df$lci)
  df$uci <- ifelse((df$trialname == "SOLOIST-WHF"& df$outcome == "HospHeartFailure") ,
                    NA, 
                   df$uci)
  df$trialname <- ifelse((df$trialname == "SOLOIST-WHF"& df$outcome == "HospHeartFailure"),
                   NA, 
                   df$trialname)
  df$trialname <- ifelse(is.na(df$hr), NA, 
                         df$trialname)
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
df$outcomen <- ifelse(df$outcome == "allcauseMort", "All-cause mortality",
                      df$outcomen)
df$outcomen <- ifelse(df$outcome == "CVMort", "Cardiovascular mortality",
                      df$outcomen)
df$outcomen <- ifelse(df$outcome == "HospHeartFailure", "Hospitalization for heart failure",
                      df$outcomen)
df$outcomen <- ifelse(df$outcome == "stroke", "Stroke",
                      df$outcomen)
df$outcomen <- ifelse(df$outcome == "MI", "Myocardial infarction",
                      df$outcomen)
df$outcomen <- ifelse(df$outcome == "sustGFRdecl", "Composite renal outcome",
                      df$outcomen)
outcomes <- unique(df$outcome) #vector of outcomes

df$Class <- ifelse(df$type == "GLP1", "GLP-1RAs", "SGLT2i")

# Meta analysis by outcome/drugclass (loghr)
l.mg <- list()
l.ms <- list()
lnames.g <- list()
lnames.s <- list()
for(i in 1:length(outcomes)){
  # I: Subset data per outcome/drugclass
  d.o <- subset(df, subset = (df$outcome == outcomes[i]))
  d.g <- subset(d.o, subset = (Class == "GLP-1RAs"))
  d.s <- subset(d.o, subset = (Class == "SGLT2i"))
  rowsg <- nrow(d.g)
  rowss <- nrow(d.s)
 if(rowsg<13){
  d.g[seq(rowsg+1, 13),] <- NA
 } else{
   print("")
 }
  if(rowss<13){
    d.s[seq(rowss+1, 13),] <- NA
  } else{
    print("")
  }


  # Empty rows for figure proportions purposes

 
  d.g$trialname <- ifelse(is.na(d.g$trialname)==T, " ", d.g$trialname )
  d.s$trialname <- ifelse(is.na(d.s$trialname)==T, " ", d.s$trialname )
  
  # II. Metagen for each outcome/class
  # II.1  model.g for GLP-1RA
  l.mg[[i]] <- metagen(TE         = d.g$loghr,
                     seTE         = d.g$logsei,
                     data         = d.g,
                     method.ci    = "z",
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
                     method.ci    = "z",
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
   tiff(paste0("~/Documents/_meta/_forest/meta_loghr_glp1_", outcomes[i],".tif"), 
      width = 6, height = 4, units = 'in', res = 1000)  
  par(mfrow = c(1,1), oma= c(1,1,1,1)*0)
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
         sortvar = d.g$ordertrial,
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
  tiff(paste0("~/Documents/_meta/_forest/meta_loghr_sglt2i_", outcomes[i],".tif"), 
       width = 6, height = 4, units = 'in', res = 1000)  

   

  forest(l.ms[[i]], 
         atransf  = log, 
         overall = T,
         sortvar = d.s$ordertrial,
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
  lnames.g[[i]] <- paste0("~/Documents/_meta/_forest/meta_loghr_glp1_", outcomes[i],".tif")
  lnames.s[[i]] <- paste0("~/Documents/_meta/_forest/meta_loghr_sglt2i_", outcomes[i],".tif")
  
}



# Export image in a grid
p1g <- image_read(lnames.g[[5]])
p1s <- image_read(lnames.s[[5]])
p2g <- image_read(lnames.g[[3]])
p2s <- image_read(lnames.s[[3]])
p3g <- image_read(lnames.g[[2]])
p3s <- image_read(lnames.s[[2]])
p4g <- image_read(lnames.g[[6]])
p4s <- image_read(lnames.s[[6]])
      
      p1g <- ggplot() + background_image(p1g) 
      p1s <- ggplot() + background_image(p1s) 
      p2g <- ggplot() + background_image(p2g) 
      p2s <- ggplot() + background_image(p2s) 
      p3g <- ggplot() + background_image(p3g) 
      p3s <- ggplot() + background_image(p3s) 
      p4g <- ggplot() + background_image(p4g) 
      p4s <- ggplot() + background_image(p4s) 

p <- (p1g + p1s) / (p2g + p2s) /(p3g + p3s)/(p4g + p4s)
ggsave("/Users/jose/Documents/_meta/_rev3/Figure1.tiff",plot = p, dpi = 1000, width =  12, height =  16)

# Same exercise for seconday outcomes (no need to have empty spaces) ====
lnames.g2 <- list()
lnames.s2 <- list()
l.mg2 <- list()
l.ms2 <- list()
outcomes2 <- c("allcauseMort", "MI", "stroke")
for(i in 1:length(outcomes2)){
  # I: Subset data per outcome/drugclass
  d.o <- subset(df, subset = (df$outcome == outcomes2[i]))
  d.g <- subset(d.o, subset = (Class == "GLP-1RAs"))
  d.s <- subset(d.o, subset = (Class == "SGLT2i"))

  
  # Empty rows for figure proportions purposes
  
  
  # II. Metagen for each outcome/class
  # II.1  model.g for GLP-1RA
  l.mg2[[i]] <- metagen(TE         = d.g$loghr,
                       seTE         = d.g$logsei,
                       data         = d.g,
                       method.ci    = "z",
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
  l.ms2[[i]] <- metagen(TE           = d.s$loghr,
                       seTE         = d.s$logsei,
                       data         = d.s,
                       method.ci    = "z",
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
  par(mfrow = c(1,1), oma= c(1,1,1,1)*0, mar= c(1,1,1,1)*4)
  png(paste0("~/Documents/_meta/_forest/meta_loghr_glp1_", outcomes2[i],".png"), 
      width = 6, height = 4, units = 'in', res = 300)  

  forest(l.mg2[[i]], 
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
         sortvar = d.g$ordertrial,
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
  png(paste0("~/Documents/_meta/_forest/meta_loghr_sglt2i_", outcomes2[i],".png"), 
      width = 6, height = 4, units = 'in', res = 300)  
  
  
  
  forest(l.ms2[[i]], 
         atransf  = log, 
         overall = T,
         sortvar = d.s$ordertrial,
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
  lnames.g2[[i]] <- paste0("~/Documents/_meta/_forest/meta_loghr_glp1_", outcomes2[i],".png")
  lnames.s2[[i]] <- paste0("~/Documents/_meta/_forest/meta_loghr_sglt2i_", outcomes2[i],".png")
  
}

p1g2 <- image_read(lnames.g2[[1]])
p1s2 <- image_read(lnames.s2[[1]])
p2g2 <- image_read(lnames.g2[[3]])
p2s2 <- image_read(lnames.s2[[3]])
p3g2 <- image_read(lnames.g2[[2]])
p3s2 <- image_read(lnames.s2[[2]])


p1g2 <- ggplot() + background_image(p1g2) + theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))
p1s2 <- ggplot() + background_image(p1s2) + theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))
p2g2 <- ggplot() + background_image(p2g2) + theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))
p2s2 <- ggplot() + background_image(p2s2) + theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))
p3g2 <- ggplot() + background_image(p3g2) + theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))
p3s2 <- ggplot() + background_image(p3s2) + theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))


p2.1 <- (p1g2 / p1s2) +theme(plot.margin = unit(c(0, 0,2, 0), "cm"))
p2.2 <- (p2g2 / p2s2) +theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))
p2.3 <- (p3g2 / p3s2) +theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))

ggsave("/Users/jose/Documents/_meta/_rev3/metahrall.png",plot = p2.1, dpi = 300, width =  12, height =  16)
ggsave("/Users/jose/Documents/_meta/_rev3/metahrstr.png",plot = p2.2, dpi = 300, width =  12, height =  16)
ggsave("/Users/jose/Documents/_meta/_rev3/metahrmi.png",plot = p2.3, dpi = 300, width =  12, height =  16)




