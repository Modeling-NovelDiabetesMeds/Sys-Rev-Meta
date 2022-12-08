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
library(ggplot2)
library(patchwork)
library(magick)
library(ggpubr)
# Set working directory
wd <- "~/Documents/_meta"


setwd(wd)

# Load data 
df.ard <- read.csv("_data/ardse.csv")
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


# trial with no info into blanks
df.ard$ard <- ifelse((df.ard$trialname == "SOLOIST-WHF"& df.ard$outcome == "HospHeartFailure") |
                      ( df.ard$trialname == "SCORED" & df.ard$outcome == "HospHeartFailure"),
                      NA,
                     df.ard$ard
                        )
                 
df.ard$trialname <- ifelse(is.na(df.ard$ard), " ", 
                       df.ard$trialname)

df.ard$ordertrial <- ifelse(is.na(df.ard$ard), 50, 
                           df.ard$ordertrial)


# Pooled result for leader/sustain
df.ard <- subset(df.ard, subset = !(df$trialname == "LEADER/SUSTAIN-6"))


# Meta analysis by outcome/drugclass (ard)

l.mg <- list()
l.ms <- list()
lnames.g <- list()
lnames.s <- list()


for(i in 1:length(outcomes)){
  # I: Subset data per outcome/drugclass
  d.o <- subset(df.ard, subset = (df.ard$outcome == outcomes[i]))
  d.g <- subset(d.o, subset = (Class == "GLP-1RAs"))
  d.s <- subset(d.o, subset = (Class == "SGLT2i"))
  # Empty rows for figure proportions purposes
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
  model.s<- metagen(TE           = 100*d.s$ard,
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
                     title        = "study weights",  
         verbose=TRUE, digits=5)
  # Forest plots
  # III.a GLP-1Ra
  png(paste0("~/Documents/_meta/_forest/meta_ard_glp1_", outcomes[i],".png"), 
       width = 6, height = 4, units = 'in', res = 300)    
  forest.meta(model.g, 
          sortvar =  d.g$ordertrial,
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
  dev.off()
  # III.b SGLT2i  
  png(paste0("~/Documents/_meta/_forest/meta_ard_sglt2_", outcomes[i],".png"), 
      width = 6, height = 4, units = 'in', res = 300) 
  forest(model.s, 
         #atransf  = log,
         sortvar = d.s$ordertrial,
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
  dev.off()
  
  lnames.g[[i]] <- paste0("~/Documents/_meta/_forest/meta_ard_glp1_", outcomes[i],".png")
  lnames.s[[i]] <- paste0("~/Documents/_meta/_forest/meta_ard_sglt2_", outcomes[i],".png")
  
  }

#### EXPORTING Figure 2 and Online supplement




# Export image in a grid
p1g <- image_read(lnames.g[[2]])
p1s <- image_read(lnames.s[[2]])
p2g <- image_read(lnames.g[[5]])
p2s <- image_read(lnames.s[[5]])
p3g <- image_read(lnames.g[[4]])
p3s <- image_read(lnames.s[[4]])
p4g <- image_read(lnames.g[[3]])
p4s <- image_read(lnames.s[[3]])

p1g <- ggplot() + background_image(p1g) 
p1s <- ggplot() + background_image(p1s) 
p2g <- ggplot() + background_image(p2g) 
p2s <- ggplot() + background_image(p2s) 
p3g <- ggplot() + background_image(p3g) 
p3s <- ggplot() + background_image(p3s) 
p4g <- ggplot() + background_image(p4g) 
p4s <- ggplot() + background_image(p4s) 

p <- (p1g + p1s) / (p2g + p2s) /(p3g + p3s)/(p4g + p4s)
ggsave("/Users/josemariarodriguezvaladez/Documents/_meta/_figures/_revision2/Figure2.png",plot = p, dpi = 300, width =  12, height =  16)






#################################################################################
#        Same exercise for seconday outcomes (no need to have empty spaces) ====
#################################################################################

lnames.g2 <- list()
lnames.s2 <- list()
l.mg2 <- list()
l.ms2 <- list()
outcomes2 <- c("All-cause Mortality", "Myocardial Infarction", "Stroke")

for(i in 1:length(outcomes2)){
  # I: Subset data per outcome/drugclass
  d.o <- subset(df.ard, subset = (df.ard$outcome == outcomes2[i]))
  d.o <- subset(d.o, subset = !(is.na(d.o$ard)))
  d.g <- subset(d.o, subset = (Class == "GLP-1RAs"))
  d.s <- subset(d.o, subset = (Class == "SGLT2i"))
  
  
  # Empty rows for figure proportions purposes
  
  
  # II. Metagen for each outcome/class
  # II.1  model.g for GLP-1RA
  l.mg2[[i]] <- metagen(TE           = 100*d.g$ard,
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
  l.ms2[[i]]<- metagen(TE           = 100*d.s$ard,
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
                        title        = "study weights",  
                        verbose=TRUE, digits=5)
  # III. Forest plot
  # III.a GLP-1Ra
  png(paste0("~/Documents/_meta/_forest/meta_ard_glp1_", outcomes2[i],".png"), 
      width = 6, height = 4, units = 'in', res = 300)    
  forest.meta(l.mg2[[i]], 
              sortvar =  d.g$ordertrial,
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
  dev.off()
  # III.b SGLT2i  
  png(paste0("~/Documents/_meta/_forest/meta_ard_sglt2_", outcomes2[i],".png"), 
      width = 6, height = 4, units = 'in', res = 300) 
  forest(l.ms2[[i]], 
         #atransf  = log,
         sortvar = d.s$ordertrial,
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
  dev.off()
  
  lnames.g2[[i]] <- paste0("~/Documents/_meta/_forest/meta_ard_glp1_", outcomes2[i],".png")
  lnames.s2[[i]] <- paste0("~/Documents/_meta/_forest/meta_ard_sglt2_", outcomes2[i],".png")
  
}

# Transoform images to be able to combine them
    p1g2 <- image_read(lnames.g2[[1]])
    p1s2 <- image_read(lnames.s2[[1]])
    p2g2 <- image_read(lnames.g2[[2]])
    p2s2 <- image_read(lnames.s2[[2]])
    p3g2 <- image_read(lnames.g2[[3]])
    p3s2 <- image_read(lnames.s2[[3]])
    
    
    p1g2 <- ggplot() + background_image(p1g2) + theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))
    p1s2 <- ggplot() + background_image(p1s2) + theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))
    p2g2 <- ggplot() + background_image(p2g2) + theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))
    p2s2 <- ggplot() + background_image(p2s2) + theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))
    p3g2 <- ggplot() + background_image(p3g2) + theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))
    p3s2 <- ggplot() + background_image(p3s2) + theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))

 #  create individual panels for three secondary outcomes
p2.1 <- (p1g2 / p1s2) +theme(plot.margin = unit(c(0, 0,2, 0), "cm"))
p2.2 <- (p2g2 / p2s2) +theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))
p2.3 <- (p3g2 / p3s2) +theme(plot.margin = unit(c(0, 0,0.25, 0), "cm"))

ggsave("/Users/josemariarodriguezvaladez/Documents/_meta/_figures/_revision2/metaardall.png",plot = p2.1, dpi = 300, width =  12, height =  16)
ggsave("/Users/josemariarodriguezvaladez/Documents/_meta/_figures/_revision2/metaardmi.png",plot = p2.2, dpi = 300, width =  12, height =  16)
ggsave("/Users/josemariarodriguezvaladez/Documents/_meta/_figures/_revision2/metaardstr.png",plot = p2.3, dpi = 300, width =  12, height =  16)




