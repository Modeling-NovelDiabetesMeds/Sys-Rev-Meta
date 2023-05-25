# Forest mediators
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


df <- read.csv("~/Documents/Sys-Rev-Meta/Potential_mediators/_data/mediators_md.csv")
order <- read.csv("~/Documents/Sys-Rev-Meta/Potential_mediators/_data/order.csv")

df <- merge(df, order, by = "trialname", all.x = TRUE)

df$Class <- ifelse(df$class == "GLP1", "GLP-1RA", "SGLT2i")

df$outcometitle  <- NA
df$outcometitle <- ifelse(df$outcome == "hba1c", "HbA1c ",df$outcometitle)
df$outcometitle <- ifelse(df$outcome == "bodyweight", "Bodyweight",df$outcometitle)
df$outcometitle <- ifelse(df$outcome == "sbp", "Systolic blood pressure",df$outcometitle)
df$outcometitle <- ifelse(df$outcome == "hematocrit", "Hematocrit",df$outcometitle)
df$outcometitle <- ifelse(df$outcome == "uacr", "uACR",df$outcometitle)

tc <- unique(df[c("trialname","Class")])



df <-  df %>% complete(trialname, outcome) %>%
    select(-Class) %>%
    left_join( tc, by = "trialname")

df$trialname <- ifelse(is.na(df$diff.se), " ",  df$trialname)


df$diff <- ifelse(is.na(df$diff.se), NA,  df$diff)
df$diff.se <- ifelse(is.na(df$diff.se), NA,  df$diff.se)
df$ordertrial <- ifelse(is.na(df$diff.se), NA,  df$ordertrial)
df <- df[order(df$ordertrial),]
#hba1c
d.gh <- subset(df, subset = (Class == "GLP-1RA" & outcome == "hba1c"))
        rows <- nrow(d.gh)
        if(rows<13){
          d.gh[seq(rows+1, 13),] <- NA
        } else{
          print("")
        }
d.sh <- subset(df, subset = (Class == "SGLT2i" & outcome == "hba1c"))
        rows <- nrow(d.sh)
        if(rows<13){
          d.sh[seq(rows+1, 13),] <- NA
        } else{
          print("")
        }

#bw
d.gb <- subset(df, subset = (Class == "GLP-1RA" & outcome == "bodyweight" ))
          rows <- nrow(d.gb)
          if(rows<13){
            d.gb[seq(rows+1, 13),] <- NA
          } else{
            print("")
          }
d.sb <- subset(df, subset = (Class == "SGLT2i" & outcome == "bodyweight" ))
          rows <- nrow(d.sb)
          if(rows<13){
            d.sb[seq(rows+1, 13),] <- NA
          } else{
            print("")
          }
#sbp
d.gs <- subset(df, subset = (Class == "GLP-1RA" & outcome == "sbp" ))
        rows <- nrow(d.gs)
        if(rows<13){
          d.gs[seq(rows+1, 13),] <- NA
        } else{
          print("")
        }
d.ss <- subset(df, subset = (Class == "SGLT2i" & outcome == "sbp" ))
        rows <- nrow(d.ss)
        if(rows<13){
          d.ss[seq(rows+1, 13),] <- NA
        } else{
          print("")
        }
#hem
d.shem <- subset(df, subset = (Class == "SGLT2i" & outcome == "hematocrit"))
        rows <- nrow(d.shem)
        if(rows<13){
          d.shem[seq(rows+1, 13),] <- NA
        } else{
          print("")
        }
#uacr
d.gu <- subset(df, subset = (Class == "GLP-1RA" & outcome == "uacr" ))
        rows <- nrow(d.gu)
        if(rows<13){
          d.gu[seq(rows+1, 13),] <- NA
        } else{
          print("")
        }



#hba1c ======

model.gh <- metagen(TE         = d.gh$diff,
                   seTE         = d.gh$diff.se,
                   data         = d.gh,
                   sm           = "MD",
                   studlab      = d.gh$trialname,
                   fixed        = FALSE,
                   random       = TRUE,
                   method.tau   = "REML",
                   overall      = FALSE,
                   hakn         = TRUE,
                   title        = "study weights"
)
# II.2 model.s for SGLT2i
model.sh <- metagen(TE           = d.sh$diff,
                   seTE         = d.sh$diff.se,
                   data         = d.sh,
                   sm           = "MD",
                   studlab      = d.sh$trialname,
                   fixed        = FALSE,
                   random       = TRUE,
                   method.tau   = "REML",
                   overall      = TRUE,
                   hakn         = TRUE,
                   title        = "study weights"
)
# III. Forest plot
# III.a GLP-1Ra
png(paste0("~/Documents/_mediation/forest_hba1c_glp1.png"),width = 6, height = 4, units = 'in', res = 300) 

forest(model.gh,
     #  atransf  = log, 
       overall = T,
       sortvar = ordertrial,
       studlab  = TRUE,
       fontsize = 9,
       xlim     = c(-6,2),
       at       = c(-2,-1,0,1,2),
       fs.axis  = 7,
       ref      = 0,
       level    = 0.95,
       weight.subgroup = "weight",
       smlab = "",
       col.diamond = "gray40",
       col.square  = "gold2",
       leftcols = c("studlab"),
       label = T,
       print.Q = T, print.pval.Q = T,
       print.I2 = T,
       print.tau2 = F,
       fs.hetstat = 7,
       leftlabs = paste0(unique(d.gh$outcometitle)[1], "\nGLP-1RA"),
       subgroup=F
)
dev.off()
# III.b SGLT2i
png(paste0("~/Documents/_mediation/forest_hba1c_sglt2.png"),width = 6, height = 4, units = 'in', res = 300) 
forest(model.sh,
       #  atransf  = log, 
       overall = T,
       #sortvar = ordertrial,
       studlab  = TRUE,
       fontsize = 9,
       xlim     = c(-6,2),
       at       = c(-2,-1,0,1,2),
       fs.axis  = 7,
       ref      = 0,
       level    = 0.95,
       weight.subgroup = "weight",
       smlab = "",
       col.diamond = "gold2",
       col.square  = "gray40",
       leftcols = c("studlab"),
       label = T,
       print.Q = T, print.pval.Q = T,
       print.I2 = T,
       print.tau2 = F,
       fs.hetstat = 7,
       leftlabs ="\nSGLT2i",
       subgroup=F
)
dev.off()

# BW========

model.gb <- metagen(TE         = d.gb$diff,
                    seTE         = d.gb$diff.se,
                    data         = d.gb,
                    sm           = "MD",
                    studlab      = d.gb$trialname,
                    fixed        = FALSE,
                    random       = TRUE,
                    method.tau   = "REML",
                    overall      = FALSE,
                    hakn         = TRUE,
                    title        = "study weights"
)
# II.2 model.s for SGLT2i
model.sb <- metagen(TE           = d.sb$diff,
                   seTE         = d.sb$diff.se,
                   data         = d.sb,
                   sm           = "MD",
                   studlab      = d.sb$trialname,
                   fixed        = FALSE,
                   random       = TRUE,
                   method.tau   = "REML",
                   overall      = TRUE,
                   hakn         = TRUE,
                   title        = "study weights"
)
# III. Forest plot
# III.a GLP-1Ra
png(paste0("~/Documents/_mediation/forest_bw_glp1.png"), width = 6, height = 4, units = 'in', res = 300) 

forest(model.gb,
       #  atransf  = log, 
       overall = T,
       studlab  = TRUE,
       sortvar = ordertrial,
       fontsize = 9,
       xlim     = c(-9,2),
       at       = c(-4,-2,0,2),
       fs.axis  = 7,
       ref      = 0,
       level    = 0.95,
       weight.subgroup = "weight",
       smlab = "",
       col.diamond = "gray40",
       col.square  = "gold2",
       leftcols = c("studlab"),
       label = T,
       print.Q = T, print.pval.Q = T,
       print.I2 = T,
       print.tau2 = F,
       fs.hetstat = 7,
       leftlabs =paste0(unique(d.gb$outcometitle)[1], "\nGLP-1RA"),
       subgroup=F
)
dev.off()
# III.b SGLT2i
png(paste0("~/Documents/_mediation/forest_bw_sglt2.png"), width = 6, height = 4, units = 'in', res = 300) 

forest(model.sb,
       #  atransf  = log, 
       overall = T,
       sortvar = ordertrial,
       studlab  = TRUE,
       fontsize = 9,
       xlim     = c(-9,2),
       at       = c(-4,-2,0,2),
       fs.axis  = 7,
       ref      = 0,
       level    = 0.95,
       weight.subgroup = "weight",
       smlab = "",
       col.diamond = "gold2",
       col.square  = "gray40",
       leftcols = c("studlab"),
       label = T,
       print.Q = T, print.pval.Q = T,
       print.I2 = T,
       print.tau2 = F,
       fs.hetstat = 7,
       leftlabs ="\nSGLT2i",
       subgroup=F
)
dev.off()

# SBP ==============

model.gs <- metagen(TE         = d.gs$diff,
                    seTE         = d.gs$diff.se,
                    data         = d.gs,
                    sm           = "MD",
                    studlab      = d.gs$trialname,
                    fixed        = FALSE,
                    random       = TRUE,
                    method.tau   = "REML",
                    overall      = FALSE,
                    hakn         = TRUE,
                    title        = "study weights"
)
# II.2 model.s for SGLT2i
model.ss<- metagen(TE           = d.ss$diff,
                    seTE         = d.ss$diff.se,
                    data         = d.ss,
                    sm           = "MD",
                    studlab      = d.ss$trialname,
                    fixed        = FALSE,
                    random       = TRUE,
                    method.tau   = "REML",
                    overall      = TRUE,
                    hakn         = TRUE,
                    title        = "study weights"
)
# III. Forest plot
# III.a GLP-1Ra
png(paste0("~/Documents/_mediation/forest_sbp_glp1.png"), width = 6, height = 4, units = 'in', res = 300) 

forest(model.gs,
       #  atransf  = log, 
       overall = T,
       sortvar = ordertrial,
       studlab  = TRUE,
       fontsize = 9,
       xlim     = c(-10,2),
       at       = c(-6,-4,-2,0,2),
       fs.axis  = 7,
       ref      = 0,
       level    = 0.95,
       weight.subgroup = "weight",
       smlab = "",
       col.diamond = "gray40",
       col.square  = "gold2",
       leftcols = c("studlab"),
       label = T,
       print.Q = T, print.pval.Q = T,
       print.I2 = T,
       print.tau2 = F,
       fs.hetstat = 7,
       leftlabs = paste0(unique(d.gs$outcometitle)[1], "\nGLP-1RA"),
       subgroup=F
)
dev.off()
# III.b SGLT2i
png(paste0("~/Documents/_mediation/forest_sbp_sglt2.png"), width = 6, height = 4, units = 'in', res = 300) 
forest(model.ss,
       #  atransf  = log, 
       overall = T,
       studlab  = TRUE,
       fontsize = 9,
       sortvar = ordertrial,
       xlim     = c(-10,2),
       at       = c(-6,-4,-2,0,2),
       fs.axis  = 7,
       ref      = 0,
       level    = 0.95,
       weight.subgroup = "weight",
       smlab = "",
       col.diamond = "gold2",
       col.square  = "gray40",
       leftcols = c("studlab"),
       label = T,
       print.Q = T, print.pval.Q = T,
       print.I2 = T,
       print.tau2 = F,
       fs.hetstat = 7,
       leftlabs ="\nSGLT2i",
       subgroup=F
)
dev.off()

# Hematocrit ==============


# II.2 model.s for SGLT2i
model.shem<- metagen(TE           = d.shem$diff,
                   seTE         = d.shem$diff.se,
                   data         = d.shem,
                   sm           = "MD",
                   studlab      = d.shem$trialname,
                   fixed        = FALSE,
                   random       = TRUE,
                   method.tau   = "REML",
                   overall      = TRUE,
                   hakn         = TRUE,
                   title        = "study weights"
)
# III. Forest plot
# III.a SGLT2
png(paste0("~/Documents/_mediation/forest_hem_sglt2.png"),  width = 6, height = 4, units = 'in', res = 300) 
forest(model.shem,
       #  atransf  = log, 
       overall = T,
       sortvar = ordertrial,
       studlab  = TRUE,
       fontsize = 9,
       xlim     = c(0,4),
       at       = c(1.5, 2.5,3,3.5),
       fs.axis  = 7,
       ref      = 0,
       level    = 0.95,
       weight.subgroup = "weight",
       smlab = "",
       col.diamond = "gold2",
       col.square  = "gray40",
       leftcols = c("studlab"),
       label = T,
       print.Q = T, print.pval.Q = T,
       print.I2 = T,
       print.tau2 = F,
       fs.hetstat = 7,
       leftlabs = paste0(unique(d.shem$outcometitle)[1], "\nSGLT2i"),
       subgroup=F,
)
dev.off()




# UACR ====

model.gu <- metagen(TE         = d.gu$diff,
                    seTE         = d.gu$diff.se,
                    data         = d.gu,
                    sm           = "MD",
                    studlab      = d.gu$trialname,
                    fixed        = FALSE,
                    random       = TRUE,
                    method.tau   = "REML",
                    overall      = FALSE,
                    hakn         = TRUE,
                    title        = "study weights"
)

model.su <- metagen(TE         = d.su$diff,
                    seTE         = d.su$diff.se,
                    data         = d.su,
                    sm           = "MD",
                    studlab      = d.su$trialname,
                    fixed        = TRUE,
                    random       = FALSE,
                    method.tau   = "REML",
                    overall      = FALSE,
                    hakn         = TRUE,
                    title        = "study weights"
)

# III. Forest plot
# III.a GLP-1Ra
png(paste0("~/Documents/_mediation/forest_uacr_glp1.png"), width = 6, height = 4, units = 'in', res = 300) 

forest(model.gu,
       #  atransf  = log, 
       overall = T,
       studlab  = TRUE,
       sortvar = ordertrial,
       fontsize = 9,
       xlim     = c(-15,2),
      at       = c(-10,-5,0),
       fs.axis  = 7,
       ref      = 0,
       level    = 0.95,
       weight.subgroup = "weight",
       smlab = "",
       col.diamond = "gray40",
       col.square  = "gold2",
       leftcols = c("studlab"),
       label = T,
       print.Q = T, print.pval.Q = T,
       print.I2 = T,
       print.tau2 = F,
       fs.hetstat = 7,
       leftlabs = paste0(unique(d.gu$outcometitle)[1], "\nGLP-1RA"),
       subgroup=F
)
dev.off()

png(paste0("~/Documents/_mediation/forest_uacr_sglt2.png"), width = 6, height = 4, units = 'in', res = 300) 

forest(model.su,
       #  atransf  = log, 
       overall = T,
       studlab  = TRUE,
       sortvar = ordertrial,
       fontsize = 9,
       #xlim     = c(-15,2),
       #   at       = c(0.5,0.75,1,1.5),
       fs.axis  = 7,
       ref      = 0,
       level    = 0.95,
       weight.subgroup = "weight",
       smlab = "",
       col.diamond = "gray40",
       col.square  = "gold2",
       leftcols = c("studlab"),
       label = T,
       print.Q = T, print.pval.Q = T,
       print.I2 = T,
       print.tau2 = F,
       fs.hetstat = 7,
       leftlabs = paste0(unique(d.gu$outcometitle)[1], "\nGLP-1RA"),
       subgroup=F
)
dev.off()




# Export  Figure 1: Panel forest plot of meta analysis of change in risk factors ====

#### Create panel of figuers
p1g <- image_read("~/Documents/_mediation/forest_hba1c_glp1.png")
p1s <- image_read("~/Documents/_mediation/forest_hba1c_sglt2.png")
p2g <- image_read("~/Documents/_mediation/forest_bw_glp1.png")
p2s <- image_read("~/Documents/_mediation/forest_bw_sglt2.png")
p3g <- image_read("~/Documents/_mediation/forest_sbp_glp1.png")
p3s <- image_read("~/Documents/_mediation/forest_sbp_sglt2.png")
p4g <- image_read("~/Documents/_mediation/forest_uacr_glp1.png")
p5s <- image_read("~/Documents/_mediation/forest_hem_sglt2.png")

p1g <- ggplot() + background_image(p1g) 
p1s <- ggplot() + background_image(p1s) 
p2g <- ggplot() + background_image(p2g) 
p2s <- ggplot() + background_image(p2s) 
p3g <- ggplot() + background_image(p3g) 
p3s <- ggplot() + background_image(p3s) 
p4g <- ggplot() + background_image(p4g) 
p5s <- ggplot() + background_image(p5s) 


p <- (p1g + p1s) / (p2g + p2s) /(p3g + p3s)/(p4g + p5s)


ggsave("~/Documents/_mediation/_output/_figures/figure_forest.png",plot = p, dpi = 300, width =  12, height =  16)



## Forest plot for UACR SGLT2 ====


d.su <- subset(df, subset = (Class == "SGLT2i" & outcome == "uacr" & trialname != " "))

model.su <- metagen(TE         = d.su$diff,
                    seTE         = d.su$diff.se,
                    data         = d.su,
                    sm           = "MD",
                    studlab      = d.su$trialname,
                    fixed        = TRUE,
                    random       = FALSE,
                    method.tau   = "REML",
                    overall      = FALSE,
                    hakn         = TRUE,
                    title        = "study weights"
)


png(paste0("~/Documents/_mediation/_output/_figures/figure_S_forest_uacr_sglt2.png"), width = 6, height = 4, units = 'in', res = 300) 

forest(model.su)
forest(model.su,
       #  atransf  = log, 
       overall = F,
       studlab  = TRUE,
       sortvar = ordertrial,
       fontsize = 9,
       #xlim     = c(-15,2),
       #   at       = c(0.5,0.75,1,1.5),
       fs.axis  = 7,
       ref      = 0,
       level    = 0.95,
       #weight.subgroup = "weight",
       #smlab = "",
       col.diamond = "gray40",
       col.square  = "gold2",
       leftcols = c("studlab"),
       rightcols = c("effect", "ci"),
       label = T,
       print.Q = F, print.pval.Q = F,
       print.I2 = F,
       print.tau2 = F,
       fs.hetstat = 7,
       leftlabs = paste0(unique(d.gu$outcometitle)[1], "\n SGLT2i"),
       subgroup=F
)
dev.off()








