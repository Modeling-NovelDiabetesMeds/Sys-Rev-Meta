#         #        #         #        #         #        #         #          # 
#                       Meta-regression analysis (loghr)
# Inputs  : df.csv 
# Outputs : - meta regression coefficients tables, and figures for loghr
#           and baseline cardiovascular risk
#           
#         #        #         #        #         #        #         #          #

  
library(data.table)
library(tidyr)
library(dplyr)
library(scales)
library(stringr)
library(metafor)
library(meta)
library(grid)
library(lme4)
library(stargazer)
library(knitr)

set.seed(10029)
wd <- "~/Documents/Sys-Rev-Meta"
setwd(wd)

# read complete meta dataset
meta <- fread("data/df.csv")
meta <- subset(meta, subset = (outcome != "macroalbuminuria"))
meta <- subset(meta, subset = (!(meta$outcome == "sustGFRdecl" & 
                                   (meta$trialname == "LEADER" | meta$trialname =="SUSTAIN-6"))))

# Create CVD baseline rate variable (events per 100 patient year)
meta$cvdepy <- NA
meta$cvdepy <- ifelse(meta$outcomen == "CVMort", 
                      meta$p.rate, 
                      meta$cvdepy)
# Other variables needed
meta$loghr <- log(meta$hr)
meta$logvi <- ((log(meta$uci) - log(meta$lci)) / (2 * qnorm(.975) ))^2
meta$logsei <- sqrt(meta$logvi)
meta$vi <- exp(meta$logvi)
meta$sei <- sqrt(meta$vi)

# Rename outcomes

meta$outcomelab <- meta$outcomen
meta$outcomelab <- ifelse(meta$outcomen == "allcauseMort", 
                          "All-cause mortality", meta$outcomelab )
meta$outcomelab <- ifelse(meta$outcomen == "CVMort", 
                          "Cardiovascular mortality", meta$outcomelab )
meta$outcomelab <- ifelse(meta$outcomen == "HospHF",
                          "Hospitalization for Heart Failure", meta$outcomelab )
meta$outcomelab <- ifelse(meta$outcomen == "sustGFRdecl", 
                          "Composite renal outcome", meta$outcomelab )
meta$outcomelab <- ifelse(meta$outcomen == "MACE", 
                          "MACE", meta$outcomelab )
meta$outcomelab <- ifelse(meta$outcomen == "MI", 
                          "Myocardial infarction", meta$outcomelab )
meta$outcomelab <- ifelse(meta$outcomen == "stroke", 
                          "Stroke", meta$outcomelab )


# Colors
meta$colnum2 <- ifelse(meta$type == "SGLT2", 
                       "darkorchid2", "cadetblue2")
meta$colnum2 <- alpha(meta$colnum2, 
                      alpha = 0.6)
colfit <- alpha(c("gray60", "darkorchid2", "cadetblue2"), 
                alpha = 0.6)


# complete rate for cvd for all outcomes
meta <- meta %>% 
  group_by(trialname) %>% 
  mutate(cvdepy = ifelse(is.na(cvdepy), 
                         median(cvdepy, na.rm = TRUE), 
                         cvdepy))
# Ranges
ylimm <- range(meta$loghr, na.rm =T) +  c(-1.1, 1.1)*sd(meta$loghr, na.rm =T)
xlimm <- range(meta$cvdepy, na.rm =T) + c(-0.5, 0)*sd(meta$cvdepy, na.rm =T)


meta$lograte <- log(meta$cvdepy) 
xlimm <- range(meta$lograte, na.rm =T) +   c(-0.5, 0)*sd(meta$lograte, na.rm =T)
xlimm0 <- range(meta$lograte[meta$type == "GLP1"], na.rm = T) + c(-0.75, 0.75)*sd(meta$lograte[meta$type == "GLP1"], na.rm =T)

# vector of outcome names to run loop
v.outcome <- unique(meta$outcome)
v.outcome <- c("allcauseMort", "CVMort", "MACE", "MI", "stroke", "HospHeartFailure", "sustGFRdecl")
coeftab <- as.data.frame(matrix(ncol = 6, nrow = 2*length(v.outcome)))


for( i in 1:length(v.outcome)){
  #subset dataset by outcome
  mt <- subset(meta, subset = (outcome == v.outcome[i]) )
  ran <-range((1/mt$logvi), na.rm = T)[2] - range((1/mt$logvi), na.rm = T)[1]
  mt$wsize <- (1/mt$logvi)/ran
  # metaregression, with baseline cvd rate as mediator
  r <- rma(loghr, logvi, 
           mods = ~lograte, 
           data = mt, 
           method="REML", 
           slab =trialname)
  # summary to recover quantities of interest (coefficient and pvalue)
  sr <- summary(r)
  names(coeftab) <- c("Outcome", "Class",  "Slope","CI.lb", "CI.ub", "P-value")
  beta <- round(sr$beta[2],4) # slope coefficient
  pp   <- round(sr$pval[2],3) # pvalue 
  
  # Separate analyses by sglt and glp types
  m.g <- mt[mt$type == "GLP1",]
  m.s <- mt[mt$type == "SGLT2",]
  # GLP1
  rg <- rma(loghr, logvi, 
            mods   = ~lograte, 
            data   = m.g, 
            method ="REML", 
            slab   =trialname
            )
  # summary to recover quantities of interest (coefficient and pvalue)
  sr <- summary(rg)
  coeftab[i,1] <- v.outcome[i]
  coeftab[i,2] <- "GLP1"
  coeftab[i,3] <- as.numeric(sr$beta[2])
  coeftab[i,4] <- as.numeric(sr$ci.lb[2])
  coeftab[i,5] <- as.numeric(sr$ci.ub[2])
  coeftab[i,6] <- as.numeric(sr$pval[2])

  # Plot
  logg <- log(seq(0.4, 1.6, 0.2)) # for Y axis jump in a not weird way
  par(oma = c(3,1,1,1), mfrow = c(1,2) )
  regplot(rg, pch = 19, col = m.g$colnum2, bg = "white",label =  "all", labsize = 0.75,
          lcol = c(colfit[3], "gray95", "gray95", "gray95", "gray95"), shade = c("gray97"), 
          xlab = " ", ylab = "",
          cex.axis = 1.0, las = 1,
          psize = m.g$wsize*9,
          xlim = xlimm0,
          ylim = ylimm,
          atransf = exp, at = logg
          )
  title(     xlab = "cardiovascular death events per 100 patient-yr, control group",
             ylab = paste0("Log of HR for ", v.outcome[i]),
             cex.lab = 1.0, line = 2.2)
  title(main =  "GLP- 1RAs trials", cex.main = 1)
  
  
  
  # SGLT
  rs <- rma(loghr, logvi, 
            mods = ~lograte, 
            data = m.s, 
            method="REML", 
            slab =trialname,
            # weights = weights
  )
  # summary to recover quantities of interest (coefficient and pvalue)
  sr <- summary(rs)
  
  coeftab[i+7,1] <- v.outcome[i]
  coeftab[i+7,2] <- "SGLT2i"
  coeftab[i+7,3] <- as.numeric(sr$beta[2])
  coeftab[i+7,4] <- as.numeric(sr$ci.lb[2])
  coeftab[i+7,5] <- as.numeric(sr$ci.ub[2])
  coeftab[i+7,6] <- as.numeric(sr$pval[2])
  
  regplot(rs, pch = 19, col = m.s$colnum2, bg = "white",label =  "all", labsize = 0.75,
          lcol = c(colfit[2], "gray95", "gray95", "gray95", "gray95"), shade = c("gray97"), 
          xlab = " ", ylab = "",
          cex.axis = 1.0, las = 1,
          psize = m.s$wsize*9,
          xlim = xlimm,
          ylim = ylimm,
          atransf = exp,
          axes= F, at = logg
  )
  
  title(     xlab = "cardiovascular death events per 100 patient-yr, control group",
             ylab = paste0("Log of HR for ", v.outcome[i]),
             cex.lab = 1.0, line = 2)
  title(main = "SGLT2i", cex.main = 1)
  
  title(main = paste0(unique(mt$outcomelab) ),
        cex.main = 1.2, outer = TRUE)
}






 stargazer(coeftab, out = "output/metareg_loghr.txt",
          summary = F,type = "text", 
          title = "Meta-regression coefficients, by drug class",
          notes = "Log hazard ratio and baseline cardiovascular mortality rate")



 
 
 
 
 
 
