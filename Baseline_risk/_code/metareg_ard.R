#         #        #         #        #         #        #         #          # 
#                       Meta-regression analysis
# Inputs  : df.csv 
# Outputs : - meta regression coefficients tables, and figures
#           for main outcomes, secondary outcomes and sensitivity analyses.
#           - ardse.csv contains per trial and outcome the simulated standard
#             errors for ARD point estimate after using bootstrapping
#         #        #         #        #         #        #         #          #

    rm(list = ls())
# preamble, wd and load data ----
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

    set.seed(10029)
    wd <- "~/Documents/Sys-Rev-Meta"
   
    
      setwd(wd)
    
    meta <- read.csv("_data/df.csv")

# A.   Data preparation (labels and ids for plots) -----
## generate id per trial/group for regression plots
    meta <- meta[order(meta$trialname),]
    meta <- meta %>% 
      group_by(type, outcome) %>%
      mutate(id = row_number())  
    meta$id <- ifelse(meta$trialname == "LEADER/SUSTAIN-6", 10, meta$id)
    id0 <- median(meta$id[meta$trialname =="SUSTAIN-6"])
    meta$id <- ifelse(meta$trialname == "SUSTAIN-6", id0, meta$id)

## Outcome labels, used for later in plot
    meta$outcomelab <- meta$outcomen
    meta$outcomelab <- ifelse(meta$outcomen == "CVMort", 
                              "Cardiovascular Mortality", meta$outcomelab)
    meta$outcomelab <- ifelse(meta$outcomen == "HospHF",
                              "Hospitalization for Heart Failure", 
                              meta$outcomelab )
    meta$outcomelab <- ifelse(meta$outcomen == "sustGFRdecl", 
                              "Composite Renal Outcome", meta$outcomelab)
    meta$outcomelab <- ifelse(meta$outcomen == "MACE", 
                              "MACE", meta$outcomelab )
    meta$outcomelab <- ifelse(meta$outcomen == "allcauseMort", 
                              "All-cause Mortality", meta$outcomelab)
    meta$outcomelab <- ifelse(meta$outcomen == "MI", 
                              "Myocardial Infarction", meta$outcomelab)
    meta$outcomelab <- ifelse(meta$outcomen == "stroke", 
                              "Stroke", meta$outcomelab )

# B.   Computing variables of interest (loghr)       ----
## Create CVD baseline rate variable (events per 100 patient year)
##       cvdepy: cardiovascular death events pero 100 patient-years
##       from p.rate (placebo group rate per outcome)  
    meta$cvdepy <- NA
    meta$cvdepy <- ifelse(meta$outcomen == "CVMort", 
                          meta$p.rate, 
                          meta$cvdepy)
  # assign value of cardiovascular rate to all outcomes
      meta <- meta%>%
        group_by(trialname) %>% 
        mutate(cvdepy = ifelse(is.na(cvdepy), 
                               median(cvdepy, na.rm = TRUE), 
                               cvdepy))

## Hazard Ratios: log, variance and standard errors
  meta$loghr  <- log(meta$hr)
  meta$logvi  <- ((log(meta$uci) - log(meta$lci)) / (2 * qnorm(.975) ))^2
  meta$logsei <- sqrt(meta$logvi)
  meta$vi     <- (exp(meta$logvi)-1)*exp(2*meta$loghr + meta$logvi) # log normal variance
  meta$sei    <- sqrt(meta$vi)
  
## For EMPA KIDNEY we can compute t.rate and p.rate from HFCV rates  
    hcv.t.rate <- 2.04      # HHFCV/100py treatment 
    hcv.p.rate <- 2.37      # HHFCV/100py placebo 
    hcv.t.n <- 131          # HHFCV No. treatment
    hcv.p.n <- 152          # HHFCV No. placebo
    hf.t.n <- 88            # HHF No. treatment
    hf.p.n <- 107           # HHF No. placebo
    ### HF rates for EMPA KIDNEY
    hf.t.rate <- hcv.t.rate * (hf.t.n / hcv.t.n )
    hf.p.rate <- hcv.p.rate * (hf.p.n / hcv.p.n )
    meta$t.rate <- ifelse(meta$trialname == "EMPA-KIDNEY" & 
                          meta$outcome == "HospHeartFailure",
                          hf.t.rate,
                          meta$t.rate
                          )
    meta$p.rate <- ifelse(meta$trialname == "EMPA-KIDNEY" & 
                          meta$outcome == "HospHeartFailure",
                          hf.p.rate,
                          meta$p.rate
                          )
    
# Average risk difference (point estimates)
    # 5 year projection
     Ty <- 5
    # Rates - 5 years
     meta$baserateT    <- (1-exp(-meta$p.rate/100*Ty))
     meta$rateratio    <- meta$t.rate / meta$p.rate
     meta$treatrateT   <- (1-exp(-meta$rateratio*meta$p.rate/100*Ty))
    # Average risk difference
     meta$ard          <- meta$baserateT - meta$treatrateT
    # Number needed to treat
     meta$NNT          <- 1/meta$ard

# C.   Simulate variance/se for ARD point estimates -----
# Variables needed for simulation: log rate ratio and patient-years per arm
    meta$rateratio    <- meta$t.rate / meta$p.rate
    meta$se.rateratio <- sqrt(1/meta$n0 + 1/meta$n1) #se of ln of rate ratio
    meta$lograteratio <- log(meta$rateratio)
    meta$vi.rateratio <- meta$se.rateratio^2 

    meta$npy1 <- 100*(meta$n1 / meta$t.rate) # Total patient years treat
    meta$npy0 <- 100*(meta$n0 / meta$p.rate) # Total patient years control

    # List to store results
    lse <- list()
    # Vector of outcomes interest
    v.outcomes <- c("Cardiovascular Mortality", "MACE",
                    "Hospitalization for Heart Failure", 
                    "Composite Renal Outcome", "All-cause Mortality", 
                    "Myocardial Infarction", "Stroke") 
    # Number of simulations S
    S <- 100000
# Run simulation exercise for each outcome
      for(j in 1:length(v.outcomes)){
        m0 <- subset(meta, subset = (outcomelab == v.outcomes[j]))
        # Generate variance for ARD
        # Simulate a data frame vectors size S for p.rate and Ty = 5 years to get a 
        # simulated vector of arr and compute its variance
        # for each study generate p.rate hr and arr matrix to compute variance of ard
        lsim <- list()
        v.trial <- unique(m0$trialname)
        vector <- rep(NA, length(v.trial))
        se <- as.data.frame(cbind(vector,vector,vector))
        names(se) <- c("trialname", "ardse", "ardse2")
        # Simulations    
        for(i in 1:length(v.trial)){
          m <- subset(m0, subset = (trialname == v.trial[i]))
          # simulate random variables
          # simulate t.rate vector (Poisson)
          v.rt <- rpois(S, lambda = m$n1) / m$npy1 
          # simulate p.rate vector  (Poisson)
          v.rc <- rpois(S, lambda = m$n0) / m$npy0 
          # simulate rate ratio vector (log normal)
          v.rr <- rlnorm(S,  mean= m$lograteratio, sd = m$se.rateratio) 
          x <- list(vrc = v.rc, vrt = v.rt, vrr = v.rr)
          lapply(x, mean)
          sim.df <- as.data.frame(cbind(v.rt,  v.rc, v.rr))
          # Method 1 (using rates)
          sim.df$rate.t.5   <- 1-exp(-Ty * sim.df$v.rt) 
          # Method 2, using rate ratio
          sim.df$rate.t.5.2 <- 1-exp(-Ty * exp(sim.df$v.rr)*sim.df$v.rc) 
          # Rate 5 years for control
          sim.df$rate.c.5   <- 1-exp(-Ty * sim.df$v.rc) 
          
          # Compute ARD
          sim.df$ard   <- sim.df$rate.c.5  - sim.df$rate.t.5     # Method 1
          sim.df$ard.2 <- sim.df$rate.c.5  - sim.df$rate.t.5.2   # Method 2
          # List to compare means
          x <- list(vrc   = v.rc, 
                    vrt   = v.rt, 
                    rt5   = sim.df$rate.t.5, 
                    rt5.2 = sim.df$rate.t.5.2 , 
                    rc5   = sim.df$rate.c.5, 
                    ard   = sim.df$ard, 
                    ard.2 = sim.df$ard.2
          )
          lapply(x, mean)
          # Recover quantities of interest 
          sdard   <- sd(sim.df$ard)
          sdard.2 <- sd(sim.df$ard.2)
          se[i,1] <- m$trialname
          se[i,2] <- sd(sim.df$ard)
          se[i,3] <- sd(sim.df$ard.2)
        }
        
        se$ardvi <- se$ardse^2
        se$ardvi.2  <- se$ardse2^2
        se$outcomelab <- v.outcomes[j]
        lse[[j]] <- se
      }
# Merge variances back to main dataframe
# Dataframe of ard variances per trial/outcome
    mse  <- rbind(lse[[1]],lse[[2]],lse[[3]], lse[[4]], lse[[5]], lse[[6]], lse[[7]])
    # merge into main dataframe: 
    meta <- merge(meta, mse, by = c("trialname", "outcomelab"), 
                  all.x = TRUE, 
                  all.y= TRUE)
    arddf <- meta[, c("trialname","type", "outcomelab","ordertrial", "ard", "ardse",
                      "ardse2", "ardvi", "ardvi.2", "cvdepy", "outcome", "outcomen", "id")]
    arddf$ardse2 <- ifelse(arddf$trialname == "DELIVER",  arddf$ardse,  arddf$ardse2)
    arddf$ardse2 <- ifelse(arddf$trialname == "EMPA-KIDNEY",  arddf$ardse,  arddf$ardse2)
    
    # Export ardse.csv dataframe (used in forest plots)  
    write.csv(arddf, "_data/ardse.csv")


# -----------------------------------------------------------------------------#    
#                     Meta regression analyses                                 #    
# -----------------------------------------------------------------------------#        
    
meta <- read.csv("_data/ardse.csv")
# Wsize/ Wsize2: Ad-hoc Graphical Parameter used for bubble sizes, 
#                function proportional of variances
    meta <- meta %>%
        group_by(outcome) %>%
        mutate(ran = range(1/(10000*ardvi.2), na.rm = T)[2] -
               range(1/(10000*ardvi.2), na.rm = T)[1])
    meta$wsize   <- 1/(100*meta$ardvi.2)/meta$ran
    meta$wsize2  <- 8*log(1+meta$wsize^(1/3) ) 

# I.1  Meta regression: Models with cvdepy as mediator: -----
  meta$lograte <- log(meta$cvdepy) 
#  7 outcomes, model by drug class
  # lists to store models
    mod.glp1  <- list() 
    mod.sglt2 <- list()
    v.type <- c("GLP1","SGLT2")

for(o in 1:length(v.outcomes)){ #loop runs over each outcome 'o' and each drug class 't'
  
  for(t in 1:length(v.type)){
    df.0 <- subset(meta, subset = (outcomelab == v.outcomes[o] &
                                   type == v.type[t]))
  
  if(v.outcomes[o] == "Hospitalization for Heart Failure"){ #leave out two trials with inconsistent definition for HospHF
    df <- subset(df.0, 
                 subset = (trialname != "SOLOIST-WHF" & trialname != "SCORED") 
                 )
  } else {
    df <- df.0
  }

  # Model (random effects) for outcome o class t
  r <- rma(-100*ard, (100^2)*ardvi.2, # times 100 to re-scale to percentage points
           mods = ~ cvdepy, 
           data = df, 
           method="REML", 
           slab = trialname)
  if( t == 1){ #Store in corresponding type-list (GLP1-RA/SGLT2i)
    mod.glp1[[o]]  <- r
    } else {
    mod.sglt2[[o]] <- r
    }
  }
  }

  # summary to recover quantities of interest (coefficient and pvalue)
    s.s <- lapply(mod.sglt2, summary)
    s.g <- lapply(mod.glp1, summary)
#primary outcomes

    # empty data frame to store
    cc <- as.data.frame(cbind(rep(v.outcomes[1:7],2), 
                          c(rep("GLP-1RA",7), 
                            rep("SGLT2i",7)), 
                          rep(NA,14), 
                          rep(NA,14),
                          rep(NA,14),
                          rep(NA,14)))

names(cc)   <- c("Outcome", "Class", "Slope","CI.lb", "CI.ub", "P-value")

    # extract quantities of interest summary lists
      # Primary outcomes
      for(i in 1:7){
        cc[i,3] <- round(s.g[[i]]$beta[2],2)
        cc[i,4] <- round(s.g[[i]]$ci.lb[2],2)
        cc[i,5] <- round(s.g[[i]]$ci.ub[2],2)
        cc[i,6] <- round(s.g[[i]]$pval[2],2)
        cc[i+7,3] <- round(s.s[[i]]$beta[2],2)
        cc[i+7,4] <- round(s.s[[i]]$ci.lb[2],2)
        cc[i+7,5] <- round(s.s[[i]]$ci.ub[2],2)
        cc[i+7,6] <- round(s.s[[i]]$pval[2],2)
      }
  # Sort in order of outcomes for tables
  ccorder <- c(5,1,2,6,7,3,4,12,8,9,13,14,10,11)
  cc <- cc[ccorder,]
  cc
  cc[,c(3,4,5,6)] <- round(cc[,c(3,4,5,6)],2)
  # Export into text tables
    stargazer(cc, out = "_output/metareg_ard.txt",
              summary = F,type = "text", 
              title = "Meta-regression coefficients, by drug class",
              notes = "Absolute risk difference and baseline cardiovascular mortality rate")

    
# I.2  Meta regression: Figure ------
    # Colors by drug class 
    meta$border <- ifelse(meta$type == "SGLT2", 
                       "darkslategray3", 
                       "indianred3")
    meta$colnum <- ifelse(meta$type == "SGLT2", 
                           "darkslategray3", 
                           "transparent")
    meta$colnum2 <- ifelse(meta$type == "SGLT2", alpha(meta$colnum,alpha = 0.4),
                           "transparent")
  # A list to store data frames per outcome, to run plot over a loop.
    # Data frames contain data points for bubbles (x = cvdepy, y = ard, weight)
    l.do <- list()
  # 1. Df for each outcome, 
  l.do[[1]]  <- meta[meta$outcomen == "CVMort",]
  l.do[[2]]  <- meta[meta$outcomen == "MACE",]
  l.do[[3]]  <- meta[meta$outcomen == "HospHF"& 
                       (meta$trialname != "SCORED" & meta$trialname != "SOLOIST-WHF"),]
        ###  Remove LEADER and SUSTAIN from Kidney
  l.do[[4]]  <- meta[meta$outcomen == "sustGFRdecl" & 
                       (meta$trialname != "LEADER" & meta$trialname != "SUSTAIN 6"),]
  # Object to title subplots
  v.titles <- c("Cardiovascular mortality",
                "MACE", 
                "Hospitalization for heart failure",
                "Composite renal outcome"
                )


  
#  png("_output/Figure3.png", width = 9, height = 9, units = 'in', res = 300)  
  par(oma = c(2,1,1,1), mfrow = c(2,2), mar = c(5,5,3,2)*0.75)
  for(i in 1:4){
    # Plot elements: 1. white canvas, 2. confint (sglt2 + glp1) + 
    #                3. bubbles 4. Model fit (line) 5. Trial number (text)
    #                6. titles (outcome)
  # 1. White canvas 
    plot(x = 0, y = 0, type ='n', 
         xlim = range(meta$cvdepy, na.rm = T), 
         ylim = c(-14,4), 
         xlab= "", 
         ylab ="", 
         axes = F)
  # 2. Confidence interval: predict confint for both glp-1ras and sglt2i models
    # create domain for independent variable
    newx.glp <- seq(min(l.do[[i]]$cvdepy[l.do[[i]]$type == "GLP1" & 
                                           !is.na(l.do[[i]]$ard)], na.rm = T),
                    max(l.do[[i]]$cvdepy[l.do[[i]]$type == "GLP1" & 
                                           !is.na(l.do[[i]]$ard)], na.rm = T),
                length.out = 1000)
    # predict glp1 ard for given x1 ( ard =f(x1) )
    preds.glp <- predict(mod.glp1[[i]], 
                     newmods = newx.glp)
    # create domain for independent variable
    newx.sglt <- seq(min(l.do[[i]]$cvdepy[l.do[[i]]$type == "SGLT2" & 
                                            !is.na(l.do[[i]]$ard)], na.rm = T),
                     max(l.do[[i]]$cvdepy[l.do[[i]]$type == "SGLT2" & 
                                            !is.na(l.do[[i]]$ard)], na.rm = T),
                    length.out = 1000)
    # predict sglt2 ard for given x2 ( ard =f(x2) )
    preds.sglt <- predict(mod.sglt2[[i]], 
                         newmods = newx.sglt)
      # Create CI polygons with computed coordinates
    polygon(x = c(newx.sglt, rev(newx.sglt)), 
            y = c(preds.sglt$ci.lb, rev(preds.sglt$ci.ub)),
            col = alpha("darkslategray3", alpha = 0.4),
            border = "transparent")
    polygon(x = c(newx.glp, rev(newx.glp)), 
            y = c(preds.glp$ci.lb, rev(preds.glp$ci.ub)),
            col = alpha("indianred3", alpha = 0.4),
            border = "transparent")
  # 3. Point data (bubbles) with size = f(1/var)
  points(x =   l.do[[i]]$cvdepy, 
         y =   -100*l.do[[i]]$ard, # 100 to make it percentaje points
         col = l.do[[i]]$border,
         bg =  l.do[[i]]$colnum2,
         pch = 21, 
         cex = ((0.25)*l.do[[i]]$wsize2)^(1.5))
  # Model for GLP1-RA
  segments(x0 = min(newx.glp), 
           x1 = max(newx.glp),
           y0 = mod.glp1[[i]]$beta[1] + 
             min(newx.glp)*mod.glp1[[i]]$beta[2], 
           y1 = mod.glp1[[i]]$beta[1] + 
             max(newx.glp)*mod.glp1[[i]]$beta[2],
           col =  "indianred3",
           lwd = 2)
  # 4. Model for SGLT2i
  segments(x0 = min(newx.sglt), 
           x1 = max(newx.sglt),
           y0 = mod.sglt2[[i]]$beta[1] + 
             min(newx.sglt)*mod.sglt2[[i]]$beta[2], 
           y1 = mod.sglt2[[i]]$beta[1] + 
             max(newx.sglt)*mod.sglt2[[i]]$beta[2],
           col =  "darkslategray4",
           lwd = 2)
  # 4. Add id to identify trial
  text(x = jitter(l.do[[i]]$cvdepy[!is.na(l.do[[i]]$ard)],20), 
       y = jitter(-100*l.do[[i]]$ard[!is.na(l.do[[i]]$ard)],20),  
       l.do[[i]]$id[!is.na(l.do[[i]]$ard)], cex = 0.65 )
  # 5. title label per outcome
  title(main = v.titles[i], cex.main = 0.9, line = 0.7)
  # Axes titles

  title(ylab = "Absolute risk difference (%)",
        xlab = "Cardiovascular mortality rate (/100 py) in control group", 
        cex.lab = 0.8, line = 2.3)

  # Axes
  axis(1, 
       at = pretty(range(meta$cvdepy, na.rm = T)),
       cex.axis = 0.75,
       col.axis = "gray40", 
       col.ticks = "gray40", 
       col = "gray40", 
       las = 1)
  axis(2, 
       #at = pretty(range(100*meta$ard, na.rm = T), n = 5),
       at  = c(4, 0,-4, -8, -10,-14),
       labels = c(4, 0,-4, -8,-10, -14),
       cex.axis = 0.75, 
       col.axis = "gray40", 
       col.ticks = "gray40", 
       col = "gray40", 
       las = 1)  
  if(i == 2){
    legend("bottomright", c("GLP-1RA", "SGLT2i"), 
           cex = 0.85, 
           pch = 21, 
           bg = "transparent",
           pt.bg = c("transparent", "darkslategray3" ),
           col = c("indianred3", "darkslategray3" ), 
           box.col = "transparent")
  }else{
    print("hi")

  }
  abline(h = 0, lty = 2, col = "gray60")
  text( x = 13.8, y = 3.2, labels = "Favors \n control", cex = 0.8, srt = 90)
  text( x = 13.8, y = -3.7, labels = "Favors \n treatment", cex = 0.8, srt = 90)
  
}
#  dev.off()  


# I.3  Meta regression: Figure for secondary outcomes -----

  # A list to store data frames per secondary outcome, to run plot over a loop.
  # Data frames contain data points for bubbles (x = cvdepy, y = ard, weight)
  l.do.2 <- list()
  # 1. Df for each outcome, 
  l.do.2[[1]]  <- meta[meta$outcomen == "allcauseMort",]
  l.do.2[[2]]  <- meta[meta$outcomen == "MI",]
  l.do.2[[3]]  <- meta[meta$outcomen == "stroke",]

  # Object to title subplots
  v.titles <- v.outcomes[5:7]
  
  #png("_output/metareg_ard_panel_secondary.png", width = 15, height = 5, units = 'in', res = 300)  
  par(oma = c(2,1,0,0), mfrow = c(1,3), mar = c(5,5,3,2)*0.8)
    for(i in 5:7){
    # Plot elements: 1. white canvas, 2. confint (sglt2 + glp1) + 
    #                3. bubbles 4. Model fit (line) 5. Trial number (text)
    #                6. titles (outcome)
    # 1. White canvas 
    plot(x = 0, y = 0, type ='n', 
         xlim = range(meta$cvdepy, na.rm = T), 
         ylim = c(-14,4), 
         xlab= "", 
         ylab ="", 
         axes = F)
    # 2. Confidence interval: predict confint for both glp-1ras and sglt2i models
    # create domain for independent variable
    newx.glp <- seq(min(l.do.2[[i-4]]$cvdepy[l.do.2[[i-4]]$type == "GLP1" & 
                                           !is.na(l.do.2[[i-4]]$ard)], na.rm = T),
                    max(l.do.2[[i-4]]$cvdepy[l.do.2[[i-4]]$type == "GLP1" & 
                                           !is.na(l.do.2[[i-4]]$ard)], na.rm = T),
                    length.out = 1000)
    # predict glp1 ard for given x1 ( ard =f(x1) )
    preds.glp <- predict(mod.glp1[[i]], 
                         newmods = newx.glp)
    # create domain for independent variable
    newx.sglt <- seq(min(l.do.2[[i-4]]$cvdepy[l.do.2[[i-4]]$type == "SGLT2" & 
                                            !is.na(l.do.2[[i-4]]$ard)], na.rm = T),
                     max(l.do.2[[i-4]]$cvdepy[l.do.2[[i-4]]$type == "SGLT2" & 
                                            !is.na(l.do.2[[i-4]]$ard)], na.rm = T),
                     length.out = 1000)
    # predict sglt2 ard for given x2 ( ard =f(x2) )
    preds.sglt <- predict(mod.sglt2[[i]], 
                          newmods = newx.sglt)
    # Create CI polygons with computed coordinates
    polygon(x = c(newx.sglt, rev(newx.sglt)), 
            y = c(preds.sglt$ci.lb, rev(preds.sglt$ci.ub)),
            col = alpha("darkslategray3", alpha = 0.4),
            border = "transparent")
    polygon(x = c(newx.glp, rev(newx.glp)), 
            y = c(preds.glp$ci.lb, rev(preds.glp$ci.ub)),
            col = alpha("indianred3", alpha = 0.4),
            border = "transparent")
    # 3. Point data (bubbles) with size = f(1/var)
    points(x =   l.do.2[[i-4]]$cvdepy, 
           y =   -100*l.do.2[[i-4]]$ard, # 100 to make it percentaje points
           col = l.do.2[[i-4]]$border,
           bg =  l.do.2[[i-4]]$colnum2,
           pch = 21, 
           cex = ((0.25)*l.do.2[[i-4]]$wsize2)^(1.5))
    # Model for GLP1-RA
    segments(x0 = min(l.do.2[[i-4]]$cvdepy[l.do.2[[i-4]]$type == "GLP1"], na.rm = T), 
             x1 = max(l.do.2[[i-4]]$cvdepy[l.do.2[[i-4]]$type == "GLP1"], na.rm = T),
             y0 = mod.glp1[[i]]$beta[1] + 
               min(l.do.2[[i-4]]$cvdepy[l.do.2[[i-4]]$type == "GLP1"], na.rm = T)*mod.glp1[[i]]$beta[2], 
             y1 = mod.glp1[[i]]$beta[1] + 
               max(l.do.2[[i-4]]$cvdepy[l.do.2[[i-4]]$type == "GLP1"], na.rm = T)*mod.glp1[[i]]$beta[2],
             col =  "indianred3",
             lwd = 2)
    # 4. Model for SGLT2i
    segments(x0 = min(l.do.2[[i-4]]$cvdepy[l.do.2[[i-4]]$type == "SGLT2"& !is.na(l.do.2[[i-4]]$ard)], na.rm = T), 
             x1 = max(l.do.2[[i-4]]$cvdepy[l.do.2[[i-4]]$type == "SGLT2"& !is.na(l.do.2[[i-4]]$ard)], na.rm = T),
             y0 = mod.sglt2[[i]]$beta[1] + 
               min(l.do.2[[i-4]]$cvdepy[l.do.2[[i-4]]$type == "SGLT2"& !is.na(l.do.2[[i-4]]$ard)], na.rm = T)*mod.sglt2[[i]]$beta[2], 
             y1 = mod.sglt2[[i]]$beta[1] + 
               max(l.do.2[[i-4]]$cvdepy[l.do.2[[i-4]]$type == "SGLT2"& !is.na(l.do.2[[i-4]]$ard)], na.rm = T)*mod.sglt2[[i]]$beta[2],
             col =  "darkslategray4",
             lwd = 2)
    # 4. Add id to identify trial
    text(x = jitter(l.do.2[[i-4]]$cvdepy[!is.na(l.do.2[[i-4]]$ard)],20), 
         y = jitter(-100*l.do.2[[i-4]]$ard[!is.na(l.do.2[[i-4]]$ard)],20),  
         l.do.2[[i-4]]$id[!is.na(l.do.2[[i-4]]$ard)], cex = 0.65 )
    # 5. title label per outcome
    title(main = v.titles[i-4], cex.main = 1.5, line = 0.1)
    # Axes titles

      title(ylab = "Absolute risk difference (%)", cex.lab = 1, line = 2)
      title(xlab = "Cardiovascular mortality rate (/100 py) in control group", cex.lab = 1, line = 2.5)
    # Axes
    axis(1, 
         at = pretty(range(meta$cvdepy, na.rm = T)),
         cex.axis = 0.75,
         col.axis = "gray40", 
         col.ticks = "gray40", 
         col = "gray40", 
         las = 1)
    axis(2, 
         #at = pretty(range(100*meta$ard, na.rm = T), n = 5),
         at  = c(4, 0,-4, -8, -10,-14),
         labels = c(4, 0,-4, -8, -10,-14),
         cex.axis = 0.75, 
         col.axis = "gray40", 
         col.ticks = "gray40", 
         col = "gray40", 
         las = 1)  
    if(i == 7){
      legend("bottomright", c("GLP-1RA", "SGLT2i"), 
             cex = 0.85, 
             pch = 21, 
             bg = "transparent",
             pt.bg = c("transparent", "darkslategray3" ),
             col = c("indianred3", "darkslategray3" ), 
             box.col = "transparent")
    }else{
      print("hi")
    }
  }
 # dev.off() 
# =========================================================================#
#                   Sensitivity Analyses                                   #  
# =========================================================================#
  
# S.1  Sensitivity Analysis CVD: Remove SCORE and SOLOIST --------
v.titles <- v.titles[1:4]
df.sens <- subset(meta, subset = (outcomelab == v.outcomes[o] & type == v.type[t] & 
                                    !(trialname == "SOLOIST-WHF"| trialname == "SCORED"))
                  )
                    
df.sens  <- meta[meta$outcomen == "CVMort" & 
                     !(meta$trialname == "SOLOIST-WHF" |
                       meta$trialname == "SCORED"),]
  
msens.g <- rma(100*ard, (100^2)*ardvi.2, 
                mods = ~ cvdepy, 
                data = df.sens[df.sens$type == "GLP1",], 
                method="REML", 
                slab = trialname)
msens.s <- rma(100*ard, (100^2)*ardvi.2, 
                mods = ~ cvdepy, 
                 data = df.sens[df.sens$type == "SGLT2",], 
                 method="REML", 
                 slab = trialname)


# S.1 Sensitivity Analysis all outcomes: Remove Score and soloist ------
    modsens.glp1  <- list() 
    modsens.sglt2 <- list()
    l.do.sens <- list()
    v.type <- c("GLP1","SGLT2")
            
    for(o in 1:length(v.outcomes)){ #loop runs over each outcome 'o' and each drug class 't'
          
      for(t in 1:length(v.type)){
            
        df.sens <- subset(meta, subset = (outcomelab == v.outcomes[o] & type == v.type[t] & 
                                        !(trialname == "SOLOIST-WHF"| trialname == "SCORED"))
                          )
        l.do.sens[[o]] <- subset(meta, subset = (outcomelab == v.outcomes[o] & 
                                                   !(trialname == "SCORED"|
                                                       trialname == "SOLOIST-WHF"  ))
        )
        
       # Model (random effects) for outcome o class t
            r <- rma(-100*ard, (100^2)*ardvi.2, # times 100 to re-scale to percentage points
                     mods = ~ cvdepy, 
                     data = df.sens, 
                     method="REML", 
                     slab = trialname)
            if( t == 1){ #Store in corresponding type-list (GLP1-RA/SGLT2i)
              modsens.glp1[[o]]  <- r
            } else {
              modsens.sglt2[[o]] <- r
            }
          }
        }
        
        # summary to recover quantities of interest (coefficient and pvalue)
        s.s <- lapply(modsens.sglt2, summary)
        s.g <- lapply(modsens.glp1, summary)
        #primary outcomes
        
        # empty data frame to store
        cc.s1 <- as.data.frame(cbind(rep(v.outcomes[1:7],2), 
                                    c(rep("GLP-1RA",7), 
                                      rep("SGLT2i",7)), 
                                    rep(NA,14), 
                                    rep(NA,14),
                                    rep(NA,14),
                                    rep(NA,14)))
        
        names(cc.s1)   <- c("Outcome", "Class", "Slope","CI.lb", "CI.ub", "P-value")
        
        # extract quantities of interest summary lists
        # Primary outcomes
        for(i in 1:7){
          cc.s1[i,3] <- round(s.g[[i]]$beta[2],2)
          cc.s1[i,4] <- round(s.g[[i]]$ci.lb[2],2)
          cc.s1[i,5] <- round(s.g[[i]]$ci.ub[2],2)
          cc.s1[i,6] <- round(s.g[[i]]$pval[2],2)
          cc.s1[i+7,3] <- round(s.s[[i]]$beta[2],2)
          cc.s1[i+7,4] <- round(s.s[[i]]$ci.lb[2],2)
          cc.s1[i+7,5] <- round(s.s[[i]]$ci.ub[2],2)
          cc.s1[i+7,6] <- round(s.s[[i]]$pval[2],2)
        }
        cc.s1 <- cc.s1[ccorder,]
        cc.s1
        # Export into text tables
        stargazer(cc.s1, 
                  out = "_output/metareg_ard_s1.txt",
                  summary = F,type = "text", 
                  title = "Meta-regression coefficients, by drug class, S1 - SGLTi1/2",
                  notes = "Absolute risk difference and baseline cardiovascular mortality rate. Sensitivity")
        
      png("metareg_ard_cdv_1.png", width = 9, height = 9, units = 'in', res = 300) 
      par(oma = c(2,1,1,1), mfrow = c(2,2), mar = c(5,5,3,2)*0.75)
      for(i in 1:4){
        # Plot elements: 1. white canvas, 2. confint (sglt2 + glp1) + 
        #                3. bubbles 4. Model fit (line) 5. Trial number (text)
        #                6. titles (outcome)
        # 1. White canvas 
        
        plot(x = 0, y = 0, type ='n', 
             xlim = c(0,14), 
             ylim = c(-14,4), 
             xlab= "", 
             ylab ="", 
             axes = F)
        abline(h = 0, lty = 2, col = "gray60")
        text( x = 13.5, y = 2.5, labels = "Favors \n control", cex = 1.0, srt = 90)
        text( x = 13.5, y = -2.5, labels = "Favors \n treatment", cex = 1.0, srt = 90)
        
        # 2. Confidence interval: predict confint for both glp-1ras and sglt2i models
        # create domain for independent variable
        newx.glp <- seq(min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1" & 
                                                    !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                        max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1" & 
                                                    !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                        length.out = 1000)
        # predict glp1 ard for given x1 ( ard =f(x1) )
        preds.glp <- predict(modsens.glp1[[i]], 
                             newmods = newx.glp)
        # create domain for independent variable
        newx.sglt <- seq(min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2" & 
                                                     !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                         max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2" & 
                                                     !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                         length.out = 1000)
        # predict sglt2 ard for given x2 ( ard =f(x2) )
        preds.sglt <- predict(modsens.sglt2[[i]], 
                              newmods = newx.sglt)
        # Create CI polygons with computed coordinates
        polygon(x = c(newx.sglt, rev(newx.sglt)), 
                y = c(preds.sglt$ci.lb, rev(preds.sglt$ci.ub)),
                col = alpha("darkslategray3", alpha = 0.4),
                border = "transparent")
        polygon(x = c(newx.glp, rev(newx.glp)), 
                y = c(preds.glp$ci.lb, rev(preds.glp$ci.ub)),
                col = alpha("indianred3", alpha = 0.4),
                border = "transparent")
        # 3. Point data (bubbles) with size = f(1/var)
        points(x =   l.do.sens[[i]]$cvdepy, 
               y =   -100*l.do.sens[[i]]$ard, # 100 to make it percentaje points
               col = l.do.sens[[i]]$border,
               bg =  l.do.sens[[i]]$colnum2,
               pch = 21, 
               cex = ((0.25)*l.do.sens[[i]]$wsize2)^(1.5))
        # Model for GLP1-RA
        segments(x0 = min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1"], na.rm = T), 
                 x1 = max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1"], na.rm = T),
                 y0 = modsens.glp1[[i]]$beta[1] + 
                   min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1"], na.rm = T)*modsens.glp1[[i]]$beta[2], 
                 y1 = modsens.glp1[[i]]$beta[1] + 
                   max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1"], na.rm = T)*modsens.glp1[[i]]$beta[2],
                 col =  "indianred3",
                 lwd = 2)
        # 4. Model for SGLT2i
        segments(x0 = min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T), 
                 x1 = max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                 y0 = modsens.sglt2[[i]]$beta[1] + 
                   min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T)*modsens.sglt2[[i]]$beta[2], 
                 y1 = modsens.sglt2[[i]]$beta[1] + 
                   max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T)*modsens.sglt2[[i]]$beta[2],
                 col =  "darkslategray4",
                 lwd = 2)
        # 4. Add id to identify trial
        text(x = jitter(l.do.sens[[i]]$cvdepy[!is.na(l.do.sens[[i]]$ard)],20), 
             y = jitter(-100*l.do.sens[[i]]$ard[!is.na(l.do.sens[[i]]$ard)],20),  
             l.do.sens[[i]]$id[!is.na(l.do.sens[[i]]$ard)], cex = 0.65 )
        # 5. title label per outcome
        title(main = v.titles[i], cex.main = 0.9, line = 0.7)
        # Axes titles
        if(i  %% 2 != 0){
          title(ylab = "Absolute risk difference (%)", cex.lab = 0.8, line = 2)
        } else {
          title(ylab = " Absolute risk difference (%)",  cex.lab = 0.8, line = 1.5)
        }
        if(i  >2){
          title(xlab = "Cardiovascular mortality rate (/100 py) in control group", cex.lab = 0.8, line = 2.5)
        } else {
          title(xlab = " Cardiovascular mortality rate (/100 py) in control group ", cex.lab = 0.8, line = 2.5)
        }
        # Axes
        axis(1, 
             at = c(-1,seq(0,14,2)),
             labels = c("",seq(0,12,2), ""),
             cex.axis = 0.75,
             col.axis = "gray40", 
             col.ticks = "gray40", 
             col = "gray40", 
             las = 1)
        axis(2, 
             #at = pretty(range(100*meta$ard, na.rm = T), n = 5),
             at  = seq(4,-16,-2),
             labels = c(seq(4,-14,-2),""),
             cex.axis = 0.75, 
             col.axis = "gray40", 
             col.ticks = "gray40", 
             col = "gray40", 
             las = 1)  
        if(i == 4){
          legend("bottomright", c("GLP-1RA", "SGLT2i"), 
                 cex = 0.85, 
                 pch = 21, 
                 bg = "transparent",
                 pt.bg = c("transparent", "darkslategray3" ),
                 col = c("indianred3", "darkslategray3" ), 
                 box.col = "transparent")
        }else{
          print("hi")
        }
      }
      
      dev.off()  
      
      
# S.2  Sensitivity Analysis CVD: Remove ELIXA and SOLOIST --------
  df.sens  <- meta[meta$outcomen == "CVMort" & 
                         !(meta$trialname == "ELIXA" |
                             meta$trialname == "SOLOIST"),]
  
  msens.g <- rma(-100*ard, (100^2)*ardvi.2, 
                  mods = ~ cvdepy, 
                  data = df.sens[df.sens$type == "GLP1",], 
                  method="REML", 
                  slab = trialname)
  msens.s <- rma(1-00*ard, (100^2)*ardvi.2, 
                  mods = ~ cvdepy, 
                  data = df.sens[df.sens$type == "SGLT2",], 
                  method="REML", 
                  slab = trialname)
  
  l.do.sens <- list()
  modsens.glp1  <- list() 
  modsens.sglt2 <- list()
  v.type <- c("GLP1","SGLT2")
  
  for(o in 1:length(v.outcomes)){ #loop runs over each outcome 'o' and each drug class 't'
    
    for(t in 1:length(v.type)){
      
      df.sens <- subset(meta, subset = (outcomelab == v.outcomes[o] & type == v.type[t] & 
                                          !(trialname == "ELIXA"| trialname == "SOLOIST-WHF"))
                        )

      l.do.sens[[o]] <- subset(meta, subset = (outcomelab == v.outcomes[o] & 
                                                 !(trialname == "ELIXA"|
                                                     trialname == "SOLOIST-WHF"  ))
                               )
                        
    
      # Model (random effects) for outcome o class t
      r <- rma(-100*ard, (100^2)*ardvi.2, # times 100 to re-scale to percentage points
               mods = ~ cvdepy, 
               data = df.sens, 
               method="REML", 
               slab = trialname)
      if( t == 1){ #Store in corresponding type-list (GLP1-RA/SGLT2i)
        modsens.glp1[[o]]  <- r
      } else {
        modsens.sglt2[[o]] <- r
      }
    }
  }
  
  # summary to recover quantities of interest (coefficient and pvalue)
  s.s <- lapply(modsens.sglt2, summary)
  s.g <- lapply(modsens.glp1, summary)
  #primary outcomes
  
  # empty data frame to store
  cc.s2 <- as.data.frame(cbind(rep(v.outcomes[1:7],2), 
                              c(rep("GLP-1RA",7), 
                                rep("SGLT2i",7)), 
                              rep(NA,14), 
                              rep(NA,14),
                              rep(NA,14),
                              rep(NA,14)))
  
  names(cc.s2)   <- c("Outcome", "Class", "Slope","CI.lb", "CI.ub", "P-value")

  
  # extract quantities of interest summary lists
  # Primary outcomes
  for(i in 1:7){
    cc.s2[i,3] <- round(s.g[[i]]$beta[2],2)
    cc.s2[i,4] <- round(s.g[[i]]$ci.lb[2],2)
    cc.s2[i,5] <- round(s.g[[i]]$ci.ub[2],2)
    cc.s2[i,6] <- round(s.g[[i]]$pval[2],2)
    cc.s2[i+7,3] <- round(s.s[[i]]$beta[2],2)
    cc.s2[i+7,4] <- round(s.s[[i]]$ci.lb[2],2)
    cc.s2[i+7,5] <- round(s.s[[i]]$ci.ub[2],2)
    cc.s2[i+7,6] <- round(s.s[[i]]$pval[2],2)
  }
  cc.s2 <- cc.s2[ccorder,]
  cc.s2
  # Export into text tables
  stargazer(cc.s2, 
            out = "_output/metareg_ard_s2.txt",
            summary = F,type = "text", 
            title = "Meta-regression coefficients, by drug class, S2 - accute",
            notes = "Absolute risk difference and baseline cardiovascular mortality rate. Sensitivity")
  
  
  png("metareg_ard_cdv_s2.png", width = 9, height = 9, units = 'in', res = 300) 
  
  par(oma = c(2,1,1,1), mfrow = c(2,2), mar = c(5,5,3,2)*0.75)
  for(i in 1:4){
    # Plot elements: 1. white canvas, 2. confint (sglt2 + glp1) + 
    #                3. bubbles 4. Model fit (line) 5. Trial number (text)
    #                6. titles (outcome)
    # 1. White canvas 
    
    plot(x = 0, y = 0, type ='n', 
         xlim = c(0,14), 
         ylim = c(-14,4), 
         xlab= "", 
         ylab ="", 
         axes = F)
    abline(h = 0, lty = 2, col = "gray60")
    text( x = 13.5, y = 2.5, labels = "Favors \n control", cex = 1.0, srt = 90)
    text( x = 13.5, y = -2.5, labels = "Favors \n treatment", cex = 1.0, srt = 90)
    
    # 2. Confidence interval: predict confint for both glp-1ras and sglt2i models
    # create domain for independent variable
    newx.glp <- seq(min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1" & 
                                                !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                    max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1" & 
                                                !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                    length.out = 1000)
    # predict glp1 ard for given x1 ( ard =f(x1) )
    preds.glp <- predict(modsens.glp1[[i]], 
                         newmods = newx.glp)
    # create domain for independent variable
    newx.sglt <- seq(min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2" & 
                                                 !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                     max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2" & 
                                                 !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                     length.out = 1000)
    # predict sglt2 ard for given x2 ( ard =f(x2) )
    preds.sglt <- predict(modsens.sglt2[[i]], 
                          newmods = newx.sglt)
    # Create CI polygons with computed coordinates
    polygon(x = c(newx.sglt, rev(newx.sglt)), 
            y = c(preds.sglt$ci.lb, rev(preds.sglt$ci.ub)),
            col = alpha("darkslategray3", alpha = 0.4),
            border = "transparent")
    polygon(x = c(newx.glp, rev(newx.glp)), 
            y = c(preds.glp$ci.lb, rev(preds.glp$ci.ub)),
            col = alpha("indianred3", alpha = 0.4),
            border = "transparent")
    # 3. Point data (bubbles) with size = f(1/var)
    points(x =   l.do.sens[[i]]$cvdepy, 
           y =   -100*l.do.sens[[i]]$ard, # 100 to make it percentaje points
           col = l.do.sens[[i]]$border,
           bg =  l.do.sens[[i]]$colnum2,
           pch = 21, 
           cex = ((0.25)*l.do.sens[[i]]$wsize2)^(1.5))
    # Model for GLP1-RA
    segments(x0 = min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1"], na.rm = T), 
             x1 = max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1"], na.rm = T),
             y0 = modsens.glp1[[i]]$beta[1] + 
               min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1"], na.rm = T)*modsens.glp1[[i]]$beta[2], 
             y1 = modsens.glp1[[i]]$beta[1] + 
               max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1"], na.rm = T)*modsens.glp1[[i]]$beta[2],
             col =  "indianred3",
             lwd = 2)
    # 4. Model for SGLT2i
    segments(x0 = min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T), 
             x1 = max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T),
             y0 = modsens.sglt2[[i]]$beta[1] + 
               min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T)*modsens.sglt2[[i]]$beta[2], 
             y1 = modsens.sglt2[[i]]$beta[1] + 
               max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T)*modsens.sglt2[[i]]$beta[2],
             col =  "darkslategray4",
             lwd = 2)
    # 4. Add id to identify trial
    text(x = jitter(l.do.sens[[i]]$cvdepy[!is.na(l.do.sens[[i]]$ard)],20), 
         y = jitter(-100*l.do.sens[[i]]$ard[!is.na(l.do.sens[[i]]$ard)],20),  
         l.do.sens[[i]]$id[!is.na(l.do.sens[[i]]$ard)], cex = 0.65 )
    # 5. title label per outcome
    title(main = v.titles[i], cex.main = 0.9, line = 0.7)
    # Axes titles
    if(i  %% 2 != 0){
      title(ylab = "Absolute risk difference (%)", cex.lab = 0.8, line = 2)
    } else {
      title(ylab = " Absolute risk difference (%)",  cex.lab = 0.8, line = 1.5)
    }
    if(i  >2){
      title(xlab = "Cardiovascular mortality rate (/100 py) in control group", cex.lab = 0.8, line = 2.5)
    } else {
      title(xlab = " Cardiovascular mortality rate (/100 py) in control group ", cex.lab = 0.8, line = 2.5)
    }
    # Axes
    axis(1, 
         at = c(-1,seq(0,14,2)),
         labels = c("",seq(0,12,2), ""),
         cex.axis = 0.75,
         col.axis = "gray40", 
         col.ticks = "gray40", 
         col = "gray40", 
         las = 1)
    axis(2, 
         #at = pretty(range(100*meta$ard, na.rm = T), n = 5),
         at  = seq(4,-16,-2),
         labels = c(seq(4,-14,-2),""),
         cex.axis = 0.75, 
         col.axis = "gray40", 
         col.ticks = "gray40", 
         col = "gray40", 
         las = 1)  
    if(i == 2){
      legend("bottomright", c("GLP-1RA", "SGLT2i"), 
             cex = 0.85, 
             pch = 21, 
             bg = "transparent",
             pt.bg = c("transparent", "darkslategray3" ),
             col = c("indianred3", "darkslategray3" ), 
             box.col = "transparent")
    }else{
      print("hi")
    }
  }
  
  dev.off()  
  
  
# S.3  Sensitivity Analysis Renal: Remove ELIXA and CREDENCE --------
    l.do.sens <- list()
    modsens.glp1  <- list() 
    modsens.sglt2 <- list()
    v.type <- c("GLP1","SGLT2")
    
   
    for(o in 1:length(v.outcomes)){ #loop runs over each outcome 'o' and each drug class 't'
    for(t in 1:length(v.type)){
        
        df.sens <- subset(meta, subset = (outcomelab == "Composite Renal Outcome" & 
                                            type == v.type[t] & 
                                            !(trialname == "ELIXA"|
                                              trialname == "CREDENCE"|
                                              trialname == "DAPA-CKD"|
                                              trialname == "DAPA-HF"|
                                              trialname == "SCORED" |
                                              trialname == "SOLOIST-WHF"  ))
                          )
        l.do.sens[[1]] <- subset(meta, subset = (outcomelab == "Composite Renal Outcome" & 
                                            !(trialname == "ELIXA"|
                                                trialname == "CREDENCE"|
                                                trialname == "DAPA-CKD"|
                                                trialname == "DAPA-HF"|
                                                trialname == "SCORED" |
                                                trialname == "SOLOIST-WHF"  ))
        )
        # Model (random effects) for outcome o class t
        r <- rma(-100*ard, (100^2)*ardvi.2, # times 100 to re-scale to percentage points
                 mods = ~ cvdepy, 
                 data = df.sens, 
                 method="REML", 
                 slab = trialname)
        if( t == 1){ #Store in corresponding type-list (GLP1-RA/SGLT2i)
          modsens.glp1[[1]]  <- r
        } else {
          modsens.sglt2[[1]] <- r
        }
      }
    }
    
    # summary to recover quantities of interest (coefficient and pvalue)
    s.s <- lapply(modsens.sglt2, summary)
    s.g <- lapply(modsens.glp1, summary)
    
    
    # empty data frame to store
    cc.s3 <- as.data.frame(cbind(rep(v.outcomes[4],1), 
                                c(rep("GLP-1RA",1), 
                                  rep("SGLT2i",1)), 
                                rep(NA,1), 
                                rep(NA,1),
                                rep(NA,1),
                                rep(NA,1)))
    
    names(cc.s3)   <- c("Outcome", "Class", "Slope","CI.lb", "CI.ub", "P-value")
    
    # extract quantities of interest summary lists
    # Primary outcome
    for(i in 1){
      cc.s3[i,3] <- round(s.g[[i]]$beta[2],2)
      cc.s3[i,4] <- round(s.g[[i]]$ci.lb[2],2)
      cc.s3[i,5] <- round(s.g[[i]]$ci.ub[2],2)
      cc.s3[i,6] <- round(s.g[[i]]$pval[2],2)
      cc.s3[i+1,3] <- round(s.s[[i]]$beta[2],2)
      cc.s3[i+1,4] <- round(s.s[[i]]$ci.lb[2],2)
      cc.s3[i+1,5] <- round(s.s[[i]]$ci.ub[2],2)
      cc.s3[i+1,6] <- round(s.s[[i]]$pval[2],2)
    }
    
    cc.s3 
    # Export into text tables
    stargazer(cc.s3, 
              out = "_output/metareg_ard_s3.txt",
              summary = F,type = "text", 
              title = "Meta-regression coefficients, by drug class. S3 -renal outcome variable",
              notes = "Absolute risk difference and baseline cardiovascular mortality rate. Sensitivity")
    
    

    png("metareg_ard_cdv_s3.png", width = 9, height = 9, units = 'in', res = 300) 
    par(oma = c(3,2,0,0), mfrow = c(1,1), mar = c(5,4,3,0)*0.75)
    for(i in 1){
      # Plot elements: 1. white canvas, 2. confint (sglt2 + glp1) + 
      #                3. bubbles 4. Model fit (line) 5. Trial number (text)
      #                6. titles (outcome)
      # 1. White canvas 
      
      plot(x = 0, y = 0, type ='n', 
           xlim = c(0,14), 
           ylim = c(-14,4), 
           xlab= "", 
           ylab ="", 
           axes = F)
      abline(h = 0, lty = 2, col = "gray60")
      text( x = 13.5, y = 2.5, labels = "Favors \n control", cex = 1.0, srt = 90)
      text( x = 13.5, y = -2.5, labels = "Favors \n treatment", cex = 1.0, srt = 90)
      
      # 2. Confidence interval: predict confint for both glp-1ras and sglt2i models
      # create domain for independent variable
      newx.glp <- seq(min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1" & 
                                             !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                      max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1" & 
                                             !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                      length.out = 1000)
      # predict glp1 ard for given x1 ( ard =f(x1) )
      preds.glp <- predict(modsens.glp1[[i]], 
                           newmods = newx.glp)
      # create domain for independent variable
      newx.sglt <- seq(min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2" & 
                                              !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                       max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2" & 
                                              !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                       length.out = 1000)
      # predict sglt2 ard for given x2 ( ard =f(x2) )
      preds.sglt <- predict(modsens.sglt2[[i]], 
                            newmods = newx.sglt)
      # Create CI polygons with computed coordinates
      polygon(x = c(newx.sglt, rev(newx.sglt)), 
              y = c(preds.sglt$ci.lb, rev(preds.sglt$ci.ub)),
              col = alpha("darkslategray3", alpha = 0.4),
              border = "transparent")
      polygon(x = c(newx.glp, rev(newx.glp)), 
              y = c(preds.glp$ci.lb, rev(preds.glp$ci.ub)),
              col = alpha("indianred3", alpha = 0.4),
              border = "transparent")
      # 3. Point data (bubbles) with size = f(1/var)
      points(x =   l.do.sens[[i]]$cvdepy, 
             y =   -100*l.do.sens[[i]]$ard, # 100 to make it percentaje points
             col = l.do.sens[[i]]$border,
             bg =  l.do.sens[[i]]$colnum2,
             pch = 21, 
             cex = ((0.25)*l.do.sens[[i]]$wsize2)^(1.5))
      # Model for GLP1-RA
      segments(x0 = min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1"], na.rm = T), 
               x1 = max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1"], na.rm = T),
               y0 = modsens.glp1[[i]]$beta[1] + 
                 min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1"], na.rm = T)*modsens.glp1[[i]]$beta[2], 
               y1 = modsens.glp1[[i]]$beta[1] + 
                 max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1"], na.rm = T)*modsens.glp1[[i]]$beta[2],
               col =  "indianred3",
               lwd = 2)
      # 4. Model for SGLT2i
      segments(x0 = min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T), 
               x1 = max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T),
               y0 = modsens.sglt2[[i]]$beta[1] + 
                 min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T)*modsens.sglt2[[i]]$beta[2], 
               y1 = modsens.sglt2[[i]]$beta[1] + 
                 max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T)*modsens.sglt2[[i]]$beta[2],
               col =  "darkslategray4",
               lwd = 2)
      # 4. Add id to identify trial
      text(x = jitter(l.do.sens[[i]]$cvdepy[!is.na(l.do.sens[[i]]$ard)],20), 
           y = jitter(-100*l.do.sens[[i]]$ard[!is.na(l.do.sens[[i]]$ard)],20),  
           l.do.sens[[i]]$id[!is.na(l.do.sens[[i]]$ard)], cex = 0.65 )
      # 5. title label per outcome
      title(main = "Composite renal outcome", cex.main = 0.9, line = 0.7)
      # Axes titles
      if(i  %% 2 != 0){
        title(ylab = "Absolute risk difference (%)", cex.lab = 0.8, line = 2)
      } else {
        title(ylab = " Absolute risk difference (%)",  cex.lab = 0.8, line = 1.5)
      }
      if(i  >2){
        title(xlab = "Cardiovascular mortality rate (/100 py) in control group", cex.lab = 0.8, line = 2.5)
      } else {
        title(xlab = " Cardiovascular mortality rate (/100 py) in control group ", cex.lab = 0.8, line = 2.5)
      }
      # Axes
      axis(1, 
           at = c(-1,seq(0,14,2)),
           labels = c("",seq(0,12,2), ""),
           cex.axis = 0.75,
           col.axis = "gray40", 
           col.ticks = "gray40", 
           col = "gray40", 
           las = 1)
      axis(2, 
           #at = pretty(range(100*meta$ard, na.rm = T), n = 5),
           at  = seq(4,-16,-2),
           labels = c(seq(4,-14,-2),""),
           cex.axis = 0.75, 
           col.axis = "gray40", 
           col.ticks = "gray40", 
           col = "gray40", 
           las = 1)  
      if(i == 2){
        legend("bottomright", c("GLP-1RA", "SGLT2i"), 
               cex = 0.85, 
               pch = 21, 
               bg = "transparent",
               pt.bg = c("transparent", "darkslategray3" ),
               col = c("indianred3", "darkslategray3" ), 
               box.col = "transparent")
      }else{
        print("hi")
      }
    }
  
dev.off()    
  

# S4. Sensitivity Analysis: No approval =========
modsens.glp1  <- list() 
modsens.sglt2 <- list()
v.type <- c("GLP1","SGLT2")

for(o in 1:length(v.outcomes)){ #loop runs over each outcome 'o' and each drug class 't'
  
  for(t in 1:length(v.type)){
    
    df.sens <- subset(meta, subset = (outcomelab == v.outcomes[o] & type == v.type[t] & 
                                        !(trialname == "AMPLITUDE-O"| trialname == "FREEDOM-CVO"|
                                            trialname == "Harmony Outcomes" |trialname == "SCORED"|
                                            trialname == "SOLOIST-WHF"))
    )
    l.do.sens[[o]] <- subset(meta, subset = (outcomelab == v.outcomes[o] & 
                                               !(trialname == "AMPLITUDE-O"| trialname == "FREEDOM-CVO"|
                                                   trialname == "Harmony Outcomes" |trialname == "SCORED"|
                                                   trialname == "SOLOIST-WHF"))
    )
    
    # Model (random effects) for outcome o class t
    r <- rma(-100*ard, (100^2)*ardvi.2, # times 100 to re-scale to percentage points
             mods = ~ cvdepy, 
             data = df.sens, 
             method="REML", 
             slab = trialname)
    if( t == 1){ #Store in corresponding type-list (GLP1-RA/SGLT2i)
      modsens.glp1[[o]]  <- r
    } else {
      modsens.sglt2[[o]] <- r
    }
  }
}

# summary to recover quantities of interest (coefficient and pvalue)
s.s <- lapply(modsens.sglt2, summary)
s.g <- lapply(modsens.glp1, summary)
#primary outcomes

# empty data frame to store
cc.s4 <- as.data.frame(cbind(rep(v.outcomes[1:7],2), 
                            c(rep("GLP-1RA",7), 
                              rep("SGLT2i",7)), 
                            rep(NA,14), 
                            rep(NA,14),
                            rep(NA,14),
                            rep(NA,14)))

names(cc.s4)   <- c("Outcome", "Class", "Slope","CI.lb", "CI.ub", "P-value")

# extract quantities of interest summary lists
# Primary outcomes
for(i in 1:7){
  cc.s4[i,3] <- round(s.g[[i]]$beta[2],2)
  cc.s4[i,4] <- round(s.g[[i]]$ci.lb[2],2)
  cc.s4[i,5] <- round(s.g[[i]]$ci.ub[2],2)
  cc.s4[i,6] <- round(s.g[[i]]$pval[2],2)
  cc.s4[i+7,3] <- round(s.s[[i]]$beta[2],2)
  cc.s4[i+7,4] <- round(s.s[[i]]$ci.lb[2],2)
  cc.s4[i+7,5] <- round(s.s[[i]]$ci.ub[2],2)
  cc.s4[i+7,6] <- round(s.s[[i]]$pval[2],2)
}
cc.s4 <- cc.s4[ccorder,]
cc.s4
# Export into text tables
stargazer(cc.s4, 
          out = "_output/metareg_ard_s4.txt",
          summary = F,type = "text", 
          title = "Meta-regression coefficients, by drug class, S4 - FDA approval/ not in market",
          notes = "Absolute risk difference and baseline cardiovascular mortality rate. Sensitivity")

png("metareg_ard_cdv_s4.png", width = 9, height = 9, units = 'in', res = 300) 
par(oma = c(2,1,1,1), mfrow = c(2,2), mar = c(5,5,3,2)*0.75)
for(i in 1:4){
  # Plot elements: 1. white canvas, 2. confint (sglt2 + glp1) + 
  #                3. bubbles 4. Model fit (line) 5. Trial number (text)
  #                6. titles (outcome)
  # 1. White canvas 
  
  plot(x = 0, y = 0, type ='n', 
       xlim = c(-0,14), 
       ylim = c(-14,4), 
       xlab= "", 
       ylab ="", 
       axes = F)
  abline(h = 0, lty = 2, col = "gray60")
  text( x = 13.5, y = 2.5, labels = "Favors \n control", cex = 1.0, srt = 90)
  text( x = 13.5, y = -2.5, labels = "Favors \n treatment", cex = 1.0, srt = 90)
  
  # 2. Confidence interval: predict confint for both glp-1ras and sglt2i models
  # create domain for independent variable
  newx.glp <- seq(min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1" & 
                                              !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                  max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1" & 
                                              !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                  length.out = 1000)
  # predict glp1 ard for given x1 ( ard =f(x1) )
  preds.glp <- predict(modsens.glp1[[i]], 
                       newmods = newx.glp)
  # create domain for independent variable
  newx.sglt <- seq(min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2" & 
                                               !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                   max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2" & 
                                               !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                   length.out = 1000)
  # predict sglt2 ard for given x2 ( ard =f(x2) )
  preds.sglt <- predict(modsens.sglt2[[i]], 
                        newmods = newx.sglt)
  # Create CI polygons with computed coordinates
  polygon(x = c(newx.sglt, rev(newx.sglt)), 
          y = c(preds.sglt$ci.lb, rev(preds.sglt$ci.ub)),
          col = alpha("darkslategray3", alpha = 0.4),
          border = "transparent")
  polygon(x = c(newx.glp, rev(newx.glp)), 
          y = c(preds.glp$ci.lb, rev(preds.glp$ci.ub)),
          col = alpha("indianred3", alpha = 0.4),
          border = "transparent")
  # 3. Point data (bubbles) with size = f(1/var)
  points(x =   l.do.sens[[i]]$cvdepy, 
         y =   -100*l.do.sens[[i]]$ard, # 100 to make it percentaje points
         col = l.do.sens[[i]]$border,
         bg =  l.do.sens[[i]]$colnum2,
         pch = 21, 
         cex = ((0.25)*l.do.sens[[i]]$wsize2)^(1.5))
  # Model for GLP1-RA
  segments(x0 = min(newx.glp), 
           x1 = max(newx.glp),
           y0 = modsens.glp1[[i]]$beta[1] + 
             min(newx.glp)*modsens.glp1[[i]]$beta[2], 
           y1 = modsens.glp1[[i]]$beta[1] + 
             max(newx.glp)*modsens.glp1[[i]]$beta[2],
           col =  "indianred3",
           lwd = 2)
  # 4. Model for SGLT2i
  segments(x0 = min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T), 
           x1 = max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T),
           y0 = modsens.sglt2[[i]]$beta[1] + 
             min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T)*modsens.sglt2[[i]]$beta[2], 
           y1 = modsens.sglt2[[i]]$beta[1] + 
             max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2"& !is.na(l.do.sens[[i]]$ard)], na.rm = T)*modsens.sglt2[[i]]$beta[2],
           col =  "darkslategray4",
           lwd = 2)
  # 4. Add id to identify trial
  text(x = jitter(l.do.sens[[i]]$cvdepy[!is.na(l.do.sens[[i]]$ard)],20), 
       y = jitter(-100*l.do.sens[[i]]$ard[!is.na(l.do.sens[[i]]$ard)],20),  
       l.do.sens[[i]]$id[!is.na(l.do.sens[[i]]$ard)], cex = 0.65 )
  # 5. title label per outcome
  title(main = v.titles[i], cex.main = 0.9, line = 0.7)
  # Axes titles
  if(i  %% 2 != 0){
    title(ylab = "Absolute risk difference (%)", cex.lab = 0.8, line = 2)
  } else {
    title(ylab = " Absolute risk difference (%)",  cex.lab = 0.8, line = 1.5)
  }
  if(i  >2){
    title(xlab = "Cardiovascular mortality rate (/100 py) in control group", cex.lab = 1, line = 2.5)
  } else {
    title(xlab = " Cardiovascular mortality rate (/100 py) in control group ", cex.lab = 1, line = 2.5)
  }
  # Axes
  axis(1, 
       at = c(-1,seq(0,14,2)),
       labels = c("",seq(0,12,2), ""),
       cex.axis = 0.75,
       col.axis = "gray40", 
       col.ticks = "gray40", 
       col = "gray40", 
       las = 1)
  axis(2, 
       #at = pretty(range(100*meta$ard, na.rm = T), n = 5),
       at  = seq(4,-16,-2),
       labels = c(seq(4,-14,-2),""),
       cex.axis = 0.75, 
       col.axis = "gray40", 
       col.ticks = "gray40", 
       col = "gray40", 
       las = 1)  
  if(i == 2){
    legend("bottomright", c("GLP-1RA", "SGLT2i"), 
           cex = 0.85, 
           pch = 21, 
           bg = "transparent",
           pt.bg = c("transparent", "darkslategray3" ),
           col = c("indianred3", "darkslategray3" ), 
           box.col = "transparent")
  }else{
    print("hi")
  }
}

dev.off()  


# S5. Hist HF =======
modsens.glp1  <- list() 
modsens.sglt2 <- list()
v.type <- c("GLP1","SGLT2")

for(o in 1:length(v.outcomes)){ #loop runs over each outcome 'o' and each drug class 't'
  
  for(t in 1:length(v.type)){
    
    df.sens <- subset(meta, subset = (outcomelab == v.outcomes[o] & type == v.type[t] & 
                                        !(trialname == "EMPEROR-Reduced"| trialname == "EMPEROR-Preserved"|
                                            trialname == "DAPA-HF"|
                                            trialname == "DELIVER"|
                                            trialname == "SOLOIST-WHF"))
    )
    l.do.sens[[o]] <- subset(meta, subset = (outcomelab == v.outcomes[o] & 
                                               !(trialname == "EMPEROR-Reduced"| trialname == "EMPEROR-Preserved"|
                                                   trialname == "DAPA-HF"|
                                                   trialname == "DELIVER"|
                                                   trialname == "SOLOIST-WHF"))
    )
    
    # Model (random effects) for outcome o class t
    r <- rma(-100*ard, (100^2)*ardvi.2, # times 100 to re-scale to percentage points
             mods = ~ cvdepy, 
             data = df.sens, 
             method="REML", 
             slab = trialname)
    if( t == 1){ #Store in corresponding type-list (GLP1-RA/SGLT2i)
      modsens.glp1[[o]]  <- r
    } else {
      modsens.sglt2[[o]] <- r
    }
  }
}

# summary to recover quantities of interest (coefficient and pvalue)
s.s <- lapply(modsens.sglt2, summary)
s.g <- lapply(modsens.glp1, summary)
#primary outcomes

# empty data frame to store
cc.s5 <- as.data.frame(cbind(rep(v.outcomes[1:7],2), 
                            c(rep("GLP-1RA",7), 
                              rep("SGLT2i",7)), 
                            rep(NA,14), 
                            rep(NA,14),
                            rep(NA,14),
                            rep(NA,14)))

names(cc.s5)   <- c("Outcome", "Class", "Slope","CI.lb", "CI.ub", "P-value")

# extract quantities of interest summary lists
# Primary outcomes
for(i in 1:7){
  cc.s5[i,3] <- round(s.g[[i]]$beta[2],2)
  cc.s5[i,4] <- round(s.g[[i]]$ci.lb[2],2)
  cc.s5[i,5] <- round(s.g[[i]]$ci.ub[2],2)
  cc.s5[i,6] <- round(s.g[[i]]$pval[2],2)
  cc.s5[i+7,3] <- round(s.s[[i]]$beta[2],2)
  cc.s5[i+7,4] <- round(s.s[[i]]$ci.lb[2],2)
  cc.s5[i+7,5] <- round(s.s[[i]]$ci.ub[2],2)
  cc.s5[i+7,6] <- round(s.s[[i]]$pval[2],2)
}
cc.s5 <- cc.s5[ccorder,]
cc.s5
# Export into text tables
stargazer(cc.s5, 
          out = "_output/metareg_ard_s5.txt",
          summary = F,type = "text", 
          title = "Meta-regression coefficients, by drug class, S5 - History of HF",
          notes = "Absolute risk difference and baseline cardiovascular mortality rate. Sensitivity")

png("metareg_ard_cdv_s5.png", width = 9, height = 9, units = 'in', res = 300) 
par(oma = c(2,1,1,1), mfrow = c(2,2), mar = c(5,5,3,2)*0.75)
for(i in 1:4){
  # Plot elements: 1. white canvas, 2. confint (sglt2 + glp1) + 
  #                3. bubbles 4. Model fit (line) 5. Trial number (text)
  #                6. titles (outcome)
  # 1. White canvas 
  
  plot(x = 0, y = 0, type ='n', 
       xlim = c(-0.2,3), 
       ylim = c(-12,4), 
       xlab= "", 
       ylab ="", 
       axes = F)
  abline(h = 0, lty = 2, col = "gray60")
  text( x = 3.0, y = 2.5, labels = "Favors \n control", cex = 1.0, srt = 90)
  text( x = 3.0, y = -2.5, labels = "Favors \n treatment", cex = 1.0, srt = 90)
  
  # 2. Confidence interval: predict confint for both glp-1ras and sglt2i models
  # create domain for independent variable
  newx.glp <- seq(min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1" & 
                                              !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                  max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1" & 
                                              !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                  length.out = 1000)
  # predict glp1 ard for given x1 ( ard =f(x1) )
  preds.glp <- predict(modsens.glp1[[i]], 
                       newmods = newx.glp)
  # create domain for independent variable
  newx.sglt <- seq(min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2" & 
                                               !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                   max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2" & 
                                               !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                   length.out = 1000)
  # predict sglt2 ard for given x2 ( ard =f(x2) )
  preds.sglt <- predict(modsens.sglt2[[i]], 
                        newmods = newx.sglt)
  # Create CI polygons with computed coordinates
  polygon(x = c(newx.sglt, rev(newx.sglt)), 
          y = c(preds.sglt$ci.lb, rev(preds.sglt$ci.ub)),
          col = alpha("darkslategray3", alpha = 0.6),
          border = "transparent")
  polygon(x = c(newx.glp, rev(newx.glp)), 
          y = c(preds.glp$ci.lb, rev(preds.glp$ci.ub)),
          col = alpha("indianred3", alpha = 0.4),
          border = "transparent")
  # 3. Point data (bubbles) with size = f(1/var)
  points(x =   l.do.sens[[i]]$cvdepy, 
         y =   -100*l.do.sens[[i]]$ard, # 100 to make it percentaje points
         col = l.do.sens[[i]]$border,
         bg =  l.do.sens[[i]]$colnum2,
         pch = 21, 
         cex = ((0.25)*l.do.sens[[i]]$wsize2)^(1.5))
  # Model for GLP1-RA
  segments(x0 = min(newx.glp), 
           x1 = max(newx.glp),
           y0 = modsens.glp1[[i]]$beta[1] + 
             min(newx.glp)*modsens.glp1[[i]]$beta[2], 
           y1 = modsens.glp1[[i]]$beta[1] + 
             max(newx.glp)*modsens.glp1[[i]]$beta[2],
           col =  "indianred3",
           lwd = 2)
  # 4. Model for SGLT2i
  segments(x0 = min(newx.sglt), 
           x1 = max(newx.sglt),
           y0 = modsens.sglt2[[i]]$beta[1] + 
             min(newx.sglt)*modsens.sglt2[[i]]$beta[2], 
           y1 = modsens.sglt2[[i]]$beta[1] + 
             max(newx.sglt)*modsens.sglt2[[i]]$beta[2],
           col =  "darkslategray4",
           lwd = 2)
  # 4. Add id to identify trial
  text(x = jitter(l.do.sens[[i]]$cvdepy[!is.na(l.do.sens[[i]]$ard)],20), 
       y = jitter(-100*l.do.sens[[i]]$ard[!is.na(l.do.sens[[i]]$ard)],20),  
       l.do.sens[[i]]$id[!is.na(l.do.sens[[i]]$ard)], cex = 0.65 )
  # 5. title label per outcome
  title(main = v.titles[i], cex.main = 0.9, line = 0.7)
  # Axes titles
  if(i  %% 2 != 0){
    title(ylab = "Absolute risk difference (%)", cex.lab = 1, line = 2.5)
  } else {
    title(ylab = " Absolute risk difference (%)",  cex.lab = 1, line = 2.5)
  }
  if(i  >2){
    title(xlab = "Cardiovascular mortality rate (/100 py) in control group", cex.lab = 1, line = 2.5)
  } else {
    title(xlab = " Cardiovascular mortality rate (/100 py) in control group ", cex.lab = 1, line = 2.5)
  }
  # Axes
  axis(1, 
       at = c(-0.2,seq(0,4,0.5)),
       labels = c("",seq(0,2.5,0.5), "","",""),
       cex.axis = 0.75,
       col.axis = "gray40", 
       col.ticks = "gray40", 
       col = "gray40", 
       las = 1)
  axis(2, 
       #at = pretty(range(100*meta$ard, na.rm = T), n = 5),
       at  = seq(4,-14,-2),
       labels = c(seq(4,-12,-2),""),
       cex.axis = 0.75, 
       col.axis = "gray40", 
       col.ticks = "gray40", 
       col = "gray40", 
       las = 1)  
  if(i == 4){
    legend("bottomright", c("GLP-1RA", "SGLT2i"), 
           cex = 0.85, 
           pch = 21, 
           bg = "transparent",
           pt.bg = c("transparent", "darkslategray3" ),
           col = c("indianred3", "darkslategray3" ), 
           box.col = "transparent")
  }else{
    print("hi")
  }
}


dev.off() 


# S.6 MACE w/o Declare ========

modsens.glp1  <- list() 
modsens.sglt2 <- list()
v.type <- c("GLP1","SGLT2")


for(o in 1){ #loop runs over each outcome 'o' and each drug class 't'
  
  for(t in 1:length(v.type)){
    
    df.sens <- subset(meta, subset = (outcomelab == "MACE" & type == v.type[t] & 
                                        !(trialname == "DECLARE-TIMI 58"))
    )
    l.do.sens[[o]] <- subset(meta, subset = (outcomelab == "MACE" & 
                                               !(trialname == "DECLARE-TIMI 58"))
    )
    
    # Model (random effects) for outcome o class t
    r <- rma(-100*ard, (100^2)*ardvi.2, # times 100 to re-scale to percentage points
             mods = ~ cvdepy, 
             data = df.sens, 
             method="REML", 
             slab = trialname)
    if( t == 1){ #Store in corresponding type-list (GLP1-RA/SGLT2i)
      modsens.glp1[[o]]  <- r
    } else {
      modsens.sglt2[[o]] <- r
    }
  }
}

# summary to recover quantities of interest (coefficient and pvalue)
s.s <- lapply(modsens.sglt2, summary)
s.g <- lapply(modsens.glp1, summary)

# empty data frame to store
cc.s6 <- as.data.frame(cbind(rep("MACE",2), 
                             c(rep("GLP-1RA",1), 
                               rep("SGLT2i",1)), 
                             rep(NA,2), 
                             rep(NA,2),
                             rep(NA,2),
                             rep(NA,2)))

names(cc.s6)   <- c("Outcome", "Class", "Slope","CI.lb", "CI.ub", "P-value")

# extract quantities of interest summary lists
# Primary outcomes
for(i in 1){
  cc.s6[i,3] <- round(s.g[[i]]$beta[2],2)
  cc.s6[i,4] <- round(s.g[[i]]$ci.lb[2],2)
  cc.s6[i,5] <- round(s.g[[i]]$ci.ub[2],2)
  cc.s6[i,6] <- round(s.g[[i]]$pval[2],2)
  cc.s6[i+1,3] <- round(s.s[[i]]$beta[2],2)
  cc.s6[i+1,4] <- round(s.s[[i]]$ci.lb[2],2)
  cc.s6[i+1,5] <- round(s.s[[i]]$ci.ub[2],2)
  cc.s6[i+1,6] <- round(s.s[[i]]$pval[2],2)
}

# Export into text tables
stargazer(cc.s6, 
          out = "_output/metareg_ard_s6.txt",
          summary = F,type = "text", 
          title = "Meta-regression coefficients, by drug class, S6 - Inconsistent MACE",
          notes = "Absolute risk difference and baseline cardiovascular mortality rate. Sensitivity")

png("metareg_ard_cdv_s6.png", width = 9, height = 9, units = 'in', res = 300) 
par(oma = c(2,1,1,1), mfrow = c(1,1), mar = c(5,5,3,0)*0.75)
for(i in 1:1){
  # Plot elements: 1. white canvas, 2. confint (sglt2 + glp1) + 
  #                3. bubbles 4. Model fit (line) 5. Trial number (text)
  #                6. titles (outcome)
  # 1. White canvas 
  
  plot(x = 0, y = 0, type ='n', 
       xlim = c(-0.2,3), 
       ylim = c(-14,4), 
       xlab= "", 
       ylab ="", 
       axes = F)
  abline(h = 0, lty = 2, col = "gray60")
  text( x = 3.0, y = 2.5, labels = "Favors \n control", cex = 1.0, srt = 90)
  text( x = 3.0, y = -2.5, labels = "Favors \n treatment", cex = 1.0, srt = 90)
  
  # 2. Confidence interval: predict confint for both glp-1ras and sglt2i models
  # create domain for independent variable
  newx.glp <- seq(min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1" & 
                                              !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                  max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "GLP1" & 
                                              !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                  length.out = 1000)
  # predict glp1 ard for given x1 ( ard =f(x1) )
  preds.glp <- predict(modsens.glp1[[i]], 
                       newmods = newx.glp)
  # create domain for independent variable
  newx.sglt <- seq(min(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2" & 
                                               !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                   max(l.do.sens[[i]]$cvdepy[l.do.sens[[i]]$type == "SGLT2" & 
                                               !is.na(l.do.sens[[i]]$ard)], na.rm = T),
                   length.out = 1000)
  # predict sglt2 ard for given x2 ( ard =f(x2) )
  preds.sglt <- predict(modsens.sglt2[[i]], 
                        newmods = newx.sglt)
  # Create CI polygons with computed coordinates
  polygon(x = c(newx.sglt, rev(newx.sglt)), 
          y = c(preds.sglt$ci.lb, rev(preds.sglt$ci.ub)),
          col = alpha("darkslategray3", alpha = 0.6),
          border = "transparent")
  polygon(x = c(newx.glp, rev(newx.glp)), 
          y = c(preds.glp$ci.lb, rev(preds.glp$ci.ub)),
          col = alpha("indianred3", alpha = 0.4),
          border = "transparent")
  # 3. Point data (bubbles) with size = f(1/var)
  points(x =   l.do.sens[[i]]$cvdepy, 
         y =   -100*l.do.sens[[i]]$ard, # 100 to make it percentaje points
         col = l.do.sens[[i]]$border,
         bg =  l.do.sens[[i]]$colnum2,
         pch = 21, 
         cex = ((0.25)*l.do.sens[[i]]$wsize2)^(1.5))
  # Model for GLP1-RA
  segments(x0 = min(newx.glp), 
           x1 = max(newx.glp),
           y0 = modsens.glp1[[i]]$beta[1] + 
             min(newx.glp)*modsens.glp1[[i]]$beta[2], 
           y1 = modsens.glp1[[i]]$beta[1] + 
             max(newx.glp)*modsens.glp1[[i]]$beta[2],
           col =  "indianred3",
           lwd = 2)
  # 4. Model for SGLT2i
  segments(x0 = min(newx.sglt), 
           x1 = max(newx.sglt),
           y0 = modsens.sglt2[[i]]$beta[1] + 
             min(newx.sglt)*modsens.sglt2[[i]]$beta[2], 
           y1 = modsens.sglt2[[i]]$beta[1] + 
             max(newx.sglt)*modsens.sglt2[[i]]$beta[2],
           col =  "darkslategray4",
           lwd = 2)
  # 4. Add id to identify trial
  text(x = jitter(l.do.sens[[i]]$cvdepy[!is.na(l.do.sens[[i]]$ard)],20), 
       y = jitter(-100*l.do.sens[[i]]$ard[!is.na(l.do.sens[[i]]$ard)],20),  
       l.do.sens[[i]]$id[!is.na(l.do.sens[[i]]$ard)], cex = 0.65 )
  # 5. title label per outcome
  title(main = "MACE", cex.main = 0.9, line = 0.7)
  # Axes titles
  if(i  %% 2 != 0){
    title(ylab = "Absolute risk difference (%)", cex.lab = 0.8, line = 2)
  } else {
    title(ylab = " Absolute risk difference (%)",  cex.lab = 0.8, line = 1.5)
  }
  if(i  >2){
    title(xlab = "Cardiovascular mortality rate (/100 py) in control group", cex.lab = 0.8, line = 1.6)
  } else {
    title(xlab = " Cardiovascular mortality rate (/100 py) in control group ", cex.lab = 0.8, line = 1.6)
  }
  # Axes
  axis(1, 
       at = c(-0.3, seq(0,3,0.5)),
       labels = c("", seq(0,2.5,0.5),""),
       cex.axis = 0.75,
       col.axis = "gray40", 
       col.ticks = "gray40", 
       col = "gray40", 
       las = 1)
  axis(2, 
       #at = pretty(range(100*meta$ard, na.rm = T), n = 5),
       at  = seq(4,-16,-2),
       labels = c(seq(4,-14,-2),""),
       cex.axis = 0.75, 
       col.axis = "gray40", 
       col.ticks = "gray40", 
       col = "gray40", 
       las = 1)  
  if(i == 1){
    legend("bottomright", c("GLP-1RA", "SGLT2i"), 
           cex = 0.85, 
           pch = 21, 
           bg = "transparent",
           pt.bg = c("transparent", "darkslategray3" ),
           col = c("indianred3", "darkslategray3" ), 
           box.col = "transparent")
  }else{
    print("hi")
  }
}

dev.off() 


