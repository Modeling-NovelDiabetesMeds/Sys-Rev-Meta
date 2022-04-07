#         #        #         #        #         #        #         #          # 
#                       Meta-regression analysis
# Inputs  : df.csv 
# Outputs : - meta regression coefficients tables, and figures
#           for main outcomes, secondary outcomes and sensitivity analyses.
#           - ardse.csv contains per trial and outcome the simulated standard
#             errors for ARD point estimate after using bootstrapping
#         #        #         #        #         #        #         #          #

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
    
    meta <- fread("data/df.csv")

# A.   Data preparation (labels and ids for plots) -----
# generate id per trial/group for regression plots
      meta <- meta %>% 
        group_by(type, outcome) %>%
        mutate(id = row_number())  
      meta$id <- ifelse(meta$trialname == "LEADER/SUSTAIN-6", 10, meta$id)
        id0 <- median(meta$id[meta$trialname =="SUSTAIN-6"])
        meta$id <- ifelse(meta$trialname == "SUSTAIN-6", id0, meta$id)

# Outcome labels, used for later in plot
        
        meta$outcomelab <- meta$outcomen
        meta$outcomelab <- ifelse(meta$outcomen == "CVMort", 
                                  "Cardiovascular Mortality", meta$outcomelab )
        meta$outcomelab <- ifelse(meta$outcomen == "HospHF",
                                  "Hospitalisation for Heart Failure", 
                                  meta$outcomelab )
        meta$outcomelab <- ifelse(meta$outcomen == "sustGFRdecl", 
                                  "Composite Renal Outcome", meta$outcomelab )
        meta$outcomelab <- ifelse(meta$outcomen == "MACE", 
                                  "MACE", meta$outcomelab )
        meta$outcomelab <- ifelse(meta$outcomen == "allcauseMort", 
                                  "All-cause Mortality", meta$outcomelab )
        meta$outcomelab <- ifelse(meta$outcomen == "MI", 
                                  "Myocardial Infarction", meta$outcomelab )
        meta$outcomelab <- ifelse(meta$outcomen == "stroke", 
                                  "Stroke", meta$outcomelab )

# B.   Computing variables of interest (loghr)       ----
# Create CVD baseline incidence rate variable (events per 100 patient year)
#       cvdepy: cardiovascular death events pero 100 patient-years
#       from p.rate (placebo group incidence rate per outcome)  
    meta$cvdepy <- NA
    meta$cvdepy <- ifelse(meta$outcomen == "CVMort", 
                          meta$p.rate, 
                          meta$cvdepy)
  # spread value of cardiovascular incidence rate to all outcomes
      meta <- meta%>%
        group_by(trialname) %>% 
        mutate(cvdepy = ifelse(is.na(cvdepy), 
                               median(cvdepy, na.rm = TRUE), 
                               cvdepy))

# Hazard Ratios: log, variance and standard errors
  
  meta$loghr  <- log(meta$hr)
  meta$logvi  <- ((log(meta$uci) - log(meta$lci)) / (2 * qnorm(.975) ))^2
  meta$logsei <- sqrt(meta$logvi)
  meta$vi     <- (exp(meta$logvi)-1)*exp(2*meta$loghr + meta$logvi) # log normal variance
  meta$sei    <- sqrt(meta$vi)

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
    meta$rateratio <- meta$t.rate / meta$p.rate
    meta$se.rateratio <- sqrt(1/meta$n0 + 1/meta$n1) #se of ln of rate ratio
    meta$lograteratio <- log(meta$rateratio)
    meta$vi.rateratio <- meta$se.rateratio^2 

    meta$npy1 <- 100*(meta$n1 / meta$t.rate) # Total patient years treat
    meta$npy0 <- 100*(meta$n0 / meta$p.rate) # Total patient years control

    # List to store results
    lse <- list()
    # Vector of outcomes interest
    v.outcomes <- c("Cardiovascular Mortality", "MACE",
                    "Hospitalisation for Heart Failure", 
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
    arddf <- meta[, c("trialname","type", "outcomelab", "ard", "ardse",
                      "ardse2", "ardvi", "ardvi.2")]
# Export ardse.csv dataframe (used in forest plots)  
      write.csv(arddf, "data/ardse.csv")

# Wsize/ Wsize2: Parameter used for bubble sizes, function proportional of variances
     meta <- meta %>%
        group_by(outcome) %>%
        mutate(ran = range(1/(10000*ardvi.2), na.rm = T)[2] -
               range(1/(10000*ardvi.2), na.rm = T)[1])
    meta$wsize   <- 1/(100*meta$ardvi.2)/meta$ran
    meta$wsize2  <- 8*log(1+meta$wsize^(1/3) ) 

# I.1  Meta regression: Models with cvdepy as mediator: -----
  meta$lograte <- log(meta$cvdepy) 
#  7 outcomes, model by drugclass
  # lists to store models
    mod.glp1  <- list() 
    mod.sglt2 <- list()
    v.type <- c("GLP1","SGLT2")

for(o in 1:length(v.outcomes)){ #loop runs over each outcome 'o' and each drug class 't'
  
  for(t in 1:length(v.type)){
    df.0 <- subset(meta, subset = (outcomelab == v.outcomes[o] &
                                   type == v.type[t]))
  
  if(v.outcomes[o] == "HospHF"){ #leave out two trials with inconsistent definition for HospHF
    df <- subset(df.0, 
                 subset = (trialname != "SOLOIST-WHF" & trialname != "SCORED") 
                 )
  } else {
    df <- df.0
  }

  # Model (random effects) for outcome o class t
  r <- rma(100*ard, (100^2)*ardvi.2, # times 100 to re-scale to percentage points
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
    cc <- as.data.frame(cbind(rep(v.outcomes[1:4],2), 
                          c(rep("GLP-1RA",4), 
                            rep("SGLT2i",4)), 
                          rep(NA,8), 
                          rep(NA,8),
                          rep(NA,8),
                          rep(NA,8)))
# secondary outcomes
    cc.2 <- as.data.frame(cbind(rep(v.outcomes[5:7],2), 
                              c(rep("GLP-1RA",3), 
                                rep("SGLT2i",3)), 
                              rep(NA,6), 
                              rep(NA,6),
                              rep(NA,6),
                              rep(NA,6)))
names(cc)   <- c("Outcome", "Class", "Slope","CI.lb", "CI.ub", "P-value")
names(cc.2) <- c("Outcome", "Class", "Slope","CI.lb", "CI.ub", "P-value")

    # extract quantities of interest summary lists
      # Primary outcomes
      for(i in 1:4){
        cc[i,3] <- round(s.g[[i]]$beta[2],3)
        cc[i,4] <- round(s.g[[i]]$ci.lb[2],3)
        cc[i,5] <- round(s.g[[i]]$ci.ub[2],3)
        cc[i,6] <- round(s.g[[i]]$pval[2],3)
        cc[i+4,3] <- round(s.s[[i]]$beta[2],3)
        cc[i+4,4] <- round(s.s[[i]]$ci.lb[2],3)
        cc[i+4,5] <- round(s.s[[i]]$ci.ub[2],3)
        cc[i+4,6] <- round(s.s[[i]]$pval[2],3)
        }
      # Secondary outcomes
      for(i in 5:7){
        cc.2[i-4,3] <- round(s.g[[i]]$beta[2],3)
        cc.2[i-4,4] <- round(s.g[[i]]$ci.lb[2],3)
        cc.2[i-4,5] <- round(s.g[[i]]$ci.ub[2],3)
        cc.2[i-4,6] <- round(s.g[[i]]$pval[2],3)
        cc.2[i-1,3] <- round(s.s[[i]]$beta[2],3)
        cc.2[i-1,4] <- round(s.s[[i]]$ci.lb[2],3)
        cc.2[i-1,5] <- round(s.s[[i]]$ci.ub[2],3)
        cc.2[i-1,6] <- round(s.s[[i]]$pval[2],3)
      }
  # Export into text tables
    stargazer(cc, out = "output/metareg_primary.txt",
              summary = F,type = "text", 
              title = "Meta-regression coefficients for primary outcomes, by drug class",
              notes = "Log hazard ratio (DV) and baseline cardiovascular risk (IV)")
    stargazer(cc.2, out = "output/metareg_secondary.txt",
              summary = F,type = "text", 
              title = "Meta-regression coefficients for primary outcomes, by drug class",
              notes = "Log hazard ratio (DV) and baseline cardiovascular risk (IV)")
    
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
  v.titles <- c("Cardiovascular Mortality",
                "MACE", 
                "Hospitalisation for Heart Failure",
                "Composite Renal Outcome"
                )



  #png("plots/metareg_ard_panel.png", width = 9, height = 9, units = 'in', res = 300)  
  par(oma = c(5,3,1,1), mfrow = c(2,2), mar = c(4,4,3,2)*0.75)
  for(i in 1:4){
    # Plot elements: 1. white canvas, 2. confint (sglt2 + glp1) + 
    #                3. bubbles 4. Model fit (line) 5. Trial number (text)
    #                6. titles (outcome)
  # 1. White canvas 
    plot(x = 0, y = 0, type ='n', 
         xlim = range(meta$cvdepy, na.rm = T), 
         ylim = c(-3,14), 
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
            col = alpha("darkslategray3", alpha = 0.2),
            border = "transparent")
    polygon(x = c(newx.glp, rev(newx.glp)), 
            y = c(preds.glp$ci.lb, rev(preds.glp$ci.ub)),
            col = alpha("indianred3", alpha = 0.3),
            border = "transparent")
  # 3. Point data (bubbles) with size = f(1/var)
  points(x =   l.do[[i]]$cvdepy, 
         y =   100*l.do[[i]]$ard, # 100 to make it percentaje points
         col = l.do[[i]]$border,
         bg =  l.do[[i]]$colnum2,
         pch = 21, 
         cex = ((0.25)*l.do[[i]]$wsize2)^(1.5))
  # Model for GLP1-RA
  segments(x0 = min(l.do[[i]]$cvdepy[l.do[[i]]$type == "GLP1"], na.rm = T), 
           x1 = max(l.do[[i]]$cvdepy[l.do[[i]]$type == "GLP1"], na.rm = T),
           y0 = mod.glp1[[i]]$beta[1] + 
             min(l.do[[i]]$cvdepy[l.do[[i]]$type == "GLP1"], na.rm = T)*mod.glp1[[i]]$beta[2], 
           y1 = mod.glp1[[i]]$beta[1] + 
             max(l.do[[i]]$cvdepy[l.do[[i]]$type == "GLP1"], na.rm = T)*mod.glp1[[i]]$beta[2],
           col =  "indianred3",
           lwd = 2)
  # 4. Model for SGLT2i
  segments(x0 = min(l.do[[i]]$cvdepy[l.do[[i]]$type == "SGLT2"& !is.na(l.do[[i]]$ard)], na.rm = T), 
           x1 = max(l.do[[i]]$cvdepy[l.do[[i]]$type == "SGLT2"& !is.na(l.do[[i]]$ard)], na.rm = T),
           y0 = mod.sglt2[[i]]$beta[1] + 
             min(l.do[[i]]$cvdepy[l.do[[i]]$type == "SGLT2"& !is.na(l.do[[i]]$ard)], na.rm = T)*mod.sglt2[[i]]$beta[2], 
           y1 = mod.sglt2[[i]]$beta[1] + 
             max(l.do[[i]]$cvdepy[l.do[[i]]$type == "SGLT2"& !is.na(l.do[[i]]$ard)], na.rm = T)*mod.sglt2[[i]]$beta[2],
           col =  "darkslategray4",
           lwd = 2)
  # 4. Add id to identify trial
  text(x = jitter(l.do[[i]]$cvdepy[!is.na(l.do[[i]]$ard)],20), 
       y = jitter(100*l.do[[i]]$ard[!is.na(l.do[[i]]$ard)],20),  
       l.do[[i]]$id[!is.na(l.do[[i]]$ard)], cex = 0.65 )
  # 5. title label per outcome
  title(main = v.titles[i], cex.main = 0.9, line = 0.7)
  # Axes titles
  if(i  %% 2 != 0){
  title(ylab = "Absolute Risk Difference (%)", cex.lab = 0.8, line = 2)
    } else {
      title(ylab = "", cex.lab = 0.6, line = 1)
    }
  if(i  >2){
  title(xlab = "Cardiovascular Mortality Rate (/100py) in Control Group", cex.lab = 0.8, line = 2)
   } else {
      title(xlab = "", cex.lab = 1, line = 2)
    }
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
       at  = c(-4,-2, 0,2,4,6, 8,10,12),
       labels = c("",-2, 0,2,4,6, 8,10,12),
       cex.axis = 0.75, 
       col.axis = "gray40", 
       col.ticks = "gray40", 
       col = "gray40", 
       las = 1)  
  if(i == 2){
    legend("topright", c("GLP-1RA", "SGLT2i"), 
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
  #dev.off()  


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
  
  # png("plots/metareg_ard_panel_secondary.png", width = 15, height = 5, units = 'in', res = 300)  
  par(oma = c(3,2,0,0), mfrow = c(1,3), mar = c(4,4,3,2)*0.5)
    for(i in 5:7){
    # Plot elements: 1. white canvas, 2. confint (sglt2 + glp1) + 
    #                3. bubbles 4. Model fit (line) 5. Trial number (text)
    #                6. titles (outcome)
    # 1. White canvas 
    plot(x = 0, y = 0, type ='n', 
         xlim = range(meta$cvdepy, na.rm = T), 
         ylim = c(-3,14), 
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
            col = alpha("darkslategray3", alpha = 0.2),
            border = "transparent")
    polygon(x = c(newx.glp, rev(newx.glp)), 
            y = c(preds.glp$ci.lb, rev(preds.glp$ci.ub)),
            col = alpha("indianred3", alpha = 0.3),
            border = "transparent")
    # 3. Point data (bubbles) with size = f(1/var)
    points(x =   l.do.2[[i-4]]$cvdepy, 
           y =   100*l.do.2[[i-4]]$ard, # 100 to make it percentaje points
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
         y = jitter(100*l.do.2[[i-4]]$ard[!is.na(l.do.2[[i-4]]$ard)],20),  
         l.do.2[[i-4]]$id[!is.na(l.do.2[[i-4]]$ard)], cex = 0.65 )
    # 5. title label per outcome
    title(main = v.titles[i-4], cex.main = 1.5, line = 0.1)
    # Axes titles
    if(i  %% 2 != 0){
      title(ylab = "Absolute Risk Difference (%)", cex.lab = 0.8, line = 2)
    } else {
      title(ylab = "", cex.lab = 0.6, line = 1)
    }
    if(i  >2){
      title(xlab = "Cardiovascular Mortality Rate (/100py) in Control Group", cex.lab = 0.8, line = 2)
    } else {
      title(xlab = "", cex.lab = 1, line = 2)
    }
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
         at  = c(-4,-2, 0,2,4,6, 8,10,12),
         labels = c("",-2, 0,2,4,6, 8,10,12),
         cex.axis = 0.75, 
         col.axis = "gray40", 
         col.ticks = "gray40", 
         col = "gray40", 
         las = 1)  
    if(i == 2){
      legend("topright", c("GLP-1RA", "SGLT2i"), 
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
  #dev.off() 
# S.1  Sensitivity Analysis CVD: Remove SCORE and SOLOIST --------
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

#png("plots/metareg_ard_cdv_1.png", width = 6, height = 6, units = 'in', res = 300)  
  par(oma = c(3,2,1,1), mfrow = c(1,1), mar = c(4,4,3,2)*0.75)
    # 1. White canvas 
    plot(x = 0, y = 0, type ='n', 
         xlim = range(meta$cvdepy, na.rm = T), 
         ylim = c(-3,14), 
         xlab= "", 
         ylab ="", 
         axes = F)
          # 2. Confidence interval: predict confint for both glp-1ras and sglt2i models
          # create domain for independent variable
          newx.glp <- seq(min(df.sens$cvdepy[df.sens$type == "GLP1"], na.rm = T),
                          max(df.sens$cvdepy[df.sens$type == "GLP1"], na.rm = T),
                          length.out = 1000)
          # predict glp1 ard for given x1 ( ard =f(x1) )
          preds.glp <- predict(msens.g, 
                               newmods = newx.glp)
          # create domain for independent variable
          newx.sglt <- seq(min(df.sens$cvdepy[df.sens$type == "SGLT2"], na.rm = T),
                           max(df.sens$cvdepy[df.sens$type == "SGLT2"], na.rm = T),
                           length.out = 1000)
          # predict sglt2 ard for given x2 ( ard =f(x2) )
          preds.sglt <- predict(msens.s, 
                                newmods = newx.sglt)
          # Create CI polygons with computed coordinates
          polygon(x = c(newx.sglt, rev(newx.sglt)), 
                  y = c(preds.sglt$ci.lb, rev(preds.sglt$ci.ub)),
                  col = alpha("darkslategray3", alpha = 0.1),
                  border = "transparent")
          polygon(x = c(newx.glp, rev(newx.glp)), 
                  y = c(preds.glp$ci.lb, rev(preds.glp$ci.ub)),
                  col = alpha("indianred3", alpha = 0.3),
                  border = "transparent")
          # 3. Point data (bubbles) with size = f(1/var)
          points(x =   df.sens$cvdepy, 
                 y =   100*df.sens$ard, # 100 to make it percentaje points
                 col = df.sens$border,
                 bg =  df.sens$colnum2,
                 pch = 21, 
                 cex = ((0.25)*df.sens$wsize2)^(1.5))
          # Model for GLP1-RA
          segments(x0 = min(df.sens$cvdepy[df.sens$type == "GLP1"], na.rm = T), 
                   x1 = max(df.sens$cvdepy[df.sens$type == "GLP1"], na.rm = T),
                   y0 = msens.g$beta[1] + 
                     min(df.sens$cvdepy[df.sens$type == "GLP1"], na.rm = T)*msens.g$beta[2], 
                   y1 = msens.g$beta[1] + 
                     max(df.sens$cvdepy[df.sens$type == "GLP1"], na.rm = T)*msens.g$beta[2],
                   col =  "indianred2",
                   lwd = 2)
          # 4. Model for SGLT2i
          segments(x0 = min(df.sens$cvdepy[df.sens$type == "SGLT2"], na.rm = T), 
                   x1 = max(df.sens$cvdepy[df.sens$type == "SGLT2"], na.rm = T),
                   y0 = msens.s$beta[1] + 
                     min(df.sens$cvdepy[df.sens$type == "SGLT2"], na.rm = T)*msens.s$beta[2], 
                   y1 = msens.s$beta[1] + 
                     max(df.sens$cvdepy[df.sens$type == "SGLT2"], na.rm = T)*msens.s$beta[2],
                   col =  "darkslategray3",
                   lwd = 2)
          # 4. Add id to identify trial
          text(x = jitter(df.sens$cvdepy[!is.na(df.sens$ard)],20), 
               y = jitter(100*df.sens$ard[!is.na(df.sens$ard)],20),  
               df.sens$id[!is.na(df.sens$ard)], cex = 0.65 )
          # 5. title label per outcome
          title(main = "Cardiovascular Mortality", cex.main = 0.9, line = 1)
          # Axes titles
            title(ylab = "Absolute Risk Difference (%)", cex.lab = 0.8, line = 2)
            title(xlab = "Cardiovascular Mortality Rate (/100py) in Control Group", cex.lab = 0.8, line = 2)
                    # Axes
          axis(1, 
               at = pretty(range(meta$cvdepy, na.rm = T)),
               cex.axis = 0.65,
               col.axis = "gray40", 
               col.ticks = "gray40", 
               col = "gray40", 
               las = 1)
          axis(2, 
               at  = c(-4,-2, 0,2,4,6, 8,10,12),
               labels = c("",-2, 0,2,4,6, 8,10,12),
               cex.axis = 0.65, 
               col.axis = "gray40", 
               col.ticks = "gray40", 
               col = "gray40", 
               las = 1)  

            legend("topright", c("GLP-1RA", "SGLT2i"), 
                   cex = 0.85, 
                   pch = 21, 
                   bg = "transparent",
                   pt.bg = c("transparent", "darkslategray3" ),
                   col = c("indianred3", "darkslategray3" ), 
                   box.col = "transparent")
   
      #  dev.off()  

# S.2  Sensitivity Analysis CVD: Remove ELIXA and SOLOIST --------
  df.sens  <- meta[meta$outcomen == "CVMort" & 
                         !(meta$trialname == "ELIXA" |
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
  
#  png("plots/metareg_ard_cdv_2.png", width = 6, height = 6, units = 'in', res = 300)  
  par(oma = c(3,2,1,1), mfrow = c(1,1), mar = c(4,4,3,2)*0.75)
  # 1. White canvas 
  plot(x = 0, y = 0, type ='n', 
       xlim = range(meta$cvdepy, na.rm = T), 
       ylim = c(-5,14), 
       xlab= "", 
       ylab ="", 
       axes = F)
  # 2. Confidence interval: predict confint for both glp-1ras and sglt2i models
  # create domain for independent variable
  newx.glp <- seq(min(df.sens$cvdepy[df.sens$type == "GLP1"], na.rm = T),
                  max(df.sens$cvdepy[df.sens$type == "GLP1"], na.rm = T),
                  length.out = 1000)
  # predict glp1 ard for given x1 ( ard =f(x1) )
  preds.glp <- predict(msens.g, 
                       newmods = newx.glp)
  # create domain for independent variable
  newx.sglt <- seq(min(df.sens$cvdepy[df.sens$type == "SGLT2"], na.rm = T),
                   max(df.sens$cvdepy[df.sens$type == "SGLT2"], na.rm = T),
                   length.out = 1000)
  # predict sglt2 ard for given x2 ( ard =f(x2) )
  preds.sglt <- predict(msens.s, 
                        newmods = newx.sglt)
  # Create CI polygons with computed coordinates
  polygon(x = c(newx.sglt, rev(newx.sglt)), 
          y = c(preds.sglt$ci.lb, rev(preds.sglt$ci.ub)),
          col = alpha("darkslategray3", alpha = 0.1),
          border = "transparent")
  polygon(x = c(newx.glp, rev(newx.glp)), 
          y = c(preds.glp$ci.lb, rev(preds.glp$ci.ub)),
          col = alpha("indianred3", alpha = 0.3),
          border = "transparent")
  # 3. Point data (bubbles) with size = f(1/var)
  points(x =   df.sens$cvdepy, 
         y =   100*df.sens$ard, # 100 to make it percentaje points
         col = df.sens$border,
         bg =  df.sens$colnum2,
         pch = 21, 
         cex = ((0.25)*df.sens$wsize2)^(1.5))
  # Model for GLP1-RA
  segments(x0 = min(df.sens$cvdepy[df.sens$type == "GLP1"], na.rm = T), 
           x1 = max(df.sens$cvdepy[df.sens$type == "GLP1"], na.rm = T),
           y0 = msens.g$beta[1] + 
             min(df.sens$cvdepy[df.sens$type == "GLP1"], na.rm = T)*msens.g$beta[2], 
           y1 = msens.g$beta[1] + 
             max(df.sens$cvdepy[df.sens$type == "GLP1"], na.rm = T)*msens.g$beta[2],
           col =  "indianred2")
  # 4. Model for SGLT2i
  segments(x0 = min(df.sens$cvdepy[df.sens$type == "SGLT2"], na.rm = T), 
           x1 = max(df.sens$cvdepy[df.sens$type == "SGLT2"], na.rm = T),
           y0 = msens.s$beta[1] + 
             min(df.sens$cvdepy[df.sens$type == "SGLT2"], na.rm = T)*msens.s$beta[2], 
           y1 = msens.s$beta[1] + 
             max(df.sens$cvdepy[df.sens$type == "SGLT2"], na.rm = T)*msens.s$beta[2],
           col =  "darkslategray3")
  # 4. Add id to identify trial
  text(x = jitter(df.sens$cvdepy[!is.na(df.sens$ard)],20), 
       y = jitter(100*df.sens$ard[!is.na(df.sens$ard)],20),  
       df.sens$id[!is.na(df.sens$ard)], cex = 0.65 )
  
  
  # 5. title label per outcome
  title(main = "Cardiovascular Mortality", cex.main = 0.9, line = 1)
  # Axes titles
  title(ylab = "Absolute Risk Difference (%)", cex.lab = 0.8, line = 2)
  title(xlab = "Cardiovascular Mortality Rate (/100py) in Control Group", cex.lab = 0.8, line = 2)
  # Axes
  axis(1, 
       at = pretty(range(meta$cvdepy, na.rm = T)),
       cex.axis = 0.65,
       col.axis = "gray40", 
       col.ticks = "gray40", 
       col = "gray40", 
       las = 1)
  axis(2, 
       at  = c(-5,-2, 0,2,4,6, 8,10,12),
       labels = c(-5,-2, 0,2,4,6, 8,10,12),
       cex.axis = 0.65, 
       col.axis = "gray40", 
       col.ticks = "gray40", 
       col = "gray40", 
       las = 1)  
  
  legend("topright", c("GLP1-RAs", "SGLT2i"), 
         cex = 0.85, 
         pch = 21, 
         bg = "transparent",
         pt.bg = c("transparent", "darkslategray3" ),
         col = c("indianred3", "darkslategray3" ), 
         box.col = "transparent")
  
#  dev.off()              
        
# S.3  Sensitivity Analysis Renal: Remove ELIXA and CREDENCE --------
  df.sens  <- meta[meta$outcomen == "sustGFRdecl" & 
                       !(meta$trialname == "ELIXA" |
                           meta$trialname == "CREDENCE" |
                           meta$trialname == "DAPA-CKD" |
                           meta$trialname == "DAPA-HF" |
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
  
#  png("plots/metareg_ard_ren.png", width = 6, height = 6, units = 'in', res = 300)  
    par(oma = c(3,2,1,1), mfrow = c(1,1), mar = c(4,4,3,2)*0.75)
    # 1. White canvas 
    plot(x = 0, y = 0, type ='n', 
         xlim = range(meta$cvdepy, na.rm = T), 
        ylim = c(-3.93,14),  
         xlab= "", 
         ylab ="", 
         axes = F)
    # 2. Confidence interval: predict confint for both glp-1ras and sglt2i models
    # create domain for independent variable
    newx.glp <- seq(min(df.sens$cvdepy[df.sens$type == "GLP1" & !is.na(df.sens$ard) ], na.rm = T),
                    max(df.sens$cvdepy[df.sens$type == "GLP1"&  !is.na(df.sens$ard)], na.rm = T),
                    length.out = 1000)
    # predict glp1 ard for given x1 ( ard =f(x1) )
    preds.glp <- predict(msens.g, 
                         newmods = newx.glp)
    # create domain for independent variable
    newx.sglt <- seq(min(df.sens$cvdepy[df.sens$type == "SGLT2"& !is.na(df.sens$ard)], na.rm = T),
                     max(df.sens$cvdepy[df.sens$type == "SGLT2"& !is.na(df.sens$ard)], na.rm = T),
                     length.out = 1000)
    # predict sglt2 ard for given x2 ( ard =f(x2) )
    preds.sglt <- predict(msens.s, 
                          newmods = newx.sglt)
    # Create CI polygons with computed coordinates
    polygon(x = c(newx.sglt, rev(newx.sglt)), 
            y = c(preds.sglt$ci.lb, rev(preds.sglt$ci.ub)),
            col = alpha("darkslategray3", alpha = 0.1),
            border = "transparent")
    polygon(x = c(newx.glp, rev(newx.glp)), 
            y = c(preds.glp$ci.lb, rev(preds.glp$ci.ub)),
            col = alpha("indianred3", alpha = 0.3),
            border = "transparent")
    # 3. Point data (bubbles) with size = f(1/var)
    points(x =   df.sens$cvdepy, 
           y =   100*df.sens$ard, # 100 to make it percentaje points
           col = df.sens$border,
           bg =  df.sens$colnum2,
           pch = 21, 
           cex = ((0.25)*df.sens$wsize2)^(1.5))
    # Model for GLP1-RA
    segments(x0 = min(df.sens$cvdepy[df.sens$type == "GLP1"& !is.na(df.sens$ard)], na.rm = T), 
             x1 = max(df.sens$cvdepy[df.sens$type == "GLP1"& !is.na(df.sens$ard)], na.rm = T),
             y0 = msens.g$beta[1] + 
               min(df.sens$cvdepy[df.sens$type == "GLP1"& !is.na(df.sens$ard)], na.rm = T)*msens.g$beta[2], 
             y1 = msens.g$beta[1] + 
               max(df.sens$cvdepy[df.sens$type == "GLP1"& !is.na(df.sens$ard)], na.rm = T)*msens.g$beta[2],
             col =  "indianred2")
    # 4. Model for SGLT2i
    segments(x0 = min(df.sens$cvdepy[df.sens$type == "SGLT2"& !is.na(df.sens$ard)], na.rm = T), 
             x1 = max(df.sens$cvdepy[df.sens$type == "SGLT2"& !is.na(df.sens$ard)], na.rm = T),
             y0 = msens.s$beta[1] + 
               min(df.sens$cvdepy[df.sens$type == "SGLT2"& !is.na(df.sens$ard)], na.rm = T)*msens.s$beta[2], 
             y1 = msens.s$beta[1] + 
               max(df.sens$cvdepy[df.sens$type == "SGLT2"& !is.na(df.sens$ard)], na.rm = T)*msens.s$beta[2],
             col =  "darkslategray3")
    # 4. Add id to identify trial
    text(x = jitter(df.sens$cvdepy[!is.na(df.sens$ard)],20), 
         y = jitter(100*df.sens$ard[!is.na(df.sens$ard)],20),  
         df.sens$id[!is.na(df.sens$ard)], cex = 0.65 )
    
   
    # 5. title label per outcome
    title(main = "Composite Renal Outcome", cex.main = 0.9, line = 1)
    # Axes titles
    title(ylab = "Absolute Risk Difference (%)", cex.lab = 0.8, line = 2)
    title(xlab = "Cardiovascular Mortality Rate (/100py) in Control Group", cex.lab = 0.8, line = 2)
    # Axes
    axis(1, 
         at = pretty(range(meta$cvdepy, na.rm = T)),
         cex.axis = 0.65,
         col.axis = "gray40", 
         col.ticks = "gray40", 
         col = "gray40", 
         las = 1)
    axis(2, 
         #at = pretty(range(100*df.sens$ard, na.rm = T), n = 5),
         at  = c(-4,-2, 0,2,4,6, 8,10,12),
         labels = c(-4,-2, 0,2,4,6, 8,10,12),
         cex.axis = 0.65,
         col.axis = "gray40", 
         col.ticks = "gray40", 
         col = "gray40", 
         las = 1)  
    
    legend("topright", c("GLP1-RAs", "SGLT2i"), 
           cex = 0.85, 
           pch = 21, 
           bg = "transparent",
           pt.bg = c("transparent", "darkslategray3" ),
           col = c("indianred3", "darkslategray3" ), 
           box.col = "transparent")
    
#    dev.off()  
        
  
    
#                     End of script               #       
    
    
