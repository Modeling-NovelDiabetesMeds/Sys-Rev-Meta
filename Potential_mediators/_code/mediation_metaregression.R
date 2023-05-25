rm(list = ls())
library(stargazer)
library(data.table)
library(tidyr)
library(dplyr)
library(scales)
library(stringr)
library(metafor)
library(meta)
library(lme4)
library(tidyr)
library(data.table)
library(modelsummary)
set.seed(85040)
##  Set working directory
wd <- "~/Documents/Sys-Rev-Meta/Potential_mediators"
setwd(wd)
# Load data ------
  med <- read.csv("_data/mediators_md.csv")
  hr <- read.csv("_data/df.csv")

# Reshape and merge data
    medw <- dcast(med,   trialname ~ outcome, value.var = "diff")
    medt <- dcast(med,   trialname ~ outcome, value.var = "time")
    medw <- medw[,
                 c("trialname", "bodyweight", 
                   "hba1c", "sbp", "hematocrit", "uacr")]
    names(medw) <- c("trialname", "diff.bodyweight", "diff.hba1c", 
                     "diff.sbp", "diff.hematocrit", "diff.uacr")
    medt <- medt[,c("trialname", "bodyweight", "hba1c", 
                    "sbp", "hematocrit", "uacr")]
    names(medt) <- c("trialname", "t.bodyweight", 
                     "t.hba1c", "t.sbp", "t.hematocrit", "t.uacr")
    
    med <- merge(medw, medt, by = "trialname")
    med <- med[,c("trialname",
                  "diff.hba1c","t.hba1c",
                  "diff.bodyweight","t.bodyweight",
                  "diff.sbp","t.sbp",
                  "diff.hematocrit","t.hematocrit",
                  "diff.uacr","t.uacr")]
  # Merge
  df <- merge(hr, med, by = "trialname", all.x = TRUE)
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

# OUTCOME LABELS
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
    df$col <- ifelse(df$Class == "SGLT2i", "slateblue4", "maroon4")
    df$col <- alpha(df$col, alpha = 0.5)
    
# Model for each outcome(logHR), mediator (iv), and drugclass (type or subset)
  # Vectors to run over models
      v.outcome <- unique(df$outcome)
      v.class <- unique(df$Class)
  # Data frames to collect results    
      coef.hba1c <- data.frame(Class = NA,
                               Mediator =  rep("hba1c", 14),
                               Outcome = NA, 
                               Slope =NA,
                               Pval = NA,
                               LCI = NA,
                               UCI = NA)
      coef.sbp <- data.frame(Class = NA,
                               Mediator =  rep("sbp", 14),
                               Outcome = NA, 
                               Slope =NA,
                               Pval = NA,
                               LCI = NA,
                               UCI = NA)
      coef.bw<- data.frame(Class = NA,
                               Mediator =  rep("bw", 14),
                               Outcome = NA, 
                               Slope =NA,
                               Pval = NA,
                               LCI = NA,
                               UCI = NA)
      coef.hem<- data.frame(Class = NA,
                           Mediator =  rep("hematocrit", 14),
                           Outcome = NA, 
                           Slope =NA,
                           Pval = NA,
                           LCI = NA,
                           UCI = NA)
      coef.uacr<- data.frame(Class = NA,
                            Mediator =  rep("uacr", 14),
                            Outcome = NA, 
                            Slope =NA,
                            Pval = NA,
                            LCI = NA,
                            UCI = NA)
# List to store models for plotting
    l.bw   <- list()
    l.uacr <- list()
    l.sbp  <- list()
    l.hema <- list()
    l.hba1c <- list()
    listnames <- c(paste0("glp.",v.outcome), paste0("sglt.", v.outcome))
    
    l.hba1c.int   <- list()
    l.bw.int   <- list()
    l.sbp.int  <- list()
    l.hema.int   <- list()
    l.uacr.int  <- list()
    
    l.bw.anova   <- list()
    l.hba1c.anova <- list()
    l.sbp.anova  <- list()
    l.hema.anova  <- list()
    l.uacr.anova  <- list()

# Mediator: Hba1c====
    for(t in 1:length(v.class)){
    for(o in 1:length(v.outcome)){
          d0 <- subset(df, 
                       subset = (df$outcome == v.outcome[o] & df$Class == v.class[t]))  
          m0 <- rma(yi = loghr, 
                    vi = logvi, 
                    mods = ~ diff.hba1c,
                    data   = d0, 
                    method ="REML", 
                    slab = trialname
                    )
          s0 <- summary(m0)
          l.hba1c[[o + 7*(t-1)]] <- m0
          coef.hba1c[o + 7*(t-1),1] <- v.class[t]
          coef.hba1c[o + 7*(t-1),3] <- v.outcome[o]
          coef.hba1c[o + 7*(t-1),4] <- s0$beta[2] 
          coef.hba1c[o + 7*(t-1),5] <- s0$pval[2] 
          coef.hba1c[o + 7*(t-1),6] <- s0$ci.lb[2] 
          coef.hba1c[o + 7*(t-1),7] <- s0$ci.ub[2] 
    }
    }
            # interaction
    for(o in 1:length(v.outcome)){
      for(j in 1:2){
      d0 <- subset(df, 
                   subset = (df$outcome == v.outcome[o]))  
      m0 <- rma(yi = loghr, 
                vi = logvi, 
                mods = ~ diff.hba1c* as.factor(Class),
                data   = d0, 
                method ="REML", 
                slab = trialname
      )
      m1 <- rma(yi = loghr, 
                vi = logvi, 
                mods = ~ diff.hba1c+ as.factor(Class),
                data   = d0, 
                method ="REML", 
                slab = trialname
      )
      l.hba1c.anova[[o]] <- anova(m1, m0)
      l.hba1c.int[[o]] <- m0
      s0 <- summary(m0)
     
      }
    }

# Mediator: Systolic blood pressure=====
for(t in 1:length(v.class)){
for(o in 1:length(v.outcome)){

    d0 <- subset(df, subset = (df$outcome == v.outcome[o] & df$Class == v.class[t]))  
    m0 <- rma(loghr, logvi, 
              mods = ~ diff.sbp,
              data   = d0, 
              method ="REML", 
              slab = trialname
    )
    l.sbp[[o + 7*(t-1)]] <- m0
    s0 <- summary(m0)
    coef.sbp[o + 7*(t-1),1] <- v.class[t]
    coef.sbp[o + 7*(t-1),3] <- v.outcome[o]
    coef.sbp[o + 7*(t-1),4] <- s0$beta[2] 
    coef.sbp[o + 7*(t-1),5] <- s0$pval[2] 
    coef.sbp[o + 7*(t-1),6] <- s0$ci.lb[2] 
    coef.sbp[o + 7*(t-1),7] <- s0$ci.ub[2] 
    print(unique(d0$outcome))
  }
}

  # Interaction
  for(o in 1:length(v.outcome)){
    
    d0 <- subset(df, subset = (df$outcome == v.outcome[o] ))  
    m0 <- rma(loghr, logvi, 
              mods = ~ diff.sbp*as.factor(Class),
              data   = d0, 
              method ="REML", 
              slab = trialname
    )
    m1 <- rma(loghr, logvi, 
              mods = ~ diff.sbp+as.factor(Class),
              data   = d0, 
              method ="REML", 
              slab = trialname
    )
    l.sbp.anova[[o]] <- anova(m1, m0)
    l.sbp.int[[o]] <- m0
    s0 <- summary(m0)
  }
# Mediator: Bodyweight=======
for(t in 1:length(v.class)){
for(o in 1:length(v.outcome)){

    d0 <- subset(df, subset = (df$outcome == v.outcome[o] & df$Class == v.class[t]))  
    m0 <- rma(loghr, logvi, 
              mods = ~ diff.bodyweight,
              data   = d0, 
              method ="REML", 
              slab = trialname
    )
    l.bw[[o + 7*(t-1)]] <- m0
    s0 <- summary(m0)
    coef.bw[o + 7*(t-1),1] <- v.class[t]
    coef.bw[o + 7*(t-1),3] <- v.outcome[o]
    coef.bw[o + 7*(t-1),4] <- s0$beta[2] 
    coef.bw[o + 7*(t-1),5] <- s0$pval[2] 
    coef.bw[o + 7*(t-1),6] <- s0$ci.lb[2] 
    coef.bw[o + 7*(t-1),7] <- s0$ci.ub[2] 
  }
}
    # interaction
  for(o in 1:length(v.outcome)){
    
    d0 <- subset(df, subset = (df$outcome == v.outcome[o]))  
    m0 <- rma(loghr, logvi, 
              mods = ~ diff.bodyweight*as.factor(Class),
              data   = d0, 
              method ="ML", 
              slab = trialname
    )
    m1 <- rma(loghr, logvi, 
              mods = ~ diff.bodyweight+as.factor(Class),
              data   = d0, 
              method ="ML", 
              slab = trialname
    )
    l.bw.anova[[o]] <- anova(m1, m0)
    l.bw.int[[o]] <- m0
    s0 <- summary(m0)
  }

          
# Mediator: Hematrocrit ====

for(t in 1:length(v.class)){
  for(o in 1:length(v.outcome)){
    
    d0 <- subset(df, subset = (df$outcome == v.outcome[o] & df$Class == v.class[t]))  
    tryCatch(m0 <- rma(loghr, logvi, 
                       mods = ~ diff.hematocrit,
                       data   = d0, 
                       method ="REML", 
                       slab = trialname),
             error = function(e){
               message(paste("ERROR",v.outcome[o], v.class[t]))
             }
    )
    l.hema[[o + 7*(t-1)]] <- m0
    s0 <- summary(m0)
    coef.hem[o + 7*(t-1),1] <- v.class[t]
    coef.hem[o + 7*(t-1),3] <- v.outcome[o]
    coef.hem[o + 7*(t-1),4] <- s0$beta[2] 
    coef.hem[o + 7*(t-1),5] <- s0$pval[2] 
    coef.hem[o + 7*(t-1),6] <- s0$ci.lb[2] 
    coef.hem[o + 7*(t-1),7] <- s0$ci.ub[2] 
  }
}


for(o in 1:length(v.outcome)){
  
  d0 <- subset(df, subset = (df$outcome == v.outcome[o]))  

  tryCatch(m0 <- rma(loghr, logvi, 
            mods = ~ diff.hematocrit*as.factor(Class),
            data   = d0, 
            method ="REML", 
            slab = trialname),
  error = function(e){ message(paste("ERROR",v.outcome[o], v.class[t]))}
           )
  tryCatch(m1 <- rma(loghr, logvi, 
                     mods = ~ diff.hematocrit+as.factor(Class),
                     data   = d0, 
                     method ="REML", 
                     slab = trialname),
           error = function(e){ message(paste("ERRORm1",v.outcome[o], v.class[t]))}
           )

  tryCatch(  l.hema.int[[o]] <- m0,
           error = function(e){ message(paste("ERRORint",v.outcome[o], v.class[t]))}
  )

}

# Mediator: uACR=====
for(t in 1:length(v.class)){
  for(o in 1:length(v.outcome)){
    
    d0 <- subset(df, subset = (df$outcome == v.outcome[o] & df$Class == v.class[t]))  
    tryCatch(m0 <- rma(loghr, logvi, 
              mods = ~ diff.uacr,
              data   = d0, 
              method ="REML", 
              slab = trialname),
             error = function(e){
               message(paste("ERROR",v.outcome[o], v.class[t]))
             }
    )
    l.uacr[[o + 7*(t-1)]] <- m0
    s0 <- summary(m0)
    coef.uacr[o + 7*(t-1),1] <- v.class[t]
    coef.uacr[o + 7*(t-1),3] <- v.outcome[o]
    coef.uacr[o + 7*(t-1),4] <- s0$beta[2] 
    coef.uacr[o + 7*(t-1),5] <- s0$pval[2] 
    coef.uacr[o + 7*(t-1),6] <- s0$ci.lb[2] 
    coef.uacr[o + 7*(t-1),7] <- s0$ci.ub[2] 
  }
}
# interaction
for(o in 1:length(v.outcome)){
  
  d0 <- subset(df, subset = (df$outcome == v.outcome[o]))  
  
  tryCatch(m0 <- rma(loghr, logvi, 
                     mods = ~ diff.uacr*as.factor(Class),
                     data   = d0, 
                     method ="REML", 
                     slab = trialname),
           error = function(e){ message(paste("ERROR",v.outcome[o], v.class[t]))}
  )
  tryCatch(m1 <- rma(loghr, logvi, 
                     mods = ~ diff.uacr+as.factor(Class),
                     data   = d0, 
                     method ="REML", 
                     slab = trialname),
           error = function(e){ message(paste("ERRORm1",v.outcome[o], v.class[t]))}
  )
  
  tryCatch(  l.uacr.int[[o]] <- m0,
             error = function(e){ message(paste("ERRORint",v.outcome[o], v.class[t]))}
  )
  
}


# Output: Coefficient tables =====
      l.bw    <- setNames(l.bw, listnames)
      l.sbp   <- setNames(l.sbp, listnames)
      l.hba1c <- setNames(l.hba1c, listnames)
      l.hema  <- setNames(l.hema, listnames)
      listnamesuacr <- c(paste0("glp.",v.outcome))
      l.uacr  <- setNames(l.uacr, listnames)
      
      
      listnames.int <- c(paste("int",v.outcome))
      l.bw.int    <- setNames(l.bw.int, listnames.int)
      l.sbp.int   <- setNames(l.sbp.int, listnames.int)
      l.hba1c.int <- setNames(l.hba1c.int, listnames.int)
      
      l.bw.anova <- setNames(l.bw.anova, v.outcome)
      l.sbp.anova <- setNames(l.sbp.anova, v.outcome)
      l.hba1c.anova <- setNames(l.hba1c.anova, v.outcome)
      
      
      coef.sbp.int <- data.frame(Class = "Interaction",
                             Mediator =  rep("sbp", 7),
                             Outcome = NA, 
                             Slope.intercept = NA,
                             Pval.intercept  = NA,
                             LCI.intercet = NA,
                             UCI.intercept = NA,
                             Slope.med =NA,
                             Pval.med = NA,
                             LCI.med = NA,
                             UCI.med = NA,
                             Slope.class =NA,
                             Pval.class = NA,
                             LCI.class = NA,
                             UCI.class = NA,
                             Slope.int =NA,
                             Pval.int = NA,
                             LCI.int = NA,
                             UCI.int = NA,
                             Pval.LRT = NA
                             )
      coef.hba1c.int <- data.frame(Class = "Interaction",
                                   Mediator =  rep("hba1c", 7),
                                   Outcome = NA, 
                                   Slope.intercept = NA,
                                   Pval.intercept  = NA,
                                   LCI.intercet = NA,
                                   UCI.intercept = NA,
                                   Slope.med =NA,
                                   Pval.med = NA,
                                   LCI.med = NA,
                                   UCI.med = NA,
                                   Slope.class =NA,
                                   Pval.class = NA,
                                   LCI.class = NA,
                                   UCI.class = NA,
                                   Slope.int =NA,
                                   Pval.int = NA,
                                   LCI.int = NA,
                                   UCI.int = NA,
                                   Pval.LRT = NA
                                   )
      coef.bw.int <- data.frame(Class = 'Interaction',
                                Mediator =  rep("bw", 7),
                                Outcome = NA,
                                Slope.intercept = NA,
                                Pval.intercept  = NA,
                                LCI.intercet = NA,
                                UCI.intercept = NA,
                                Slope.med =NA,
                                Pval.med = NA,
                                LCI.med = NA,
                                UCI.med = NA,
                                Slope.class =NA,
                                Pval.class = NA,
                                LCI.class = NA,
                                UCI.class = NA,
                                Slope.int =NA,
                                Pval.int = NA,
                                LCI.int = NA,
                                UCI.int = NA,
                                Pval.LRT = NA
                                )
      
      coef.hema.int <- data.frame(Class = 'Interaction',
                                Mediator =  rep("hematocrit", 7),
                                Outcome = NA,
                                Slope.intercept = NA,
                                Pval.intercept  = NA,
                                LCI.intercet = NA,
                                UCI.intercept = NA,
                                Slope.med =NA,
                                Pval.med = NA,
                                LCI.med = NA,
                                UCI.med = NA,
                                Slope.class =NA,
                                Pval.class = NA,
                                LCI.class = NA,
                                UCI.class = NA,
                                Slope.int =NA,
                                Pval.int = NA,
                                LCI.int = NA,
                                UCI.int = NA,
                                Pval.LRT = NA
      )
      
      coef.uacr.int <- data.frame(Class = 'Interaction',
                                Mediator =  rep("uacr", 7),
                                Outcome = NA,
                                Slope.intercept = NA,
                                Pval.intercept  = NA,
                                LCI.intercet = NA,
                                UCI.intercept = NA,
                                Slope.med =NA,
                                Pval.med = NA,
                                LCI.med = NA,
                                UCI.med = NA,
                                Slope.class =NA,
                                Pval.class = NA,
                                LCI.class = NA,
                                UCI.class = NA,
                                Slope.int =NA,
                                Pval.int = NA,
                                LCI.int = NA,
                                UCI.int = NA,
                                Pval.LRT = NA
      )
      
      
      for(i in 1:length(v.outcome)){
        coef.hba1c.int[i,3] <- v.outcome[i]
        coef.hba1c.int[i,4] <- l.hba1c.int[[i]]$beta[1]
        coef.hba1c.int[i,5] <- l.hba1c.int[[i]]$pval[1]
        coef.hba1c.int[i,6] <- l.hba1c.int[[i]]$ci.lb[1]
        coef.hba1c.int[i,7] <- l.hba1c.int[[i]]$ci.ub[1]
        coef.hba1c.int[i,8] <- l.hba1c.int[[i]]$beta[2]
        coef.hba1c.int[i,9] <- l.hba1c.int[[i]]$pval[2]
        coef.hba1c.int[i,10] <- l.hba1c.int[[i]]$ci.lb[2]
        coef.hba1c.int[i,11] <- l.hba1c.int[[i]]$ci.ub[2]
        coef.hba1c.int[i,12] <- l.hba1c.int[[i]]$beta[3]
        coef.hba1c.int[i,13] <- l.hba1c.int[[i]]$pval[3]
        coef.hba1c.int[i,14] <- l.hba1c.int[[i]]$ci.lb[3]
        coef.hba1c.int[i,15] <- l.hba1c.int[[i]]$ci.ub[3]
        coef.hba1c.int[i,16] <- l.hba1c.int[[i]]$beta[4]
        coef.hba1c.int[i,17] <- l.hba1c.int[[i]]$pval[4]
        coef.hba1c.int[i,18] <- l.hba1c.int[[i]]$ci.lb[4]
        coef.hba1c.int[i,19] <- l.hba1c.int[[i]]$ci.ub[4]
        coef.hba1c.int[i,20] <- l.hba1c.anova[[i]]$pval
        
        coef.bw.int[i,3] <- v.outcome[i]
        coef.bw.int[i,4] <- l.bw.int[[i]]$beta[1]
        coef.bw.int[i,5] <- l.bw.int[[i]]$pval[1]
        coef.bw.int[i,6] <- l.bw.int[[i]]$ci.lb[1]
        coef.bw.int[i,7] <- l.bw.int[[i]]$ci.ub[1]
        coef.bw.int[i,8] <- l.bw.int[[i]]$beta[2]
        coef.bw.int[i,9] <- l.bw.int[[i]]$pval[2]
        coef.bw.int[i,10] <- l.bw.int[[i]]$ci.lb[2]
        coef.bw.int[i,11] <- l.bw.int[[i]]$ci.ub[2]
        coef.bw.int[i,12] <- l.bw.int[[i]]$beta[3]
        coef.bw.int[i,13] <- l.bw.int[[i]]$pval[3]
        coef.bw.int[i,14] <- l.bw.int[[i]]$ci.lb[3]
        coef.bw.int[i,15] <- l.bw.int[[i]]$ci.ub[3]
        coef.bw.int[i,16] <- l.bw.int[[i]]$beta[4]
        coef.bw.int[i,17] <- l.bw.int[[i]]$pval[4]
        coef.bw.int[i,18] <- l.bw.int[[i]]$ci.lb[4]
        coef.bw.int[i,19] <- l.bw.int[[i]]$ci.ub[4]
        coef.bw.int[i,20] <- l.bw.anova[[i]]$pval
        
        
        coef.sbp.int[i,3] <- v.outcome[i]
        coef.sbp.int[i,4] <- l.sbp.int[[i]]$beta[1]
        coef.sbp.int[i,5] <- l.sbp.int[[i]]$pval[1]
        coef.sbp.int[i,6] <- l.sbp.int[[i]]$ci.lb[1]
        coef.sbp.int[i,7] <- l.sbp.int[[i]]$ci.ub[1]
        coef.sbp.int[i,8] <- l.sbp.int[[i]]$beta[2]
        coef.sbp.int[i,9] <- l.sbp.int[[i]]$pval[2]
        coef.sbp.int[i,10] <- l.sbp.int[[i]]$ci.lb[2]
        coef.sbp.int[i,11] <- l.sbp.int[[i]]$ci.ub[2]
        coef.sbp.int[i,12] <- l.sbp.int[[i]]$beta[3]
        coef.sbp.int[i,13] <- l.sbp.int[[i]]$pval[3]
        coef.sbp.int[i,14] <- l.sbp.int[[i]]$ci.lb[3]
        coef.sbp.int[i,15] <- l.sbp.int[[i]]$ci.ub[3]
        coef.sbp.int[i,16] <- l.sbp.int[[i]]$beta[4]
        coef.sbp.int[i,17] <- l.sbp.int[[i]]$pval[4]
        coef.sbp.int[i,18] <- l.sbp.int[[i]]$ci.lb[4]
        coef.sbp.int[i,19] <- l.sbp.int[[i]]$ci.ub[4]
        coef.sbp.int[i,20] <- l.sbp.anova[[i]]$pval
        
        coef.uacr.int[i,3] <- v.outcome[i]
        coef.uacr.int[i,4] <- l.uacr.int[[i]]$beta[1]
        coef.uacr.int[i,5] <- l.uacr.int[[i]]$pval[1]
        coef.uacr.int[i,6] <- l.uacr.int[[i]]$ci.lb[1]
        coef.uacr.int[i,7] <- l.uacr.int[[i]]$ci.ub[1]
        coef.uacr.int[i,8] <- l.uacr.int[[i]]$beta[2]
        coef.uacr.int[i,9] <- l.uacr.int[[i]]$pval[2]
        coef.uacr.int[i,10] <- l.uacr.int[[i]]$ci.lb[2]
        coef.uacr.int[i,11] <- l.uacr.int[[i]]$ci.ub[2]
        coef.uacr.int[i,12] <- l.uacr.int[[i]]$beta[3]
        coef.uacr.int[i,13] <- l.uacr.int[[i]]$pval[3]
        coef.uacr.int[i,14] <- l.uacr.int[[i]]$ci.lb[3]
        coef.uacr.int[i,15] <- l.uacr.int[[i]]$ci.ub[3]
        coef.uacr.int[i,16] <- l.uacr.int[[i]]$beta[4]
        coef.uacr.int[i,17] <- l.uacr.int[[i]]$pval[4]
        coef.uacr.int[i,18] <- l.uacr.int[[i]]$ci.lb[4]
        coef.uacr.int[i,19] <- l.uacr.int[[i]]$ci.ub[4]
        #coef.uacr.int[i,20] <- l.uacr.anova[[i]]$pval
        
        coef.hema.int[i,3] <- v.outcome[i]
        coef.hema.int[i,4] <- l.hema.int[[i]]$beta[1]
        coef.hema.int[i,5] <- l.hema.int[[i]]$pval[1]
        coef.hema.int[i,6] <- l.hema.int[[i]]$ci.lb[1]
        coef.hema.int[i,7] <- l.hema.int[[i]]$ci.ub[1]
        coef.hema.int[i,8] <- l.hema.int[[i]]$beta[2]
        coef.hema.int[i,9] <- l.hema.int[[i]]$pval[2]
        coef.hema.int[i,10] <- l.hema.int[[i]]$ci.lb[2]
        coef.hema.int[i,11] <- l.hema.int[[i]]$ci.ub[2]
        coef.hema.int[i,12] <- l.hema.int[[i]]$beta[3]
        coef.hema.int[i,13] <- l.hema.int[[i]]$pval[3]
        coef.hema.int[i,14] <- l.hema.int[[i]]$ci.lb[3]
        coef.hema.int[i,15] <- l.hema.int[[i]]$ci.ub[3]
        coef.hema.int[i,16] <- l.hema.int[[i]]$beta[4]
        coef.hema.int[i,17] <- l.hema.int[[i]]$pval[4]
        coef.hema.int[i,18] <- l.hema.int[[i]]$ci.lb[4]
        coef.hema.int[i,19] <- l.hema.int[[i]]$ci.ub[4]
        #coef.hema.int[i,20] <- l.hema.anova[[i]]$pval
      }  
      coef.ints <- rbind(coef.hba1c.int,coef.bw.int,coef.sbp.int, coef.hema.int, coef.uacr.int)
      coefs <- rbind(coef.hba1c,coef.bw,coef.sbp, coef.hem, coef.uacr)

### PLOTS ======

# Panel of six plots: by hbac for MACE MI stroke

df$col0 <- ifelse(df$type =="GLP1", "gold", "gray70")
df$col <- alpha(df$col0, 0.8)
df$cola <- alpha(df$col, 0.5)
cola <- alpha(c( "gold", "gray70"), 0.5)
# Subsetting and ranges

  # subsett by outcome
  dmace   <- df[df$outcome == "MACE",]
  dmi     <- df[df$outcome == "MI",]
  dstroke <- df[df$outcome == "stroke",]
  
  
  dacm <- df[df$outcome == "allcauseMort",]
  dcvm <- df[df$outcome == "CVMort",]
  dhf <- df[df$outcome == "HospHeartFailure",]
  dcro <- df[df$outcome == "sustGFRdecl",]
  



# I. Models predictions ===========
## I.A Domain of prediction (x)
  # I.A.a HBA1c
    # Primary outcomes
        # I.A.a.1 MACE
        newx.hba1c.mace.glp <- seq(min(dmace$diff.hba1c[dmace$type == "GLP1" & is.na(dmace$loghr) == F], na.rm = T),
                        max(dmace$diff.hba1c[dmace$type == "GLP1"& is.na(dmace$loghr) == F], na.rm = T),
                        length.out = 50)
        newx.hba1c.mace.sglt <- seq(min(dmace$diff.hba1c[dmace$type == "SGLT2"& is.na(dmace$loghr) == F], na.rm = T),
                         max(dmace$diff.hba1c[dmace$type == "SGLT2"& is.na(dmace$loghr) == F], na.rm = T),
                         length.out = 50)
        # I.A.a.2 MI
        newx.hba1c.mi.glp <- seq(min(dmi$diff.hba1c[dmi$type == "GLP1"& is.na(dmi$loghr) == F], na.rm = T),
                              max(dmi$diff.hba1c[dmi$type == "GLP1"& is.na(dmi$loghr) == F], na.rm = T),
                              length.out = 50)
        newx.hba1c.mi.sglt <- seq(min(dmi$diff.hba1c[dmi$type == "SGLT2"& is.na(dmi$loghr) == F], na.rm = T),
                               max(dmi$diff.hba1c[dmi$type == "SGLT2"& is.na(dmi$loghr) == F], na.rm = T),
                               length.out = 50)
        # I.A.a.3 stroke
        newx.hba1c.stroke.glp <- seq(min(dstroke$diff.hba1c[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                              max(dstroke$diff.hba1c[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                              length.out = 50)
        newx.hba1c.stroke.sglt <- seq(min(dstroke$diff.hba1c[dstroke$type == "SGLT2"& is.na(dstroke$loghr) == F], na.rm = T),
                               max(dstroke$diff.hba1c[dstroke$type == "SGLT2"& is.na(dstroke$loghr) == F], na.rm = T),
                               length.out = 50)
    # Secondary outcomes
        # I.A.a.4 all cause mortality
        newx.hba1c.acm.glp <- seq(min(dacm$diff.hba1c[dacm$type == "GLP1" & is.na(dacm$loghr) == F], na.rm = T),
                                   max(dacm$diff.hba1c[dacm$type == "GLP1"& is.na(dacm$loghr) == F], na.rm = T),
                                   length.out = 50)
        newx.hba1c.acm.sglt <- seq(min(dacm$diff.hba1c[dacm$type == "SGLT2"& is.na(dacm$loghr) == F], na.rm = T),
                                    max(dacm$diff.hba1c[dacm$type == "SGLT2"& is.na(dacm$loghr) == F], na.rm = T),
                                    length.out = 50)
        # I.A.a.5 CV mortality
        newx.hba1c.cvm.glp <- seq(min(dcvm$diff.hba1c[dcvm$type == "GLP1"& is.na(dcvm$loghr) == F], na.rm = T),
                                 max(dcvm$diff.hba1c[dcvm$type == "GLP1"& is.na(dcvm$loghr) == F], na.rm = T),
                                 length.out = 50)
        newx.hba1c.cvm.sglt <- seq(min(dcvm$diff.hba1c[dcvm$type == "SGLT2"& is.na(dcvm$loghr) == F], na.rm = T),
                                  max(dcvm$diff.hba1c[dcvm$type == "SGLT2"& is.na(dcvm$loghr) == F], na.rm = T),
                                  length.out = 50)
        # I.A.a.6 Hospitalization for heart failure
        newx.hba1c.hf.glp <- seq(min(dhf$diff.hba1c[dhf$type == "GLP1"& is.na(dhf$loghr) == F], na.rm = T),
                                     max(dhf$diff.hba1c[dhf$type == "GLP1"& is.na(dhf$loghr) == F], na.rm = T),
                                     length.out = 100)
        newx.hba1c.hf.sglt <- seq(min(dhf$diff.hba1c[dhf$type == "SGLT2"& is.na(dhf$loghr) == F], na.rm = T),
                                      max(dhf$diff.hba1c[dhf$type == "SGLT2"& is.na(dhf$loghr) == F], na.rm = T),
                                      length.out = 50)
        # I.A.a.7 Composite renal outcome
        newx.hba1c.cro.glp <- seq(min(dstroke$diff.hba1c[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                                     max(dstroke$diff.hba1c[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                                     length.out = 50)
        newx.hba1c.cro.sglt <- seq(min(dcro$diff.hba1c[dcro$type == "SGLT2"& is.na(dcro$loghr) == F], na.rm = T),
                                      max(dcro$diff.hba1c[dcro$type == "SGLT2"& is.na(dcro$loghr) == F], na.rm = T),
                                      length.out = 50)
    
  # I.A.b Bodyweight
        # I.A.b.1 MACE
        newx.bw.mace.glp <- seq(min(dmace$diff.bodyweight[dmace$type == "GLP1"& is.na(dmace$loghr) == F], na.rm = T),
                                    max(dmace$diff.bodyweight[dmace$type == "GLP1"& is.na(dmace$loghr) == F], na.rm = T),
                                    length.out = 50)
        newx.bw.mace.sglt <- seq(min(dmace$diff.bodyweight[dmace$type == "SGLT2"& is.na(dmace$loghr) == F], na.rm = T),
                                    max(dmace$diff.bodyweight[dmace$type == "SGLT2"& is.na(dmace$loghr) == F], na.rm = T),
                                    length.out = 50)
        # I.A.b.1 MI
        newx.bw.mi.glp <- seq(min(dmi$diff.bodyweight[dmi$type == "GLP1"& is.na(dmi$loghr) == F], na.rm = T),
                                 max(dmi$diff.bodyweight[dmi$type == "GLP1"& is.na(dmi$loghr) == F], na.rm = T),
                                 length.out = 50)
        newx.bw.mi.sglt <- seq(min(dmi$diff.bodyweight[dmi$type == "SGLT2"& is.na(dmi$loghr) == F], na.rm = T),
                                  max(dmi$diff.bodyweight[dmi$type == "SGLT2"& is.na(dmi$loghr) == F], na.rm = T),
                                  length.out = 50)
        # I.A.b.1 stroke
        newx.bw.stroke.glp <- seq(min(dstroke$diff.bodyweight[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                                     max(dstroke$diff.bodyweight[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                                     length.out = 50)
        newx.bw.stroke.sglt <- seq(min(dstroke$diff.bodyweight[dstroke$type == "SGLT2"& is.na(dstroke$loghr) == F], na.rm = T),
                                      max(dstroke$diff.bodyweight[dstroke$type == "SGLT2"& is.na(dstroke$loghr) == F], na.rm = T),
                                      length.out = 50)
        # Secondary outcomes
        # I.A.b.4 all cause mortality
        newx.bw.acm.glp <- seq(min(dacm$diff.bodyweight[dacm$type == "GLP1" & is.na(dacm$loghr) == F], na.rm = T),
                                  max(dacm$diff.bodyweight[dacm$type == "GLP1"& is.na(dacm$loghr) == F], na.rm = T),
                                  length.out = 50)
        newx.bw.acm.sglt <- seq(min(dacm$diff.bodyweight[dacm$type == "SGLT2"& is.na(dacm$loghr) == F], na.rm = T),
                                   max(dacm$diff.bodyweight[dacm$type == "SGLT2"& is.na(dacm$loghr) == F], na.rm = T),
                                   length.out = 50)
        # I.A.b.5 CV mortality
        newx.bw.cvm.glp <- seq(min(dcvm$diff.bodyweight[dcvm$type == "GLP1"& is.na(dcvm$loghr) == F], na.rm = T),
                                  max(dcvm$diff.bodyweight[dcvm$type == "GLP1"& is.na(dcvm$loghr) == F], na.rm = T),
                                  length.out = 50)
        newx.bw.cvm.sglt <- seq(min(dcvm$diff.bodyweight[dcvm$type == "SGLT2"& is.na(dcvm$loghr) == F], na.rm = T),
                                   max(dcvm$diff.bodyweight[dcvm$type == "SGLT2"& is.na(dcvm$loghr) == F], na.rm = T),
                                   length.out = 50)
        # I.A.b.6 Hospitalization for heart failure
        newx.bw.hf.glp <- seq(min(dhf$diff.bodyweight[dhf$type == "GLP1"& is.na(dhf$loghr) == F], na.rm = T),
                                 max(dhf$diff.bodyweight[dhf$type == "GLP1"& is.na(dhf$loghr) == F], na.rm = T),
                                 length.out = 50)
        newx.bw.hf.sglt <- seq(min(dhf$diff.bodyweight[dhf$type == "SGLT2"& is.na(dhf$loghr) == F], na.rm = T),
                                  max(dhf$diff.bodyweight[dhf$type == "SGLT2"& is.na(dhf$loghr) == F], na.rm = T),
                                  length.out = 50)
        # I.A.b.7 Composite renal outcome
        newx.bw.cro.glp <- seq(min(dstroke$diff.bodyweight[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                                  max(dstroke$diff.bodyweight[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                                  length.out = 50)
        newx.bw.cro.sglt <- seq(min(dcro$diff.bodyweight[dcro$type == "SGLT2"& is.na(dcro$loghr) == F], na.rm = T),
                                   max(dcro$diff.bodyweight[dcro$type == "SGLT2"& is.na(dcro$loghr) == F], na.rm = T),
                                   length.out = 50)
        
    # I.A.c Systolic blood pressure
        # I.A.c.1 MACE
        newx.sbp.mace.glp <- seq(min(dmace$diff.sbp[dmace$type == "GLP1"& is.na(dmace$loghr) == F], na.rm = T),
                                max(dmace$diff.sbp[dmace$type == "GLP1"& is.na(dmace$loghr) == F], na.rm = T),
                                length.out = 50)
        newx.sbp.mace.sglt <- seq(min(dmace$diff.sbp[dmace$type == "SGLT2"& is.na(dmace$loghr) == F], na.rm = T),
                                 max(dmace$diff.sbp[dmace$type == "SGLT2"& is.na(dmace$loghr) == F], na.rm = T),
                                 length.out = 50)
        # I.A.c.2 MI
        newx.sbp.mi.glp <- seq(min(dmi$diff.sbp[dmi$type == "GLP1"& is.na(dmi$loghr) == F], na.rm = T),
                              max(dmi$diff.sbp[dmi$type == "GLP1"& is.na(dmi$loghr) == F], na.rm = T),
                              length.out = 50)
        newx.sbp.mi.sglt <- seq(min(dmi$diff.sbp[dmi$type == "SGLT2"& is.na(dmi$loghr) == F], na.rm = T),
                               max(dmi$diff.sbp[dmi$type == "SGLT2"& is.na(dmi$loghr) == F], na.rm = T),
                               length.out = 50)
        # I.A.c.3  stroke
        newx.sbp.stroke.glp <- seq(min(dstroke$diff.sbp[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                                  max(dstroke$diff.sbp[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                                  length.out = 50)
        newx.sbp.stroke.sglt <- seq(min(dstroke$diff.sbp[dstroke$type == "SGLT2"& is.na(dstroke$loghr) == F], na.rm = T),
                                   max(dstroke$diff.sbp[dstroke$type == "SGLT2"& is.na(dstroke$loghr) == F], na.rm = T),
                                   length.out = 50)
        
        # Secondary outcomes
        # I.A.c.4 all cause mortality
        newx.sbp.acm.glp <- seq(min(dacm$diff.sbp[dacm$type == "GLP1" & is.na(dacm$loghr) == F], na.rm = T),
                               max(dacm$diff.sbp[dacm$type == "GLP1"& is.na(dacm$loghr) == F], na.rm = T),
                               length.out = 50)
        newx.sbp.acm.sglt <- seq(min(dacm$diff.sbp[dacm$type == "SGLT2"& is.na(dacm$loghr) == F], na.rm = T),
                                max(dacm$diff.sbp[dacm$type == "SGLT2"& is.na(dacm$loghr) == F], na.rm = T),
                                length.out = 50)
        # I.A.c.5 CV mortality
        newx.sbp.cvm.glp <- seq(min(dcvm$diff.sbp[dcvm$type == "GLP1"& is.na(dcvm$loghr) == F], na.rm = T),
                               max(dcvm$diff.sbp[dcvm$type == "GLP1"& is.na(dcvm$loghr) == F], na.rm = T),
                               length.out = 50)
        newx.sbp.cvm.sglt <- seq(min(dcvm$diff.sbp[dcvm$type == "SGLT2"& is.na(dcvm$loghr) == F], na.rm = T),
                                max(dcvm$diff.sbp[dcvm$type == "SGLT2"& is.na(dcvm$loghr) == F], na.rm = T),
                                length.out = 50)
        # I.A.c.6 Hospitalization for heart failure
        newx.sbp.hf.glp <- seq(min(dhf$diff.sbp[dhf$type == "GLP1"& is.na(dhf$loghr) == F], na.rm = T),
                              max(dhf$diff.sbp[dhf$type == "GLP1"& is.na(dhf$loghr) == F], na.rm = T),
                              length.out = 50)
        newx.sbp.hf.sglt <- seq(min(dhf$diff.sbp[dhf$type == "SGLT2"& is.na(dhf$loghr) == F], na.rm = T),
                               max(dhf$diff.sbp[dhf$type == "SGLT2"& is.na(dhf$loghr) == F], na.rm = T),
                               length.out = 50)
        # I.A.c.7 Composite renal outcome
        newx.sbp.cro.glp <- seq(min(dstroke$diff.sbp[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                               max(dstroke$diff.sbp[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                               length.out = 50)
        newx.sbp.cro.sglt <- seq(min(dcro$diff.sbp[dcro$type == "SGLT2"& is.na(dcro$loghr) == F], na.rm = T),
                                max(dcro$diff.sbp[dcro$type == "SGLT2"& is.na(dcro$loghr) == F], na.rm = T),
                                length.out = 50)
        
# I.A.d Hematocrit      
        # I.A.d.1 MACE 
        newx.hema.mace.sglt <- seq(min(dmace$diff.hematocrit[dmace$type == "SGLT2"& is.na(dmace$loghr) == F], na.rm = T),
                                   max(dmace$diff.hematocrit[dmace$type == "SGLT2"& is.na(dmace$loghr) == F], na.rm = T),
                                    length.out = 50)
        # I.A.d.2 MI 
        newx.hema.mi.sglt <- seq(min(dmi$diff.hematocrit[dmi$type == "SGLT2"& is.na(dmi$loghr) == F], na.rm = T),
                                 max(dmi$diff.hematocrit[dmi$type == "SGLT2"& is.na(dmi$loghr) == F], na.rm = T),
                                  length.out = 50)
        # I.A.d.3 stroke 
        newx.hema.stroke.sglt <- seq(min(dstroke$diff.hematocrit[dstroke$type == "SGLT2"& is.na(dstroke$loghr) == F], na.rm = T),
                                     max(dstroke$diff.hematocrit[dstroke$type == "SGLT2"& is.na(dstroke$loghr) == F], na.rm = T),
                                      length.out = 50)
        # I.A.d.4 All cause mortality 
        newx.hema.acm.sglt <- seq(min(dacm$diff.hematocrit[dacm$type == "SGLT2"& is.na(dacm$loghr) == F], na.rm = T),
                                  max(dacm$diff.hematocrit[dacm$type == "SGLT2"& is.na(dacm$loghr) == F], na.rm = T),
                                 length.out = 50)
        # I.A.d.5 CV mortality 
        newx.hema.cvm.sglt <- seq(min(dcvm$diff.hematocrit[dcvm$type == "SGLT2"& is.na(dcvm$loghr) == F], na.rm = T),
                                  max(dcvm$diff.hematocrit[dcvm$type == "SGLT2"& is.na(dcvm$loghr) == F], na.rm = T),
                                 length.out = 50)
        # I.A.d.6 Heart failure 
        newx.hema.hf.sglt <- seq(min(dhf$diff.hematocrit[dhf$type == "SGLT2"& is.na(dhf$loghr) == F], na.rm = T),
                                 max(dhf$diff.hematocrit[dhf$type == "SGLT2"& is.na(dhf$loghr) == F], na.rm = T),
                                length.out = 50)
        # I.A.d.1 Comosite renal outcome 
        newx.hema.cro.sglt <- seq(min(dcro$diff.hematocrit[dcro$type == "SGLT2"& is.na(dcro$loghr) == F], na.rm = T),
                                  max(dcro$diff.hematocrit[dcro$type == "SGLT2"& is.na(dcro$loghr) == F], na.rm = T),
                                 length.out = 50)
  # I.A.e UACR  
        # I.A.e.1 MACE 
        newx.uacr.mace.glp <- seq(min(dmace$diff.uacr[dmace$type == "GLP1"& is.na(dmace$loghr) == F], na.rm = T),
                                  max(dmace$diff.uacr[dmace$type == "GLP1"& is.na(dmace$loghr) == F], na.rm = T),
                                   length.out = 50)
        # I.A.e.1 MI
        newx.uacr.mi.glp <- seq(min(dmi$diff.uacr[dmi$type == "GLP1"& is.na(dmi$loghr) == F], na.rm = T),
                                 max(dmi$diff.uacr[dmi$type == "GLP1"& is.na(dmi$loghr) == F], na.rm = T),
                                 length.out = 50)
        # I.A.e.1 stroke
        newx.uacr.stroke.glp <- seq(min(dstroke$diff.uacr[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                                    max(dstroke$diff.uacr[dstroke$type == "GLP1"& is.na(dstroke$loghr) == F], na.rm = T),
                                     length.out = 50)
        # I.A.e.4 All-cause mortality
        newx.uacr.acm.glp <- seq(min(dacm$diff.uacr[dacm$type == "GLP1" & is.na(dacm$loghr) == F], na.rm = T),
                                 max(dacm$diff.uacr[dacm$type == "GLP1"& is.na(dacm$loghr) == F], na.rm = T),
                                length.out = 50)
        # I.A.e.5 Cardiovascular mortality
        newx.uacr.cvm.glp <- seq(min(dcvm$diff.uacr[dcvm$type == "GLP1"& is.na(dcvm$loghr) == F], na.rm = T),
                                 max(dcvm$diff.uacr[dcvm$type == "GLP1"& is.na(dcvm$loghr) == F], na.rm = T),
                                length.out = 50)
        # I.A.e.6 Heart failure
        newx.uacr.hf.glp <- seq(min(dhf$diff.uacr[dhf$type == "GLP1"& is.na(dhf$loghr) == F], na.rm = T),
                                max(dhf$diff.uacr[dhf$type == "GLP1"& is.na(dhf$loghr) == F], na.rm = T),
                               length.out = 50)
        # I.A.e.7 Composite renal outcome
        newx.uacr.cro.glp <- seq(min(dcro$diff.uacr[dcro$type == "GLP1"& is.na(dcro$loghr) == F], na.rm = T),
                                 max(dcro$diff.uacr[dcro$type == "GLP1"& is.na(dcro$loghr) == F], na.rm = T),
                                length.out = 50)
        
  # Prediction given newX vectors
        #hba1c
          preds.hba1c.mace.glp        <- predict(l.hba1c$glp.MACE,  newmods = newx.hba1c.mace.glp)
          preds.hba1c.mace.sglt       <- predict(l.hba1c$sglt.MACE,  newmods = newx.hba1c.mace.sglt)
          preds.hba1c.mi.glp          <- predict(l.hba1c$glp.MI,  newmods = newx.hba1c.mi.glp)
          preds.hba1c.mi.sglt         <- predict(l.hba1c$sglt.MI,  newmods = newx.hba1c.mi.sglt)
          preds.hba1c.stroke.glp      <- predict(l.hba1c$glp.stroke,  newmods = newx.hba1c.stroke.glp)
          preds.hba1c.stroke.sglt     <- predict(l.hba1c$sglt.stroke,  newmods = newx.hba1c.stroke.sglt)
        
            preds.hba1c.acm.glp        <- predict(l.hba1c$glp.allcauseMort,  newmods = newx.hba1c.acm.glp)
            preds.hba1c.acm.sglt       <- predict(l.hba1c$sglt.allcauseMort,  newmods = newx.hba1c.acm.sglt)
            preds.hba1c.cvm.glp        <- predict(l.hba1c$glp.CVMort,  newmods = newx.hba1c.cvm.glp)
            preds.hba1c.cvm.sglt       <- predict(l.hba1c$sglt.CVMort,  newmods = newx.hba1c.cvm.sglt)
            preds.hba1c.hf.glp         <- predict(l.hba1c$glp.HospHeartFailure,  newmods = newx.hba1c.hf.glp)
            preds.hba1c.hf.sglt        <- predict(l.hba1c$sglt.HospHeartFailure,  newmods = newx.hba1c.hf.sglt)
            preds.hba1c.cro.glp        <- predict(l.hba1c$glp.sustGFRdecl,  newmods = newx.hba1c.cro.glp)
            preds.hba1c.cro.sglt       <- predict(l.hba1c$sglt.sustGFRdecl,  newmods = newx.hba1c.cro.sglt)
       # bodyweight
            preds.bw.mace.glp           <- predict(l.bw$glp.MACE,  newmods = newx.bw.mace.glp)
            preds.bw.mace.sglt          <- predict(l.bw$sglt.MACE,  newmods = newx.bw.mace.sglt)
            preds.bw.mi.glp             <- predict(l.bw$glp.MI,  newmods = newx.bw.mi.glp)
            preds.bw.mi.sglt            <- predict(l.bw$sglt.MI,  newmods = newx.bw.mi.sglt)
            preds.bw.stroke.glp         <- predict(l.bw$glp.stroke,  newmods = newx.bw.stroke.glp)
            preds.bw.stroke.sglt        <- predict(l.bw$sglt.stroke,  newmods = newx.bw.stroke.sglt)
          
            preds.bw.acm.glp        <- predict(l.bw$glp.allcauseMort,  newmods = newx.bw.acm.glp)
            preds.bw.acm.sglt       <- predict(l.bw$sglt.allcauseMort,  newmods = newx.bw.acm.sglt)
            preds.bw.cvm.glp        <- predict(l.bw$glp.CVMort,  newmods = newx.bw.cvm.glp)
            preds.bw.cvm.sglt       <- predict(l.bw$sglt.CVMort,  newmods = newx.bw.cvm.sglt)
            preds.bw.hf.glp         <- predict(l.bw$glp.HospHeartFailure,  newmods = newx.bw.hf.glp)
            preds.bw.hf.sglt        <- predict(l.bw$sglt.HospHeartFailure,  newmods = newx.bw.hf.sglt)
            preds.bw.cro.glp        <- predict(l.bw$glp.sustGFRdecl,  newmods = newx.bw.cro.glp)
            preds.bw.cro.sglt       <- predict(l.bw$sglt.sustGFRdecl,  newmods = newx.bw.cro.sglt)
        # systolic blood pressure
            preds.sbp.mace.glp           <- predict(l.sbp$glp.MACE,  newmods = newx.sbp.mace.glp)
            preds.sbp.mace.sglt          <- predict(l.sbp$sglt.MACE,  newmods = newx.sbp.mace.sglt)
            preds.sbp.mi.glp             <- predict(l.sbp$glp.MI,  newmods = newx.sbp.mi.glp)
            preds.sbp.mi.sglt            <- predict(l.sbp$sglt.MI,  newmods = newx.sbp.mi.sglt)
            preds.sbp.stroke.glp         <- predict(l.sbp$glp.stroke,  newmods = newx.sbp.stroke.glp)
            preds.sbp.stroke.sglt        <- predict(l.sbp$sglt.stroke,  newmods = newx.sbp.stroke.sglt)
            
            
            preds.sbp.acm.glp        <- predict(l.sbp$glp.allcauseMort,  newmods = newx.sbp.acm.glp)
            preds.sbp.acm.sglt       <- predict(l.sbp$sglt.allcauseMort,  newmods = newx.sbp.acm.sglt)
            preds.sbp.cvm.glp        <- predict(l.sbp$glp.CVMort,  newmods = newx.sbp.cvm.glp)
            preds.sbp.cvm.sglt       <- predict(l.sbp$sglt.CVMort,  newmods = newx.sbp.cvm.sglt)
            preds.sbp.hf.glp         <- predict(l.sbp$glp.HospHeartFailure,  newmods = newx.sbp.hf.glp)
            preds.sbp.hf.sglt        <- predict(l.sbp$sglt.HospHeartFailure,  newmods = newx.sbp.hf.sglt)
            preds.sbp.cro.glp        <- predict(l.sbp$glp.sustGFRdecl,  newmods = newx.sbp.cro.glp)
            preds.sbp.cro.sglt       <- predict(l.sbp$sglt.sustGFRdecl,  newmods = newx.sbp.cro.sglt)
      # hematocrit
            preds.hema.mace.sglt          <- predict(l.hema$sglt.MACE,  newmods = newx.hema.mace.sglt)
            preds.hema.mi.sglt            <- predict(l.hema$sglt.MI,  newmods = newx.hema.mi.sglt)
            preds.hema.stroke.sglt        <- predict(l.hema$sglt.stroke,  newmods = newx.hema.stroke.sglt)
            preds.hema.acm.sglt            <- predict(l.hema$sglt.allcauseMort,  newmods = newx.hema.mi.sglt)
            preds.hema.cvm.sglt            <- predict(l.hema$sglt.CVMort,  newmods = newx.hema.mi.sglt)
            preds.hema.hf.sglt            <- predict(l.hema$sglt.HospHeartFailure,  newmods = newx.hema.mi.sglt)
            preds.hema.cro.sglt            <- predict(l.hema$sglt.sustGFRdecl,  newmods = newx.hema.mi.sglt)
          
      # uacr
            preds.uacr.mace.glp          <- predict(l.uacr$glp.MACE,  newmods = newx.uacr.mace.glp)
            preds.uacr.mi.glp            <- predict(l.uacr$glp.MI,  newmods = newx.uacr.mi.glp)
            preds.uacr.stroke.glp        <- predict(l.uacr$glp.stroke,  newmods = newx.uacr.stroke.glp)
            
            preds.uacr.acm.glp            <- predict(l.uacr$glp.allcauseMort,  newmods = newx.uacr.acm.glp)
            preds.uacr.cvm.glp            <- predict(l.uacr$glp.CVMort,  newmods = newx.uacr.cvm.glp)
            preds.uacr.hf.glp            <- predict(l.uacr$glp.HospHeartFailure,  newmods = newx.uacr.hf.glp)
            preds.uacr.cro.glp            <- predict(l.uacr$glp.sustGFRdecl,  newmods = newx.uacr.cro.glp)
            
            
 # Store x and ys in lists (may be useful for plots)   
        l.newx <- list()
        l.preds <- list()
            l.preds[[1]]   <- preds.hba1c.mace.glp
            l.preds[[2]]   <- preds.hba1c.mace.sglt
            l.preds[[3]]   <- preds.hba1c.mi.glp
            l.preds[[4]]   <- preds.hba1c.mi.sglt
            l.preds[[5]]   <- preds.hba1c.stroke.glp
            l.preds[[6]]   <- preds.hba1c.stroke.sglt
            l.preds[[7]]   <- preds.bw.mace.glp
            l.preds[[8]]   <- preds.bw.mace.sglt
            l.preds[[9]]   <- preds.bw.mi.glp
            l.preds[[10]]  <- preds.bw.mi.sglt
            l.preds[[11]]  <- preds.bw.stroke.glp
            l.preds[[12]]  <- preds.bw.stroke.sglt
            l.preds[[13]]  <- preds.sbp.mace.glp
            l.preds[[14]]  <- preds.sbp.mace.sglt
            l.preds[[15]]  <- preds.sbp.mi.glp
            l.preds[[16]]  <- preds.sbp.mi.sglt
            l.preds[[17]]  <- preds.sbp.stroke.glp
            l.preds[[18]]  <- preds.sbp.stroke.sglt
            l.preds[[19]]  <- preds.hema.mace.sglt
            l.preds[[20]]  <- preds.hema.mi.sglt
            l.preds[[21]]  <- preds.hema.stroke.sglt
            l.preds[[22]]  <- preds.uacr.mace.glp
            l.preds[[23]]  <- preds.uacr.mi.glp
            l.preds[[24]]  <- preds.uacr.stroke.glp
            
            
            l.newx[[1]]   <- newx.hba1c.mace.glp
            l.newx[[2]]   <- newx.hba1c.mace.sglt
            l.newx[[3]]   <- newx.hba1c.mi.glp
            l.newx[[4]]   <- newx.hba1c.mi.sglt
            l.newx[[5]]   <- newx.hba1c.stroke.glp
            l.newx[[6]]   <- newx.hba1c.stroke.sglt
            l.newx[[7]]   <- newx.bw.mace.glp
            l.newx[[8]]   <- newx.bw.mace.sglt
            l.newx[[9]]   <- newx.bw.mi.glp
            l.newx[[10]]  <- newx.bw.mi.sglt
            l.newx[[11]]  <- newx.bw.stroke.glp
            l.newx[[12]]  <- newx.bw.stroke.sglt
            l.newx[[13]]  <- newx.sbp.mace.glp
            l.newx[[14]]  <- newx.sbp.mace.sglt
            l.newx[[15]]  <- newx.sbp.mi.glp
            l.newx[[16]]  <- newx.sbp.mi.sglt
            l.newx[[17]]  <- newx.sbp.stroke.glp
            l.newx[[18]]  <- newx.sbp.stroke.sglt
            l.newx[[19]]  <- newx.hema.mace.sglt
            l.newx[[20]]  <- newx.hema.mi.sglt
            l.newx[[21]]  <- newx.hema.stroke.sglt
            l.newx[[22]]  <- newx.uacr.mace.glp
            l.newx[[23]]  <- newx.uacr.mi.glp
            l.newx[[24]]  <- newx.uacr.stroke.glp

        
  # Domains for axes
  xlim.hba1c <- c(-1.4,0)
  xlim.bw <- c(-5.0,0)
  xlims.bw <- seq(-5.0, 0, 1)
  xlim.sbp       <- c(-5,1)
  xlims.sbp       <- seq(-5,1,1)
  xlim.hema <-  seq(2.0,3.2, .3)
  xlim.uacr     <- range(df$diff.uacr[df$type == "GLP1"], na.rm = TRUE)
  
  m <- matrix(NA, nrow = 21, ncol = 2)
  for(i in 1:21){
  m[i,1] <- min(l.preds[[i]]$ci.lb)  
  m[i,2] <- max(l.preds[[i]]$ci.ub) 
  }
  xlim.uacr <- c(-10,-1)
  xlims.uacr <- seq(-10,2,2)
  
  yseq <- c(0.4, 0.5, 0.6, 0.7, 0.8,0.9, 1, 1.1,1.2,1.3, 1.4 )
  yseq <- c(0.5, 0.6, 0.8, 1,1.2, 1.4)
  yseq <- c(0.4, 0.5, 0.6, 0.8, 1, 1.2,1.3, 1.4 )
  
  ylim.mace     <- range(log(yseq))
  ylim.stroke   <- range(log(yseq))
  ylim.mi       <- range(log(yseq))
  ylim.hema     <- range(log(yseq))
  
  yaxis.mace <- log(yseq)
  yaxis.mace.lab <- yseq
  yaxis.mi <- log(yseq)
  yaxis.mi.lab <-yseq
  yaxis.stroke <- log(yseq)
  yaxis.stroke.lab <- yseq
  yaxis.hema <- log(yseq)
  yaxis.hema.lab <- yseq
  
  
  
  yseq <- c(0.4, 0.6, 0.7, 0.8,0.9, 1, 1.1,1.2,1.3, 1.4 )
  yseq <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.7)
  
  yaxis   <- (log(yseq))
  ylim   <- range((log(yseq)))
  ylim.acm   <- range(log(yseq))
  ylim.cvm  <- range(log(yseq))
  ylim.cro   <- range(log(yseq))
  ylim.hf   <- range(log(yseq))
  ylim.mace   <- range(log(yseq))
  ylim.stroke <- range(log(yseq))
  ylim.mi     <- range(log(yseq))
  ylim.stroke <- range(log(yseq))
  ylim.mi     <- range(log(yseq))
  
  yaxis.mace <- log(yseq)
  yaxis.mace.lab <- yseq
  yaxis.mi <- log(yseq)
  yaxis.mi.lab <- yseq
  yaxis.stroke <- log(yseq)
  yaxis.stroke.lab <- yseq
  
  yaxis.stroke <- log(yseq)
  yaxis.stroke.lab <- yseq
  
  yaxis.acm <- log(yseq)
  yaxis.acm.lab <- yseq
  yaxis.cvm <- log(yseq)
  yaxis.cvm.lab <- yseq
  yaxis.hf <- log(yseq)
  yaxis.hf.lab <- yseq
  ylim.cro <- log(c(0.4,1.7))
  ycroseq <- c(0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.7)
  yaxis.cro.lab <- ycroseq
  yaxis.cro <- log(ycroseq)
  

############### Figure A and B with all risk factors =========
# Figure 1: 'primary' outcomes: MACE, MI, stroke, CVD (4x5)
  # 1.x Hba1C, 2.x SBP, 3.x BW 4.x Hematocrit 5.x uACR
  # y.1 MACE y.2 MI y.3 stroke y.4 CVD
# Figure 2: 'secondary' outcomes: HHF, CRO, ACM (3x5)
  # 1.x Hba1C, 2.x SBP, 3.x BW 4.x Hematocrit 5.x uACR
  # y.1 HHF y.2 CRO y.3 ACM

png("~/Documents/_mediation/_output/_figures/figure_mediator_figureA.png", width = 12, height = 9.6, units = 'in', res = 300)  
# Figure A # ======
      par(mfcol = c(4,5), oma = c(1,1,1,1)*1.3)  
      par(mar = c(1,2,0.5,0)*2)
# Col 1: Hba1C          #             #
      # 1.1 HBA1c MACE
        plot(x = 0, y = 0, type = 'n', xlim = xlim.hba1c, ylim = ylim.mace,
             axes = F, xlab = "", ylab = "")
        title(main = expression(paste("HbA1c")), line = 0.5)
        #CI, points, predictions
        polygon(x = c(newx.hba1c.mace.glp, rev(newx.hba1c.mace.glp)), 
                y = c(preds.hba1c.mace.glp$ci.lb, rev(preds.hba1c.mace.glp$ci.ub)),
                col = cola[1], border = "transparent")
        polygon(x = c(newx.hba1c.mace.sglt, rev(newx.hba1c.mace.sglt)), 
                y = c(preds.hba1c.mace.sglt$ci.lb, rev(preds.hba1c.mace.sglt$ci.ub)),
                col = cola[2], border = "transparent")
        points(x = dmace$diff.hba1c, y = dmace$loghr, col = "black", bg = dmace$col, 
               cex = (dmace$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
        lines(x = newx.hba1c.mace.glp, y = preds.hba1c.mace.glp$pred, 
              col = "black", lwd = 1.2)
        lines(x = newx.hba1c.mace.sglt, y = preds.hba1c.mace.sglt$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.mace, labels = yaxis.mace.lab, las = 2)
        axis(1, las = 1, labels = F)
        title(ylab = "Hazard ratio", line = 2.7, cex = 1.2)
      
      # 2.1 HBA1c MI
        plot(x = 0, y = 0, type = 'n', xlim = xlim.hba1c, ylim = ylim.mi,
             axes = F, xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.hba1c.mi.glp, rev(newx.hba1c.mi.glp)), 
                y = c(preds.hba1c.mi.glp$ci.lb, rev(preds.hba1c.mi.glp$ci.ub)),
                col = cola[1], border = "transparent")
        polygon(x = c(newx.hba1c.mi.sglt, rev(newx.hba1c.mi.sglt)), 
                y = c(preds.hba1c.mi.sglt$ci.lb, rev(preds.hba1c.mi.sglt$ci.ub)),
                col = cola[2], border = "transparent")
        points(x = dmi$diff.hba1c, y = dmi$loghr, col = "black", bg = dmi$col, 
               cex = (dmi$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
        lines(x = newx.hba1c.mi.glp, y = preds.hba1c.mi.glp$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        lines(x = newx.hba1c.mi.sglt, y = preds.hba1c.mi.sglt$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.mi, labels = yaxis.mi.lab, las = 2)
        axis(1, las = 1, labels = F) 
        title(ylab = "Hazard ratio", line = 2.7, cex = 1.2)
        
      # 3.1 HBA1c stroke
        plot(x = 0, y = 0, type = 'n', xlim = xlim.hba1c, ylim = ylim.stroke,
             axes = F, xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.hba1c.stroke.glp, rev(newx.hba1c.stroke.glp)), 
                y = c(preds.hba1c.stroke.glp$ci.lb, rev(preds.hba1c.stroke.glp$ci.ub)),
                col = cola[1], border = "transparent")
        polygon(x = c(newx.hba1c.stroke.sglt, rev(newx.hba1c.stroke.sglt)), 
                y = c(preds.hba1c.stroke.sglt$ci.lb, rev(preds.hba1c.stroke.sglt$ci.ub)),
                col = cola[2], border = "transparent")
        points(x = dstroke$diff.hba1c, y = dstroke$loghr, col = "black", bg = dstroke$col, 
               cex = (dstroke$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
        lines(x = newx.hba1c.stroke.glp, y = preds.hba1c.stroke.glp$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        lines(x = newx.hba1c.stroke.sglt, y = preds.hba1c.stroke.sglt$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.stroke, labels = yaxis.stroke.lab, las = 2)
        axis(1, las = 1, labels = F) 
        title(ylab = "Hazard ratio", line = 2.7, cex = 1.2)
        
      # 3.1 HBA1c CVM
        par(mar = c(2,2,0.5,0)*2)
        plot(x = 0, y = 0, type = 'n', xlim = xlim.hba1c, ylim = ylim.cvm,
             axes = F, xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.hba1c.cvm.glp, rev(newx.hba1c.cvm.glp)), 
                y = c(preds.hba1c.cvm.glp$ci.lb, rev(preds.hba1c.cvm.glp$ci.ub)),
                col = cola[1], border = "transparent")
        polygon(x = c(newx.hba1c.cvm.sglt, rev(newx.hba1c.cvm.sglt)), 
                y = c(preds.hba1c.cvm.sglt$ci.lb, rev(preds.hba1c.cvm.sglt$ci.ub)),
                col = cola[2], border = "transparent")
        points(x = dcvm$diff.hba1c, y = dcvm$loghr, col = "black", bg = dcvm$col, 
               cex = (dcvm$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
        lines(x = newx.hba1c.cvm.glp, y = preds.hba1c.cvm.glp$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        lines(x = newx.hba1c.cvm.sglt, y = preds.hba1c.cvm.sglt$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = yaxis.cvm.lab, las = 2)
        axis(1, las = 1, labels = T) 
        title(ylab = "Hazard ratio", line = 2.7, cex = 1.2)
        title(xlab = expression(paste(Delta, "(treatment - control), %")), line = 2.5)
        
#### Bodyweight        
        par(mar = c(1,2,0.5,0)*2)
      # 1.2 Bodyweight MACE
        plot(x = 0, y = 0, type = 'n', xlim = xlim.bw, ylim = ylim.mace,
             axes = F, xlab = "", ylab = "")
        title(main = expression(paste("Bodyweight")), line = 0.5)
        #CI, points, predictions
        polygon(x = c(newx.bw.mace.glp, rev(newx.bw.mace.glp)), 
                y = c(preds.bw.mace.glp$ci.lb, rev(preds.bw.mace.glp$ci.ub)),
                col = cola[1], border = "transparent")
        polygon(x = c(newx.bw.mace.sglt, rev(newx.bw.mace.sglt)), 
                y = c(preds.bw.mace.sglt$ci.lb, rev(preds.bw.mace.sglt$ci.ub)),
                col = cola[2], border = "transparent")
        points(x = dmace$diff.bodyweight, y = dmace$loghr, col = "black", bg = dmace$col, 
               cex = (dmace$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
        lines(x = newx.bw.mace.glp, y = preds.bw.mace.glp$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        lines(x = newx.bw.mace.sglt, y = preds.bw.mace.sglt$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, labels = F)

     # 2.2 Bodyweight MI
        plot(x = 0, y = 0, type = 'n', xlim = xlim.bw, ylim = ylim.mi,
             axes = F, xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.bw.mi.glp, rev(newx.bw.mi.glp)), 
                y = c(preds.bw.mi.glp$ci.lb, rev(preds.bw.mi.glp$ci.ub)),
                col = cola[1], border = "transparent")
        polygon(x = c(newx.bw.mi.sglt, rev(newx.bw.mi.sglt)), 
                y = c(preds.bw.mi.sglt$ci.lb, rev(preds.bw.mi.sglt$ci.ub)),
                col = cola[2], border = "transparent")
        points(x = dmi$diff.bodyweight, y = dmi$loghr, col = "black", bg = dmi$col, 
               cex = (dmi$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
        lines(x = newx.bw.mi.glp, y = preds.bw.mi.glp$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        lines(x = newx.bw.mi.sglt, y = preds.bw.mi.sglt$pred, 
              col = "black", lwd = 1.2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, labels = F) 

      # 3.2 Bodyweight stroke
        plot(x = 0, y = 0, type = 'n', xlim = xlim.bw, ylim = ylim.stroke,
             axes = F, xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.bw.stroke.glp, rev(newx.bw.stroke.glp)), 
                y = c(preds.bw.stroke.glp$ci.lb, rev(preds.bw.stroke.glp$ci.ub)),
                col = cola[1], border = "transparent")
        polygon(x = c(newx.bw.stroke.sglt, rev(newx.bw.stroke.sglt)), 
                y = c(preds.bw.stroke.sglt$ci.lb, rev(preds.bw.stroke.sglt$ci.ub)),
                col = cola[2], border = "transparent")
        points(x = dstroke$diff.bodyweight, y = dstroke$loghr, col = "black", bg = dstroke$col, 
               cex = (dstroke$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
        lines(x = newx.bw.stroke.glp, y = preds.bw.stroke.glp$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        lines(x = newx.bw.stroke.sglt, y = preds.bw.stroke.sglt$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, labels = F) 

      #  4.2 Bodyweight CVM
        par(mar = c(2,2,0.5,0)*2)
        plot(x = 0, y = 0, type = 'n', xlim = xlim.bw, ylim = ylim.cvm,
             axes = F, xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.bw.cvm.glp, rev(newx.bw.cvm.glp)), 
                y = c(preds.bw.cvm.glp$ci.lb, rev(preds.bw.cvm.glp$ci.ub)),
                col = cola[1], border = "transparent")
        polygon(x = c(newx.bw.cvm.sglt, rev(newx.bw.cvm.sglt)), 
                y = c(preds.bw.cvm.sglt$ci.lb, rev(preds.bw.cvm.sglt$ci.ub)),
                col = cola[2], border = "transparent")
        points(x = dcvm$diff.bodyweight, y = dcvm$loghr, col = "black", bg = dcvm$col, 
               cex = (dcvm$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
        lines(x = newx.bw.cvm.glp, y = preds.bw.cvm.glp$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        lines(x = newx.bw.cvm.sglt, y = preds.bw.cvm.sglt$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, labels = T) 
        title(xlab = expression(paste(Delta, "(treatment - control), kg.")), line = 2.5)
        
        
#### Systolic blood pressure        
        par(mar = c(1,2,0.5,0)*2)
        # 1.2 Bodyweight MACE
        plot(x = 0, y = 0, type = 'n', xlim = xlim.sbp, ylim = ylim.mace,
             axes = F, xlab = "", ylab = "")
        title(main = expression(paste("Systolic blood pressure")), line = 0.5)
        #CI, points, predictions
        polygon(x = c(newx.sbp.mace.glp, rev(newx.sbp.mace.glp)), 
                y = c(preds.sbp.mace.glp$ci.lb, rev(preds.sbp.mace.glp$ci.ub)),
                col = cola[1], border = "transparent")
        polygon(x = c(newx.sbp.mace.sglt, rev(newx.sbp.mace.sglt)), 
                y = c(preds.sbp.mace.sglt$ci.lb, rev(preds.sbp.mace.sglt$ci.ub)),
                col = cola[2], border = "transparent")
        points(x = dmace$diff.bodyweight, y = dmace$loghr, col = "black", bg = dmace$col, 
               cex = (dmace$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
        lines(x = newx.sbp.mace.glp, y = preds.sbp.mace.glp$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        lines(x = newx.sbp.mace.sglt, y = preds.sbp.mace.sglt$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, labels = F)
        
        # 2.2 Bodyweight MI
        plot(x = 0, y = 0, type = 'n', xlim = xlim.sbp, ylim = ylim.mi,
             axes = F, xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.sbp.mi.glp, rev(newx.sbp.mi.glp)), 
                y = c(preds.sbp.mi.glp$ci.lb, rev(preds.sbp.mi.glp$ci.ub)),
                col = cola[1], border = "transparent")
        polygon(x = c(newx.sbp.mi.sglt, rev(newx.sbp.mi.sglt)), 
                y = c(preds.sbp.mi.sglt$ci.lb, rev(preds.sbp.mi.sglt$ci.ub)),
                col = cola[2], border = "transparent")
        points(x = dmi$diff.bodyweight, y = dmi$loghr, col = "black", bg = dmi$col, 
               cex = (dmi$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
        lines(x = newx.sbp.mi.glp, y = preds.sbp.mi.glp$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        lines(x = newx.sbp.mi.sglt, y = preds.sbp.mi.sglt$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, labels = F) 
        
        # 3.2 Bodyweight stroke
        plot(x = 0, y = 0, type = 'n', xlim = xlim.sbp, ylim = ylim.stroke,
             axes = F, xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.sbp.stroke.glp, rev(newx.sbp.stroke.glp)), 
                y = c(preds.sbp.stroke.glp$ci.lb, rev(preds.sbp.stroke.glp$ci.ub)),
                col = cola[1], border = "transparent")
        polygon(x = c(newx.sbp.stroke.sglt, rev(newx.sbp.stroke.sglt)), 
                y = c(preds.sbp.stroke.sglt$ci.lb, rev(preds.sbp.stroke.sglt$ci.ub)),
                col = cola[2], border = "transparent")
        points(x = dstroke$diff.sbp, y = dstroke$loghr, col = "black", bg = dstroke$col, 
               cex = (dstroke$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
        lines(x = newx.sbp.stroke.glp, y = preds.sbp.stroke.glp$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        lines(x = newx.sbp.stroke.sglt, y = preds.sbp.stroke.sglt$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, labels = F) 
        
        #  4.2 Bodyweight CVM
        par(mar = c(2,2,0.5,0)*2)
        plot(x = 0, y = 0, type = 'n', xlim = xlim.sbp, ylim = ylim.cvm,
             axes = F, xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.sbp.cvm.glp, rev(newx.sbp.cvm.glp)), 
                y = c(preds.sbp.cvm.glp$ci.lb, rev(preds.sbp.cvm.glp$ci.ub)),
                col = cola[1], border = "transparent")
        polygon(x = c(newx.sbp.cvm.sglt, rev(newx.sbp.cvm.sglt)), 
                y = c(preds.sbp.cvm.sglt$ci.lb, rev(preds.sbp.cvm.sglt$ci.ub)),
                col = cola[2], border = "transparent")
        points(x = dcvm$diff.sbp, y = dcvm$loghr, col = "black", bg = dcvm$col, 
               cex = (dcvm$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
        lines(x = newx.sbp.cvm.glp, y = preds.sbp.cvm.glp$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        lines(x = newx.sbp.cvm.sglt, y = preds.sbp.cvm.sglt$pred, 
              col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, labels = T) 
        title(xlab = expression(paste(Delta, "(treatment - control), mmHg")), line = 2.5)
      
#### Hematocrit
    # Plot 4.1 hematocrit MACE
        par(mar = c(1,2,0.5,0)*2)
        plot(x = 0, y = 0, type = 'n', xlim = range(xlim.hema), ylim = ylim.mace ,
             axes = F, xlab = "", ylab = "")
        title(main = expression(paste("Hematocrit")), line = 0.5)
        #CI, points, predictions
        polygon(x = c(newx.hema.mace.sglt, rev(newx.hema.mace.sglt)), 
                y = c(preds.hema.mace.sglt$ci.lb, rev(preds.hema.mace.sglt$ci.ub)),
                col = cola[2],border = "transparent")
        points(x = dmace$diff.hematocrit, y = dmace$loghr, col = "black", bg = dmace$col,
               cex = (dmace$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
        lines(newx.hema.mace.sglt, preds.hema.mace.sglt$pred, col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, labels = F, at = xlim.hema) 
        
    # Plot 4.2 hematocrit MI
        plot(x = 0, y = 0, type = 'n', xlim = range(xlim.hema), ylim = ylim.mi ,
             axes = F, xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.hema.mi.sglt, rev(newx.hema.mi.sglt)), 
                y = c(preds.hema.mi.sglt$ci.lb, rev(preds.hema.mi.sglt$ci.ub)),
                col = cola[2],border = "transparent")
        points(x = dmi$diff.hematocrit, y = dmi$loghr, col = "black", bg = dmi$col,
               pch = 21 , lwd = 0.4 )
        lines(newx.hema.mi.sglt, preds.hema.mi.sglt$pred, col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, labels = F, at = xlim.hema) 
        
    # Plot 4.3 hematocrit stroke
        plot(x = 0, y = 0, type = 'n', xlim = range(xlim.hema), ylim = ylim.stroke ,
             axes = F,xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.hema.stroke.sglt, rev(newx.hema.stroke.sglt)), 
                y = c(preds.hema.stroke.sglt$ci.lb, rev(preds.hema.stroke.sglt$ci.ub)),
                col = cola[2],
                border = "transparent")
        points(x = dstroke$diff.hematocrit, y = dstroke$loghr, 
               col = "black", bg = dmi$col,
               cex = (dstroke$logvi*500)^(0.5),
               pch = 21 , lwd = 0.4 )
        lines(newx.hema.stroke.sglt, preds.hema.stroke.sglt$pred, col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, labels = F, at = xlim.hema) 
        # Plot 4.4 hematocrit CVM
        par(mar = c(2,2,0.5,0)*2)
        plot(x = 0, y = 0, type = 'n', xlim = range(xlim.hema), ylim = ylim.cvm ,
             axes = F, xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.hema.cvm.sglt, rev(newx.hema.cvm.sglt)), 
                y = c(preds.hema.cvm.sglt$ci.lb, rev(preds.hema.cvm.sglt$ci.ub)),
                col = cola[2],
                border = "transparent")
        points(x = dcvm$diff.hematocrit[dcvm$type == "SGLT2"], y = dcvm$loghr[dcvm$type == "SGLT2"], 
               col = "black", bg = dcvm$col[dcvm$type == "SGLT2"],
               cex = (dcvm$logvi[dcvm$type == "SGLT2"]*500)^(0.5),
               pch = 21 , lwd = 0.4 )
        lines(x = newx.hema.cvm.sglt, y = preds.hema.cvm.sglt$pred, col = "gray40", lwd = 0.85, lty = 2)
        # reference line
        abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, labels = T, at = xlim.hema) 
        title(xlab = expression(paste(Delta, "(treatment - control), %")), line = 2.5)
        
        
    # Plot 5.1 UACR MACE
        par(mar = c(1,2,0.5,0)*2)
        plot(x = 0, y = 0, type = 'n', xlim = range(xlim.uacr)+c(0,3) , ylim = ylim.mace ,
             axes = F, xlab = "", ylab = "")
        legend("topright", 
               pch =c(19,19,NA), lwd =c(NA,NA,2), 
               col = c("gold","gray70", "black"), c("GLP-1RA", "SGLT2i", "P<0.05"), 
               cex = 0.75, box.col = "transparent", ncol = 1)
        #CI, points, predictions
        polygon(x = c(newx.uacr.mace.glp, rev(newx.uacr.mace.glp)), 
                y = c(preds.uacr.mace.glp$ci.lb, rev(preds.uacr.mace.glp$ci.ub)),
                col = cola[1],
                border = "transparent")
        points(x = dmace$diff.uacr[dmace$type == "GLP1"], y = dmace$loghr[dmace$type == "GLP1"], 
               col = "black", bg = dmace$col[dmace$type == "GLP1"],
               cex = (dmace$logvi[dmace$type == "GLP1"]*500)^(0.5),
               pch = 21 , lwd = 0.4 )
        lines(newx.uacr.mace.glp, preds.uacr.mace.glp$pred, col = "gray40", lwd = 0.8, lty = 3)
        title(main = expression(paste("uACR")), line = 0.5)
        # reference line
        segments(x0 = -10, x1 = 2, y0 = log(1), col = "gray20", lwd = 1, lty = 3)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, at = c(xlims.uacr), labels = F)
        text( x = 1.9, y = -.25, labels = "MACE", cex = 1.2, srt = 270)
        text( x = -0.6, y = log(1.2), labels = "Favors \n control", cex = 0.9, srt = 0)
        text( x = -0.6, y = log(0.65), labels = "Favors \n treatment", cex = 0.9, srt = 0)
        
        
    # Plot 5.2 UACR MI
  
        plot(x = 0, y = 0, type = 'n', xlim = range(xlim.uacr)+c(0,3), ylim = ylim.mi ,
             axes = F,
             xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.uacr.mi.glp, rev(newx.uacr.mi.glp)), 
                y = c(preds.uacr.mi.glp$ci.lb, rev(preds.uacr.mi.glp$ci.ub)),
                col = cola[1],
                border = "transparent")
        points(x = dmi$diff.uacr[dmi$type == "GLP1"], y = dmi$loghr[dmi$type == "GLP1"], 
               col = "black", bg = dmi$col[dmi$type == "GLP1"],
               cex = (dmi$logvi[dmi$type == "GLP1"]*500)^(0.5),
               pch = 21 , lwd = 0.4 )
        lines(newx.uacr.mi.glp, preds.uacr.mi.glp$pred, col = "gray40", lwd = 0.8, lty = 3)
        # reference line
        segments(x0 = -10, x1 = 2, y0 = log(1), col = "gray20", lwd = 1, lty = 3)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, at = xlims.uacr, labels = F)
        text( x = 1.9, y = -.25, labels = "Myocardial infarction", cex = 1.2, srt = 270)
       text( x = -0.6, y = log(1.2), labels = "Favors \n control", cex = 0.9, srt = 0)
       text( x = -0.6, y = log(0.65), labels = "Favors \n treatment", cex = 0.9, srt = 0)
        
    # Plot 5.3 uacr stroke

            plot(x = 0, y = 0, type = 'n', xlim = range(xlim.uacr)+c(0,3), ylim = ylim.stroke ,
                 axes = F,
                 xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.uacr.stroke.glp, rev(newx.uacr.stroke.glp)), 
                y = c(preds.uacr.stroke.glp$ci.lb, rev(preds.uacr.stroke.glp$ci.ub)),
                col = cola[1],
                border = "transparent")
        points(x = dstroke$diff.uacr[dstroke$type == "GLP1"], y = dstroke$loghr[dstroke$type == "GLP1"], 
               col = "black", bg = dstroke$col[dstroke$type == "GLP1"],
               cex = (dstroke$logvi[dstroke$type == "GLP1"]*500)^(0.5),
               pch = 21 , lwd = 0.4 )
        lines(newx.uacr.stroke.glp, preds.uacr.stroke.glp$pred, col = "gray40", lwd = 0.8, lty = 3)
        # reference line
        segments(x0 = -10, x1 = 2, y0 = log(1), col = "gray20", lwd = 1, lty = 3)
        # axes
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, at = xlims.uacr, labels = F)
        text( x = 1.9, y = 0, labels = "Stroke", cex = 1.2, srt = 270)
        text( x = -0.6, y = log(1.2), labels = "Favors \n control", cex = 0.9, srt = 0)
        text( x = -0.6, y = log(0.65), labels = "Favors \n treatment", cex = 0.9, srt = 0)
        
      # 5.4 UACR CVM
        par(mar = c(2,2,0.5,0)*2)
        plot(x = 0, y = 0, type = 'n', xlim = xlim.uacr+c(0,3), ylim = ylim.cvm ,
             axes = F,
             xlab = "", ylab = "")
        #CI, points, predictions
        polygon(x = c(newx.uacr.cvm.glp, rev(newx.uacr.cvm.glp)), 
                y = c(preds.uacr.cvm.glp$ci.lb, rev(preds.uacr.cvm.glp$ci.ub)),
                col = cola[1],
                border = "transparent")
        points(x = dcvm$diff.uacr[dcvm$type == "GLP1"], y = dcvm$loghr[dcvm$type == "GLP1"], 
               col = "black", bg = dcvm$col[dcvm$type == "GLP1"],
               cex = (dcvm$logvi[dcvm$type == "GLP1"]*500)^(0.5),
               pch = 21 , lwd = 0.4 )
        lines(newx.uacr.cvm.glp, preds.uacr.cvm.glp$pred,  col = "gray40", lwd = 0.8, lty = 3)
        # reference line
        segments(x0 = -10, x1 = 2, y0 = log(1), col = "gray20", lwd = 1, lty = 3)
        # Axes properties
        axis(2, at = yaxis.cvm, labels = F)
        axis(1, las = 1, at = xlims.uacr, labels = T)
        title(xlab = expression(paste(Delta, "(treatment - control), mg/g")), line = 2.5)
        text( x = 1.9, y = -.25, labels = "Cardiovacular mortality", cex = 1.2, srt = 270)
       text( x = -0.6, y = log(1.2), labels = "Favors \n control", cex = 0.9, srt = 0)
       text( x = -0.6, y = log(0.65), labels = "Favors \n treatment", cex = 0.9, srt = 0)
        
        
dev.off()        
      

##############################





##### Figure B

png("~/Documents/_mediation/_output/_figures/figure_mediator_figureB.png", width = 12, height = 7.2, units = 'in', res = 300)  


par(mfcol = c(3,5), oma = c(1,1,1,4)*1.3)  
par(mar = c(1,2,0.5,0)*2)
# Col 1: Hba1C          #             #
      # 1.1 HBA1c HHF
          plot(x = 0, y = 0, type = 'n', xlim = xlim.hba1c, ylim = ylim.hf,
               axes = F, xlab = "", ylab = "")
          title(main = expression(paste("HbA1c")), line = 0.5)
          #CI, points, predictions
          polygon(x = c(newx.hba1c.hf.glp, rev(newx.hba1c.hf.glp)), 
                  y = c(preds.hba1c.hf.glp$ci.lb, rev(preds.hba1c.hf.glp$ci.ub)),
                  col = cola[1], border = "transparent")
          polygon(x = c(newx.hba1c.hf.sglt, rev(newx.hba1c.hf.sglt)), 
                  y = c(preds.hba1c.hf.sglt$ci.lb, rev(preds.hba1c.hf.sglt$ci.ub)),
                  col = cola[2], border = "transparent")
          points(x = dhf$diff.hba1c, y = dhf$loghr, col = "black", bg = dhf$col, 
                 cex = (dhf$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
          lines(x = newx.hba1c.hf.glp, y = preds.hba1c.hf.glp$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          lines(x = newx.hba1c.hf.sglt, y = preds.hba1c.hf.sglt$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          
          # reference line
          abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
          # Axes properties
          axis(2, at = yaxis.hf, labels = yaxis.hf.lab, las = 2)
          axis(1, las = 1, labels = F)
          title(ylab = "Hazard ratio", line = 2.7, cex = 1.2)

    # 2.1 HBA1c CRO
          plot(x = 0, y = 0, type = 'n', xlim = xlim.hba1c, ylim = ylim.mi,
               axes = F, xlab = "", ylab = "")
          #CI, points, predictions
          polygon(x = c(newx.hba1c.cro.glp, rev(newx.hba1c.cro.glp)), 
                  y = c(preds.hba1c.cro.glp$ci.lb, rev(preds.hba1c.cro.glp$ci.ub)),
                  col = cola[1], border = "transparent")
          polygon(x = c(newx.hba1c.cro.sglt, rev(newx.hba1c.cro.sglt)), 
                  y = c(preds.hba1c.cro.sglt$ci.lb, rev(preds.hba1c.cro.sglt$ci.ub)),
                  col = cola[2], border = "transparent")
          points(x = dmi$diff.hba1c, y = dmi$loghr, col = "black", bg = dmi$col, 
                 cex = (dmi$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
          lines(x = newx.hba1c.cro.glp, y = preds.hba1c.cro.glp$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          lines(x = newx.hba1c.cro.sglt, y = preds.hba1c.cro.sglt$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          # reference line
          abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
          # Axes properties
          axis(2, at = yaxis.cro, labels = yaxis.cro.lab, las = 2)
          axis(1, las = 1, labels = F) 
          title(ylab = "Hazard ratio", line = 2.7, cex = 1.2)

    # 3.1 HBA1c ACM
          par(mar = c(2,2,0.5,0)*2)
          plot(x = 0, y = 0, type = 'n', xlim = xlim.hba1c, ylim = ylim.cvm,
               axes = F, xlab = "", ylab = "")
          #CI, points, predictions
          polygon(x = c(newx.hba1c.acm.glp, rev(newx.hba1c.acm.glp)), 
                  y = c(preds.hba1c.acm.glp$ci.lb, rev(preds.hba1c.acm.glp$ci.ub)),
                  col = cola[1], border = "transparent")
          polygon(x = c(newx.hba1c.acm.sglt, rev(newx.hba1c.acm.sglt)), 
                  y = c(preds.hba1c.acm.sglt$ci.lb, rev(preds.hba1c.acm.sglt$ci.ub)),
                  col = cola[2], border = "transparent")
          points(x = dacm$diff.hba1c, y = dacm$loghr, col = "black", bg = dacm$col, 
                 cex = (dacm$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
          lines(x = newx.hba1c.acm.glp, y = preds.hba1c.acm.glp$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          lines(x = newx.hba1c.acm.sglt, y = preds.hba1c.acm.sglt$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          # reference line
          abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
          # Axes properties
          axis(2, at = yaxis.acm, labels = yaxis.acm.lab, las = 2)
          axis(1, las = 1, labels = T) 
          title(ylab = "Hazard ratio", line = 2.7, cex = 1.2)
          title(xlab = expression(paste(Delta, "(treatment - control), %")), line = 2.5)

#### Bodyweight        
     
      # 1.2 Bodyweight HF
          par(mar = c(1,2,0.5,0)*2)
          plot(x = 0, y = 0, type = 'n', xlim = xlim.bw, ylim = ylim.hf,
               axes = F, xlab = "", ylab = "")
          title(main = expression(paste("Bodyweight")), line = 0.5)
          #CI, points, predictions
          polygon(x = c(newx.bw.hf.glp, rev(newx.bw.hf.glp)), 
                  y = c(preds.bw.hf.glp$ci.lb, rev(preds.bw.hf.glp$ci.ub)),
                  col = cola[1], border = "transparent")
          polygon(x = c(newx.bw.hf.sglt, rev(newx.bw.hf.sglt)), 
                  y = c(preds.bw.hf.sglt$ci.lb, rev(preds.bw.hf.sglt$ci.ub)),
                  col = cola[2], border = "transparent")
          points(x = dhf$diff.bodyweight, y = dhf$loghr, col = "black", bg = dhf$col, 
                 cex = (dhf$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
          lines(x = newx.bw.hf.glp, y = preds.bw.hf.glp$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          lines(x = newx.bw.hf.sglt, y = preds.bw.hf.sglt$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          # reference line
          abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
          # Axes properties
          axis(2, at = yaxis.acm, labels = F)
          axis(1, las = 1, labels = F)

      # 2.2 Bodyweight CRO
          plot(x = 0, y = 0, type = 'n', xlim = xlim.bw, ylim = ylim.cro,
               axes = F, xlab = "", ylab = "")
          #CI, points, predictions
          polygon(x = c(newx.bw.cro.glp, rev(newx.bw.cro.glp)), 
                  y = c(preds.bw.cro.glp$ci.lb, rev(preds.bw.cro.glp$ci.ub)),
                  col = cola[1], border = "transparent")
          polygon(x = c(newx.bw.cro.sglt, rev(newx.bw.cro.sglt)), 
                  y = c(preds.bw.cro.sglt$ci.lb, rev(preds.bw.cro.sglt$ci.ub)),
                  col = cola[2], border = "transparent")
          points(x = dmi$diff.bodyweight, y = dmi$loghr, col = "black", bg = dmi$col, 
                 cex = (dmi$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
          lines(x = newx.bw.cro.glp, y = preds.bw.cro.glp$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          lines(x = newx.bw.cro.sglt, y = preds.bw.cro.sglt$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          # reference line
          abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
          # Axes properties
          axis(2, at = yaxis.acm, labels = F)
          axis(1, las = 1, labels = F) 

      #  4.2 Bodyweight ACM
          par(mar = c(2,2,0.5,0)*2)
          plot(x = 0, y = 0, type = 'n', xlim = xlim.bw, ylim = ylim.acm,
               axes = F, xlab = "", ylab = "")
          #CI, points, predictions
          polygon(x = c(newx.bw.acm.glp, rev(newx.bw.acm.glp)), 
                  y = c(preds.bw.acm.glp$ci.lb, rev(preds.bw.acm.glp$ci.ub)),
                  col = cola[1], border = "transparent")
          polygon(x = c(newx.bw.acm.sglt, rev(newx.bw.acm.sglt)), 
                  y = c(preds.bw.acm.sglt$ci.lb, rev(preds.bw.acm.sglt$ci.ub)),
                  col = cola[2], border = "transparent")
          points(x = dacm$diff.bodyweight, y = dacm$loghr, col = "black", bg = dacm$col, 
                 cex = (dacm$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
          lines(x = newx.bw.acm.glp, y = preds.bw.acm.glp$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          lines(x = newx.bw.acm.sglt, y = preds.bw.acm.sglt$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          # reference line
          abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
          # Axes properties
          axis(2, at = yaxis.acm, labels = F)
          axis(1, las = 1, labels = T) 
          title(xlab = expression(paste(Delta, "(treatment - control), kg.")), line = 2.5)


#### Systolic blood pressure        

    # 1.2 SBP HHF
          par(mar = c(1,2,0.5,0)*2)
          plot(x = 0, y = 0, type = 'n', xlim = xlim.sbp, ylim = ylim.hf,
               axes = F, xlab = "", ylab = "")
          title(main = expression(paste("Systolic blood pressure")), line = 0.5)
          #CI, points, predictions
          polygon(x = c(newx.sbp.hf.glp, rev(newx.sbp.hf.glp)), 
                  y = c(preds.sbp.hf.glp$ci.lb, rev(preds.sbp.hf.glp$ci.ub)),
                  col = cola[1], border = "transparent")
          polygon(x = c(newx.sbp.hf.sglt, rev(newx.sbp.hf.sglt)), 
                  y = c(preds.sbp.hf.sglt$ci.lb, rev(preds.sbp.hf.sglt$ci.ub)),
                  col = cola[2], border = "transparent")
          points(x = dhf$diff.bodyweight, y = dhf$loghr, col = "black", bg = dhf$col, 
                 cex = (dhf$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
          lines(x = newx.sbp.hf.glp, y = preds.sbp.hf.glp$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          lines(x = newx.sbp.hf.sglt, y = preds.sbp.hf.sglt$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          # reference line
          abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
          # Axes properties
          axis(2, at = yaxis.acm, labels = F)
          axis(1, las = 1, labels = F)

    # 2.2 SBP CRO
          plot(x = 0, y = 0, type = 'n', xlim = xlim.sbp, ylim = ylim.cro,
               axes = F, xlab = "", ylab = "")
          #CI, points, predictions
          polygon(x = c(newx.sbp.cro.glp, rev(newx.sbp.cro.glp)), 
                  y = c(preds.sbp.cro.glp$ci.lb, rev(preds.sbp.cro.glp$ci.ub)),
                  col = cola[1], border = "transparent")
          polygon(x = c(newx.sbp.cro.sglt, rev(newx.sbp.cro.sglt)), 
                  y = c(preds.sbp.cro.sglt$ci.lb, rev(preds.sbp.cro.sglt$ci.ub)),
                  col = cola[2], border = "transparent")
          points(x = dmi$diff.bodyweight, y = dmi$loghr, col = "black", bg = dmi$col, 
                 cex = (dmi$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
          lines(x = newx.sbp.cro.glp, y = preds.sbp.cro.glp$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          lines(x = newx.sbp.cro.sglt, y = preds.sbp.cro.sglt$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          # reference line
          abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
          # Axes properties
          axis(2, at = yaxis.acm, labels = F)
          axis(1, las = 1, labels = F) 


      #  4.2 SBP ACM
          par(mar = c(2,2,0.5,0)*2)
          plot(x = 0, y = 0, type = 'n', xlim = xlim.sbp, ylim = ylim.acm,
               axes = F, xlab = "", ylab = "")
          #CI, points, predictions
          polygon(x = c(newx.sbp.acm.glp, rev(newx.sbp.acm.glp)), 
                  y = c(preds.sbp.acm.glp$ci.lb, rev(preds.sbp.acm.glp$ci.ub)),
                  col = cola[1], border = "transparent")
          polygon(x = c(newx.sbp.acm.sglt, rev(newx.sbp.acm.sglt)), 
                  y = c(preds.sbp.acm.sglt$ci.lb, rev(preds.sbp.acm.sglt$ci.ub)),
                  col = cola[2], border = "transparent")
          points(x = dacm$diff.sbp, y = dacm$loghr, col = "black", bg = dacm$col, 
                 cex = (dacm$logvi*500)^(0.5), pch = 21 , lwd = 0.4 )
          lines(x = newx.sbp.acm.glp, y = preds.sbp.acm.glp$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          lines(x = newx.sbp.acm.sglt, y = preds.sbp.acm.sglt$pred, 
                col = "gray40", lwd = 0.85, lty = 2)
          # reference line
          abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
          # Axes properties
          axis(2, at = yaxis.acm, labels = F)
          axis(1, las = 1, labels = T) 
          title(xlab = expression(paste(Delta, "(treatment - control), mmHg")), line = 2.5)

#### Hematocrit
        # Plot 1.4 hematocrit HF
          par(mar = c(1,2,0.5,0)*2)
          plot(x = 0, y = 0, type = 'n', xlim = range(xlim.hema), ylim = ylim.hf ,
               axes = F, xlab = "", ylab = "")
          title(main = expression(paste("Hematocrit")), line = 0.5)
          #CI, points, predictions
          polygon(x = c(newx.hema.hf.sglt, rev(newx.hema.hf.sglt)), 
                  y = c(preds.hema.hf.sglt$ci.lb, rev(preds.hema.hf.sglt$ci.ub)),
                  col = cola[2],border = "transparent")
          points(x = dhf$diff.hematocrit[dhf$type == "SGLT2"], y = dhf$loghr[dhf$type == "SGLT2"], 
                 col = "black", bg = dhf$col[dhf$type == "SGLT2"],
                 cex = (dhf$logvi[dhf$type == "GLP1"]*500)^(0.5),
                 pch = 21 , lwd = 0.4 )
          lines(newx.hema.hf.sglt, preds.hema.hf.sglt$pred, col = "gray40", lwd = 0.85, lty = 2)
          
          # reference line
          abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
          # Axes properties
          axis(2, at = yaxis.acm, labels = F)
          axis(1, las = 1, labels = F, at = xlim.hema) 
          
        # Plot 2.4 hematocrit cro
          plot(x = 0, y = 0, type = 'n', xlim = range(xlim.hema), ylim = ylim.cro ,
               axes = F, xlab = "", ylab = "")
          #CI, points, predictions
          polygon(x = c(newx.hema.cro.sglt, rev(newx.hema.cro.sglt)), 
                  y = c(preds.hema.cro.sglt$ci.lb, rev(preds.hema.cro.sglt$ci.ub)),
                  col = cola[2],border = "transparent")
          points(x = dcro$diff.hematocrit[dcro$type == "SGLT2"], y = dcro$loghr[dcro$type == "SGLT2"], 
                  col = "black", bg = dcro$col[dcro$type == "SGLT2"],
                  cex = (dcro$logvi[dcro$type == "SGLT2"]*500)^(0.5),
                  pch = 21 , lwd = 0.4 )
          lines(newx.hema.cro.sglt, preds.hema.cro.sglt$pred, col = "gray40", lwd = 0.85, lty = 2)
          # reference line
          abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
          # Axes properties
          axis(2, at = yaxis.cro, labels = F)
          axis(1, las = 1, labels = F) 
          
      # Plot 3.4 hematocrit ACM
          par(mar = c(2,2,0.5,0)*2)
          plot(x = 0, y = 0, type = 'n', xlim = range(xlim.hema), ylim = ylim.acm ,
               axes = F, xlab = "", ylab = "")
          #CI, points, predictions
          polygon(x = c(newx.hema.acm.sglt, rev(newx.hema.acm.sglt)), 
                  y = c(preds.hema.acm.sglt$ci.lb, rev(preds.hema.acm.sglt$ci.ub)),
                  col = cola[2],
                  border = "transparent")
          points(x = dacm$diff.hematocrit[dacm$type == "SGLT2"], y = dacm$loghr[dacm$type == "SGLT2"], 
                 col = "black", bg = dacm$col[dacm$type == "SGLT2"],
                 cex = (dacm$logvi[dacm$type == "SGLT2"]*500)^(0.5),
                 pch = 21 , lwd = 0.4 )
          lines(x = newx.hema.acm.sglt, y = preds.hema.acm.sglt$pred, col = "gray40", lwd = 0.85, lty = 2)
          # reference line
          abline(a = log(1), b = 0, lty = 3, col = "gray20", lwd = 1)
          # Axes properties
          axis(2, at = yaxis.acm, labels = F)
          axis(1, las = 1, labels = T, at = xlim.hema) 
          title(xlab = expression(paste(Delta, "(treatment - control), %")), line = 2.5)

          
      # Plot 5.1 UACR hf
          par(mar = c(1,2,0.5,0)*2)
          plot(x = 0, y = 0, type = 'n', xlim = range(xlim.uacr)+c(0,3), ylim = ylim.hf ,
               axes = F, xlab = "", ylab = "")
          legend("topright", 
                 pch =c(19,19,NA), lwd =c(NA,NA,2), 
                 col = c("gold","gray70", "black"), c("GLP-1RA", "SGLT2i", "P<0.05"), 
                 cex = 0.75, box.col = "transparent", ncol = 1)
          #CI, points, predictions
          polygon(x = c(newx.uacr.hf.glp, rev(newx.uacr.hf.glp)), 
                  y = c(preds.uacr.hf.glp$ci.lb, rev(preds.uacr.hf.glp$ci.ub)),
                  col = cola[1],
                  border = "transparent")
          points(x = dhf$diff.uacr[dhf$type == "GLP1"], y = dhf$loghr[dhf$type == "GLP1"], 
                 col = "black", bg = dhf$col[dhf$type == "GLP1"],
                 cex = (dhf$logvi[dhf$type == "GLP1"]*500)^(0.5),
                 pch = 21 , lwd = 0.4 )
          lines(newx.uacr.hf.glp, preds.uacr.hf.glp$pred, col = "gray40", lwd = 0.8, lty = 3)
          title(main = expression(paste("uACR")), line = 0.5)
          
          # reference line
          segments(x0 = -10, x1 = 2, y0 = log(1), col = "gray20", lwd = 1, lty = 3)
          # Axes properties
          axis(2, at = yaxis.acm, labels = F)
          axis(1, las = 1, at = xlims.uacr, labels = F)
          text( x = 2.1, y = -.1, labels = "Heart failure", cex = 1.1, srt = 270)
          text( x = -0.45, y = log(1.15), labels = "Favors \n control", cex = 0.9, srt = 0)
          text( x = -0.45, y = log(0.85), labels = "Favors \n treatment", cex = 0.9, srt = 0)
          
          
      # Plot 5.2 UACR MI
          plot(x = 0, y = 0, type = 'n', xlim = range(xlim.uacr)+c(0,3), ylim = ylim.cro ,
               axes = F,
               xlab = "", ylab = "")
          #CI, points, predictions
          polygon(x = c(newx.uacr.cro.glp, rev(newx.uacr.cro.glp)), 
                  y = c(preds.uacr.cro.glp$ci.lb, rev(preds.uacr.cro.glp$ci.ub)),
                  col = cola[1],
                  border = "transparent")
          points(x = dcro$diff.uacr[dcro$type == "GLP1"], y = dcro$loghr[dmi$type == "GLP1"], 
                 col = "black", bg = dmi$col[dcro$type == "GLP1"],
                 cex = (dcro$logvi[dmi$type == "GLP1"]*500)^(0.5),
                 pch = 21 , lwd = 0.4 )
          lines(newx.uacr.cro.glp, preds.uacr.cro.glp$pred, col = "gray40", lwd = 0.8, lty = 3)
          # reference line
          segments(x0 = -10, x1 = 2, y0 = log(1), col = "gray20", lwd = 1, lty = 3)
          # Axes properties
          axis(2, at = yaxis.acm, labels = F)
          axis(1, at = xlims.uacr, labels = F)
          text( x = 2.1, y = -.25, labels = "Composite renal outcome", cex = 1.1, srt = 270)
          text( x = -0.45, y = log(1.15), labels = "Favors \n control", cex = 0.9, srt = 0)
         text( x = -0.45, y = log(0.85), labels = "Favors \n treatment", cex = 0.9, srt = 0)
          
      # 5.3 UACR acm
          par(mar = c(2,2,0.5,0)*2)
          plot(x = 0, y = 0, type = 'n', xlim = xlim.uacr+c(0,3), ylim = ylim.acm ,
               axes = F,
               xlab = "", ylab = "")
          #CI, points, predictions
          polygon(x = c(newx.uacr.acm.glp, rev(newx.uacr.acm.glp)), 
                  y = c(preds.uacr.acm.glp$ci.lb, rev(preds.uacr.acm.glp$ci.ub)),
                  col = cola[1],
                  border = "transparent")
          points(x = dacm$diff.uacr[dacm$type == "GLP1"], y = dacm$loghr[dacm$type == "GLP1"], 
                 col = "black", bg = dacm$col[dacm$type == "GLP1"],
                 cex = (dacm$logvi[dacm$type == "GLP1"]*500)^(0.5),
                 pch = 21 , lwd = 0.4 )
          lines(newx.uacr.acm.glp, preds.uacr.acm.glp$pred,  col = "gray40", lwd = 0.8, lty = 3)
          # reference line
          segments(x0 = -10, x1 = 2, y0 = log(1), col = "gray20", lwd = 1, lty = 3)
          
          # Axes properties
          axis(2, at = yaxis.acm, labels = F)
          axis(1, las = 1, at = xlims.uacr, labels = T)
          title(xlab = expression(paste(Delta, "(treatment - control), mg/g")), line = 2.5)
          text( x = 2.1, y = -.25, labels = "All-cause mortality", cex = 1.1, srt = 270)
          text( x = -0.45, y = log(1.15), labels = "Favors \n control", cex = 0.9, srt = 0)
          text( x = -0.45, y = log(0.85), labels = "Favors \n treatment", cex = 0.9, srt = 0)
          
dev.off()
                


