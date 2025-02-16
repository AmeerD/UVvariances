## All sims
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(sampling)
library(survey)
library(purrr)
library(argparse)

load("WAdat.rda")
source("1_functions.R")

## -----------------------------------------
## Load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("--nreps", type = "double", default = 10,
                    help = "number of replicates for each set of params")
args <- parser$parse_args()
jobid <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(jobid)

nsim <- args$nreps
set.seed(jobid)

estfile <- paste("res/est", "_", jobid, ".txt", sep="")
varfile <- paste("res/var", "_", jobid, ".txt", sep="")

## -----------------------------------------
## Set up simulation data
## -----------------------------------------
## First set up the population data (two responses: real-scale income and log-scale income)
WAdat <- WAdat %>% arrange(puma) %>% filter(income > 0) %>% 
  mutate(log.income = log(income), idx=row_number()) %>% 
  group_by(puma) %>% mutate(nc = n()) %>% ungroup
## Define linear regression formula
form.real <- income ~ 1 + age + transit + workhrs
form.log <- log.income ~ 1 + age + transit + workhrs
## Define population totals
N <- nrow(WAdat)
Xall <- model.matrix(form.real, WAdat) ## Same for both formulas so just use real
XTXinv <- solve(t(Xall) %*% Xall) 
tx <- colSums(Xall)
yobs.real <- mean(WAdat$income)
yobs.log <- mean(WAdat$log.income)

# Define simulation parameters
nsamp <- c(30, 50, 100, 250, 500) ## This is the number of samples to be drawn from WAdat

# Set up results arrays - 2500 sims x 5 samples x 5 methods x 5 designs x 2 outcomes 
## First define the labels for each array dimension
samps <- as.character(nsamp)
methods <- c("asymptotic", "IJ", "Hdecomp", "Hdecomp.2b", "Empirical")
designs <- c("bernoulli", "poisson", "srswor", "clustered", "stratified")
outcomes <- c("real", "log")

# Start simulation:
## 1. Choose the sample size x design
## 2. For each simulation, draw a sample
## 3. Estimate the population average for each outcome
## 4. Estimate the variance using each method
for (ns in nsamp) {
  for (d in designs) {
    ## Skip small sample stratified designs (61 PUMAs not compatible with sample sizes under 100)
    if (d == "stratified" & ns < 100) {
      next
    }
    
    ## First copmute the inclusion probabilities given the sample size and design
    IP <- IPsim(ns, N, d, fulldat=WAdat) 
    
    for (i in 1:nsim) {
      ## Draw a sample
      samp <- getsamp(WAdat, IP, ns, N, d)
      
      ## Construct the 2nd order joint inclusion probability matrix (and Delta)
      ## Note that the diagonal entries are the 1st order inclusion probabilities
      pimat <- getpi(IP, samp, ns, N, d)
      Delta <- pimat - outer(samp$Prob, samp$Prob)
      
      ## Construct the sample X matrices GREG estimates
      mm <- model.matrix(form.real, samp) ## Same for both formulas so just use real
      mpi <- matrix(rep(samp$Prob, ncol(mm)), ncol=ncol(mm), byrow=F)
      tHT <- colSums(mm/mpi)
      
      ## Iterate through outcomes
      for (inc in outcomes) {
        # SETUP
        ## Set the outcome vector
        if (inc == "real") {
          yvec <- samp$income
        } else {
          yvec <- samp$log.income
        }
        
        # ESTIMATION
        ## Compute GREG estimates and residuals
        omega <- 1/samp$Prob + (tx - tHT) %*% XTXinv %*% t(mm/mpi)
        ahat <- sum(omega*yvec)/N
        lm.res <- yvec - as.vector(mm %*% XTXinv %*% t(mm) %*% (yvec/samp$Prob))
        
        ## Save point estimate
        write(c(i, ns, nrow(samp), d, inc, ahat), file = estfile, append=TRUE, ncolumns = 6)
        
        # VARIANCES 
        ## Asymptotic GREG variance estimate
        ASYtime <- system.time(ASYvar <- (sum((Delta/pimat)*outer(lm.res/samp$Prob,lm.res/samp$Prob)))/(N^2))
        write(c(i, ns, "asymptotic", d, inc, ASYtime[3], ASYvar), file = varfile, append=TRUE, ncolumns = 7)
        
        ## Infinitessimal jackknife
        IJtime <- system.time(IJvar <- GREGIJ(yvec, mm, tx, samp$Prob, ahat, XTXinv))
        write(c(i, ns, "IJ", d, inc, IJtime[3], IJvar), file = varfile, append=TRUE, ncolumns = 7)
        
        ## H-decomposition
        ### No tau2b
        HDtime <- system.time(HDvar <- GREGHdecomp(yvec, mm, tx, samp$Prob, ahat, XTXinv, t2b=F))
        write(c(i, ns, "Hdecomp", d, inc, HDtime[3], HDvar), file = varfile, append=TRUE, ncolumns = 7)
        ### With tau2b
        HDtime <- system.time(HDvar <- GREGHdecomp(yvec, mm, tx, samp$Prob, ahat, XTXinv, t2b=T))
        write(c(i, ns, "Hdecomp.2b", d, inc, HDtime[3], HDvar), file = varfile, append=TRUE, ncolumns = 7)
      }
    }
  }
  
  print(paste0("Done sample size ", ns))
}

