IPsim <- function(n, N, design, fulldat) {
  if (design == "bernoulli") {
    return(inclusionprobabilities(rep(1, N), n))
  } else if (design == "poisson") {
    baseweight <- 0.5*fulldat$pweight + 0.5*mean(fulldat$pweight)
    return(inclusionprobabilities(1/baseweight, n))
  } else if (design == "srswor") {
    return(rep(n/N, N))
  } else if (design == "clustered") { 
    npuma <- fulldat$puma %>% unique %>% length
    if (n <= 100) {
      return(fulldat %>% mutate(IP = ((n/2)/npuma)*(2/nc)) %>% pull(IP))
    } else {
      return(fulldat %>% mutate(IP = (50/npuma)*((n/50)/nc)) %>% pull(IP))
    }
  } else if (design == "stratified") {
    nstrata <- length(unique(fulldat$puma))
    nperstrata <- max(ceiling(n/nstrata), 2)
    return(inclusionprobastrata(as.numeric(as.factor(fulldat$puma)), rep(nperstrata, nstrata)))
  }
}

getsamp <- function(dat, IP, n, N, design) {
  if (design == "bernoulli") {
    return(getdata(dat, UPpoisson(IP)) %>% 
             mutate(idxsamp = row_number(), Prob=IP[idx], pop=N) )
  } else if (design == "poisson") {
    return(getdata(dat, UPpoisson(IP)) %>% 
             mutate(idxsamp = row_number(), Prob=IP[idx], pop=N))
  } else if (design == "srswor") {
    return(getdata(dat, srswor(n, N)) %>% 
             mutate(idxsamp = row_number(), Prob=n/N, pop=N))
  } else if (design == "clustered") {
    if (n < 100) {
      return(getdata(dat, mstage(dat, stage=c("cluster", ""), varnames=c("puma", "idx"), 
                                 size=list(n/2, rep(2, n/2)), 
                                 method=c("srswor", "srswor")))[[2]] %>% 
               mutate(idxsamp = row_number(), pop=N))
    } else {
      return(getdata(dat, mstage(dat, stage=c("cluster", ""), varnames=c("puma", "idx"), 
                                 size=list(50, rep(n/50, 50)), 
                                 method=c("srswor", "srswor")))[[2]] %>% 
               mutate(idxsamp = row_number(), pop=N))
    }
  } else if (design == "stratified") {
    nstrata <- length(unique(dat$puma))
    nperstrata <- max(ceiling(n/nstrata), 2)
    return(getdata(dat, sampling::strata(dat, "puma", size=rep(nperstrata, nstrata), 
                                         method = "srswor")) %>%
             mutate(idxsamp = row_number(), pop=N))
  }
}

# Function to compute the joint inclusion probability matrix
# (first-order inclusion probabilities given on the diagonal)
# ncluster is the number of clusters sampled
# nper is the number of units sampled in each cluster/strata
getpi <- function(IP, samp, n, N, design) {
  if (design %in% c("bernoulli", "poisson")) {
    temp <- outer(IP[samp$idx], IP[samp$idx])
    diag(temp) <- IP[samp$idx]
    return(temp)
  } else if (design == "srswor") {
    temp <- matrix((n*(n-1))/(N*(N-1)), nrow=n, ncol=n)
    diag(temp) <- IP[samp$idx]
    return(temp)
  } else if (design == "clustered") { 
    ncluster <- min(50, n/2)
    nper <- n/ncluster
    temp <- diag(IP[samp$idx])
    
    for (i in 2:n) {
      for (j in 1:(i-1)) {
        if (samp$puma[i] == samp$puma[j]) {
          temp[i,j] <- temp[j,i] <- (ncluster/61)*((nper*(nper-1))/(samp$nc[i]*(samp$nc[i]-1)))
        } else {
          temp[i,j] <- temp[j,i] <- ((ncluster*(ncluster-1))/(61*60))*(nper/samp$nc[i])*(nper/samp$nc[j])
        }
      }
    }
    return(temp)
  } else if (design == "stratified") {
    nper <- max(ceiling(n/61), 2)
    temp <- outer(IP[samp$idx], IP[samp$idx])
    diag(temp) <- IP[samp$idx]
    
    for (i in 2:nrow(samp)) {
      for (j in 1:(i-1)) {
        if (samp$puma[i] == samp$puma[j]) {
          temp[i,j] <- temp[j,i] <- (nper*(nper-1))/(samp$nc[i]*(samp$nc[i]-1))
        }
      }
    }
    return(temp)
  }
}

makedes <- function(dat, design) {
  if (design == "bernoulli") {
    return(svydesign(~1, probs=~Prob, pps=poisson_sampling(dat$Prob), data=dat))
  } else if (design == "poisson") {
    return(svydesign(~1, probs=~Prob, pps=poisson_sampling(dat$Prob), data=dat))
  } else if (design == "srswor") {
    return(svydesign(~1, probs=~Prob, fpc=~pop, data=dat))
  } else if (grepl("clustered", design)) {
    return(svydesign(~puma + idx, probs=~Prob, fpc=~pop + nc, data=dat))
  } else if (design == "stratified") {
    return(svydesign(~1, probs=~Prob, strata=~puma, fpc=~pop, data=dat))
  }
}

## H-decomposition based variance estimation for GREG
GREGHdecomp <- function(ysamp, xsamp, tx, sIP, ahat, Sinv, t2b=F) {
  N <- tx[1]
  n <- nrow(xsamp)
  piratio <- (1-sIP)/sIP
  pi2 <- outer(sIP, sIP)
  Nht <- N - 1/sIP
  
  # Prep V-statistic kernels
  h.unsym <- (1 + N*((matrix(rep(tx/N, n), nrow=n, byrow=T) - xsamp/sIP) %*% Sinv %*% t(xsamp)))*matrix(rep(ysamp/sIP, n), nrow=n, byrow=T)
  hV.full <- 0.5*(h.unsym + t(h.unsym))
  hV.partial <- 0.5*(1 + tx %*% Sinv %*% t(xsamp)) * ysamp/sIP
  
  # Prep U-statistic kernels
  ## nxn symmetric matrix of fully observed kernels
  hU.full <- ((N-1)/(2*N))*(2*hV.full + outer(diag(hV.full), diag(hV.full), FUN="+")/(N-1))
  diag(hU.full) <- 0
  ## length n vector of partially observed kernels
  hU.partial <- as.vector(((N-1)/(2*N))*(2*hV.partial + diag(hV.full)/(N-1)))
  
  # Prep observed theta' values
  ## Compute fully observed theta values
  theta <- (hU.full-outer(hU.partial, hU.partial, FUN="+"))*pi2 + outer(hU.partial*sIP, hU.partial*sIP, FUN="+")
  diag(theta) <- 0
  ## Compute weighted estimates of theta' 
  theta.prime <- colSums(theta/matrix(rep(sIP,n), nrow=n, byrow=F))/Nht #(colSums(1/matrix(rep(sIP,n), nrow=n, byrow=F))-1)
  
  # Prep observed phi' terms
  ## Start with unobserved phi terms (phi(Ii=0) = hU(0,1)*(pi.j-pi.ij)/(1-pi.i)=hU(0,1)*pi.j under Poisson sampling)
  ## Then, when we take weighted estimates, it becomes sum_{j\in S}(hU(0,1))
  ### Observed
  hUp.mat <- matrix(rep(hU.partial, n), nrow=n, byrow=T)
  diag(hUp.mat) <- 0
  # phi0p.observed <- rowSums(hUp.mat)/(N-1) ## This never gets used
  ### Unobserved
  phi0p.unobserved <- sum(hU.partial)/(N-1)
  # phi1.observed <- colSums(hU.full-hUp.mat)/(N-1) + hU.partial
  ## Next do the observed phi terms, phi'(Ii=1)
  phi1mat <- (hU.full-hUp.mat)*matrix(rep(sIP,n), nrow=n, byrow=F) + hUp.mat
  phi1p.observed <- colSums(phi1mat/matrix(rep(sIP,n), nrow=n, byrow=F))/Nht
  
  # Estimate tau1 and the bias correction
  ## Start with tau1
  tau1 <- sum((phi1p.observed - theta.prime)^2/(1-sIP))/(N^2)
  ## For the bias correction, construct (phi_{1;i,j}(I_i=1) - theta_{i,j})/pi.j for all observed pairs
  bcmat <- (phi1mat - theta)/matrix(rep(sIP,n), nrow=n, byrow=F)

  phi.observed <- colSums(phi1mat/matrix(rep(sIP,n), nrow=n, byrow=F))/(N-1) ## SAME AS phi1p.observed
  
  ## Compute correction term (HT variance of h_{1;i}(I_i=1)) for each observed unit
  bc.tau1.vec <- rep(NA, n)
  for (i in 1:n) {
    bc.tau1.vec[i] <- sum((1-sIP[-i])*(bcmat[-i,i]^2))/(1-sIP[i])/(Nht[i]^2)
  }
  bc.tau1 <- sum(bc.tau1.vec)/(N^2)
  tau1.BCF <- max(tau1 - bc.tau1, 0)
  
  # Estimate tau2a, tau2b, and the bias correction
  ## Compute hU-phi.i-phi.j+a terms for tau2a. There are three cases:
  ## 1. Fully observed kernels
  t2a.1 <- hU.full - outer(phi1p.observed, phi1p.observed, FUN="+") + ahat
  ## 2. Partially observed kernels
  t2a.2 <- hU.partial - phi1p.observed - phi0p.unobserved + ahat
  ## 3. Fully unobserved kernels
  t2a.3 <- 0 - phi0p.unobserved - phi0p.unobserved + ahat
  ## Compute tau2a
  tau2a <- sum(t2a.1[lower.tri(t2a.1)]^2) + (N-n)*sum(t2a.2^2) + choose(N-n,2)*t2a.3^2
  
  ## Compute the bias correction
  bc.tau2.mat <- matrix(NA, nrow=n, ncol=n)
  for (i in 2:n) {
    for (j in 1:i) {
      bc.tau2.mat[i,j] <- (((Nht[i]+Nht[j])/(Nht[i]*Nht[j]))^2)*sum((1-sIP)*(bcmat[,i]+bcmat[,j])*(bcmat[,i]+bcmat[,j]))
    }
  }
  bc.tau2 <- sum((bc.tau2.mat/pi2)[lower.tri(bc.tau2.mat)])
  
  
  if (t2b) {
    ## Compute tau2b
    tau2b.temp <- ((theta - outer(theta.prime, theta.prime, FUN="+") + ahat)^2)/pi2
    tau2b <- sum(tau2b.temp[lower.tri(tau2b.temp)])
    
    tau2.BCF2b <- max(4*(tau2a - tau2b - bc.tau2)/((N*(N-1))^2), 0)
    vhat2b <- 4*tau1.BCF + tau2.BCF2b
    return(vhat2b)
  } else {
    tau2.BCF <- max(4*(tau2a - bc.tau2)/((N*(N-1))^2), 0)
    vhat <- 4*tau1.BCF + tau2.BCF
    return(vhat)
  }
}

## Infinitessimal Jackknife based variance estimation for GREG
GREGIJ <- function(ysamp, xsamp, tx, sIP, ahat, Sinv) {
  N <- tx[1]
  n <- nrow(xsamp)
  piratio <- (1-sIP)/sIP
  pi2 <- outer(sIP, sIP)
  
  # Prep V-statistic kernels
  h.unsym <- (1 + N*((matrix(rep(tx/N, n), nrow=n, byrow=T) - xsamp/sIP) %*% Sinv %*% t(xsamp)))*matrix(rep(ysamp/sIP, n), nrow=n, byrow=T) ##PROBLEM!!
  hV.full <- 0.5*(h.unsym + t(h.unsym))
  hV.partial <- 0.5*(1 + tx %*% Sinv %*% t(xsamp)) * ysamp/sIP
  
  # Prep U-statistic kernels
  ## nxn symmetric matrix of fully observed kernels
  hU.full <- ((N-1)/(2*N))*(2*hV.full + outer(diag(hV.full), diag(hV.full), FUN="+")/(N-1))
  diag(hU.full) <- 0
  ## length n vector of partially observed kernels
  hU.partial <- as.vector(((N-1)/(2*N))*(2*hV.partial + diag(hV.full)/(N-1)))
  
  # Prep observed phi' terms
  ## Start with unobserved phi terms (phi(Ii=0) = hU(0,1)*(pi.j-pi.ij)/(1-pi.i)=hU(0,1)*pi.j under Poisson sampling)
  ## Then, when we take weighted estimates, it becomes sum_{j\in S}(hU(0,1))
  ### Observed
  hUp.mat <- matrix(rep(hU.partial, n), nrow=n, byrow=T)
  diag(hUp.mat) <- 0
  ### Unobserved
  phi0p.unobserved <- sum(hU.partial)/(N-1)
  # phi1.observed <- colSums(hU.full-hUp.mat)/(N-1) + hU.partial
  ## Next do the observed phi terms, phi'(Ii=1)
  phi1mat <- (hU.full-hUp.mat)*matrix(rep(sIP,n), nrow=n, byrow=F) + hUp.mat
  phi1p.observed <- colSums(phi1mat/matrix(rep(sIP,n), nrow=n, byrow=F))/(N-1)
  
  return( 4*((N-1)^2)*(sum((phi1p.observed - ahat)^2) + (N-n)*(phi0p.unobserved-ahat)^2)/((N^2)*(N-2)^2) )
}






