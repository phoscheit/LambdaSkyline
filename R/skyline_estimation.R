# Main estimation function : takes as input a coalescentIntervals object, as well as the alpha parameter of the Beta coalescent, and the grouping parameter
# epsilon, and returns a skyline object, along with its log-likelihood

Skyline.multi.coalescentIntervals <- function(Inter,alpha,epsilon=0.0){
  N <- Inter$interval.count # Total number of intervals in the phylogeny
  Factors <- sapply(Inter$lineages,function(x) TotalRate(x,alpha)) # Total rate of coalescence on each interval
  Estimate <- numeric(N) # Logs the skyline estimate
  Coal <- Inter$coalescences # Records if the intervals end with a coalescence or not (regardless of samplings)
  Involved <- Inter$NumInvolved
  
  NumSampling <- 0
  CurrEstimate <- 0
  for(i in 1:N){
    if(Coal[i]){
      CurrEstimate <- CurrEstimate + Factors[i]*Inter$interval.length[i]
      for(j in (i-NumSampling):i){
        Estimate[j] <- CurrEstimate
      }
      NumSampling <- 0
      CurrEstimate <- 0
    }
    else{
      CurrEstimate <- CurrEstimate + Factors[i]*Inter$interval.length[i]
      NumSampling <- NumSampling + 1
    }
  }
  
  NInter <- N
  
  if(epsilon > 0.0){
    CollInter <- collapsed.intervals(Inter,epsilon)
    NColl <- CollInter$collapsed.interval.count
    EstimateColl <- numeric(NColl)
    for(i in 1:NColl){
      EstimateColl[i] <- mean(Estimate[CollInter$collapsed.interval==i])
    }
    for(i in 1:N){
      Estimate[i] <- EstimateColl[CollInter$collapsed.interval[i]]
    }
    NInter <- NColl
  }
  
  Comb <- numeric(N)
  Lambdas <- numeric(N)
  Contrib <- numeric(N)
  for(i in 1:N){
    if(Coal[i]){
      # Comb[i] <- lchoose(Inter$lineages[i],Involved[i])
      Lambdas[i] <- lbeta(Involved[i]-alpha,Inter$lineages[i]-Involved[i]+alpha) - lbeta(2-alpha,alpha)
      Contrib[i] <- Comb[i]+Lambdas[i]-log(Estimate[i])-Factors[i]*(Inter$interval.length[i])/Estimate[i]
    }
    else{
      Contrib[i] <- -Factors[i]*(Inter$interval.length[i])/Estimate[i]
    }
  }
  obj <- list(
    time=cumsum(Inter$interval.length),
    interval.length=Inter$interval.length,
    population.size = Estimate,
    parameter.count=NInter,
    epsilon=epsilon,
    logL = sum(Contrib)  
  )
  class(obj)<- "skyline"
  return(obj)
}

Skyline.multi.phylo <- function(phylo,alpha,epsilon=0.0){
  Inter <- coalescent.intervals.multi(phylo)
  Skyline.multi.coalescentIntervals(Inter,alpha,epsilon)
}



# Performs a maximum-likelihood estimation of the alpha parameter, with the multifurcating skyline
# as demographic function.
# epsilon is the grouping parameter, as defined originally in (Strimmer, Pybus 2001)

BetaCoal.maxlik <- function(phylo,epsilon=0.0) {
  Inter <- coalescent.intervals.multi(phylo)
  optimx(par=1.5, fn=function(x) -Skyline.multi.coalescentIntervals(Inter,x,epsilon)$logL,lower=0.001,upper=1.999, method="bobyqa")
}

# Computes and plots the likelihood landscape of a given phylogeny as a function of alpha

FitnessLandscape <- function(phylo,epsilon=0.0){
  Inter <- coalescent.intervals.multi(x = phylo)
  Lik <- sapply(seq(from=1,to=1.99,length.out = 100),function(x) Skyline.multi.coalescentIntervals(Inter,alpha=x,epsilon)$logL)
  graphics::plot(seq(from=1,to=1.99,length.out = 100),Lik)
  return(Lik)
}
