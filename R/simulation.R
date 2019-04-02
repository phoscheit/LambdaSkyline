#' This function simulates the ancestral process of a sample of n lineages in a larger population which has size of order
#' N and fluctuates according to pop_size. The simulation is based on the results of (Schweinsberg 2003) and relies
#' on supercritical branching processes with heavy-tailed offspring distribution. 
#'
#' @title Simulate Beta-coalescents with varying population size
#' @param n Number of sampled lineages at present time
#' @param alpha The parameter of the Beta(2-alpha,alpha)-coalescent
#' @param N Scaling parameter of the population size
#' @param pop_size Vector containing the population size factor (true population size at time i is pop_size[i]*N)
#' 
#' @export

simcoal <- function(n,alpha,N,pop_size){
  T_F <- length(pop_size) # Number of generations to be simulated
  t <- T_F # Current time index
  
  b <- n # Initial number of blocks
  V <- as.list(1:n) # Initial partition (dust)
  
  Lengths <- vector("numeric",length=n) 
  Branches <- as.list(vector(length=n,mode="numeric")) # Vector containing the current branch lengths
  
  Labels <- sort(sample(1:(ceiling(pop_size[t]*N)),size=n,replace=FALSE)) # Will contain the labels (lineage numbers) of ancestral lineages
  Lines <- as.list(1:n) # Will contain an index (in 1:n) of one representative lineage for each block
  Trees <- vector("list",n) # Will contain the current string describing the tree
  for (i in 1:n){
    Trees[[i]]<- ape::read.tree(text=paste("(",as.character(Labels[i]),":0);",sep=""))
    Trees[[i]]$root.edge <- 0
  }
  
  Reseau <- distr::Lattice(pivot=2,width=1,Length=10000)
  Poids <- vector(mode="numeric",length=10000)
  for(i in 1:9998) {
    Poids[i+1]<-i^(-alpha-1)
  }
  Poids[10000] <- 10000^(-alpha)
  Somme <- sum(Poids)
  for(i in 1:9999) {
    Poids[i+1]<- Poids[i+1]*2/(3*Somme)
  }
  Poids[1]<-1/3
  LoiAlpha <- distr::LatticeDistribution(lattice=Reseau,prob=Poids)
  rm(Poids,Somme,Reseau)
  
  while(b!=1 & t!=1){
    X <- distr::r(LoiAlpha)(ceiling(pop_size[t-1]*N))
    NumDes <- sum(X)
    DesPot <- vector(mode="numeric",length=NumDes)
    j<-1
    for(k in 1:(ceiling(pop_size[t-1]*N))){
      for(l in 1:X[k]){
        DesPot[j]<-k
        j<-j+1
      }
    }
    Ancestors <- sample(DesPot,size=b,replace=FALSE) # We sample the ancestors of the current lineages among potential descendants
    
    for(j in 1:b){ 
      Trees[[j]]$root.edge <- (Trees[[j]]$root.edge)+1
    }
    
    if(length(unique(Ancestors))!=b){# If there are coalescences
      ToRemove <- vector()
      for(lab in unique(Ancestors)){ # Enumerate all lineages at time i
        indices <- which(Ancestors==lab) 
        if(length(indices)!=1){ # More than one block is coalescing into lineage lab
          i0 <- 0
          for(k in indices){
            if(i0 ==0){ # First lineage becomes basis for this coalescence
              i0 <- k
            }
            else{
              Trees[[i0]] <- ape::collapse.singles(Trees[[i0]]+Trees[[k]])
              b <- b-1
            }
          }
          ToRemove <- c(ToRemove,setdiff(indices,i0))
        }
      }
      Trees <- Trees[-ToRemove]
      Lines <- Lines[-ToRemove]
    }
    print(c(t,b))
    t<-t-1
  }
  Trees[[1]]
}