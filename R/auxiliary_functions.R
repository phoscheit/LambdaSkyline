# Computes the total rate of coalescence for n lineages under the Beta(2-alpha,alpha)-coalescent

total_coal_rate <- function(n,alpha) 
{ 
  if(n<2) return(0)
  Comb <- sapply(2:n,function(x) lchoose(n,x))
  Lambda <- sapply(2:n,function(x) lbeta(x-alpha,n-x+alpha))
  Lambda <- Lambda-lbeta(2-alpha,alpha)
  sum(sapply(Comb+Lambda,exp))
}

# Extension of the coalescent.intervals function in ape. Takes a phylogeny (class phylo) as input and returns a coalescentIntervals object,
# which lists the interval lengths, in backwards time, along with the number of lineages during each of those intervals. Contrary to the ape
# function, this can deal with multifurcations, as well as serially-sampled (non-ultrametric) trees.

coalescent.intervals.multi <- function(x) 
{
  if (class(x) != "phylo") stop("object \"x\" is not of class \"phylo\"")
  
  # ordered branching times
  t <- sort(branching.times.samp(x))
  
  t <- signif(t,digits=10) # Required to eliminate numerical instabilities
  Dist <- unique(t) # We group all identical times together
  lt <- length(Dist) # Total number of intervals
  
  # interval widths
  w <- numeric(lt)
  w[1] <- Dist[1]
  for (i in 2:lt) w[i] <- Dist[i] - Dist[i - 1]
  
  # Logs if a coalescence happens or not
  Coal <- logical(lt)
  
  # Logs the number of lineages involved in coalescences (breaks down in the case of simultaneous coalescences)
  NumInvolved <- integer(lt)
  
  l<-numeric(lt)
  l[1] <- 0 # Number of lineages at present time
  for(i in 2:lt){
    Nodes <- as.integer(names(t[t==Dist[i-1]])) # Identifies the nodes (leaves or internal) with timestamp Dist[i-1]
    
    CoalNodes <- Nodes[Nodes>length(x$tip.label)] # Finds all internal nodes, corresponding to a coalescence event
    SampNodes <- Nodes[Nodes<=length(x$tip.label)] # Finds all leaves, corresponding to a sampling event
    Involved <- 0 # Will record the total number of lineages coalescing at that time
    
    if(length(CoalNodes) >0){
      # if(length(CoalNodes)>1) print("Simultaneous Coalescences detected! Results meaningless!")
      Coal[i-1] <- TRUE
      for(j in 1:length(CoalNodes)){
        Edges <- x$edge[x$edge[,1]==CoalNodes[j],] # Finds how many lineages coalesce at node CoalNodes[j]
        Involved <- Involved + length(Edges[,1]) # Updates the count of coalescing lineages
        NumInvolved[i-1] <- length(Edges[,1])
      }
    }
    l[i]<-l[i-1]-Involved+length(CoalNodes)+length(SampNodes)
  }
  NumInvolved[lt] <- l[lt]
  Coal[lt] <- TRUE
  obj <- list(
    lineages=l,
    interval.length=w,
    interval.count=lt,
    total.depth =sum(w),
    coalescences = Coal,
    NumInvolved=NumInvolved)
  class(obj) <- "coalescentIntervals"
  return(obj)
}

# Extension of the branching.times function in ape. Takes a phylogeny (class phylo) and returns a vector with all event times in backwards time
# including the coalescences and the sampling events.

branching.times.samp <- function(phy)
{
  if (!inherits(phy, "phylo"))
    stop('object "phy" is not of class "phylo"')
  DistToRoot <- dist.nodes(phy)[length(phy$tip.label)+1,] # length(phy$tip.label)+1 is the label of the root, so we get all distances to the root
  TotalDepth <- max(DistToRoot)
  TotalDepth - DistToRoot # We want time backwards from the last sampling event
}

# Transforms a binary phylogeny into a multifurcating one, collapsing all branches shorter than tol

di2multi.cons <- function(phylo, tol = 1e-8)
{
  if (is.null(phylo$edge.length)) stop("the tree has no branch length")
  ## We select only the internal branches which are
  ## significantly small:
  ind <- which(phylo$edge.length < tol & phylo$edge[, 2] > length(phylo$tip.label))
  n <- length(ind)
  if (!n) return(phylo)
  ## recursive function to `propagate' node #'s in case
  ## there is a series of consecutive edges to remove
  foo <- function(ancestor, des2del,length) {
    wh <- which(phylo$edge[, 1] == des2del)
    for (k in wh) {
      if (phylo$edge[k, 2] %in% node2del) foo(ancestor, phylo$edge[k, 2],length+phylo$edge.length[k])
      else{
        phylo$edge[k, 1] <<- ancestor
        phylo$edge.length[k] <<- phylo$edge.length[k]+length
      } 
    }
  }
  node2del <- phylo$edge[ind, 2]
  anc <- phylo$edge[ind, 1]
  long <- phylo$edge.length[ind]
  for (i in 1:n) {
    if (anc[i] %in% node2del) next
    foo(anc[i], node2del[i],long[i])
  }
  phylo$edge <- phylo$edge[-ind, ]
  phylo$edge.length <- phylo$edge.length[-ind]
  phylo$Nnode <- phylo$Nnode - n
  ## Now we renumber the nodes that need to be:
  sel <- phylo$edge > min(node2del)
  for (i in which(sel))
    phylo$edge[i] <- phylo$edge[i] - sum(node2del < phylo$edge[i])
  if (!is.null(phylo$node.label))
    phylo$node.label <- phylo$node.label[-(node2del - length(phylo$tip.label))]
  phylo
}