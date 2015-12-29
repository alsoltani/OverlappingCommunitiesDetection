# Statistical Analysis of Network Data Project, ENSAE 2015.
# Robin Vogel, Clement Puppo, Alain Soltani.
# ---------------------

library(expm)
library(igraph)
library(Matrix)
#setwd("/home/Saved_Documents/Documents/Cours/M2_ENSAE/Network_stat/project/OverlappingCommunitiesDetection/")
setwd("/Users/alain/Documents/ENSAE/Statistical\ Analysis\ of\ Network\ Data/Project")

# 1. DER Algorithm.
# ---------------------

divergence <- function(nu, mu){
  
  # Divergence function.
  # ::::::::::::::::::::
  
  idx_mu = which(mu != 0)
  if(length(idx_mu)==0){
    sum(nu[idx_mu] * log(mu[idx_mu]))
  }
  else{
    if(sum(nu[-idx_mu]) != 0){
      -9000
    }
    else{
      sum(nu[idx_mu] * log(mu[idx_mu]))
    }
  }
}

der_algorithm <- function(g, L, k){
  
  # DER Algorithm, for non-overlapping communities' detection.
  # ::::::::::::::::::::
  # Returns a (1, n) matrix, containing communities for all vertices,
  # and a (k, n) containing mu's for each community.
  
  # a : adjacency matrix.
  a = get.adjacency(g, sparse=T)
  m = dim(a)[1]
  n = dim(a)[2]
  
  # di : inverse degree matrix.
<<<<<<< HEAD
  #di = diag(1 / degree(g))
  di = .sparseDiagonal(n = n,x = 1/degree(g)) # no NAs
=======
  di = diag(1 / degree(g))
  
>>>>>>> 37d39e2dda20f5b6726b84adb20ad6960b19ccf0
  # t : lists of transition matrices, from 1 to l.
  t = Matrix(di %*% a, sparse = T)
  w = Matrix(0, nrow=n, ncol=n, sparse = T)
  temp = Matrix(.sparseDiagonal(n), sparse = T)
  
  for (i in 1:L){
    temp = temp %*% t
    w = w + temp
  }
  
  # w_i : distribution corresponding to the average of
  # the empirical measures of sequences x that start at i.
  w = w / L
  
  # partition : first row contains the old partition,
  # second one the newly formed one.
  
  partition = matrix(0, nrow=2, ncol=n)
  partition[1, ] = 0
  partition[2, ] = sample(c(rep(seq(k),  n %/% k), seq(n %% k)))
  
  # mu_s : distribution of the random walk started from pi_s.
  mu = matrix(0, nrow = k, ncol = n)
  t = 0

  # While the two partition differs :
  while(Reduce("|", partition[1, ] != partition[2, ])){
    t = t+1
    
    # Step 1. Form mu_s.
    for(i in 1:k){
      s = which(partition[2, ] == i)
      mu[i,] = apply(w[s, ], 2, function(r) {weighted.mean(r, w=degree(g)[s])})
    }
    
    # Update the old partition.
    partition[1, ] = partition[2, ]
    
    # Step 2a. Construct the divergence matrix.
    divergence_matrix = matrix(0, nrow = n, ncol = k)
    for (i in 1:n){
      for(j in 1:k){

        divergence_matrix[i,j] = divergence(w[i, ], mu[j, ])
      }
      #divergence_matrix[i, ] = apply(mu, 1, function(r) {divergence(w[i, ], r)})
    }
    
    #print(w[n, ])
    #print(mu[,1])
    
    # Step 2b. Form the new partition.
    for (i in 1:n){
      partition[2, i] = which(divergence_matrix[i, ] == max(divergence_matrix[i, ]))[1]
    }
    cat("DER Algorithm : Iteration", t, "\n")
  }
  return(list(partition[2, ], mu))
}

# 2. DER Algorithm for Overlapping Communities.
# ---------------------

der_overlapping <- function(g, L, k){
  
  # Modified DER Algorithm, for overlapping communities' detection.
  # ::::::::::::::::::::
  # Returns a (k, n) matrix, containing overlapping communities C1, ... , Ck as rows.
  
  # Number of vertices.
  n = vcount(g)
  
  # Output of the classical DER algorithm : partitions, mu.
  res = der_algorithm(g, L, k)
  partition = res[1][[1]]
  mu = res[2][[1]]
  
  # pi : stationary measure of the random walk.
  pi = (1 / degree(g)) / norm(data.matrix(1 / degree(g)), "1")
  
  # m : (n, k) matrix of probabilities that the walk started at P_s, given that it finished in i.
  m = matrix(F, nrow=n, ncol=k)
  for(j in 1:k){
    s = which(partition == j)
    m[, j] = mu[j, ] / pi * mean(pi[s])
  }
  
  # For each i in V, find the most likely community given i and the associated value 1/2 * m_i[s_i].
  most_likely_community = vector(length=n)
  ms = vector(length=n)
  
  for (i in 1:n){
    most_likely_community[i] = which(m[i, ] == max(m[i, ]))[1]
    ms[i] = m[i, most_likely_community[i]]
  }
  
  # Create overlapping communities.
  overlapping_communities = matrix(F, nrow=k, ncol=n)
  for (t in 1:k){
    
    # Indexes verifiying the criterion.
    indexes = which(m[, t] >= 0.5 * ms)
    overlapping_communities[t, indexes] = t
  }
  return(overlapping_communities)
}

# 3. Test on Flickr data.
# ---------------------

# Read in edges information.
edgelist = read.csv("./Amazon_data/com-amazon.ungraph.txt", sep = "\t", skip = 4)[, 1:2]
                    #, nrows = 1000
colnames(edgelist) = c("Node Id", "Node Id")

# Because the vertex IDs in this dataset are numbers, 
# we make sure igraph knows these should be treated as characters. 
edgelist[, 1] = as.character(edgelist[, 1]) 
edgelist[, 2] = as.character(edgelist[, 2])

# Igraph needs the edgelist to be in matrix format.
edgelist = as.matrix(edgelist)
igraph_data = graph.edgelist(edgelist[, 1:2], directed = F)

# Simplify to remove duplications and from-self-to-self loops
igraph_data <- simplify(igraph_data, remove.multiple = T, remove.loops = T)

# DER Algorithm.
test <- der_algorithm(igraph_data, 5, 3)
test_overlapping <- der_overlapping(igraph_data, 5, 3)

###########################################################################
# Karate Club test
G <- read.graph("Benchmarks/karate/karate.gml", format = "gml")
test <- der_algorithm(G, 5, 3)

for(i in V(G)){
  V(G)[i]$color = test[[1]][i]
  V(G)[i]$community = test[[1]][i]
}

plot(G)
