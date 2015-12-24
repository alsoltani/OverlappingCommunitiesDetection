# Statistical Analysis of Network Data Project, ENSAE 2015.
# Robin Vogel, Clement Puppo, Alain Soltani.
# ---------------------

library(expm)
library(igraph)
library(Matrix)
setwd("/home/Saved_Documents/Documents/Cours/M2_ENSAE/Network_stat/project/OverlappingCommunitiesDetection/")

# 1. DER Algorithm.
# ---------------------

divergence <- function(nu, mu){
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
  
  # a : adjacency matrix.
  a = get.adjacency(g, sparse=T)
  m = dim(a)[1]
  n = dim(a)[2]
  
  # di : inverse degree matrix.
  di = diag(1 / degree(g))

  # t : lists of transition matrices, from 1 to l.
  t = Matrix(di %*% a, sparse = T)
  w = Matrix(0, nrow=n, ncol=n, sparse = T)
  temp = Matrix(diag(n), sparse = T)
  
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
  
  #print(c(partition[1, ] != partition[2, ]))
  
  
  # mu_s : distribution  of the random walk started from pi_s.
  mu = matrix(0, nrow = k, ncol = n)
  t = 0
  # While the two partition differs
  while(Reduce("|", partition[1, ] != partition[2, ])){
    t = t+1
    # Step 1.
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
        divergence_matrix[i,j] = divergence(w[i,],mu[j,])
      }
      #divergence_matrix[i, ] = apply(mu, 1, function(r) {divergence(w[i, ], r)})
    }
    
    #print(w[n, ])
    #print(mu[,1])
    
    # Step 2b.
    for (i in 1:n){
      partition[2, i] = which(divergence_matrix[i, ] == max(divergence_matrix[i, ]))[1]
    }
    print(t)
    #print(c(partition[1, ] != partition[2, ]))
    
  }
  print(t)
  return(list(partition[2, ], mu))
}

# Read in edges information.
edgelist = read.csv("./flickrEdges.txt", sep = " ", skip = 3, nrows = 1000)[, 1:2]
colnames(edgelist) = c("Node Id", "Node Id")

# Because the vertex IDs in this dataset are numbers, 
# we make sure igraph knows these should be treated as characters. 
edgelist[, 1] = as.character(edgelist[, 1]) 
edgelist[, 2] = as.character(edgelist[, 2])

# Igraph needs the edgelist to be in matrix format.
edgelist = as.matrix(edgelist)
igraph_data = graph.edgelist(edgelist[, 1:2], directed = F)

# Simplify to remove duplications and from-self-to-self loops
igraph_data <- simplify(igraph_data, remove.multiple = TRUE, remove.loops = TRUE)

# DER Algorithm.
test <- der_algorithm(igraph_data, 5, 3)

# Karate Club test

G <- read.graph("Benchmarks/karate/karate.gml",format = "gml")
test <- der_algorithm(G, 5, 3)

for(i in V(G)){
  V(G)[i]$color = test[[1]][i]
  V(G)[i]$community = test[[1]][i]
}

plot(G)

