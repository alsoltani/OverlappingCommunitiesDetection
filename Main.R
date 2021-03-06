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

divergence <- function(nu, mu, n){
  
  # Divergence function.
  # ::::::::::::::::::::
  
  idx_mu = which(mu != 0)
  if(length(idx_mu) == n){
    sum(nu[idx_mu] * log(mu[idx_mu]))
  }
  else{
    if(sum(nu[-idx_mu]) != 0){
      -10000
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
  
  # unordered_a : unordered adjacency matrix.
  # When forming the adjacency matrix, the indexes do not match the 
  # vertices' names. Hence the need for an permutation.
  unordered_a = get.adjacency(g, sparse=T)

  permutation_v <- as.numeric(V(g)$name)
  p = as(permutation_v, 'pMatrix')  # Permutation matrix
  inv_p = as(invPerm(permutation_v), 'pMatrix')
  
  # a : unordered adjacency matrix.
  a = inv_p %*% unordered_a %*% p
  n = dim(a)[2]
  
  # unordered_di : unordered inverse degree matrix.
  unordered_di = .sparseDiagonal(n = n, x = 1 / degree(g)) # no NAs
  di = inv_p %*% unordered_di %*% p
  
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
      mu[i, ] = apply(w[s, ], 2, function(r) {weighted.mean(r, w=degree(g)[s])})
    }
    
    # Update the old partition.
    partition[1, ] = partition[2, ]
    
    # Step 2a. Construct the divergence matrix.
    divergence_matrix = matrix(0, nrow = n, ncol = k)
    for (i in 1:n){
      for(j in 1:k){
        divergence_matrix[i,j] = divergence(w[i, ], mu[j, ], n)
      }
    }
    
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

der_overlapping <- function(der_algorithm, g, k){
  
  # Modified DER Algorithm, for overlapping communities' detection.
  # ::::::::::::::::::::
  # Returns a (k, n) matrix, containing overlapping communities C1, ... , Ck as rows.
  
  # Number of vertices.
  n = vcount(g)
  
  # Output of the classical DER algorithm : partitions, mu.
  res = der_algorithm
  partition = res[1][[1]]
  mu = res[2][[1]]
  
  # pi : stationary measure of the random walk.
  # As in the DER algorithm, 
  permutation_v <- as.numeric(V(g)$name)
  p = as(permutation_v, 'pMatrix')
  inv_p = as(invPerm(permutation_v), 'pMatrix')
  
  unordered_pi = degree(g) / norm(data.matrix(degree(g)), "1")
  pi = inv_p * unordered_pi
  
  # m : (n, k) matrix of probabilities that the walk started at P_s, given that it finished in i.
  m = matrix(F, nrow=n, ncol=k)
  for(j in 1:k){
    s = which(partition == j)
    m[, j] = mu[j, ] / pi * sum(pi[s])
  }
  
  # For each i in V, find the most likely community given i and the associated value 1/2 * m_i[s_i].
  most_likely_community = vector(length=n)
  ms = vector(length=n)
  
  for (i in 1:n){
    most_likely_community[i] = which(m[i, ] == max(m[i, ]))[1]
    ms[i] = m[i, most_likely_community[i]]
  }
  
  # Create overlapping communities.
  overlapping_communities = matrix(NA, nrow=k, ncol=n)
  for (t in 1:k){
    
    # Indexes verifiying the criterion.
    indexes = which(m[, t] >= 0.5 * ms)
    overlapping_communities[t, indexes] = indexes
  }
  
  return(overlapping_communities)
}

# 3. Test on artificial data.
# ---------------------

# Read community file.
community_files = list.files(pattern = "community.dat", recursive = T)
community = community_files[1]
communities = read.csv(community, sep = "\t", header = F)[, 2]

# Find the number of communities.
n_communities = max(as.numeric(
  unlist(lapply(communities, 
                function(x) {
                  unlist(strsplit(toString(x), " "))}
  ))))

# Read network file.
network_files = list.files(pattern = "network.dat", recursive = T)
network = network_files[1]
edgelist = read.csv(network, sep = "\t", header = F)[, 1:2]
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

# DER Algorithms.
test <- der_algorithm(igraph_data, 5, n_communities)
test_overlapping <- der_overlapping(test, igraph_data, n_communities)

# Save results.
network = network_files[1]
file_name = unlist(strsplit(network, "[/]"))[3]
write.table(test_overlapping, paste("Data/DERMemberships/", file_name, ".dat", sep=""), 
            na="", col.names = F, row.names = F)

# Plot results.
for(i in V(igraph_data)){
  V(igraph_data)[i]$color = test[[1]][i]
  V(igraph_data)[i]$community = test[[1]][i]
}

plot(igraph_data)
