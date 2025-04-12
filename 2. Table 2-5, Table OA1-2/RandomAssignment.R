#=========================================================================#
#===== Benchmark: Random Assignments + Difference-in-Means Estimator =====#
#=========================================================================#

#============================================================#
#===== Note: Permutation over time does not make sense. =====#
#=====       Should probably permute across units.      =====#
#============================================================#

#-----------------------------------------#
#--- blockTools requires R/4.4.0 ---------#
#--- gurobi requires R/4.0.2 -------------#
#--- note these different requirements ---#
#-----------------------------------------#

#--- Last Change: 20250330 10PM coded up iterative matching
#--- Last Change: 20250326 10PM: coded up Genetic Algorithm
#--- Last Change: 20250321 10PM: coded up the simulation with k-NN
#--- Last Change: 20250316 2PM: set the weights in the estimand (beta.vector) to be uniform weights 1/N.Regions
#--- Last Change: 20250310 11AM: coded up balanced blocking as a different Randomization #2 method
#--- Last Change: 20250218 8AM: coded up the simulation with regression adjustment
#--- Last Change: 20250216 1PM: fixed a few bugs
#--- Last Change: 20250130 10AM: coded up the simulation with Randomization #2 using K means
#--- Last Change: 20221204 5PM: coded up the simulation with Randomization #1
#--- Last Change: 20221205 11AM: made it parallel

# library(blockTools)
library(stats)
library(GA)
# library(GenSA)
library(reticulate)
library(stringr)

# repetition.RANDOM.SEED = 577
repetition.RANDOM.SEED = as.numeric(Sys.getenv("SGE_TASK_ID"))

# start_time <- Sys.time()


#=====Define variables=====#
#===Basics
#Basics --- Scalars
N.Regions = 15
T.naught = 25
T.prime = 20
T.total = 30
r.ob.covariates.dim = 7
F.unob.covariates.dim = 11

#===Model Primitives
#Population weights --- Vectors
beta.vector = c()
#Constants --- Vectors
delta.N.vector = c()
upsilon.I.vector = c()
#Covariates --- Matrices
Z.ob.covariates.matrix = matrix(NA, nrow = r.ob.covariates.dim, ncol = N.Regions)
mu.unob.covariates.matrix = matrix(NA, nrow = F.unob.covariates.dim, ncol = N.Regions)
#Coefficients --- Matrices
theta.ob.N.matrix = matrix(NA, nrow = r.ob.covariates.dim, ncol = T.total)
gamma.ob.I.matrix = matrix(NA, nrow = r.ob.covariates.dim, ncol = T.total)
lambda.unob.N.matrix = matrix(NA, nrow = F.unob.covariates.dim, ncol = T.total)
eta.unob.I.matrix = matrix(NA, nrow = F.unob.covariates.dim, ncol = T.total)
#Noises --- Matrices
epsilon.N.matrix = matrix(NA, nrow = N.Regions, ncol = T.total)
xi.I.matrix = matrix(NA, nrow = N.Regions, ncol = T.total)

#===Potential Outcomes
#Potential Outcomes in hindsight --- Matrices
Y.N.matrix = matrix(NA, nrow = N.Regions, ncol = T.total)
Y.I.matrix = matrix(NA, nrow = N.Regions, ncol = T.total)

#===Define the range of parameters
#Constants --- Vectors
#Range --- [0,20]
range.intercept.max = 20 #ChangeThis
#Covariates --- Matrices
#Range --- [0,1]
range.covariates.max = 1 #ChangeThis
#Coefficients --- Matrices
#Range --- [0,10]
range.coefficients.max = 10 #ChangeThis
#Noises --- Matrices
#Normal Distribution ~ N(0,1)
noise.variance = 1 #ChangeThis

Generate_Model_Primitives <- function(range.intercept.max_ = range.intercept.max,
                                      range.covariates.max_ = range.covariates.max,
                                      range.coefficients.max_ = range.coefficients.max,
                                      noise.variance_ = noise.variance,
                                      random.seed_)
  # random.seed.FixedPrimitives = 123456,
  # random.seed.RandomParts = 123456)
{
  set.seed(random.seed_)
  #===Initialize Model Primitives Again
  #Constants --- Vectors
  delta.N.vector_ = c()
  upsilon.I.vector_ = c()
  #Covariates --- Matrices
  Z.ob.covariates.matrix_ = matrix(NA, nrow = r.ob.covariates.dim, ncol = N.Regions)
  mu.unob.covariates.matrix_ = matrix(NA, nrow = F.unob.covariates.dim, ncol = N.Regions)
  #Coefficients --- Matrices
  theta.ob.N.matrix_ = matrix(NA, nrow = r.ob.covariates.dim, ncol = T.total)
  gamma.ob.I.matrix_ = matrix(NA, nrow = r.ob.covariates.dim, ncol = T.total)
  lambda.unob.N.matrix_ = matrix(NA, nrow = F.unob.covariates.dim, ncol = T.total)
  eta.unob.I.matrix_ = matrix(NA, nrow = F.unob.covariates.dim, ncol = T.total)
  #Noises --- Matrices
  epsilon.N.matrix_ = matrix(NA, nrow = N.Regions, ncol = T.total)
  xi.I.matrix_ = matrix(NA, nrow = N.Regions, ncol = T.total)
  #Potential Outcomes in hindsight --- Matrices
  Y.N.matrix_ = matrix(NA, nrow = N.Regions, ncol = T.total)
  Y.I.matrix_ = matrix(NA, nrow = N.Regions, ncol = T.total)
  
  #=====Generate Random Values=====#
  #===Basics
  #Population weights --- Vectors
  # beta.temp = runif(N.Regions)
  # beta.vector_ = beta.temp / sum(beta.temp)
  beta.vector_ = rep(1/N.Regions, N.Regions)
  #===Model Primitives
  #Constants --- Vectors
  #Range --- [0,20] #range.intercept.max_ = 20
  delta.N.vector_ = sort(c(range.intercept.max_ * runif(T.naught), range.intercept.max_ * runif(T.total - T.naught)))
  upsilon.I.vector_ = c(rep(NA, T.naught), sort(range.intercept.max_ * runif(T.total - T.naught)))
  #Covariates --- Matrices
  #Range --- [0,1] # range.covariates.max_ = 1
  for(j in 1:N.Regions)
  {
    z.temp = range.covariates.max_ * runif(r.ob.covariates.dim)
    Z.ob.covariates.matrix_[,j] = z.temp
  }
  colnames(Z.ob.covariates.matrix_) = c(1:N.Regions)
  for(j in 1:N.Regions)
  {
    mu.temp = range.covariates.max_ * runif(F.unob.covariates.dim)
    mu.unob.covariates.matrix_[,j] = mu.temp
  }
  colnames(mu.unob.covariates.matrix_) = c(1:N.Regions)
  #Coefficients --- Matrices
  #Range --- [0,10] # range.coefficients.max_ = 10
  for(t in 1:T.naught)
  {
    theta.temp = range.coefficients.max_ * runif(r.ob.covariates.dim)
    theta.ob.N.matrix_[,t] = theta.temp
  }
  for(t in (T.naught+1):T.total)
  {
    theta.temp = range.coefficients.max_ * runif(r.ob.covariates.dim)
    theta.ob.N.matrix_[,t] = theta.temp
  }
  colnames(theta.ob.N.matrix_) = c(1:T.total)
  
  for(t in 1:T.naught)
  {
    gamma.temp = rep(NA, r.ob.covariates.dim)
    gamma.ob.I.matrix_[,t] = gamma.temp
  }
  for(t in (T.naught+1):T.total)
  {
    gamma.temp = range.coefficients.max_ * runif(r.ob.covariates.dim)
    gamma.ob.I.matrix_[,t] = gamma.temp
  }
  colnames(gamma.ob.I.matrix_) = c(1:T.total)
  
  for(t in 1:T.naught)
  {
    lambda.temp = range.coefficients.max_ * runif(F.unob.covariates.dim)
    lambda.unob.N.matrix_[,t] = lambda.temp
  }
  for(t in (T.naught+1):T.total)
  {
    lambda.temp = range.coefficients.max_ * runif(F.unob.covariates.dim)
    lambda.unob.N.matrix_[,t] = lambda.temp
  }
  colnames(lambda.unob.N.matrix_) = c(1:T.total)
  
  for(t in 1:T.naught)
  {
    eta.temp = rep(NA, F.unob.covariates.dim)
    eta.unob.I.matrix_[,t] = eta.temp
  }
  for(t in (T.naught+1):T.total)
  {
    eta.temp = range.coefficients.max_ * runif(F.unob.covariates.dim)
    eta.unob.I.matrix_[,t] = eta.temp
  }
  colnames(eta.unob.I.matrix_) = c(1:T.total)
  
  #Noises --- Matrices
  #Normal Distribution ~ N(0,1) #noise.variance_ = 1
  for(t in 1:T.naught)
  {
    epsilon.temp = rnorm(N.Regions, 0, noise.variance_)
    epsilon.N.matrix_[,t] = epsilon.temp
  }
  for(t in (T.naught+1):T.total)
  {
    epsilon.temp = rnorm(N.Regions, 0, noise.variance_)
    epsilon.N.matrix_[,t] = epsilon.temp
  }
  colnames(epsilon.N.matrix_) = c(1:T.total)
  
  for(t in 1:T.naught)
  {
    K.temp = rep(NA, N.Regions)
    xi.I.matrix_[,t] = K.temp
  }
  for(t in (T.naught+1):T.total)
  {
    K.temp = rnorm(N.Regions, 0, noise.variance_)
    xi.I.matrix_[,t] = K.temp
  }
  colnames(xi.I.matrix_) = c(1:T.total)
  
  #=====Generate Potential Outcomes=====#
  #Potential Outcomes in hindsight --- Matrices
  #A Factor Model
  for(j in 1:N.Regions)
  {
    for(t in 1:T.total)
    {
      Y.N.matrix_[j,t] = delta.N.vector_[t] + theta.ob.N.matrix_[,t] %*% Z.ob.covariates.matrix_[,j] + lambda.unob.N.matrix_[,t] %*% mu.unob.covariates.matrix_[,j] + epsilon.N.matrix_[j,t]
    }
  }
  
  for(j in 1:N.Regions)
  {
    for(t in 1:T.total)
    {
      Y.I.matrix_[j,t] = upsilon.I.vector_[t] + gamma.ob.I.matrix_[,t] %*% Z.ob.covariates.matrix_[,j] + eta.unob.I.matrix_[,t] %*% mu.unob.covariates.matrix_[,j] + xi.I.matrix_[j,t]
    }
  }
  
  returned = list(beta.vector = beta.vector_, 
                  delta.N.vector = delta.N.vector_,
                  upsilon.I.vector = upsilon.I.vector_,
                  Z.ob.covariates.matrix = Z.ob.covariates.matrix_,
                  mu.unob.covariates.matrix = mu.unob.covariates.matrix_,
                  theta.ob.N.matrix = theta.ob.N.matrix_,
                  gamma.ob.I.matrix = gamma.ob.I.matrix_,
                  lambda.unob.N.matrix = lambda.unob.N.matrix_,
                  eta.unob.I.matrix = eta.unob.I.matrix_,
                  epsilon.N.matrix = epsilon.N.matrix_,
                  xi.I.matrix = xi.I.matrix_,
                  Y.N.matrix = Y.N.matrix_,
                  Y.I.matrix = Y.I.matrix_)
  return(returned)
}

permutation.test <- function(pre.intervention.residuals, post.intervention.residuals)
{
  permutation.test.vec = c(abs(pre.intervention.residuals), abs(post.intervention.residuals))
  test.statistic = sum(tail(permutation.test.vec, length(post.intervention.residuals)))
  combinations.matrix = combn(1:length(permutation.test.vec), length(post.intervention.residuals)) #each column is a combination
  permutation.statistics = c()
  for(permu.temp in 1:ncol(combinations.matrix))
  {
    permutation.statistics = c(permutation.statistics, sum(permutation.test.vec[combinations.matrix[,permu.temp]]))
  }
  returned = sum(permutation.statistics >= test.statistic) / ncol(combinations.matrix)
  return(returned)
}

#-------------------------------------------------------------#
#----- Randomization 1: less than a number of treatments -----#
#-------------------------------------------------------------#

random_assignments <- function(random.seed_, N.Regions_ = N.Regions, beta.vector_ = beta.vector, K.cardinality_ = -1)
{
  set.seed(random.seed_)
  
  treated.weights = rep(0, N.Regions_)
  control.weights = rep(0, N.Regions_)
  
  if(K.cardinality_ == -1)
  {
    K.cardinality = floor(N.Regions_/2)
  }
  if(K.cardinality_ >= 1)
  {
    K.cardinality = K.cardinality_
  }
  
  # candidate.partition = list()
  # for(partition.cardinality in 1:K.cardinality)
  # {
  #   candidate.partition = c(candidate.partition, combn(1:N.Regions_, partition.cardinality, simplify = FALSE))
  # }
  # total.candidate.partitions = length(candidate.partition)
  # 
  # random.number = runif(1, min = 0, max = 1)
  # random.index = ceiling(total.candidate.partitions * random.number)
  # 
  # treated.units = candidate.partition[[random.index]] #sort( sample(c(1:N.Regions), size=K.cardinality, replace=F) )
  # control.units = c(1:N.Regions)[-treated.units]
  # treated.weight.value = 1 / length(treated.units)
  # control.weight.value = 1 / length(control.units)
  # treated.weights[treated.units] = treated.weight.value
  # control.weights[control.units] = control.weight.value
  # 
  # treated.weights = treated.weights * beta.vector_ * N.Regions_
  # control.weights = control.weights * beta.vector_ * N.Regions_
  
  treated.units = sample(c(1:N.Regions), K.cardinality)
  control.units = setdiff(c(1:N.Regions), treated.units)
  
  treated.weight.value = 1 / length(treated.units)
  control.weight.value = 1 / length(control.units)
  treated.weights[treated.units] = treated.weight.value
  control.weights[control.units] = control.weight.value
  
  treated.weights = treated.weights * beta.vector_ * N.Regions_
  control.weights = control.weights * beta.vector_ * N.Regions_
  
  return(list(Treatment.weights = round(treated.weights, digits = 8), 
              Control.weights = round(control.weights, digits = 8)))
}


#-----------------------------------------------------#
#----- Randomization 2: stratified randomization -----#
#-----------------------------------------------------#

euclidean_distance <- function(point1, point2)
{
  sqrt(sum((point1 - point2)^2))
}

my_k_means <- function(random.seed_, N.Regions_ = N.Regions, feature.matrix_, min.size_ = 2, max_iterations_ = 100, tolerence_ = 0.001, K.cardinality_ = -1)
{
  if(K.cardinality_ == -1)
  {
    K.cardinality = floor(N.Regions_/2)
  }
  if(K.cardinality_ >= 1)
  {
    K.cardinality = K.cardinality_
  }
  
  feature.matrix_ = apply(feature.matrix_, 2, function(x) (x - mean(x)) / sd(x))
  
  my_data = data.frame(feature.matrix_)
  
  #-----------------------------------------#
  #----- Randomly seed the K centroids -----#
  #----- Then repeat many (100) times  -----#
  #-----------------------------------------#
  my.obj.value.vec = rep(NA, max_iterations_)
  
  for(times.temp in 1:max_iterations_)
  {
    this.seed = random.seed_ * max_iterations_ + times.temp
    set.seed(this.seed)
    
    # Step 1: Initialize centroids randomly
    centroid.labels = sample(N.Regions_, K.cardinality)
    centroids = my_data[centroid.labels, ]
    
    # Initialize variables
    iter = 0
    # Repeat until convergence or max iterations
    while(iter < max_iterations_) 
    {
      cluster_assignments = integer(N.Regions_)
      
      # Step 2: Assign each point to the nearest centroid
      count.size = rep(0, K.cardinality) #Initialize a vector to count how many are there in each cluster
      distances = matrix(data = NA, nrow = N.Regions_, ncol = K.cardinality)  # Initialize a matrix to store distances
      # First, calculate the distance of each point to each centroid
      for(i in 1:N.Regions_)
      {
        for(j in 1:K.cardinality)
        {
          distances[i,j] = euclidean_distance(my_data[i, ], centroids[j, ])
        }
      }
      define_my_max = round(max(distances), digits = 0) + 2
      
      # Second, assign each point to a centroid, make sure that each centroid gets min.size_ many points
      while(min(distances) < define_my_max)
      {
        indices = which(distances == min(distances), arr.ind = TRUE)
        first.indices = indices[1,]
        
        cluster_assignments[first.indices[1]] = first.indices[2]  # Assign to the nearest centroid
        distances[first.indices[1], ] = define_my_max
        count.size[first.indices[2]] = count.size[first.indices[2]] + 1
        if(count.size[first.indices[2]] >= min.size_)
        {
          distances[, first.indices[2]] = define_my_max
        }
      }
      
      #Third, assign all the remaining points to a centroid; need to calculate the distances matrix again
      distances = matrix(data = NA, nrow = N.Regions_, ncol = K.cardinality)
      for(i in 1:N.Regions_)
      {
        for(j in 1:K.cardinality)
        {
          distances[i,j] = euclidean_distance(my_data[i, ], centroids[j, ])
        }
      }
      for(i in which(cluster_assignments == 0))
      {
        cluster_assignments[i] = which.min(distances[i, ])
      }
      
      # Step 3: Update centroids as the mean of assigned points
      prev_centroids = centroids
      for (j in 1:K.cardinality)
      {
        cluster_points = my_data[cluster_assignments == j, ]
        if(nrow(cluster_points) > 0)
        {
          centroids[j, ] = colMeans(cluster_points)
        }
      }
      
      # Check for convergence (if centroids don't change much)
      if(all(abs(centroids - prev_centroids) < tolerence_))
      {
        break
      }
      
      iter = iter + 1
    }
    
    # Now calculate the objective function (sum of squared distances from each point to its centroid)
    for(j in 1:K.cardinality)
    {
      cluster_points = my_data[cluster_assignments == j, ]
      if(nrow(cluster_points) > 0)
      {
        centroids[j, ] = colMeans(cluster_points)
      }
    }
    obj.distances.vec = rep(NA, N.Regions_)
    for(i in 1:N.Regions_)
    {
      cluster.temp = cluster_assignments[i]
      obj.distances.vec[i] = euclidean_distance(my_data[i, ], centroids[cluster.temp, ])
    }
    my.obj.value.vec[times.temp] = obj.distances.vec %*% obj.distances.vec
  }
  
  this.times.temp = which.min(my.obj.value.vec)
  
  this.seed = random.seed_ * max_iterations_ + this.times.temp
  
  #----------------------------------------------------------#
  #----- Now, use the best seed and run the codes again -----#
  #----------------------------------------------------------#
  
  set.seed(this.seed)
  
  # Step 1: Initialize centroids randomly
  centroid.labels = sample(N.Regions_, K.cardinality)
  centroids = my_data[centroid.labels, ]
  
  # Initialize variables
  iter = 0
  # Repeat until convergence or max iterations
  while(iter < max_iterations_) 
  {
    cluster_assignments = integer(N.Regions_)
    
    # Step 2: Assign each point to the nearest centroid
    count.size = rep(0, K.cardinality) #Initialize a vector to count how many are there in each cluster
    distances = matrix(data = NA, nrow = N.Regions_, ncol = K.cardinality)  # Initialize a matrix to store distances
    # First, calculate the distance of each point to each centroid
    for(i in 1:N.Regions_)
    {
      for(j in 1:K.cardinality)
      {
        distances[i,j] = euclidean_distance(my_data[i, ], centroids[j, ])
      }
    }
    define_my_max = round(max(distances), digits = 0) + 2
    
    # Second, assign each point to a centroid, make sure that each centroid gets min.size_ many points
    while(min(distances) < define_my_max)
    {
      indices = which(distances == min(distances), arr.ind = TRUE)
      first.indices = indices[1,]
      
      cluster_assignments[first.indices[1]] = first.indices[2]  # Assign to the nearest centroid
      distances[first.indices[1], ] = define_my_max
      count.size[first.indices[2]] = count.size[first.indices[2]] + 1
      if(count.size[first.indices[2]] >= min.size_)
      {
        distances[, first.indices[2]] = define_my_max
      }
    }
    
    #Third, assign all the remaining points to a centroid; need to calculate the distances matrix again
    distances = matrix(data = NA, nrow = N.Regions_, ncol = K.cardinality)
    for(i in 1:N.Regions_)
    {
      for(j in 1:K.cardinality)
      {
        distances[i,j] = euclidean_distance(my_data[i, ], centroids[j, ])
      }
    }
    for(i in which(cluster_assignments == 0))
    {
      cluster_assignments[i] = which.min(distances[i, ])
    }
    
    # Step 3: Update centroids as the mean of assigned points
    prev_centroids = centroids
    for (j in 1:K.cardinality)
    {
      cluster_points = my_data[cluster_assignments == j, ]
      if(nrow(cluster_points) > 0)
      {
        centroids[j, ] = colMeans(cluster_points)
      }
    }
    
    # Check for convergence (if centroids don't change much)
    if(all(abs(centroids - prev_centroids) < tolerence_))
    {
      break
    }
    
    iter = iter + 1
  }
  
  returned = list(clusters = cluster_assignments,
                  centroids = centroids)
  return(returned)
}

next_power_of_2 <- function(n)
{
  return(2^ceiling(log2(n + 1)))
}

numerator_and_denominator <- function(N.Regions, K.cardinality)
{
  if(N.Regions %% K.cardinality == 0)
  {
    numerator_to_python = 1
    denominator_to_python = N.Regions / K.cardinality
  } else
  {
    N.augmented = next_power_of_2(N.Regions)
    if(N.Regions %% K.cardinality == 0)
    {
      numerator_to_python = 1
      denominator_to_python = N.augmented / K.cardinality
    } else
    {
      numerator_to_python = 1
      denominator_to_python = floor(N.augmented / K.cardinality)
    }
  }
  return(list(num = numerator_to_python,
              denom = denominator_to_python))
}

my_IterMatching_blocking <- function(random.seed_, N.Regions_ = N.Regions, feature.matrix_, K.cardinality_ = -1)
{
  set.seed(random.seed_)
  
  if(K.cardinality_ == -1)
  {
    K.cardinality = floor(N.Regions_/2)
  }
  if(K.cardinality_ >= 1)
  {
    K.cardinality = K.cardinality_
  }
  
  feature.matrix_ = apply(feature.matrix_, 2, function(x) (x - mean(x)) / sd(x))
  
  #--- the matrix is converted by column
  # mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, byrow = TRUE)
  # mat_string <- paste(as.vector(mat), collapse = ",")
  # mat_string
  string_to_python = paste(as.vector(round(feature.matrix_, digits = 4)), collapse = ",")
  
  my_num_denom = numerator_and_denominator(N.Regions_, K.cardinality)
  numerator_to_python = my_num_denom$num
  denominator_to_python = my_num_denom$denom
  
  python_output = system2("python", args = c("rerand.py", N.Regions_, string_to_python, numerator_to_python, denominator_to_python), stdout = TRUE)
  output_string = python_output[3]
  cleaned_string = gsub("np\\.int64", "", output_string)
  cleaned_string = gsub("[()]", "", cleaned_string)
  
  # Remove the outer brackets and split by the '], [' pattern
  cleaned_string_modified = substring(cleaned_string, 3, nchar(cleaned_string) - 2)
  cleaned_string_modified = str_sub(cleaned_string_modified)
  groups = str_split(cleaned_string_modified, "\\], \\[")
  
  # Convert each group to a numeric vector
  vector_list = lapply(groups[[1]], function(group) {
    as.integer(str_split(group, ",")[[1]])
  })
  
  length_vector_list = length(vector_list)
  vector_list_lengthes = c()
  this_center = list()
  if(length_vector_list > K.cardinality) # need to merge clusters
  {
    how_many_merges = length_vector_list - K.cardinality
    for(this.cluster in 1:length_vector_list)
    {
      vector_list_lengthes = c(vector_list_lengthes, length(vector_list[[this.cluster]]))
      if(length(vector_list[[this.cluster]]) == 1)
      {
        this_center[[this.cluster]] = feature.matrix_[vector_list[[this.cluster]]+1,]
      } else
      {
        this_center[[this.cluster]] = colMeans(feature.matrix_[vector_list[[this.cluster]]+1,])
      }
    }
    for(this.cluster in 1:how_many_merges)
    {
      small_to_large_index = order(vector_list_lengthes)[this.cluster]
      
      target_vector = this_center[[small_to_large_index]]
      distances_to_this = sapply(this_center, function(vec) euclidean_distance(vec, target_vector))
      
      sorted_indices = order(distances_to_this)
      sorted_indices_subtracted = setdiff(sorted_indices, order(vector_list_lengthes)[1:how_many_merges])
      closest_index = sorted_indices_subtracted[1]
      vector_list[[closest_index]] = c(vector_list[[closest_index]], vector_list[[small_to_large_index]])
      vector_list[[small_to_large_index]] = integer(0)
    }
    
    cluster_assignments = rep(NA, N.Regions_)
    for(this.cluster in 1:length_vector_list)
    {
      if(length(vector_list[[this.cluster]]) != 0)
      {
        cluster_assignments[vector_list[[this.cluster]] + 1] = this.cluster
      }
    }
    for(this.cluster in 1:K.cardinality)
    {
      if(!any(cluster_assignments == this.cluster))
      {
        max_indices = which(cluster_assignments == max(cluster_assignments))
        cluster_assignments[max_indices] = this.cluster
      }
    }
  } else
  {
    cluster_assignments = rep(NA, N.Regions_)
    for(this.cluster in 1:length_vector_list)
    {
      if(length(vector_list[[this.cluster]]) != 0)
      {
        cluster_assignments[vector_list[[this.cluster]] + 1] = this.cluster
      }
    }
  }
  
  returned = list(clusters = cluster_assignments)
  return(returned)
}

my_MinMaxDiameter_blocking <- function(random.seed_, N.Regions_ = N.Regions, feature.matrix_, K.cardinality_ = -1, min.units_ = 2)
{
  set.seed(random.seed_)
  
  if(K.cardinality_ == -1)
  {
    K.cardinality = floor(N.Regions_/2)
  }
  if(K.cardinality_ >= 1)
  {
    K.cardinality = K.cardinality_
  }
  
  feature.matrix_ = apply(feature.matrix_, 2, function(x) (x - mean(x)) / sd(x))
  my_data = data.frame(feature.matrix_)
  
  distance_matrix = as.matrix(dist(my_data, method = "euclidean"))
  
  #-------------------------#
  #--- Genetic Algorithm ---#
  #-------------------------#
  
  # Define the objective function
  objective_function_GA <- function(vec_, GA.min.units_ = min.units_)
  {
    vec_ = as.integer(vec_)
    max.obj.value = max(distance_matrix) + 100
    element_counts = table(vec_)
    num_unique = length(unique(vec_))
    
    #--- Calculate the objective function
    if(any(element_counts < GA.min.units_) | (num_unique < K.cardinality)) # Check if any element appears less than min.units_ times, and if there are less than K.cardinality blocks
    {
      returned = - max.obj.value
    } else # If constraint is satisfied
    {
      max.distance.vec = rep(NA, num_unique)
      for(this.temp in 1:num_unique)
      {
        this.label = unique(vec_)[this.temp]
        this.indices = which(vec_ == this.label)
        max.distance.vec[this.temp] = max(distance_matrix[this.indices, this.indices])
      }
      returned = - max(max.distance.vec)
    }
    
    return(returned)
  }
  
  # Define an alternative function (that allows (num_unique < K.cardinality) to happen)
  objective_function_GA_alternative <- function(vec_, GA.min.units_ = min.units_)
  {
    vec_ = as.integer(vec_)
    max.obj.value = max(distance_matrix) + 100
    element_counts = table(vec_)
    num_unique = length(unique(vec_))
    
    #--- Calculate the objective function
    if(any(element_counts < GA.min.units_)) # Check if any element appears less than min.units_ times
    {
      returned = - max.obj.value
    } else # If constraint is satisfied
    {
      max.distance.vec = rep(NA, num_unique)
      for(this.temp in 1:num_unique)
      {
        this.label = unique(vec_)[this.temp]
        this.indices = which(vec_ == this.label)
        max.distance.vec[this.temp] = max(distance_matrix[this.indices, this.indices])
      }
      returned = - max(max.distance.vec)
    }
    
    return(returned)
  }
  
  popSize_ = 100000
  
  # Run the Genetic Algorithm (GA)
  ga_result = ga(
    type = "real-valued",                        # Real-valued parameters
    fitness = objective_function_GA_alternative,             # Objective function
    lower = rep(1, N.Regions_),                  # Lower bound for each dimension
    upper = rep(K.cardinality+0.1, N.Regions_),  # Upper bound for each dimension
    popSize = popSize_,                          # Population size
    maxiter = 100,                               # Maximum number of generations
    run = 50,                                    # Number of iterations without improvement before stopping
    pmutation = 0.5,                             # Mutation probability
    # initialPop = initial_population            # Set the custom initial population
    monitor = FALSE                              # Display progress
  )
  cluster_assignments = as.integer(ga_result@solution[1,])
  cluster_fitnessValue = - ga_result@fitnessValue
  
  # # Run the Genetic Algorithm (GA)
  # ga_result_alternative = ga(
  #   type = "real-valued",                        # Real-valued parameters
  #   fitness = objective_function_GA_alternative, # Objective function
  #   lower = rep(1, N.Regions_),                  # Lower bound for each dimension
  #   upper = rep(K.cardinality+0.1, N.Regions_),  # Upper bound for each dimension
  #   popSize = popSize_,                          # Population size
  #   maxiter = 100,                               # Maximum number of generations
  #   run = 50,                                    # Number of iterations without improvement before stopping
  #   pmutation = 0.5,                             # Mutation probability
  #   # initialPop = initial_population            # Set the custom initial population
  #   monitor = FALSE                              # Display progress
  # )
  # cluster_assignments_alternative = as.integer(ga_result_alternative@solution[1,])
  # cluster_fitnessValue_alternative = - ga_result_alternative@fitnessValue
  # 
  # if(cluster_fitnessValue_alternative < cluster_fitnessValue)
  # {
  #   # Run the Genetic Algorithm (GA)
  #   ga_result = ga(
  #     type = "real-valued",                        # Real-valued parameters
  #     fitness = objective_function_GA_alternative,             # Objective function
  #     lower = rep(1, N.Regions_),                  # Lower bound for each dimension
  #     upper = rep(K.cardinality+0.1, N.Regions_),  # Upper bound for each dimension
  #     popSize = popSize_*10,                       # Population size
  #     maxiter = 100,                               # Maximum number of generations
  #     run = 50,                                    # Number of iterations without improvement before stopping
  #     pmutation = 0.5,                             # Mutation probability
  #     # initialPop = initial_population            # Set the custom initial population
  #     monitor = FALSE                              # Display progress
  #   )
  #   cluster_assignments = as.integer(ga_result@solution[1,])
  #   cluster_fitnessValue = - ga_result@fitnessValue
  # }
  
  return(list(clusters = cluster_assignments,
              fitnessValue = cluster_fitnessValue))
}

stratified_assignments <- function(random.seed_, N.Regions_ = N.Regions, beta.vector_ = beta.vector, feature.matrix_, min.units_ = 2, K.cardinality_ = -1)
{
  set.seed(random.seed_)
  
  treated.weights = rep(0, N.Regions_)
  control.weights = rep(0, N.Regions_)
  
  if(K.cardinality_ == -1)
  {
    K.cardinality = floor(N.Regions_/2)
  }
  if(K.cardinality_ >= 1)
  {
    K.cardinality = K.cardinality_
  }
  
  if(K.cardinality >= 2)
  {
    # my.blocks = my_blockTools_blocking(random.seed_ = random.seed_, N.Regions_ = N.Regions_, feature.matrix_ = feature.matrix_, K.cardinality_ = K.cardinality)
    # my.blocks = my_k_means(random.seed_ = random.seed_, N.Regions_ = N.Regions_, feature.matrix_ = feature.matrix_, K.cardinality_ = K.cardinality)
    
    my.blocks = my_MinMaxDiameter_blocking(random.seed_ = random.seed_, N.Regions_ = N.Regions_, feature.matrix_ = feature.matrix_, K.cardinality_ = K.cardinality)
    # my.blocks = my_IterMatching_blocking(random.seed_ = random.seed_, N.Regions_ = N.Regions_, feature.matrix_ = feature.matrix_, K.cardinality_ = K.cardinality)
    
    fitness_ = my.blocks$fitnessValue
    
    labels = my.blocks$clusters
    num_unique = length(unique(labels))
    
    if(num_unique < K.cardinality)
    {
      K.cardinality = num_unique
    }
    
    for(num.clusters.temp in 1:K.cardinality)
    {
      this.cluster.temp = which(labels == unique(labels)[num.clusters.temp])
      this.cluster.treated = sample(this.cluster.temp, 1)
      this.cluster.control = setdiff(this.cluster.temp, this.cluster.treated)
      treated.weights[this.cluster.treated] = beta.vector_[this.cluster.treated] * length(this.cluster.temp) / length(this.cluster.treated)
      control.weights[this.cluster.control] = beta.vector_[this.cluster.control] * length(this.cluster.temp) / length(this.cluster.control)
    }
  }
  if(K.cardinality == 1)
  {
    fitness_ = 0
    
    treated.units = sample(c(1:N.Regions), 1)
    control.units = setdiff(c(1:N.Regions), treated.units)
    
    treated.weight.value = 1 / length(treated.units)
    control.weight.value = 1 / length(control.units)
    treated.weights[treated.units] = treated.weight.value
    control.weights[control.units] = control.weight.value
    
    treated.weights = treated.weights * beta.vector_ * N.Regions_
    control.weights = control.weights * beta.vector_ * N.Regions_
  }
  
  return(list(Treatment.weights = round(treated.weights, digits = 8), 
              Control.weights = round(control.weights, digits = 8),
              fitnessValue = fitness_))
}




#=====================================#
#===== Generate Model Primitives =====#
#=====================================#
model.primitives = Generate_Model_Primitives(random.seed_ = repetition.RANDOM.SEED)
#-----Fetch the outputs-----#
#Population weights --- Vectors
beta.vector = model.primitives$beta.vector
#Constants --- Vectors
delta.N.vector = model.primitives$delta.N.vector
upsilon.I.vector = model.primitives$upsilon.I.vector
#Covariates --- Matrices
Z.ob.covariates.matrix = model.primitives$Z.ob.covariates.matrix
mu.unob.covariates.matrix = model.primitives$mu.unob.covariates.matrix
#Coefficients --- Matrices
theta.ob.N.matrix = model.primitives$theta.ob.N.matrix
gamma.ob.I.matrix = model.primitives$gamma.ob.I.matrix
lambda.unob.N.matrix = model.primitives$lambda.unob.N.matrix
eta.unob.I.matrix = model.primitives$eta.unob.I.matrix
#Noises --- Matrices
epsilon.N.matrix = model.primitives$epsilon.N.matrix
xi.I.matrix = model.primitives$xi.I.matrix
#Potential Outcomes --- Matrices
Y.N.matrix = model.primitives$Y.N.matrix
Y.I.matrix = model.primitives$Y.I.matrix


#------------------------------------------------#
#--- There are several options here:          ---#
#--- whether we use pre-treatment outcomes in ---#
#--- stratified randomization and matching    ---#
#------------------------------------------------#

my.feature.matrix_StA_Matching = cbind(t(Z.ob.covariates.matrix), Y.N.matrix[,1:T.naught])
# my.feature.matrix_StA_Matching = t(Z.ob.covariates.matrix)
my.feature.matrix_StA_Matching = apply(my.feature.matrix_StA_Matching, 2, function(x) (x - mean(x)) / sd(x))



#====================================================#
#===== Generate the basics about the population =====#
#====================================================#

#-----About the population average-----#
population.N.before = t(beta.vector) %*% Y.N.matrix[,1:T.naught]
population.N.after = t(beta.vector) %*% Y.N.matrix[,(T.naught+1):T.total]
population.I.after = t(beta.vector) %*% Y.I.matrix[,(T.naught+1):T.total]

true.ATE = population.I.after - population.N.after

#-----Under the sharp null hypothesis-----#
Y.null.hypothesis.matrix = matrix(NA, nrow = N.Regions, ncol = T.total)
for(j in 1:N.Regions)
{
  for(t in 1:T.total)
  {
    Y.null.hypothesis.matrix[j,t] = delta.N.vector[t] + theta.ob.N.matrix[,t] %*% Z.ob.covariates.matrix[,j] + lambda.unob.N.matrix[,t] %*% mu.unob.covariates.matrix[,j] + xi.I.matrix[j,t]
  }
}
population.null.hypothesis.after = t(beta.vector) %*% Y.null.hypothesis.matrix[,(T.naught+1):T.total]
null.hypothesis.true.ATE = population.null.hypothesis.after - population.N.after



#=========================================================#
#===== Generate one assignment and fetch the outputs =====#
#=========================================================#

# random.results = random_assignments(random.seed_ = repetition.RANDOM.SEED, K.cardinality_ = 1)
# T.weights = random.results$Treatment.weights
# C.weights = random.results$Control.weights

# stratified.results = stratified_assignments(random.seed_ = repetition.RANDOM.SEED, feature.matrix_ = my.feature.matrix_StA_Matching, K.cardinality_ = 7)
# T.weights = stratified.results$Treatment.weights
# C.weights = stratified.results$Control.weights

#===== Less than a number of treatments =====#
random.results.changing.K = list()
stratified.results.changing.K = list()
changing.K = c(1:floor(N.Regions/2))
for(K.temp in changing.K)
{
  random.results.changing.K[[K.temp]] = random_assignments(random.seed_ = repetition.RANDOM.SEED, K.cardinality_ = K.temp)
}
for(K.temp in changing.K)
{
  stratified.results.changing.K[[K.temp]] = stratified_assignments(random.seed_ = repetition.RANDOM.SEED, feature.matrix_ = my.feature.matrix_StA_Matching, min.units_ = 2, K.cardinality_ = K.temp)
}

# #--- Sanity check: fitness (maximum diameter) of the blocks ---#
# fitnessValue.vec = c()
# for(K.temp in changing.K)
# {
#   fitnessValue.vec = c(fitnessValue.vec, stratified.results.changing.K[[K.temp]]$fitnessValue)
# }
# 
# fitnessValue.vec = fitnessValue.vec[-1]
# file_path = as.character(paste0("output_Randomization/ ", repetition.RANDOM.SEED, "_BlockDiameters.txt"))
# if(any(diff(fitnessValue.vec) > 0)) # Check if the new vector is not monotone decreasing
# {
#   write.table(fitnessValue.vec, file = file_path, row.names = FALSE, col.names = FALSE)
# }




#--- Retrieve outputs for Random Assignment + Difference-in-Means (RADiM) ---#
T.I.after.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.N.after.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.N.after.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.I.after.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.null.hypothesis.after.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
estimated.ATE.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATET.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATEC.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.weights.RADiM.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
C.weights.RADiM.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
null.hypothesis.estimated.ATE.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))

blank.period.residuals.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.naught - T.prime))
T.fitted.N.before.RADiM.matrix = matrix(NA, nrow = 0, ncol = T.naught)
C.fitted.N.before.RADiM.matrix = matrix(NA, nrow = 0, ncol = T.naught)

estimated.p.value.RADiM.vector = c()
null.hypothesis.p.value.RADiM.vector = c()

for(K.temp in c(1:length(changing.K))) 
{
  T.weights.RADiM.matrix = rbind(T.weights.RADiM.matrix, random.results.changing.K[[K.temp]]$Treatment.weights)
  C.weights.RADiM.matrix = rbind(C.weights.RADiM.matrix, random.results.changing.K[[K.temp]]$Control.weights)
  
  T.I.after.RADiM.matrix = rbind(T.I.after.RADiM.matrix, t(random.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  C.N.after.RADiM.matrix = rbind(C.N.after.RADiM.matrix, t(random.results.changing.K[[K.temp]]$Control.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  T.N.after.RADiM.matrix = rbind(T.N.after.RADiM.matrix, t(random.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  C.I.after.RADiM.matrix = rbind(C.I.after.RADiM.matrix, t(random.results.changing.K[[K.temp]]$Control.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  T.null.hypothesis.after.RADiM.matrix = rbind(T.null.hypothesis.after.RADiM.matrix, t(random.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.null.hypothesis.matrix[,(T.naught+1):T.total])
  
  estimated.ATE.RADiM.matrix = rbind(estimated.ATE.RADiM.matrix, T.I.after.RADiM.matrix[K.temp,] - C.N.after.RADiM.matrix[K.temp,])
  ATET.RADiM.matrix = rbind(ATET.RADiM.matrix, T.I.after.RADiM.matrix[K.temp,] - T.N.after.RADiM.matrix[K.temp,])
  ATEC.RADiM.matrix = rbind(ATEC.RADiM.matrix, C.I.after.RADiM.matrix[K.temp,] - C.N.after.RADiM.matrix[K.temp,])
  null.hypothesis.estimated.ATE.RADiM.matrix = rbind(null.hypothesis.estimated.ATE.RADiM.matrix, T.null.hypothesis.after.RADiM.matrix[K.temp,] - C.N.after.RADiM.matrix[K.temp,])
  
  T.fitted.N.before.RADiM.matrix = rbind(T.fitted.N.before.RADiM.matrix, t(random.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.N.matrix[,1:T.naught])
  C.fitted.N.before.RADiM.matrix = rbind(C.fitted.N.before.RADiM.matrix, t(random.results.changing.K[[K.temp]]$Control.weights) %*% Y.N.matrix[,1:T.naught])
  blank.period.residuals.RADiM.matrix = rbind(blank.period.residuals.RADiM.matrix, T.fitted.N.before.RADiM.matrix[K.temp, (T.prime+1):T.naught] - C.fitted.N.before.RADiM.matrix[K.temp, (T.prime+1):T.naught])
  estimated.p.value.RADiM.vector = c(estimated.p.value.RADiM.vector, permutation.test(blank.period.residuals.RADiM.matrix[K.temp,], estimated.ATE.RADiM.matrix[K.temp,]))
  null.hypothesis.p.value.RADiM.vector = c(null.hypothesis.p.value.RADiM.vector, permutation.test(blank.period.residuals.RADiM.matrix[K.temp,], null.hypothesis.estimated.ATE.RADiM.matrix[K.temp,]))
}



#--- Retrieve outputs for Stratified Assignment + Difference-in-Means (StADiM) ---#
T.I.after.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.N.after.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.N.after.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.I.after.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.null.hypothesis.after.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
estimated.ATE.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATET.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATEC.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.weights.StADiM.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
C.weights.StADiM.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
null.hypothesis.estimated.ATE.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))

blank.period.residuals.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.naught - T.prime))
T.fitted.N.before.StADiM.matrix = matrix(NA, nrow = 0, ncol = T.naught)
C.fitted.N.before.StADiM.matrix = matrix(NA, nrow = 0, ncol = T.naught)

estimated.p.value.StADiM.vector = c()
null.hypothesis.p.value.StADiM.vector = c()

for(K.temp in c(1:length(changing.K)))
{
  T.weights.StADiM.matrix = rbind(T.weights.StADiM.matrix, stratified.results.changing.K[[K.temp]]$Treatment.weights)
  C.weights.StADiM.matrix = rbind(C.weights.StADiM.matrix, stratified.results.changing.K[[K.temp]]$Control.weights)
  
  T.I.after.StADiM.matrix = rbind(T.I.after.StADiM.matrix, t(stratified.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  C.N.after.StADiM.matrix = rbind(C.N.after.StADiM.matrix, t(stratified.results.changing.K[[K.temp]]$Control.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  T.N.after.StADiM.matrix = rbind(T.N.after.StADiM.matrix, t(stratified.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  C.I.after.StADiM.matrix = rbind(C.I.after.StADiM.matrix, t(stratified.results.changing.K[[K.temp]]$Control.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  T.null.hypothesis.after.StADiM.matrix = rbind(T.null.hypothesis.after.StADiM.matrix, t(stratified.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.null.hypothesis.matrix[,(T.naught+1):T.total])
  
  estimated.ATE.StADiM.matrix = rbind(estimated.ATE.StADiM.matrix, T.I.after.StADiM.matrix[K.temp,] - C.N.after.StADiM.matrix[K.temp,])
  ATET.StADiM.matrix = rbind(ATET.StADiM.matrix, T.I.after.StADiM.matrix[K.temp,] - T.N.after.StADiM.matrix[K.temp,])
  ATEC.StADiM.matrix = rbind(ATEC.StADiM.matrix, C.I.after.StADiM.matrix[K.temp,] - C.N.after.StADiM.matrix[K.temp,])
  null.hypothesis.estimated.ATE.StADiM.matrix = rbind(null.hypothesis.estimated.ATE.StADiM.matrix, T.null.hypothesis.after.StADiM.matrix[K.temp,] - C.N.after.StADiM.matrix[K.temp,])
  
  T.fitted.N.before.StADiM.matrix = rbind(T.fitted.N.before.StADiM.matrix, t(stratified.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.N.matrix[,1:T.naught])
  C.fitted.N.before.StADiM.matrix = rbind(C.fitted.N.before.StADiM.matrix, t(stratified.results.changing.K[[K.temp]]$Control.weights) %*% Y.N.matrix[,1:T.naught])
  blank.period.residuals.StADiM.matrix = rbind(blank.period.residuals.StADiM.matrix, T.fitted.N.before.StADiM.matrix[K.temp, (T.prime+1):T.naught] - C.fitted.N.before.StADiM.matrix[K.temp, (T.prime+1):T.naught])
  estimated.p.value.StADiM.vector = c(estimated.p.value.StADiM.vector, permutation.test(blank.period.residuals.StADiM.matrix[K.temp,], estimated.ATE.StADiM.matrix[K.temp,]))
  null.hypothesis.p.value.StADiM.vector = c(null.hypothesis.p.value.StADiM.vector, permutation.test(blank.period.residuals.StADiM.matrix[K.temp,], null.hypothesis.estimated.ATE.StADiM.matrix[K.temp,]))
}



#--- Retrieve outputs for Random Assignment + Regression Adjustment (RARegAdj) ---#
#--- regression adjustment can't use pre-treatment outcomes (dimension too big) ---#
estimated.ATE.RARegAdj.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))

for(K.temp in c(1:length(changing.K))) 
{
  T.units = random.results.changing.K[[K.temp]]$Treatment.weights > 0
  C.units = random.results.changing.K[[K.temp]]$Control.weights > 0
  
  estimated.ATE.RARegAdj.vec = c()
  for(t.temp in (T.naught+1):T.total)
  {
    Y.Reg.vec = rep(NA, N.Regions)
    Y.Reg.vec[T.units] = Y.I.matrix[T.units, t.temp]
    Y.Reg.vec[C.units] = Y.N.matrix[C.units, t.temp]
    
    my.df = data.frame(y = Y.Reg.vec, treated = as.numeric(T.units), t(Z.ob.covariates.matrix))
    
    my.lm = lm(y ~ ., data = my.df)
    
    estimated.ATE.RARegAdj.vec = c(estimated.ATE.RARegAdj.vec, my.lm$coefficients["treated"])
  }
  estimated.ATE.RARegAdj.matrix = rbind(estimated.ATE.RARegAdj.matrix, estimated.ATE.RARegAdj.vec)
}
rownames(estimated.ATE.RARegAdj.matrix) = seq(1:length(changing.K))+1



#--- Retrieve outputs for Random Assignment + 1 Nearest Neighbor (RA1NN) ---#
kNN.how.many = 1

estimated.ATE.RA1NN.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))

for(K.temp in c(1:length(changing.K))) 
{
  T.units = random.results.changing.K[[K.temp]]$Treatment.weights > 0
  C.units = random.results.changing.K[[K.temp]]$Control.weights > 0
  
  T.indices = which(T.units)
  C.indices = which(C.units)
  
  distance.matrix.C = matrix(NA, nrow = N.Regions, ncol = N.Regions)
  Y.estimated.matrix.temp = matrix(NA, nrow = 0, ncol = T.total - T.naught)
  for(this.treated in T.indices)
  {
    for(this.control in C.indices)
    {
      distance.matrix.C[this.treated, this.control] = euclidean_distance(my.feature.matrix_StA_Matching[this.treated,], my.feature.matrix_StA_Matching[this.control,])
    }
    smallest.indices = order(distance.matrix.C[this.treated, ])[1:min(kNN.how.many,sum(C.units))]
    
    smallest.indexed.matrix = matrix(NA, nrow = 0, ncol = T.total - T.naught)
    smallest.indexed.matrix = rbind(smallest.indexed.matrix, Y.N.matrix[smallest.indices,(T.naught+1):T.total])
    Y.this.found.control = colMeans(smallest.indexed.matrix)
    
    Y.this.estimated = Y.I.matrix[this.treated,(T.naught+1):T.total] - Y.this.found.control
    Y.estimated.matrix.temp = rbind(Y.estimated.matrix.temp, Y.this.estimated)
  }
  
  estimated.ATE.RA1NN.vec = colMeans(Y.estimated.matrix.temp)
  estimated.ATE.RA1NN.matrix = rbind(estimated.ATE.RA1NN.matrix, estimated.ATE.RA1NN.vec)
}
rownames(estimated.ATE.RA1NN.matrix) = seq(1:length(changing.K))+1



#--- Retrieve outputs for Random Assignment + 5 Nearest Neighbor (RA5NN) ---#
kNN.how.many = 5

estimated.ATE.RA5NN.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))

for(K.temp in c(1:length(changing.K))) 
{
  T.units = random.results.changing.K[[K.temp]]$Treatment.weights > 0
  C.units = random.results.changing.K[[K.temp]]$Control.weights > 0
  
  T.indices = which(T.units)
  C.indices = which(C.units)
  
  distance.matrix.C = matrix(NA, nrow = N.Regions, ncol = N.Regions)
  Y.estimated.matrix.temp = matrix(NA, nrow = 0, ncol = T.total - T.naught)
  for(this.treated in T.indices)
  {
    for(this.control in C.indices)
    {
      distance.matrix.C[this.treated, this.control] = euclidean_distance(my.feature.matrix_StA_Matching[this.treated,], my.feature.matrix_StA_Matching[this.control,])
    }
    smallest.indices = order(distance.matrix.C[this.treated, ])[1:min(kNN.how.many,sum(C.units))]
    
    smallest.indexed.matrix = matrix(NA, nrow = 0, ncol = T.total - T.naught)
    smallest.indexed.matrix = rbind(smallest.indexed.matrix, Y.N.matrix[smallest.indices,(T.naught+1):T.total])
    Y.this.found.control = colMeans(smallest.indexed.matrix)
    
    Y.this.estimated = Y.I.matrix[this.treated,(T.naught+1):T.total] - Y.this.found.control
    Y.estimated.matrix.temp = rbind(Y.estimated.matrix.temp, Y.this.estimated)
  }
  
  estimated.ATE.RA5NN.vec = colMeans(Y.estimated.matrix.temp)
  estimated.ATE.RA5NN.matrix = rbind(estimated.ATE.RA5NN.matrix, estimated.ATE.RA5NN.vec)
}
rownames(estimated.ATE.RA5NN.matrix) = seq(1:length(changing.K))+1




#======================================#
#=====Print and save to .txt files=====#
#======================================#

p.value.rejection.criteria = 0.05

#Average treatment effect Random Assignment + Difference-in-Means (RADiM) 
ATE.RADiM.report = rbind(true.ATE,
                         estimated.ATE.RADiM.matrix)
Loss.MAE.RADiM.vec = c(0)
Loss.MSE.RADiM.vec = c(0)
for(loss.temp in 2:nrow(ATE.RADiM.report))
{
  Loss.MAE.RADiM.vec = c(Loss.MAE.RADiM.vec, sum(abs(ATE.RADiM.report[loss.temp,] - ATE.RADiM.report[1,])) / (T.total - T.naught))
}
for(loss.temp in 2:nrow(ATE.RADiM.report))
{
  Loss.MSE.RADiM.vec = c(Loss.MSE.RADiM.vec, sum((ATE.RADiM.report[loss.temp,] - ATE.RADiM.report[1,]) * (ATE.RADiM.report[loss.temp,] - ATE.RADiM.report[1,])) / (T.total - T.naught))
}
estimated.p.values.RADiM.report = c(1,
                                    estimated.p.value.RADiM.vector)
rejection.RADiM.report = as.numeric(estimated.p.values.RADiM.report < p.value.rejection.criteria)
ATE.RADiM.report.print = cbind(ATE.RADiM.report, Loss.MAE.RADiM.vec, Loss.MSE.RADiM.vec)
my.file.name.copy = as.character(paste0("output11_Randomization/ ", repetition.RANDOM.SEED, "RandomAssignment_RADiM.txt"))
write.table(cbind(round(ATE.RADiM.report.print, digits = 3), estimated.p.values.RADiM.report, rejection.RADiM.report), file = my.file.name.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Average treatment effect Stratified Assignment + Difference-in-Means (StADiM)
ATE.StADiM.report = rbind(true.ATE,
                          estimated.ATE.StADiM.matrix)
Loss.MAE.StADiM.vec = c(0)
Loss.MSE.StADiM.vec = c(0)
for(loss.temp in 2:nrow(ATE.StADiM.report))
{
  Loss.MAE.StADiM.vec = c(Loss.MAE.StADiM.vec, sum(abs(ATE.StADiM.report[loss.temp,] - ATE.StADiM.report[1,])) / (T.total - T.naught))
}
for(loss.temp in 2:nrow(ATE.StADiM.report))
{
  Loss.MSE.StADiM.vec = c(Loss.MSE.StADiM.vec, sum((ATE.StADiM.report[loss.temp,] - ATE.StADiM.report[1,]) * (ATE.StADiM.report[loss.temp,] - ATE.StADiM.report[1,])) / (T.total - T.naught))
}
estimated.p.values.StADiM.report = c(1,
                                     estimated.p.value.StADiM.vector)
rejection.StADiM.report = as.numeric(estimated.p.values.StADiM.report < p.value.rejection.criteria)
ATE.StADiM.report.print = cbind(ATE.StADiM.report, Loss.MAE.StADiM.vec, Loss.MSE.StADiM.vec)
my.file.name.copy = as.character(paste0("output11_Randomization/ ", repetition.RANDOM.SEED, "StratifiedAssignment_StADiM.txt"))
write.table(cbind(round(ATE.StADiM.report.print, digits = 3), estimated.p.values.StADiM.report, rejection.StADiM.report), file = my.file.name.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Average treatment effect Random Assignment + Regression Adjustment (RARegAdj) 
ATE.RARegAdj.report = rbind(true.ATE,
                            estimated.ATE.RARegAdj.matrix)
Loss.MAE.RARegAdj.vec = c(0)
Loss.MSE.RARegAdj.vec = c(0)
for(loss.temp in 2:nrow(ATE.RARegAdj.report))
{
  Loss.MAE.RARegAdj.vec = c(Loss.MAE.RARegAdj.vec, sum(abs(ATE.RARegAdj.report[loss.temp,] - ATE.RARegAdj.report[1,])) / (T.total - T.naught))
}
for(loss.temp in 2:nrow(ATE.RARegAdj.report))
{
  Loss.MSE.RARegAdj.vec = c(Loss.MSE.RARegAdj.vec, sum((ATE.RARegAdj.report[loss.temp,] - ATE.RARegAdj.report[1,]) * (ATE.RARegAdj.report[loss.temp,] - ATE.RARegAdj.report[1,])) / (T.total - T.naught))
}
ATE.RARegAdj.report.print = cbind(ATE.RARegAdj.report, Loss.MAE.RARegAdj.vec, Loss.MSE.RARegAdj.vec)
my.file.name.copy = as.character(paste0("output11_Randomization/ ", repetition.RANDOM.SEED, "RandomAssignment_RARegAdj.txt"))
write.table(round(ATE.RARegAdj.report.print, digits = 3), file = my.file.name.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Average treatment effect Random Assignment + 1 Nearest Neighbor (RA1NN)
ATE.RA1NN.report = rbind(true.ATE,
                         estimated.ATE.RA1NN.matrix)
Loss.MAE.RA1NN.vec = c(0)
Loss.MSE.RA1NN.vec = c(0)
for(loss.temp in 2:nrow(ATE.RA1NN.report))
{
  Loss.MAE.RA1NN.vec = c(Loss.MAE.RA1NN.vec, sum(abs(ATE.RA1NN.report[loss.temp,] - ATE.RA1NN.report[1,])) / (T.total - T.naught))
}
for(loss.temp in 2:nrow(ATE.RA1NN.report))
{
  Loss.MSE.RA1NN.vec = c(Loss.MSE.RA1NN.vec, sum((ATE.RA1NN.report[loss.temp,] - ATE.RA1NN.report[1,]) * (ATE.RA1NN.report[loss.temp,] - ATE.RA1NN.report[1,])) / (T.total - T.naught))
}
ATE.RA1NN.report.print = cbind(ATE.RA1NN.report, Loss.MAE.RA1NN.vec, Loss.MSE.RA1NN.vec)
my.file.name.copy = as.character(paste0("output11_Randomization/ ", repetition.RANDOM.SEED, "RandomAssignment_RA1NN.txt"))
write.table(round(ATE.RA1NN.report.print, digits = 3), file = my.file.name.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Average treatment effect Random Assignment + 5 Nearest Neighbor (RA5NN)
ATE.RA5NN.report = rbind(true.ATE,
                         estimated.ATE.RA5NN.matrix)
Loss.MAE.RA5NN.vec = c(0)
Loss.MSE.RA5NN.vec = c(0)
for(loss.temp in 2:nrow(ATE.RA5NN.report))
{
  Loss.MAE.RA5NN.vec = c(Loss.MAE.RA5NN.vec, sum(abs(ATE.RA5NN.report[loss.temp,] - ATE.RA5NN.report[1,])) / (T.total - T.naught))
}
for(loss.temp in 2:nrow(ATE.RA5NN.report))
{
  Loss.MSE.RA5NN.vec = c(Loss.MSE.RA5NN.vec, sum((ATE.RA5NN.report[loss.temp,] - ATE.RA5NN.report[1,]) * (ATE.RA5NN.report[loss.temp,] - ATE.RA5NN.report[1,])) / (T.total - T.naught))
}
ATE.RA5NN.report.print = cbind(ATE.RA5NN.report, Loss.MAE.RA5NN.vec, Loss.MSE.RA5NN.vec)
my.file.name.copy = as.character(paste0("output11_Randomization/ ", repetition.RANDOM.SEED, "RandomAssignment_RA5NN.txt"))
write.table(round(ATE.RA5NN.report.print, digits = 3), file = my.file.name.copy, append = FALSE, sep = "\t", col.names=FALSE)


# end_time <- Sys.time()
# end_time - start_time

