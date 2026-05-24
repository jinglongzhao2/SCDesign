#==========================================================#
#========== Read and process Walmart Kaggle data ==========#
#==========================================================#

#--------------------------------------------#
#--- Use Lines 749 - 759 to run the codes ---#
#--------------------------------------------#

#--- Last Change: 20250327 1AM --- added multiple randomization methods
#--- Last Change: 20231109 3PM --- converted these codes into batch run

repetition.SEED = 1
# repetition.RANDOM.SEED = 123456
repetition.RANDOM.SEED = as.numeric(Sys.getenv("SGE_TASK_ID"))

library(tidyr)
library(dplyr)
library(slam)
library(Matrix)
library(GA)
library(reticulate)
library(stringr)

set.seed(repetition.RANDOM.SEED)

Weight.digits = 4

# Aggregate data by week (=1), by two weeks (=2), or by month (=4)
time.agg.SEED = ceiling(repetition.SEED / 3 / 3)
if(time.agg.SEED == 1)
{
  config.part3.name = "_weekly_"
}
if(time.agg.SEED == 2)
{
  config.part3.name = "_biweekly_"
}
if(time.agg.SEED == 3)
{
  config.part3.name = "_monthly_"
}

# Choose the outcomes: original, log, or normalized
outcomes.SEED = ceiling((repetition.SEED - (time.agg.SEED - 1) * 3 * 3) / 3)
if(outcomes.SEED == 1)
{
  config.part1.name = "Y_original"
}
if(outcomes.SEED == 2)
{
  config.part1.name = "Y_log"
}
if(outcomes.SEED == 3)
{
  config.part1.name = "Y_normalized"
}

# Choose the covariates: no, few, or many
covariates.SEED = (repetition.SEED - (time.agg.SEED - 1) * 3 * 3 - (outcomes.SEED - 1) * 3)
if(covariates.SEED == 1)
{
  config.part2.name = "_no_covariates"
}
if(covariates.SEED == 2)
{
  config.part2.name = "_few_covariates"
}
if(covariates.SEED == 3)
{
  config.part2.name = ",_many_covariates"
}

config.name = paste0(config.part1.name, config.part2.name, config.part3.name)





#======================#
#===== Parameters =====#
#======================#

# Aggregate data by week (=1), by two weeks (=2), or by month (=4)
if(time.agg.SEED == 1)
{
  merge.once.every = 1
}
if(time.agg.SEED == 2)
{
  merge.once.every = 2
}
if(time.agg.SEED == 3)
{
  merge.once.every = 4
}

# Length of the fitting periods, blank periods, and experimental periods
estimation.length.factor = 0.7
blank.length.factor = 0.2
experimental.length.factor = 1 - estimation.length.factor - blank.length.factor



#=====================#
#===== Read data =====#
#=====================#

data_frame = read.csv("Walmart.csv", sep = ',', header = TRUE)
Read_csv = as.matrix(data_frame)

unique.dates = unique(Read_csv[,2])
unique.date.names = format(as.Date(unique.dates, format = "%d-%m-%Y"), format = "%b %d, %Y")
num.unique = length(unique.dates)
for(temp in 1:num.unique)
{
  Read_csv[Read_csv[,2] == unique.dates[temp], 2] = temp # The dates mark the start date of a week
}

data = data.frame(
  Store_csv = as.numeric(Read_csv[,1]),
  Date_csv = as.numeric(Read_csv[,2]),
  Sales_csv = as.numeric(Read_csv[,3])
)

# Pivot the data to create the desired matrix format
Sales.panel = data %>% pivot_wider(names_from = Date_csv, values_from = Sales_csv)
Sales.panel = Sales.panel[,-1]
Sales.panel = matrix(unlist(Sales.panel), ncol = ncol(Sales.panel), byrow = FALSE)



#========================#
#===== Process data =====#
#========================#

# Merge into biweekly or monthly data
if(merge.once.every == 1)
{
  Sales.panel.merged = Sales.panel
  unique.date.names.merged = unique.date.names
} else
{
  Sales.panel.merged = matrix(NA, nrow = nrow(Sales.panel), ncol = 0)
  for(merge.temp in 1:floor(ncol(Sales.panel)/merge.once.every))
  {
    Sales.panel.merged = cbind(Sales.panel.merged, rowSums(Sales.panel[,((merge.temp*merge.once.every-1):(merge.temp*merge.once.every))]))
  }
  unique.date.names.merged = unique.date.names[seq(1, length(unique.date.names), by = merge.once.every)]
}

# Basics --- Scalars
N.Regions = nrow(Sales.panel.merged)
T.total = ncol(Sales.panel.merged)
print(T.total)
T.prime = floor(T.total * estimation.length.factor) # Fitting periods
T.naught = floor(T.total * (estimation.length.factor + blank.length.factor))
r.ob.covariates.dim = ncol(Read_csv) - 4 # Here, 3 out of 4 stands for the store id, date id, and sales; the last 1 is "holiday" which is a time-indicator

# Get the observed covariates --- Matrix
cov.list = list()
all.cov.list = list()
cov.average = list()
cov.merged.list = list()
for(cov.temp in 1:r.ob.covariates.dim)
{
  data = data.frame(
    Store_csv = as.numeric(Read_csv[,1]),
    Date_csv = as.numeric(Read_csv[,2]),
    cov_csv = as.numeric(Read_csv[,(cov.temp+4)])
  )
  
  all.cov.list[[cov.temp]] = data %>% pivot_wider(names_from = Date_csv, values_from = cov_csv)
  all.cov.list[[cov.temp]] = all.cov.list[[cov.temp]][,-1]
  cov.list[[cov.temp]] = data %>% pivot_wider(names_from = Date_csv, values_from = cov_csv)
  cov.list[[cov.temp]] = cov.list[[cov.temp]][,-1]
  cov.list[[cov.temp]] = cov.list[[cov.temp]][,1:(T.prime*merge.once.every)]
  cov.average[[cov.temp]] = rowMeans(cov.list[[cov.temp]])
  
  cov.merged.temp = matrix(NA, nrow = nrow(cov.list[[cov.temp]]), ncol = 0)
  for(merge.temp in 1:floor((ncol(cov.list[[cov.temp]])/4)))
  {
    cov.merged.temp = cbind(cov.merged.temp, rowMeans(cov.list[[cov.temp]][,((merge.temp*4-3):(merge.temp*4))]))
  }
  cov.merged.list[[cov.temp]] = cov.merged.temp
}





#-----------------------------#
#--- Choose the covariates ---#
#-----------------------------#

if(covariates.SEED == 1)
{
  #--- No covariates ---#
  Z.ob.covariates.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
  r.ob.covariates.dim = nrow(Z.ob.covariates.matrix)
}
if(covariates.SEED == 2)
{
  #--- Few covariates ---#
  Z.ob.covariates.matrix = matrix(NA, nrow = 4, ncol = N.Regions)
  for(cov.temp in 1:4)
  {
    Z.ob.covariates.matrix[cov.temp,] = cov.average[[cov.temp]]
  }
  r.ob.covariates.dim = nrow(Z.ob.covariates.matrix)
}
if(covariates.SEED == 3)
{
  #--- Many covariates ---#
  Z.ob.covariates.matrix = matrix(NA, nrow = 2, ncol = N.Regions)
  for(cov.temp in 1:2)
  {
    Z.ob.covariates.matrix[cov.temp,] = cov.average[[cov.temp]]
  }
  Z.ob.covariates.matrix = rbind(Z.ob.covariates.matrix, t(cov.merged.list[[3]]))
  Z.ob.covariates.matrix = rbind(Z.ob.covariates.matrix, t(cov.merged.list[[4]]))
  Z.ob.covariates.matrix = round(Z.ob.covariates.matrix, digits = 2)
  r.ob.covariates.dim = nrow(Z.ob.covariates.matrix) #--- Potentially using many (time-dependent) covariates ---#
}





#---------------------------#
#--- Choose the outcomes ---#
#---------------------------#

if(outcomes.SEED == 1)
{
  #--- Y original ---#
  # Get the outcomes (sales) --- Matrix
  Y.N.matrix = Sales.panel.merged
  # Population weights --- Vector
  f.vector = rep(1/N.Regions, N.Regions)
  target.outcomes = f.vector %*% Y.N.matrix
  plot(1:T.total, target.outcomes, type = 'l')
}
if(outcomes.SEED == 2)
{
  #--- Y log ---#
  # Get the outcomes (sales) --- Matrix
  Y.N.matrix = log(Sales.panel.merged+1) # use +1 to avoid zero values
  # Population weights --- Vector
  f.vector = rep(1/N.Regions, N.Regions)
  target.outcomes = f.vector %*% Y.N.matrix
  plot(1:T.total, target.outcomes, type = 'l')
}
if(outcomes.SEED == 3)
{
  #--- Y Normalized by size ---#
  # Get the outcomes (sales) --- Matrix
  Approximated.sizes = rowMeans(Sales.panel.merged[,1:T.prime]) # normalize by row means
  Y.N.matrix = Sales.panel.merged / Approximated.sizes
  # Population weights --- Vector
  f.vector = Approximated.sizes / sum(Approximated.sizes)
  target.outcomes = f.vector %*% (Y.N.matrix[,1:T.total]) * sum(Approximated.sizes)
  plot(1:T.total, target.outcomes, type = 'l')
}





#-------------------------------------------------------------#
#----- Randomization 1: less than a number of treatments -----#
#-------------------------------------------------------------#

random_assignments <- function(random.seed_, N.Regions_ = N.Regions, f.vector_ = f.vector, K.cardinality_ = -1)
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
  # treated.weights = treated.weights * f.vector_ * N.Regions_
  # control.weights = control.weights * f.vector_ * N.Regions_
  
  treated.units = sample(c(1:N.Regions), K.cardinality)
  control.units = setdiff(c(1:N.Regions), treated.units)
  
  treated.weight.value = 1 / length(treated.units)
  control.weight.value = 1 / length(control.units)
  treated.weights[treated.units] = treated.weight.value
  control.weights[control.units] = control.weight.value
  
  treated.weights = treated.weights * f.vector_ * N.Regions_
  control.weights = control.weights * f.vector_ * N.Regions_
  
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

stratified_assignments <- function(random.seed_, N.Regions_ = N.Regions, f.vector_ = f.vector, feature.matrix_, min.units_ = 2, K.cardinality_ = -1)
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
    
    # my.blocks = my_MinMaxDiameter_blocking(random.seed_ = random.seed_, N.Regions_ = N.Regions_, feature.matrix_ = feature.matrix_, K.cardinality_ = K.cardinality)
    my.blocks = my_IterMatching_blocking(random.seed_ = random.seed_, N.Regions_ = N.Regions_, feature.matrix_ = feature.matrix_, K.cardinality_ = K.cardinality)
    
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
      treated.weights[this.cluster.treated] = f.vector_[this.cluster.treated] * length(this.cluster.temp) / length(this.cluster.treated)
      control.weights[this.cluster.control] = f.vector_[this.cluster.control] * length(this.cluster.temp) / length(this.cluster.control)
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
    
    treated.weights = treated.weights * f.vector_ * N.Regions_
    control.weights = control.weights * f.vector_ * N.Regions_
  }
  
  return(list(Treatment.weights = round(treated.weights, digits = 8), 
              Control.weights = round(control.weights, digits = 8),
              fitnessValue = fitness_))
}


#---------------------------------------------------#
#--- We assume that interventions have no effect ---#
#---------------------------------------------------#
Y.I.matrix = Y.N.matrix

population.N.after = t(f.vector) %*% Y.N.matrix[,(T.naught+1):T.total]
population.I.after = t(f.vector) %*% Y.I.matrix[,(T.naught+1):T.total]

true.ATE = population.I.after - population.N.after

#===============================================#
#===== Compare with randomized assignments =====#
#===============================================#

my.feature.matrix_StA_Matching = Y.N.matrix[,1:T.naught]
my.feature.matrix_StA_Matching = apply(my.feature.matrix_StA_Matching, 2, function(x) (x - mean(x)) / sd(x))


#===== Less than a number of treatments =====#
random.results.changing.K = list()
stratified.results.changing.K = list()
changing.K = c(1:5)
for(K.temp in changing.K)
{
  random.results.changing.K[[K.temp]] = random_assignments(random.seed_ = repetition.RANDOM.SEED, K.cardinality_ = K.temp)
}
for(K.temp in changing.K)
{
  stratified.results.changing.K[[K.temp]] = stratified_assignments(random.seed_ = repetition.RANDOM.SEED, feature.matrix_ = my.feature.matrix_StA_Matching, min.units_ = 2, K.cardinality_ = K.temp)
}



#--- Retrieve outputs for Random Assignment + Difference-in-Means (RADiM) ---#
T.I.after.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.N.after.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.N.after.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.I.after.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
estimated.ATE.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATET.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATEC.RADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.weights.RADiM.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
C.weights.RADiM.matrix = matrix(NA, nrow = 0, ncol = N.Regions)

T.fitted.N.before.RADiM.matrix = matrix(NA, nrow = 0, ncol = T.naught)
C.fitted.N.before.RADiM.matrix = matrix(NA, nrow = 0, ncol = T.naught)

for(K.temp in c(1:length(changing.K))) 
{
  T.weights.RADiM.matrix = rbind(T.weights.RADiM.matrix, random.results.changing.K[[K.temp]]$Treatment.weights)
  C.weights.RADiM.matrix = rbind(C.weights.RADiM.matrix, random.results.changing.K[[K.temp]]$Control.weights)
  
  T.I.after.RADiM.matrix = rbind(T.I.after.RADiM.matrix, t(random.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  C.N.after.RADiM.matrix = rbind(C.N.after.RADiM.matrix, t(random.results.changing.K[[K.temp]]$Control.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  T.N.after.RADiM.matrix = rbind(T.N.after.RADiM.matrix, t(random.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  C.I.after.RADiM.matrix = rbind(C.I.after.RADiM.matrix, t(random.results.changing.K[[K.temp]]$Control.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  
  estimated.ATE.RADiM.matrix = rbind(estimated.ATE.RADiM.matrix, T.I.after.RADiM.matrix[K.temp,] - C.N.after.RADiM.matrix[K.temp,])
  ATET.RADiM.matrix = rbind(ATET.RADiM.matrix, T.I.after.RADiM.matrix[K.temp,] - T.N.after.RADiM.matrix[K.temp,])
  ATEC.RADiM.matrix = rbind(ATEC.RADiM.matrix, C.I.after.RADiM.matrix[K.temp,] - C.N.after.RADiM.matrix[K.temp,])
  
  T.fitted.N.before.RADiM.matrix = rbind(T.fitted.N.before.RADiM.matrix, t(random.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.N.matrix[,1:T.naught])
  C.fitted.N.before.RADiM.matrix = rbind(C.fitted.N.before.RADiM.matrix, t(random.results.changing.K[[K.temp]]$Control.weights) %*% Y.N.matrix[,1:T.naught])
}



#--- Retrieve outputs for Stratified Assignment + Difference-in-Means (StADiM) ---#
T.I.after.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.N.after.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.N.after.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.I.after.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
estimated.ATE.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATET.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATEC.StADiM.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.weights.StADiM.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
C.weights.StADiM.matrix = matrix(NA, nrow = 0, ncol = N.Regions)

T.fitted.N.before.StADiM.matrix = matrix(NA, nrow = 0, ncol = T.naught)
C.fitted.N.before.StADiM.matrix = matrix(NA, nrow = 0, ncol = T.naught)

for(K.temp in c(1:length(changing.K)))
{
  T.weights.StADiM.matrix = rbind(T.weights.StADiM.matrix, stratified.results.changing.K[[K.temp]]$Treatment.weights)
  C.weights.StADiM.matrix = rbind(C.weights.StADiM.matrix, stratified.results.changing.K[[K.temp]]$Control.weights)
  
  T.I.after.StADiM.matrix = rbind(T.I.after.StADiM.matrix, t(stratified.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  C.N.after.StADiM.matrix = rbind(C.N.after.StADiM.matrix, t(stratified.results.changing.K[[K.temp]]$Control.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  T.N.after.StADiM.matrix = rbind(T.N.after.StADiM.matrix, t(stratified.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  C.I.after.StADiM.matrix = rbind(C.I.after.StADiM.matrix, t(stratified.results.changing.K[[K.temp]]$Control.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  
  estimated.ATE.StADiM.matrix = rbind(estimated.ATE.StADiM.matrix, T.I.after.StADiM.matrix[K.temp,] - C.N.after.StADiM.matrix[K.temp,])
  ATET.StADiM.matrix = rbind(ATET.StADiM.matrix, T.I.after.StADiM.matrix[K.temp,] - T.N.after.StADiM.matrix[K.temp,])
  ATEC.StADiM.matrix = rbind(ATEC.StADiM.matrix, C.I.after.StADiM.matrix[K.temp,] - C.N.after.StADiM.matrix[K.temp,])
  
  T.fitted.N.before.StADiM.matrix = rbind(T.fitted.N.before.StADiM.matrix, t(stratified.results.changing.K[[K.temp]]$Treatment.weights) %*% Y.N.matrix[,1:T.naught])
  C.fitted.N.before.StADiM.matrix = rbind(C.fitted.N.before.StADiM.matrix, t(stratified.results.changing.K[[K.temp]]$Control.weights) %*% Y.N.matrix[,1:T.naught])
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
ATE.RADiM.report.print = cbind(ATE.RADiM.report, Loss.MAE.RADiM.vec, Loss.MSE.RADiM.vec)
my.file.name.copy = as.character(paste0("output_Randomization/ ", repetition.RANDOM.SEED, "RandomAssignment_RADiM.txt"))
write.table(round(ATE.RADiM.report.print, digits = 3), file = my.file.name.copy, append = FALSE, sep = "\t", col.names=FALSE)

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
ATE.StADiM.report.print = cbind(ATE.StADiM.report, Loss.MAE.StADiM.vec, Loss.MSE.StADiM.vec)
my.file.name.copy = as.character(paste0("output_Randomization/ ", repetition.RANDOM.SEED, "StratifiedAssignment_StADiM.txt"))
write.table(round(ATE.StADiM.report.print, digits = 3), file = my.file.name.copy, append = FALSE, sep = "\t", col.names=FALSE)

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
my.file.name.copy = as.character(paste0("output_Randomization/ ", repetition.RANDOM.SEED, "RandomAssignment_RARegAdj.txt"))
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
my.file.name.copy = as.character(paste0("output_Randomization/ ", repetition.RANDOM.SEED, "RandomAssignment_RA1NN.txt"))
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
my.file.name.copy = as.character(paste0("output_Randomization/ ", repetition.RANDOM.SEED, "RandomAssignment_RA5NN.txt"))
write.table(round(ATE.RA5NN.report.print, digits = 3), file = my.file.name.copy, append = FALSE, sep = "\t", col.names=FALSE)




