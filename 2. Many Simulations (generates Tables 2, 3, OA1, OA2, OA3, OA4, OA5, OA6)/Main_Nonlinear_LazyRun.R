#======================================================================================#
#========== This is a lazy run to generate all results under nonlinear model ==========#
#========== Two tasks: (1) compare different formulations                    ==========#
#==========            (2) compare with randomized experiments               ==========#
#======================================================================================#

#==================================================#
#===== TASK 1: Compare different formulations =====#
#==================================================#

#-------------------------------------------------------------------#
#----- Gurobi cannot be called too often in cluster computing ------#
#----- we use open source solvers 'limSolve' and 'quadprog' --------#
#----- Gurobi parameter: Threads = 1 (for cluster computing) -------#
#-------------------------------------------------------------------#

#--- Last Change: 20221106 4PM
#--- Last Change: 20221107 8PM
#--- Last Change: 20221108 3PM
#--- Last Change: 20221117 3PM
#--- Last Change: 20221117 10PM --- changed: recentering and rescaling of X.matrix
#---------------------------------- added: The Weakly Targeted formulation
#--- Last Change: 20221121 9PM  --- round the solutions to ~10^(-3) because numerical error could be up to ~10^(-4)
#--- Last Change: 20221201 10PM --- added MSE (4 occasions in the outputs), put the new formulation earlier
#--- Last Change: 20221205 12PM --- run everything one more time to double check
#--- Last Change: 20221212 10AM --- added controllable digits and trimming 
#---------------------------------- the edits on 20221121 9PM might not be correct because ~10^(-4) errors could lead to problems with rounding to ~10^(-3)
#--- Last Change: 20250316 2PM  --- set the weights in the estimand (beta.vector) to be uniform weights 1/N.Regions
#---------------------------------- in the "Analyzing_Outputs.R" file we changed the MSE into RMSE to make it comparable with MAE
#---------------------------------- for simplicity, we do not use print_model.primitives to print the model primitives
#--- Last Change: 20250901 9PM  --- fixed a bug when all the covariates take the same value, which causes scaling to divide by 0
#--- Last Change: 20250906 2PM  --- renamed some random variables to align with the names in the paper
#--- Last Change: 20250906 3PM  --- a nonlinear setup using exponential function
#---------------------------------- changed range.coefficients.max from 10 to 3 to reduce the magnitude of errors


# setwd("/project/jinglonggrp/Synthetic Experiment/")
# repetition.RANDOM.SEED = 123456
repetition.RANDOM.SEED = as.numeric(Sys.getenv("SGE_TASK_ID"))

library(slam)
library(dplyr)
library(gurobi)
library(Matrix)
library(limSolve)
library(quadprog)

start_time <- Sys.time()

#=========================================#
#=====Below we define model variables=====#
#=========================================#

#===Basics
Weight.digits = 2
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
range.coefficients.max = 3 #ChangeThis #20250906 3PM changed this from 10 to 3
#Noises --- Matrices
#Normal Distribution ~ N(0,1)
noise.variance = 1 #ChangeThis

#=======================================================#
#=====Below we define two data generating processes=====#
#=======================================================#

Generate_Model_Primitives <- function(range.intercept.max_ = range.intercept.max,
                                      range.covariates.max_ = range.covariates.max,
                                      range.coefficients.max_ = range.coefficients.max,
                                      noise.variance_ = noise.variance,
                                      random.seed_)
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
    xi.temp = rep(NA, N.Regions)
    xi.I.matrix_[,t] = xi.temp
  }
  for(t in (T.naught+1):T.total)
  {
    xi.temp = rnorm(N.Regions, 0, noise.variance_)
    xi.I.matrix_[,t] = xi.temp
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

Generate_Nonlinear_Primitives <- function(range.intercept.max_ = range.intercept.max,
                                          range.covariates.max_ = range.covariates.max,
                                          range.coefficients.max_ = range.coefficients.max,
                                          noise.variance_ = noise.variance,
                                          random.seed_)
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
  beta.temp = runif(N.Regions)
  beta.vector_ = beta.temp / sum(beta.temp)
  #===Model Primitives
  #Constants --- Vectors
  #Range --- [0,20] #range.intercept.max_ = 20
  delta.N.vector_ = sort(c(range.intercept.max_ * runif(T.naught), range.intercept.max_ * runif(T.total - T.naught)))
  upsilon.I.vector_ = c(rep(NA, T.naught), sort(range.intercept.max_ * runif(T.total - T.naught)))
  #Covariates --- Matrices
  #Range --- [0,1] # range.covariates.max_ = 1
  for(j in 1:N.Regions)
  {
    z.temp = range.covariates.max_ * (runif(r.ob.covariates.dim) - 0.5)
    Z.ob.covariates.matrix_[,j] = z.temp
  }
  colnames(Z.ob.covariates.matrix_) = c(1:N.Regions)
  for(j in 1:N.Regions)
  {
    mu.temp = range.covariates.max_ * (runif(F.unob.covariates.dim) - 0.5)
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
    xi.temp = rep(NA, N.Regions)
    xi.I.matrix_[,t] = xi.temp
  }
  for(t in (T.naught+1):T.total)
  {
    xi.temp = rnorm(N.Regions, 0, noise.variance_)
    xi.I.matrix_[,t] = xi.temp
  }
  colnames(xi.I.matrix_) = c(1:T.total)
  
  #=====Generate Potential Outcomes=====#
  #Potential Outcomes in hindsight --- Matrices
  #A Factor Model
  for(j in 1:N.Regions)
  {
    for(t in 1:T.total)
    {
      Y.N.matrix_[j,t] = delta.N.vector_[t] + exp(theta.ob.N.matrix_[,t] %*% Z.ob.covariates.matrix_[,j]) + exp(lambda.unob.N.matrix_[,t] %*% mu.unob.covariates.matrix_[,j]) + epsilon.N.matrix_[j,t]
    }
  }
  
  for(j in 1:N.Regions)
  {
    for(t in 1:T.total)
    {
      Y.I.matrix_[j,t] = upsilon.I.vector_[t] + exp(gamma.ob.I.matrix_[,t] %*% Z.ob.covariates.matrix_[,j]) + exp(eta.unob.I.matrix_[,t] %*% mu.unob.covariates.matrix_[,j]) + xi.I.matrix_[j,t]
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

#================================================#
#=====Below we define different formulations=====#
#================================================#

#--- Support functions 

norm_square <- function(x) sum(x^2)

l.zero.norm <- function(vector1, vector2)
{
  vector1.binary = as.numeric(vector1 != 0)
  vector2.binary = as.numeric(vector2 != 0)
  return(sum((vector1.binary - vector2.binary) != 0))
}

#==================================#
#=====The Original Formulation=====#
#==================================#

Synthetic_Experiment <- function(T.prime_, r.ob.covariates.dim_, N.Regions_, Y.N.matrix, Z.ob.covariates.matrix, beta.vector_, digits_ = Weight.digits)
{
  beta.vec = beta.vector_
  
  M.dim = T.prime_ + r.ob.covariates.dim_
  N.dim = N.Regions_
  
  X.matrix = matrix(NA, nrow = M.dim, ncol = N.dim)
  for(i in 1:M.dim)
  {
    for(j in 1:N.dim)
    {
      if(i <= T.prime_)
      {
        X.matrix[i,j] = Y.N.matrix[j, i]
      }
      else
      {
        X.matrix[i,j] = Z.ob.covariates.matrix[(i - T.prime_),j]
      }
    }
  }
  #--- rescale each row ---#
  row.means = apply(X.matrix, 1, mean)
  row.stdevs = apply(X.matrix, 1, sd)
  row.stdevs[row.stdevs == 0] = 1 #--- just in case there are covariates that has no variation (in this case we do not need to scale or do anything)
  X.matrix = (X.matrix - matrix(row.means, nrow = M.dim, ncol = N.dim, byrow = FALSE)) / row.stdevs
  
  center.vector = X.matrix %*% beta.vec
  
  #=====Use Gurobi to solve the weights=====#
  #===Preparations
  model = list()
  
  #=====Write the objective
  #===Quadratic terms
  Q.obj = matrix(NA, nrow = 2 * N.dim, ncol = 2 * N.dim)
  # #Diagonal elements
  # for(q.pointer in 1:N.dim)
  # {
  #   sum.temp = 0
  #   for(i in 1:M.dim)
  #   {
  #     sum.temp = sum.temp + X.matrix[i,q.pointer]^2
  #   }
  #   Q.obj[q.pointer, q.pointer] = sum.temp
  #   Q.obj[q.pointer + N.dim, q.pointer + N.dim] = sum.temp
  # }
  # #Off-diagonal elements
  #Diagonal blocks
  for(q.pointer1 in 1:N.dim)
  {
    for(q.pointer2 in 1:N.dim)
    {
      sum.temp = 0
      for(i in 1:M.dim)
      {
        sum.temp = sum.temp + X.matrix[i,q.pointer1] * X.matrix[i,q.pointer2]
      }
      Q.obj[q.pointer1, q.pointer2] = sum.temp
      Q.obj[q.pointer1 + N.dim, q.pointer2 + N.dim] = sum.temp
    }
  }
  #Finishing Q.obj
  Q.obj[is.na(Q.obj)] = 0
  #===Linear terms
  #Note there are coefficients "2"!
  C.obj = rep(NA, 2 * N.dim)
  for(c.pointer in 1:N.dim)
  {
    sum.temp = 0
    for(i in 1:M.dim)
    {
      sum.temp = sum.temp + X.matrix[i,c.pointer] * center.vector[i]
    }
    C.obj[c.pointer] = - 2 * sum.temp
    C.obj[c.pointer + N.dim] = - 2 * sum.temp
  }
  #===Constants (irrelavant to optimal solutions)
  C.const = 0
  sum.temp = 0
  for(i in 1:M.dim)
  {
    sum.temp = sum.temp + center.vector[i]^2
  }
  C.const = 2 * sum.temp #Each L2 norm contributes to one constant -- two constants in total
  
  #=====Write the constraints
  #===Linear constraints
  A.matrix = matrix(NA, nrow = 2, ncol = 2 * N.dim)
  for(a.pointer in 1:N.dim)
  {
    A.matrix[1,a.pointer] = 1
    A.matrix[2,a.pointer+N.dim] = 1
  }
  #Finishing A.matrix
  A.matrix[is.na(A.matrix)] = 0
  #===Quadratic constraints
  qc1 = list()
  Q.constr = matrix(NA, nrow = 2 * N.dim, ncol = 2 * N.dim)
  for(qc.pointer in 1:N.dim)
  {
    Q.constr[qc.pointer, qc.pointer + N.dim] = 1
    Q.constr[qc.pointer + N.dim, qc.pointer] = 1
  }
  #Finishing Q.constr
  Q.constr[is.na(Q.constr)] = 0
  qc1$Qc = Q.constr
  qc1$rhs = 0
  
  #=====Pass to Gurobi solver
  #Objective
  model$Q       <- Q.obj
  model$obj     <- C.obj
  #Variables
  model$vtype   <- 'C'
  #Linear constraints
  model$A       <- A.matrix
  model$rhs     <- c(1,1)
  model$sense   <- c('=', '=')
  #Quadratic constraints
  model$quadcon <- list(qc1)
  #Constants in the objective
  model$objcon  <- C.const
  
  #=====Solve the model
  params = list(Threads = 1, TimeLimit = 3600, NonConvex=2) #2 stands for non-convex optimization
  result = gurobi(model, params = params)
  
  #=====Treatment regions are the smaller set of regions
  units.temp1 = round(result$x[1:N.dim], digits = digits_)
  units.temp2 = round(result$x[(N.dim+1):(2*N.dim)], digits = digits_)
  #--- Trimming: reduce the near-zero weights to zero
  units.temp1[units.temp1 <= 10^(-digits_)] = 0
  units.temp2[units.temp2 <= 10^(-digits_)] = 0
  
  if(sum(units.temp1 != 0) > sum(units.temp2 != 0))
  {
    units.temp.buf = units.temp1
    units.temp1 = units.temp2
    units.temp2 = units.temp.buf
  }
  if(sum(units.temp1 != 0) == sum(units.temp2 != 0) && min(which(units.temp1>0)) > min(which(units.temp2>0)))
  {
    units.temp.buf = units.temp1
    units.temp1 = units.temp2
    units.temp2 = units.temp.buf
  }
  
  returned = list(Treatment.weights = units.temp1, 
                  Control.weights = units.temp2, 
                  solver.output = result)
  return(returned)
}

#--- Easier Implementation (no linear terms)

Synthetic_Experiment_Easier_Formulation <- function(T.prime_, r.ob.covariates.dim_, N.Regions_, Y.N.matrix_, Z.ob.covariates.matrix_, beta.vector_, digits_ = Weight.digits)
{
  beta.vec = beta.vector_
  
  M.dim = T.prime_ + r.ob.covariates.dim_
  N.dim = N.Regions_
  
  X.matrix = matrix(NA, nrow = M.dim, ncol = N.dim)
  for(i in 1:M.dim)
  {
    for(j in 1:N.dim)
    {
      if(i <= T.prime_)
      {
        X.matrix[i,j] = Y.N.matrix_[j, i]
      }
      else
      {
        X.matrix[i,j] = Z.ob.covariates.matrix_[(i - T.prime_),j]
      }
    }
  }
  #--- rescale each row ---#
  row.means = apply(X.matrix, 1, mean)
  row.stdevs = apply(X.matrix, 1, sd)
  row.stdevs[row.stdevs == 0] = 1
  X.matrix = (X.matrix - matrix(row.means, nrow = M.dim, ncol = N.dim, byrow = FALSE)) / row.stdevs
  
  center.vector = X.matrix %*% beta.vec
  X.effective = X.matrix - matrix(data = center.vector, nrow = nrow(X.matrix), ncol = ncol(X.matrix), byrow = FALSE)
  
  #=====Use Gurobi to solve the weights=====#
  #===Preparations
  model = list()
  
  #=====Write the objective
  #===Quadratic terms
  Q.obj = bdiag(t(X.effective) %*% X.effective, t(X.effective) %*% X.effective)
  
  #=====Write the constraints
  #===Linear constraints
  A.matrix = matrix(NA, nrow = 2, ncol = 2 * N.dim)
  for(a.pointer in 1:N.dim)
  {
    A.matrix[1,a.pointer] = 1
    A.matrix[2,a.pointer+N.dim] = 1
  }
  #Finishing A.matrix
  A.matrix[is.na(A.matrix)] = 0
  #===Quadratic constraints
  qc1 = list()
  Q.constr = matrix(NA, nrow = 2 * N.dim, ncol = 2 * N.dim)
  for(qc.pointer in 1:N.dim)
  {
    Q.constr[qc.pointer, qc.pointer + N.dim] = 1
    Q.constr[qc.pointer + N.dim, qc.pointer] = 1
  }
  #Finishing Q.constr
  Q.constr[is.na(Q.constr)] = 0
  qc1$Qc = Q.constr
  qc1$rhs = 0
  
  #=====Pass to Gurobi solver
  #Objective
  model$Q       <- Q.obj
  #Variables
  model$vtype   <- 'C'
  #Linear constraints
  model$A       <- A.matrix
  model$rhs     <- c(1,1)
  model$sense   <- c('=', '=')
  #Quadratic constraints
  model$quadcon <- list(qc1)
  
  #=====Solve the model
  params = list(Threads = 1, TimeLimit = 3600, NonConvex=2) #2 stands for non-convex optimization
  result = gurobi(model, params = params)
  
  #=====Treatment regions are the smaller set of regions
  units.temp1 = round(result$x[1:N.dim], digits = digits_)
  units.temp2 = round(result$x[(N.dim+1):(2*N.dim)], digits = digits_)
  #--- Trimming: reduce the near-zero weights to zero
  units.temp1[units.temp1 <= 10^(-digits_)] = 0
  units.temp2[units.temp2 <= 10^(-digits_)] = 0
  
  if(sum(units.temp1 != 0) > sum(units.temp2 != 0))
  {
    units.temp.buf = units.temp1
    units.temp1 = units.temp2
    units.temp2 = units.temp.buf
  }
  if(sum(units.temp1 != 0) == sum(units.temp2 != 0) && min(which(units.temp1>0)) > min(which(units.temp2>0)))
  {
    units.temp.buf = units.temp1
    units.temp1 = units.temp2
    units.temp2 = units.temp.buf
  }
  
  returned = list(Treatment.weights = units.temp1, 
                  Control.weights = units.temp2, 
                  solver.output = result)
  return(returned)
}

#=========================================#
#=====The Weakly Targeted Formulation=====#
#=========================================#

Synthetic_Experiment_Weakly_Targeted <- function(T.prime_, r.ob.covariates.dim_, N.Regions_, Y.N.matrix, Z.ob.covariates.matrix, beta.vector_, xi_ = 1, digits_ = Weight.digits)
{
  beta.vec = beta.vector_
  
  M.dim = T.prime_ + r.ob.covariates.dim_
  N.dim = N.Regions_
  
  X.matrix = matrix(NA, nrow = M.dim, ncol = N.dim)
  for(i in 1:M.dim)
  {
    for(j in 1:N.dim)
    {
      if(i <= T.prime_)
      {
        X.matrix[i,j] = Y.N.matrix[j, i]
      }
      else
      {
        X.matrix[i,j] = Z.ob.covariates.matrix[(i - T.prime_),j]
      }
    }
  }
  #--- rescale each row ---#
  row.means = apply(X.matrix, 1, mean)
  row.stdevs = apply(X.matrix, 1, sd)
  row.stdevs[row.stdevs == 0] = 1
  X.matrix = (X.matrix - matrix(row.means, nrow = M.dim, ncol = N.dim, byrow = FALSE)) / row.stdevs
  
  center.vector = X.matrix %*% beta.vec
  
  #=====Use Gurobi to solve the weights=====#
  #===Preparations
  model = list()
  
  #=====Write the objective
  #===Quadratic terms
  Q.obj = matrix(NA, nrow = 2 * N.dim, ncol = 2 * N.dim)
  # #Diagonal elements
  # for(q.pointer in 1:N.dim)
  # {
  #   sum.temp = 0
  #   for(i in 1:M.dim)
  #   {
  #     sum.temp = sum.temp + X.matrix[i,q.pointer]^2
  #   }
  #   Q.obj[q.pointer, q.pointer] = sum.temp
  #   Q.obj[q.pointer + N.dim, q.pointer + N.dim] = sum.temp
  # }
  # #Off-diagonal elements
  #Diagonal blocks - coming from (f average vs. treatment weights average)
  for(q.pointer1 in 1:N.dim)
  {
    for(q.pointer2 in 1:N.dim)
    {
      sum.temp = 0
      for(i in 1:M.dim)
      {
        sum.temp = sum.temp + X.matrix[i,q.pointer1] * X.matrix[i,q.pointer2]
      }
      Q.obj[q.pointer1, q.pointer2] = sum.temp
    }
  }
  #Finishing Q.obj
  Q.obj[is.na(Q.obj)] = 0
  #Diagonal blocks - coming from (treatment weights average vs. control weights average)
  for(q.pointer1 in 1:N.dim)
  {
    for(q.pointer2 in 1:N.dim)
    {
      sum.temp = 0
      for(i in 1:M.dim)
      {
        sum.temp = sum.temp + X.matrix[i,q.pointer1] * X.matrix[i,q.pointer2]
      }
      Q.obj[q.pointer1, q.pointer2] = Q.obj[q.pointer1, q.pointer2] + xi_ * sum.temp
      Q.obj[q.pointer1 + N.dim, q.pointer2 + N.dim] = Q.obj[q.pointer1 + N.dim, q.pointer2 + N.dim] + xi_ * sum.temp
    }
  }
  #Off-diagonal blocks - coming from (treatment weights average vs. control weights average)
  for(q.pointer1 in 1:N.dim)
  {
    for(q.pointer2 in 1:N.dim)
    {
      sum.temp = 0
      for(i in 1:M.dim)
      {
        sum.temp = sum.temp + X.matrix[i,q.pointer1] * X.matrix[i,q.pointer2]
      }
      Q.obj[q.pointer1 + N.dim, q.pointer2] = - xi_ * sum.temp
      Q.obj[q.pointer1, q.pointer2 + N.dim] = - xi_ * sum.temp
    }
  }
  #===Linear terms
  #Note there are coefficients "2"!
  #The f averaged part only comes with w weights. The (w and v) weights part are quadratic terms
  C.obj = rep(NA, 2 * N.dim)
  for(c.pointer in 1:N.dim)
  {
    sum.temp = 0
    for(i in 1:M.dim)
    {
      sum.temp = sum.temp + X.matrix[i,c.pointer] * center.vector[i]
    }
    C.obj[c.pointer] = - 2 * sum.temp
  }
  #Finishing C.obj
  C.obj[is.na(C.obj)] = 0
  #===Constants (irrelavant to optimal solutions)
  C.const = 0
  sum.temp = 0
  for(i in 1:M.dim)
  {
    sum.temp = sum.temp + center.vector[i]^2
  }
  C.const = sum.temp #The first L2 norm contributes to one constant -- the second doesn't have constant terms
  
  #=====Write the constraints
  #===Linear constraints
  A.matrix = matrix(NA, nrow = 2, ncol = 2 * N.dim)
  for(a.pointer in 1:N.dim)
  {
    A.matrix[1,a.pointer] = 1
    A.matrix[2,a.pointer+N.dim] = 1
  }
  #Finishing A.matrix
  A.matrix[is.na(A.matrix)] = 0
  #===Quadratic constraints
  qc1 = list()
  Q.constr = matrix(NA, nrow = 2 * N.dim, ncol = 2 * N.dim)
  for(qc.pointer in 1:N.dim)
  {
    Q.constr[qc.pointer, qc.pointer + N.dim] = 1
    Q.constr[qc.pointer + N.dim, qc.pointer] = 1
  }
  #Finishing Q.constr
  Q.constr[is.na(Q.constr)] = 0
  qc1$Qc = Q.constr
  qc1$rhs = 0
  
  #=====Pass to Gurobi solver
  #Objective
  model$Q       <- Q.obj
  model$obj     <- C.obj
  #Variables
  model$vtype   <- 'C'
  #Linear constraints
  model$A       <- A.matrix
  model$rhs     <- c(1,1)
  model$sense   <- c('=', '=')
  #Quadratic constraints
  model$quadcon <- list(qc1)
  #Constants in the objective
  model$objcon  <- C.const
  
  #=====Solve the model
  params = list(Threads = 1, TimeLimit = 3600, NonConvex=2) #2 stands for non-convex optimization
  result = gurobi(model, params = params)
  
  #=====Treatment regions are the smaller set of regions
  units.temp1 = round(result$x[1:N.dim], digits = digits_)
  units.temp2 = round(result$x[(N.dim+1):(2*N.dim)], digits = digits_)
  #--- Trimming: reduce the near-zero weights to zero
  units.temp1[units.temp1 <= 10^(-digits_)] = 0
  units.temp2[units.temp2 <= 10^(-digits_)] = 0
  
  returned = list(Treatment.weights = units.temp1, 
                  Control.weights = units.temp2, 
                  solver.output = result)
  return(returned)
}

#=====================================================#
#=====The Synthetic Control (Observational Study)=====#
#=====================================================#

Synthetic_Control <- function(target.vector_, X.matrix_)
{
  X.effective = X.matrix_ - matrix(data = target.vector_, nrow = nrow(X.matrix_), ncol = ncol(X.matrix_), byrow = FALSE)
  
  #======================================================================#
  #=============================== NOTICE ===============================#
  #=====Gurobi needs to check the liscense every time it is called. =====#
  #=====And Gurobi does not like to be called too frequently. ===========#
  #=====Cannot use Gurobi to solve simple optmization. ==================#
  #======================================================================#
  
  # #=====Use Gurobi to solve the weights=====#
  # #===Preparations
  # model = list()
  # 
  # #=====Write the objective
  # #===Quadratic terms
  # Q.obj = t(X.effective) %*% X.effective
  # 
  # #=====Write the constraints
  # #===Linear constraints
  # A.matrix = matrix(1, nrow = 1, ncol = ncol(X.matrix_))
  # 
  # #=====Pass to Gurobi solver
  # #Objective
  # model$Q       <- Q.obj
  # #Variables
  # model$vtype   <- 'C'
  # #Linear constraints
  # model$A       <- A.matrix
  # model$rhs     <- c(1)
  # model$sense   <- c('=')
  # 
  # #=====Solve the model
  # params = list(Threads = 1, TimeLimit = 3600, OutputFlag=0)
  # result = gurobi(model, params = params)
  # 
  # returned = list(Weights = result$x, solver.output = result)
  
  #============================================================================#
  #=====Below we use the lsei package to solve constrained least squares. =====#
  #=====https://search.r-project.org/CRAN/refmans/limSolve/html/lsei.html =====#
  #=====The input formulation is as follows: ==================================#
  #============================= min  ||Ax-b||^2 ==============================#
  #============================= s.t. Ex  = f    ==============================#
  #=============================      Gx >= h    ==============================#
  #============================================================================#
  
  AA = X.effective
  BB = rep(0, nrow(AA))
  EE = rep(1, ncol(AA))
  FF = 1
  GG = diag(1, nrow = ncol(AA))
  HH = rep(0, ncol(AA))
  result = lsei(E = EE, F = FF, A = AA, B = BB, G = GG, H = HH)
  
  returned = list(Weights = result$X, Objval = result$solutionNorm, solver.output = result)
  return(returned)
}

#=====================================#
#=====The Constrained Formulation=====#
#=====================================#

Synthetic_Experiment_Cardinality_Constraint <- function(T.prime_, r.ob.covariates.dim_, N.Regions_, Y.N.matrix_, Z.ob.covariates.matrix_, beta.vector_, K.cardinality_ = -1, digits_ = Weight.digits)
{
  if(K.cardinality_ == -1)
  {
    K.cardinality = floor(N.Regions_/2)
  }
  if(K.cardinality_ >= 1)
  {
    K.cardinality = K.cardinality_
  }
  
  beta.vec = beta.vector_
  
  M.dim = T.prime_ + r.ob.covariates.dim_
  N.dim = N.Regions_
  
  X.matrix = matrix(NA, nrow = M.dim, ncol = N.dim)
  for(i in 1:M.dim)
  {
    for(j in 1:N.dim)
    {
      if(i <= T.prime_)
      {
        X.matrix[i,j] = Y.N.matrix_[j, i]
      }
      else
      {
        X.matrix[i,j] = Z.ob.covariates.matrix_[(i - T.prime_),j]
      }
    }
  }
  #--- rescale each row ---#
  row.means = apply(X.matrix, 1, mean)
  row.stdevs = apply(X.matrix, 1, sd)
  row.stdevs[row.stdevs == 0] = 1
  X.matrix = (X.matrix - matrix(row.means, nrow = M.dim, ncol = N.dim, byrow = FALSE)) / row.stdevs
  
  center.vector = X.matrix %*% beta.vec
  
  #=====Enumerate over all possible partitions
  candidate.partition = list()
  for(partition.cardinality in 1:K.cardinality)
  {
    candidate.partition = c(candidate.partition, combn(1:N.Regions_, partition.cardinality, simplify = FALSE))
  }
  total.candidate.partitions = length(candidate.partition)
  
  loss = c()
  for(temp.candidate in 1:total.candidate.partitions)
  {
    this.candidate = candidate.partition[[temp.candidate]]
    
    if(length(this.candidate) == 1)
    {
      X.tentative.treated = matrix(X.matrix[,this.candidate], ncol = 1)
    }
    if(length(this.candidate) != 1)
    {
      X.tentative.treated = X.matrix[,this.candidate]
    }
    X.tentative.control = X.matrix[,-this.candidate]
    Synth.treated = Synthetic_Control(target.vector_ = center.vector, X.matrix_ = X.tentative.treated)
    Synth.control = Synthetic_Control(target.vector_ = center.vector, X.matrix_ = X.tentative.control)
    
    loss = c(loss, (Synth.treated$Objval + Synth.control$Objval))
  }
  found.min = which.min(loss)
  
  found.min.candidate = candidate.partition[[found.min]]
  if(length(found.min.candidate) == 1)
  {
    X.found.treated = matrix(X.matrix[,found.min.candidate], ncol = 1)
  }
  if(length(found.min.candidate) != 1)
  {
    X.found.treated = X.matrix[,found.min.candidate]
  }
  X.found.control = X.matrix[,-found.min.candidate]
  Synth.found.treated = Synthetic_Control(target.vector_ = center.vector, X.matrix_ = X.found.treated)
  Synth.found.control = Synthetic_Control(target.vector_ = center.vector, X.matrix_ = X.found.control)
  
  found.treated.weights = rep(0, N.Regions_)
  for(n.t.temp in 1:(length(found.min.candidate)))
  {
    found.treated.weights[found.min.candidate[n.t.temp]] = Synth.found.treated$Weights[n.t.temp]
  }
  this.complement = setdiff(1:N.Regions_, found.min.candidate)
  found.control.weights = rep(0, N.Regions_)
  for(n.t.temp in 1:(length(this.complement)))
  {
    found.control.weights[this.complement[n.t.temp]] = Synth.found.control$Weights[n.t.temp]
  }
  
  #=====Treatment regions are the smaller set of regions
  units.temp1 = round(found.treated.weights, digits = digits_)
  units.temp2 = round(found.control.weights, digits = digits_)
  if(sum(units.temp1 != 0) > sum(units.temp2 != 0))
  {
    units.temp.buf = units.temp1
    units.temp1 = units.temp2
    units.temp2 = units.temp.buf
  }
  if(sum(units.temp1 != 0) == sum(units.temp2 != 0) && min(which(units.temp1>0)) > min(which(units.temp2>0)))
  {
    units.temp.buf = units.temp1
    units.temp1 = units.temp2
    units.temp2 = units.temp.buf
  }
  
  returned = list(Treatment.weights = units.temp1, 
                  Control.weights = units.temp2)
  return(returned)
}

#====================================#
#=====The Unit Level Formulation=====#
#====================================#

Synthetic_Experiment_Unit_Level <- function(T.prime_, r.ob.covariates.dim_, N.Regions_, Y.N.matrix_, Z.ob.covariates.matrix_, beta.vector_, xi_ = 1, K.cardinality_ = -1, digits_ = Weight.digits)
{
  if(K.cardinality_ == -1)
  {
    K.cardinality = N.Regions_ - 1
  }
  if(K.cardinality_ >= 1)
  {
    K.cardinality = K.cardinality_
  }
  
  beta.vec = beta.vector_
  
  M.dim = T.prime_ + r.ob.covariates.dim_
  N.dim = N.Regions_
  
  X.matrix = matrix(NA, nrow = M.dim, ncol = N.dim)
  for(i in 1:M.dim)
  {
    for(j in 1:N.dim)
    {
      if(i <= T.prime_)
      {
        X.matrix[i,j] = Y.N.matrix_[j, i]
      }
      else
      {
        X.matrix[i,j] = Z.ob.covariates.matrix_[(i - T.prime_),j]
      }
    }
  }
  #--- rescale each row ---#
  row.means = apply(X.matrix, 1, mean)
  row.stdevs = apply(X.matrix, 1, sd)
  row.stdevs[row.stdevs == 0] = 1
  X.matrix = (X.matrix - matrix(row.means, nrow = M.dim, ncol = N.dim, byrow = FALSE)) / row.stdevs
  
  center.vector = X.matrix %*% beta.vec
  X.effective = X.matrix - matrix(data = center.vector, nrow = nrow(X.matrix), ncol = ncol(X.matrix), byrow = FALSE)
  
  #=====Enumerate over all possible partitions
  candidate.partition = list()
  for(partition.cardinality in 1:K.cardinality)
  {
    candidate.partition = c(candidate.partition, combn(1:N.Regions_, partition.cardinality, simplify = FALSE))
  }
  total.candidate.partitions = length(candidate.partition)
  
  loss = c()
  treated.weights = matrix(NA, nrow = 0, ncol = N.dim)
  for(temp.candidate in 1:total.candidate.partitions)
  {
    this.candidate = candidate.partition[[temp.candidate]]
    
    if(length(this.candidate) == 1)
    {
      X.tentative.treated = matrix(X.matrix[,this.candidate], ncol = 1)
      X.tentative.control = X.matrix[,-this.candidate]
    }
    if(length(this.candidate) == (N.Regions_ - 1))
    {
      X.tentative.treated = X.matrix[,this.candidate]
      X.tentative.control = matrix(X.matrix[,-this.candidate], ncol = 1)
    }
    if(length(this.candidate) != 1 && length(this.candidate) != (N.Regions_ - 1))
    {
      X.tentative.treated = X.matrix[,this.candidate]
      X.tentative.control = X.matrix[,-this.candidate]
    }
    
    #---Deal with control units first---#
    loss.each.T.vector = c()
    for(candidate.T.temp in this.candidate)
    {
      Synth.control = Synthetic_Control(target.vector_ = matrix(X.matrix[,candidate.T.temp], nrow = nrow(X.matrix), byrow = FALSE), X.matrix_ = X.tentative.control)
      loss.each.T.vector = c(loss.each.T.vector, Synth.control$Objval)
    }
    
    #---Deal with treated units then---#
    #--treat the control unit losses as linear coefficients for the treated units--#
    linear.coefficients.each.T = xi_ * loss.each.T.vector
    #--Now solve a Synthetic Control Problem with linear terms--#
    
    #======================================================================#
    #=============================== NOTICE ===============================#
    #=====Gurobi needs to check the liscense every time it is called. =====#
    #=====And Gurobi does not like to be called too frequently. ===========#
    #=====Cannot use Gurobi to solve simple optmization. ==================#
    #======================================================================#
    
    #=====Use Gurobi to solve the weights=====#
    # #===Preparations
    # model = list()
    
    #=====Write the objective
    #===Quadratic terms
    Q.obj = bdiag(t(X.effective) %*% X.effective)
    #===Linear terms
    C.obj = rep(NA, N.dim)
    for(linear.coef.temp in 1:length(this.candidate))
    {
      C.obj[this.candidate[linear.coef.temp]] = linear.coefficients.each.T[linear.coef.temp]
    }
    C.obj[is.na(C.obj)] = 0
    
    #=====Write the constraints
    #===Linear constraints
    A.matrix = matrix(NA, nrow = 1 + N.dim - length(this.candidate), ncol = ncol(X.matrix))
    for(candidate.T.temp in this.candidate)
    {
      A.matrix[1, candidate.T.temp] = 1
    }
    remaining.candidates = setdiff(1:N.Regions_, this.candidate)
    for(candidate.C.temp in 1:length(remaining.candidates))
    {
      A.matrix[candidate.C.temp+1, remaining.candidates[candidate.C.temp]] = 1 #+1 because the first equality constraint is the sum equal to one constraint
    }
    A.matrix[is.na(A.matrix)] = 0
    
    # #=====Pass to Gurobi solver
    # #Objective
    # model$Q       <- Q.obj
    # model$obj     <- C.obj
    # #Variables
    # model$vtype   <- 'C'
    # #Linear constraints
    # model$A       <- A.matrix
    # model$rhs     <- c(1, rep(0, N.dim - length(this.candidate)))
    # model$sense   <- c('=', rep("=", N.dim - length(this.candidate)))
    # 
    # #=====Solve the model
    # params = list(Threads = 1, TimeLimit = 3600, OutputFlag=0) #2 stands for non-convex optimization
    # Synth.treated = gurobi(model, params = params)
    # 
    #=====Save the objective value
    # loss = c(loss, (Synth.treated$objval))
    # treated.weights = rbind(treated.weights, round(Synth.treated$x, digits = digits_))
    
    #================================================================================#
    #=====Below we use the quadprog package to solve constrained least squares. =====#
    #=====https://cran.r-project.org/web/packages/quadprog/quadprog.pdf =============#
    #=====The input formulation is as follows: ======================================#
    #===== min  (-d'x + 1/2 x' D x)            ======================================#
    #===== s.t.                 A'x >= b       ======================================#
    #===== first 'meq' rows of A are equality constraints ===========================#
    #================================================================================#
    
    Dmat = 2 * Q.obj
    pd_Dmat = nearPD(Dmat)$mat
    dvec = - C.obj
    A.eq = A.matrix
    b.eq = c(1, rep(0, N.dim - length(this.candidate)))
    A.nonneg = matrix(NA, nrow = length(this.candidate), ncol = ncol(X.matrix))
    b.nonneg = c()
    for(A.nonneg.temp in 1:length(this.candidate))
    {
      A.nonneg[A.nonneg.temp, this.candidate[A.nonneg.temp]] = 1
      b.nonneg = c(b.nonneg, 0)
    }
    A.nonneg[is.na(A.nonneg)] = 0
    
    Synth.treated = solve.QP(Dmat = pd_Dmat,
                             dvec = dvec,
                             Amat = t(rbind(A.eq, A.nonneg)),
                             bvec = c(b.eq, b.nonneg),
                             meq = nrow(A.eq)) # this argument says the first "meq" rows of Amat are equalities
    
    loss = c(loss, (Synth.treated$value))
    treated.weights = rbind(treated.weights, round(Synth.treated$solution, digits = digits_))
    
    # Synth.treated = Synthetic_Control(target.vector_ = center.vector, X.matrix_ = X.tentative.treated)
    # #----------------# round(Synth.treated$Weights, digits = digits_)
    # loss = c(loss, (Synth.treated$solver.output$objval + xi_ * Synth.treated$Weights %*% loss.each.T.vector))
    # loss = c(loss, (Synth.treated$solver.output$objval + xi_ * sum(loss.each.T.vector)))
  }
  found.min = which.min(loss)
  found.treated.weights = treated.weights[found.min,]
  
  found.min.candidate = candidate.partition[[found.min]]
  if(length(found.min.candidate) == 1)
  {
    X.found.treated = matrix(X.matrix[,found.min.candidate], ncol = 1)
    X.found.control = X.matrix[,-found.min.candidate]
  }
  if(length(found.min.candidate) == (N.Regions_ - 1))
  {
    X.found.treated = X.matrix[,found.min.candidate]
    X.found.control = matrix(X.matrix[,-found.min.candidate], ncol = 1)
  }
  if(length(found.min.candidate) != 1 && length(found.min.candidate) != (N.Regions_ - 1))
  {
    X.found.treated = X.matrix[,found.min.candidate]
    X.found.control = X.matrix[,-found.min.candidate]
  }
  # Synth.found.treated = Synthetic_Control(target.vector_ = center.vector, X.matrix_ = X.found.treated)
  # found.treated.weights = rep(0, N.Regions_)
  # for(n.t.temp in 1:(length(found.min.candidate)))
  # {
  #   found.treated.weights[found.min.candidate[n.t.temp]] = round(Synth.found.treated$Weights[n.t.temp], digits = digits_)
  # }
  
  weights.for.each.T = matrix(NA, nrow = 0, ncol = (N.Regions_ - length(found.min.candidate)))
  for(for.T.temp in 1:length(found.min.candidate))
  {
    candidate.T.temp = found.min.candidate[for.T.temp]
    Synth.control = Synthetic_Control(target.vector_ = matrix(X.matrix[,candidate.T.temp], nrow = nrow(X.matrix), byrow = FALSE), X.matrix_ = X.found.control)
    weights.for.each.T = rbind(weights.for.each.T, round(Synth.control$Weights, digits = digits_))
  }
  this.complement = setdiff(1:N.Regions_, found.min.candidate)
  found.control.weights.for.each.T = matrix(NA, nrow = 0, ncol = N.Regions_)
  for(for.T.temp in 1:length(found.min.candidate))
  {
    found.control.weights = rep(0, N.Regions_)
    for(n.t.temp in 1:length(this.complement))
    {
      found.control.weights[this.complement[n.t.temp]] = weights.for.each.T[for.T.temp,n.t.temp]
    }
    found.control.weights.for.each.T = rbind(found.control.weights.for.each.T, found.control.weights)
  }
  found.control.weights = rep(0, N.Regions_)
  found.control.weights = found.treated.weights[found.min.candidate] %*% found.control.weights.for.each.T
  
  returned = list(Treatment.weights = round(found.treated.weights, digits = digits_), 
                  Control.weights = round(as.vector(found.control.weights), digits = digits_),
                  All.control.weights = found.control.weights.for.each.T)
  return(returned)
}

#===================================#
#=====The Penalized Formulation=====#
#===================================#

Synthetic_Experiment_Penalization <- function(T.prime_, r.ob.covariates.dim_, N.Regions_, Y.N.matrix_, Z.ob.covariates.matrix_, beta.vector_, lambda_ = 1, digits_ = Weight.digits)
{
  beta.vec = beta.vector_
  
  M.dim = T.prime_ + r.ob.covariates.dim_
  N.dim = N.Regions_
  
  X.matrix = matrix(NA, nrow = M.dim, ncol = N.dim)
  for(i in 1:M.dim)
  {
    for(j in 1:N.dim)
    {
      if(i <= T.prime_)
      {
        X.matrix[i,j] = Y.N.matrix_[j, i]
      }
      else
      {
        X.matrix[i,j] = Z.ob.covariates.matrix_[(i - T.prime_),j]
      }
    }
  }
  #--- rescale each row ---#
  row.means = apply(X.matrix, 1, mean)
  row.stdevs = apply(X.matrix, 1, sd)
  row.stdevs[row.stdevs == 0] = 1
  X.matrix = (X.matrix - matrix(row.means, nrow = M.dim, ncol = N.dim, byrow = FALSE)) / row.stdevs
  
  center.vector = X.matrix %*% beta.vec
  X.effective = X.matrix - matrix(data = center.vector, nrow = nrow(X.matrix), ncol = ncol(X.matrix), byrow = FALSE)
  
  #=====Use Gurobi to solve the weights=====#
  #===Preparations
  model = list()
  
  #=====Write the objective
  #===Quadratic terms
  Q.obj = bdiag(t(X.effective) %*% X.effective, t(X.effective) %*% X.effective)
  #===Linear terms
  C.obj = rep(NA, 2 * N.dim)
  for (linear.coef.temp in 1:N.dim)
  {
    C.obj[linear.coef.temp] = norm_square(X.effective[,linear.coef.temp])
    C.obj[linear.coef.temp + N.dim] = norm_square(X.effective[,linear.coef.temp])
  }
  C.obj = C.obj * lambda_
  
  #=====Write the constraints
  #===Linear constraints
  A.matrix = matrix(NA, nrow = 2, ncol = 2 * N.dim)
  for(a.pointer in 1:N.dim)
  {
    A.matrix[1,a.pointer] = 1
    A.matrix[2,a.pointer+N.dim] = 1
  }
  #Finishing A.matrix
  A.matrix[is.na(A.matrix)] = 0
  #===Quadratic constraints
  qc1 = list()
  Q.constr = matrix(NA, nrow = 2 * N.dim, ncol = 2 * N.dim)
  for(qc.pointer in 1:N.dim)
  {
    Q.constr[qc.pointer, qc.pointer + N.dim] = 1
    Q.constr[qc.pointer + N.dim, qc.pointer] = 1
  }
  #Finishing Q.constr
  Q.constr[is.na(Q.constr)] = 0
  qc1$Qc = Q.constr
  qc1$rhs = 0
  
  #=====Pass to Gurobi solver
  #Objective
  model$Q       <- Q.obj
  model$obj     <- C.obj
  #Variables
  model$vtype   <- 'C'
  #Linear constraints
  model$A       <- A.matrix
  model$rhs     <- c(1,1)
  model$sense   <- c('=', '=')
  #Quadratic constraints
  model$quadcon <- list(qc1)
  
  #=====Solve the model
  params = list(Threads = 1, TimeLimit = 3600, NonConvex=2) #2 stands for non-convex optimization
  result = gurobi(model, params = params)
  
  #=====Treatment regions are the smaller set of regions
  units.temp1 = round(result$x[1:N.dim], digits = digits_)
  units.temp2 = round(result$x[(N.dim+1):(2*N.dim)], digits = digits_)
  #--- Trimming: reduce the near-zero weights to zero
  units.temp1[units.temp1 <= 10^(-digits_)] = 0
  units.temp2[units.temp2 <= 10^(-digits_)] = 0
  
  if(sum(units.temp1 != 0) > sum(units.temp2 != 0))
  {
    units.temp.buf = units.temp1
    units.temp1 = units.temp2
    units.temp2 = units.temp.buf
  }
  if(sum(units.temp1 != 0) == sum(units.temp2 != 0) && min(which(units.temp1>0)) > min(which(units.temp2>0)))
  {
    units.temp.buf = units.temp1
    units.temp1 = units.temp2
    units.temp2 = units.temp.buf
  }
  
  returned = list(Treatment.weights = units.temp1, 
                  Control.weights = units.temp2, 
                  solver.output = result)
  return(returned)
}

#==========================#
#=====Permutation Test=====#
#==========================#

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

#====================================#
#===== Repeating multiple times =====#
#====================================#

repetition.output.Original.Formulation.SAVED = list()
repetition.output.Cardinality.SAVED = list()
repetition.output.Penalization.SAVED = list()
repetition.output.UnitLevel.SAVED = list()
repetition.output.WeaklyTargeted.SAVED = list()
repetition.reported.Original.Formulation.SAVED = matrix(NA, nrow = 0, ncol = (4*(T.total-T.naught) + (T.naught - T.prime) + 4))

#=====================================#
#===== Generate Model Primitives =====#
#=====================================#
# model.primitives = Generate_Model_Primitives(random.seed_ = repetition.RANDOM.SEED)
model.primitives = Generate_Nonlinear_Primitives(random.seed_ = repetition.RANDOM.SEED)
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


#======================================================#
#=====  Call Different Optimization Formulations  =====#
#===== Warning: This step easily takes many hours =====#
#======================================================#

#=====Formulation: original=====#
my.result = Synthetic_Experiment(T.prime, r.ob.covariates.dim, N.Regions, Y.N.matrix, Z.ob.covariates.matrix, beta.vector)
#-----Save and print-----#
repetition.output.Original.Formulation.SAVED[[repetition.RANDOM.SEED]] = my.result

#=====Formulation: cardinality=====#
my.results.Cardinality.changing.K = list()
changing.K = c(1:floor(N.Regions/2))
for(K.temp in changing.K)
{
  my.results.Cardinality.changing.K[[K.temp]] = Synthetic_Experiment_Cardinality_Constraint(T.prime, r.ob.covariates.dim, N.Regions, Y.N.matrix, Z.ob.covariates.matrix, beta.vector, K.cardinality_ = K.temp)
}
#-----Save and print-----#
repetition.output.Cardinality.SAVED[[repetition.RANDOM.SEED]] = my.results.Cardinality.changing.K

#=====Formulation: penalization=====#
my.results.Penalization.changing.lambda = list()
# changing.lambda = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)
changing.lambda = c(0.01, 0.1, 1, 10, 100)
for(lambda.temp in 1:length(changing.lambda))
{
  lambda.this.round = changing.lambda[lambda.temp]
  my.results.Penalization.changing.lambda[[lambda.temp]] = Synthetic_Experiment_Penalization(T.prime, r.ob.covariates.dim, N.Regions, Y.N.matrix, Z.ob.covariates.matrix, beta.vector, lambda_ = lambda.this.round)
}
#-----Save and print-----#
repetition.output.Penalization.SAVED[[repetition.RANDOM.SEED]] = my.results.Penalization.changing.lambda

#=====Formulation: unit level=====#
my.results.UnitLevel.changing.xi = list()
# changing.xi = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)
changing.xi = c(0.01, 0.1, 1, 10, 100)
for(xi.temp in 1:length(changing.xi))
{
  xi.this.round = changing.xi[xi.temp]
  my.results.UnitLevel.changing.xi[[xi.temp]] = Synthetic_Experiment_Unit_Level(T.prime, r.ob.covariates.dim, N.Regions, Y.N.matrix, Z.ob.covariates.matrix, beta.vector, xi_ = xi.this.round, K.cardinality_ = -1)
}
#-----Save and print-----#
repetition.output.UnitLevel.SAVED[[repetition.RANDOM.SEED]] = my.results.UnitLevel.changing.xi

#=====Formulation: weakly targeted=====#
my.results.WeaklyTargeted.changing.xi = list()
# changing.xi = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)
changing.xi = c(0.01, 0.1, 1, 10, 100)
for(xi.temp in 1:length(changing.xi))
{
  xi.this.round = changing.xi[xi.temp]
  my.results.WeaklyTargeted.changing.xi[[xi.temp]] = Synthetic_Experiment_Weakly_Targeted(T.prime, r.ob.covariates.dim, N.Regions, Y.N.matrix, Z.ob.covariates.matrix, beta.vector, xi_ = xi.this.round)
}
#-----Save and print-----#
repetition.output.WeaklyTargeted.SAVED[[repetition.RANDOM.SEED]] = my.results.WeaklyTargeted.changing.xi



#========================================#
#=====Fetch the Optimization Outputs=====#
#========================================#

#-----------------------------------------------------------#
#---Rounding is needed. ------------------------------------#
#---Since the QCQP is solved numerically, the unselected ---#
#---weights are not exactly zero, ~10^(-6) error. ----------#
#-----------------------------------------------------------#

#=====Formulation: original=====#
T.weights = round(my.result$Treatment.weights, digits = 3)
C.weights = round(my.result$Control.weights, digits = 3)

population.N.before = t(beta.vector) %*% Y.N.matrix[,1:T.naught]
T.fitted.N.before = t(T.weights) %*% Y.N.matrix[,1:T.naught]
C.fitted.N.before = t(C.weights) %*% Y.N.matrix[,1:T.naught]

population.N.after = t(beta.vector) %*% Y.N.matrix[,(T.naught+1):T.total]
population.I.after = t(beta.vector) %*% Y.I.matrix[,(T.naught+1):T.total]
T.I.after = t(T.weights) %*% Y.I.matrix[,(T.naught+1):T.total]
C.N.after = t(C.weights) %*% Y.N.matrix[,(T.naught+1):T.total]

true.ATE = population.I.after - population.N.after
estimated.ATE = T.I.after - C.N.after
ATET.T.weighted = t(T.weights) %*% Y.I.matrix[,(T.naught+1):T.total] - t(T.weights) %*% Y.N.matrix[,(T.naught+1):T.total]
ATEC.C.weighted = t(C.weights) %*% Y.I.matrix[,(T.naught+1):T.total] - t(C.weights) %*% Y.N.matrix[,(T.naught+1):T.total]

blank.period.residuals = T.fitted.N.before[(T.prime+1):T.naught] - C.fitted.N.before[(T.prime+1):T.naught]
estimated.p.value = permutation.test(blank.period.residuals, estimated.ATE)
MeanAbsoluteError = sum(abs(true.ATE - estimated.ATE))/length(estimated.ATE)

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

T.null.hypothesis.after = t(T.weights) %*% Y.null.hypothesis.matrix[,(T.naught+1):T.total]
null.hypothesis.estimated.ATE = T.null.hypothesis.after - C.N.after
null.hypothesis.p.value = permutation.test(blank.period.residuals, null.hypothesis.estimated.ATE)




#=====Formulation: cardinality=====#
T.I.after.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.N.after.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.N.after.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.I.after.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.null.hypothesis.after.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
estimated.ATE.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATET.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATEC.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.weights.Cardinality.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
C.weights.Cardinality.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
null.hypothesis.estimated.ATE.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))

blank.period.residuals.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.naught - T.prime))
T.fitted.N.before.Cardinality.matrix = matrix(NA, nrow = 0, ncol = T.naught)
C.fitted.N.before.Cardinality.matrix = matrix(NA, nrow = 0, ncol = T.naught)

estimated.p.value.Cardinality.vector = c()
null.hypothesis.p.value.Cardinality.vector = c()
for(K.temp in changing.K[length(changing.K):1]) 
  #--- Note: per rule of swapping treatment and control weights, we are filling the matrix from the last row
{
  tentative.T.weights = my.results.Cardinality.changing.K[[K.temp]]$Treatment.weights
  tentative.C.weights = my.results.Cardinality.changing.K[[K.temp]]$Control.weights
  
  if(nrow(T.weights.Cardinality.matrix) > 0)
  {
    last.row.T.weights = T.weights.Cardinality.matrix[nrow(T.weights.Cardinality.matrix),]
    last.row.C.weights = C.weights.Cardinality.matrix[nrow(C.weights.Cardinality.matrix),]
  } else
  {
    last.row.T.weights = T.weights
    last.row.C.weights = C.weights
  }
  l0discrepancy.no.flip = l.zero.norm(last.row.T.weights, tentative.T.weights) + l.zero.norm(last.row.C.weights, tentative.C.weights)
  l0discrepancy.flip = l.zero.norm(last.row.C.weights, tentative.T.weights) + l.zero.norm(last.row.T.weights, tentative.C.weights)
  
  if(sum(tentative.T.weights != 0) == sum(tentative.C.weights != 0) && l0discrepancy.no.flip > l0discrepancy.flip)
  { #flip treatment and control weights
    use.T.weights = tentative.C.weights
    use.C.weights = tentative.T.weights
  } else
  {
    use.T.weights = tentative.T.weights
    use.C.weights = tentative.C.weights
  }
  
  T.weights.Cardinality.matrix = rbind(T.weights.Cardinality.matrix, use.T.weights)
  C.weights.Cardinality.matrix = rbind(C.weights.Cardinality.matrix, use.C.weights)
  
  T.I.after.Cardinality.matrix = rbind(T.I.after.Cardinality.matrix, t(use.T.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  C.N.after.Cardinality.matrix = rbind(C.N.after.Cardinality.matrix, t(use.C.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  T.N.after.Cardinality.matrix = rbind(T.N.after.Cardinality.matrix, t(use.T.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  C.I.after.Cardinality.matrix = rbind(C.I.after.Cardinality.matrix, t(use.C.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  T.null.hypothesis.after.Cardinality.matrix = rbind(T.null.hypothesis.after.Cardinality.matrix, t(use.T.weights) %*% Y.null.hypothesis.matrix[,(T.naught+1):T.total])
  
  estimated.ATE.Cardinality.matrix = rbind(estimated.ATE.Cardinality.matrix, T.I.after.Cardinality.matrix[(length(changing.K)+1-K.temp),] - C.N.after.Cardinality.matrix[(length(changing.K)+1-K.temp),])
  ATET.Cardinality.matrix = rbind(ATET.Cardinality.matrix, T.I.after.Cardinality.matrix[(length(changing.K)+1-K.temp),] - T.N.after.Cardinality.matrix[(length(changing.K)+1-K.temp),])
  ATEC.Cardinality.matrix = rbind(ATEC.Cardinality.matrix, C.I.after.Cardinality.matrix[(length(changing.K)+1-K.temp),] - C.N.after.Cardinality.matrix[(length(changing.K)+1-K.temp),])
  null.hypothesis.estimated.ATE.Cardinality.matrix = rbind(null.hypothesis.estimated.ATE.Cardinality.matrix, T.null.hypothesis.after.Cardinality.matrix[(length(changing.K)+1-K.temp),] - C.N.after.Cardinality.matrix[(length(changing.K)+1-K.temp),])
  
  T.fitted.N.before.Cardinality.matrix = rbind(T.fitted.N.before.Cardinality.matrix, t(use.T.weights) %*% Y.N.matrix[,1:T.naught])
  C.fitted.N.before.Cardinality.matrix = rbind(C.fitted.N.before.Cardinality.matrix, t(use.C.weights) %*% Y.N.matrix[,1:T.naught])
  blank.period.residuals.Cardinality.matrix = rbind(blank.period.residuals.Cardinality.matrix, T.fitted.N.before.Cardinality.matrix[(length(changing.K)+1-K.temp), (T.prime+1):T.naught] - C.fitted.N.before.Cardinality.matrix[(length(changing.K)+1-K.temp), (T.prime+1):T.naught])
  estimated.p.value.Cardinality.vector = c(estimated.p.value.Cardinality.vector, permutation.test(blank.period.residuals.Cardinality.matrix[(length(changing.K)+1-K.temp),], estimated.ATE.Cardinality.matrix[(length(changing.K)+1-K.temp),]))
  null.hypothesis.p.value.Cardinality.vector = c(null.hypothesis.p.value.Cardinality.vector, permutation.test(blank.period.residuals.Cardinality.matrix[(length(changing.K)+1-K.temp),], null.hypothesis.estimated.ATE.Cardinality.matrix[(length(changing.K)+1-K.temp),]))
}

T.I.after.Cardinality.matrix = T.I.after.Cardinality.matrix[nrow(T.I.after.Cardinality.matrix):1,]
C.N.after.Cardinality.matrix = C.N.after.Cardinality.matrix[nrow(C.N.after.Cardinality.matrix):1,]
T.N.after.Cardinality.matrix = T.N.after.Cardinality.matrix[nrow(T.N.after.Cardinality.matrix):1,]
C.I.after.Cardinality.matrix = C.I.after.Cardinality.matrix[nrow(C.I.after.Cardinality.matrix):1,]
T.null.hypothesis.after.Cardinality.matrix = T.null.hypothesis.after.Cardinality.matrix[nrow(T.null.hypothesis.after.Cardinality.matrix):1,]
estimated.ATE.Cardinality.matrix = estimated.ATE.Cardinality.matrix[nrow(estimated.ATE.Cardinality.matrix):1,]
ATET.Cardinality.matrix = ATET.Cardinality.matrix[nrow(ATET.Cardinality.matrix):1,]
ATEC.Cardinality.matrix = ATEC.Cardinality.matrix[nrow(ATEC.Cardinality.matrix):1,]
T.weights.Cardinality.matrix = T.weights.Cardinality.matrix[nrow(T.weights.Cardinality.matrix):1,]
C.weights.Cardinality.matrix = C.weights.Cardinality.matrix[nrow(C.weights.Cardinality.matrix):1,]
null.hypothesis.estimated.ATE.Cardinality.matrix = null.hypothesis.estimated.ATE.Cardinality.matrix[nrow(null.hypothesis.estimated.ATE.Cardinality.matrix):1,]

blank.period.residuals.Cardinality.matrix = blank.period.residuals.Cardinality.matrix[nrow(blank.period.residuals.Cardinality.matrix):1,]
T.fitted.N.before.Cardinality.matrix = T.fitted.N.before.Cardinality.matrix[nrow(T.fitted.N.before.Cardinality.matrix):1,]
C.fitted.N.before.Cardinality.matrix = C.fitted.N.before.Cardinality.matrix[nrow(C.fitted.N.before.Cardinality.matrix):1,]

estimated.p.value.Cardinality.vector = estimated.p.value.Cardinality.vector[length(estimated.p.value.Cardinality.vector):1]
null.hypothesis.p.value.Cardinality.vector = null.hypothesis.p.value.Cardinality.vector[length(null.hypothesis.p.value.Cardinality.vector):1]

#=====Formulation: penalization=====#
T.I.after.Penalization.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.N.after.Penalization.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.N.after.Penalization.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.I.after.Penalization.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.null.hypothesis.after.Penalization.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
estimated.ATE.Penalization.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATET.Penalization.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATEC.Penalization.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.weights.Penalization.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
C.weights.Penalization.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
null.hypothesis.estimated.ATE.Penalization.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))

blank.period.residuals.Penalization.matrix = matrix(NA, nrow = 0, ncol = (T.naught - T.prime))
T.fitted.N.before.Penalization.matrix = matrix(NA, nrow = 0, ncol = T.naught)
C.fitted.N.before.Penalization.matrix = matrix(NA, nrow = 0, ncol = T.naught)

estimated.p.value.Penalization.vector = c()
null.hypothesis.p.value.Penalization.vector = c()
for(lambda.temp in c(1:length(changing.lambda)))
{
  tentative.T.weights = my.results.Penalization.changing.lambda[[lambda.temp]]$Treatment.weights
  tentative.C.weights = my.results.Penalization.changing.lambda[[lambda.temp]]$Control.weights
  
  if(nrow(T.weights.Penalization.matrix) > 0)
  {
    last.row.T.weights = T.weights.Penalization.matrix[nrow(T.weights.Penalization.matrix),]
    last.row.C.weights = C.weights.Penalization.matrix[nrow(C.weights.Penalization.matrix),]
  } else
  {
    last.row.T.weights = T.weights
    last.row.C.weights = C.weights
  }
  l0discrepancy.no.flip = l.zero.norm(last.row.T.weights, tentative.T.weights) + l.zero.norm(last.row.C.weights, tentative.C.weights)
  l0discrepancy.flip = l.zero.norm(last.row.C.weights, tentative.T.weights) + l.zero.norm(last.row.T.weights, tentative.C.weights)
  
  if(sum(tentative.T.weights != 0) == sum(tentative.C.weights != 0) && l0discrepancy.no.flip > l0discrepancy.flip)
  { #flip treatment and control weights
    use.T.weights = tentative.C.weights
    use.C.weights = tentative.T.weights
  } else
  {
    use.T.weights = tentative.T.weights
    use.C.weights = tentative.C.weights
  }
  
  T.weights.Penalization.matrix = rbind(T.weights.Penalization.matrix, use.T.weights)
  C.weights.Penalization.matrix = rbind(C.weights.Penalization.matrix, use.C.weights)
  
  T.I.after.Penalization.matrix = rbind(T.I.after.Penalization.matrix, t(use.T.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  C.N.after.Penalization.matrix = rbind(C.N.after.Penalization.matrix, t(use.C.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  T.N.after.Penalization.matrix = rbind(T.N.after.Penalization.matrix, t(use.T.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  C.I.after.Penalization.matrix = rbind(C.I.after.Penalization.matrix, t(use.C.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  T.null.hypothesis.after.Penalization.matrix = rbind(T.null.hypothesis.after.Penalization.matrix, t(use.T.weights) %*% Y.null.hypothesis.matrix[,(T.naught+1):T.total])
  
  estimated.ATE.Penalization.matrix = rbind(estimated.ATE.Penalization.matrix, T.I.after.Penalization.matrix[lambda.temp,] - C.N.after.Penalization.matrix[lambda.temp,])
  ATET.Penalization.matrix = rbind(ATET.Penalization.matrix, T.I.after.Penalization.matrix[lambda.temp,] - T.N.after.Penalization.matrix[lambda.temp,])
  ATEC.Penalization.matrix = rbind(ATEC.Penalization.matrix, C.I.after.Penalization.matrix[lambda.temp,] - C.N.after.Penalization.matrix[lambda.temp,])
  null.hypothesis.estimated.ATE.Penalization.matrix = rbind(null.hypothesis.estimated.ATE.Penalization.matrix, T.null.hypothesis.after.Penalization.matrix[lambda.temp,] - C.N.after.Penalization.matrix[lambda.temp,])
  
  T.fitted.N.before.Penalization.matrix = rbind(T.fitted.N.before.Penalization.matrix, t(use.T.weights) %*% Y.N.matrix[,1:T.naught])
  C.fitted.N.before.Penalization.matrix = rbind(C.fitted.N.before.Penalization.matrix, t(use.C.weights) %*% Y.N.matrix[,1:T.naught])
  blank.period.residuals.Penalization.matrix = rbind(blank.period.residuals.Penalization.matrix, T.fitted.N.before.Penalization.matrix[lambda.temp, (T.prime+1):T.naught] - C.fitted.N.before.Penalization.matrix[lambda.temp, (T.prime+1):T.naught])
  estimated.p.value.Penalization.vector = c(estimated.p.value.Penalization.vector, permutation.test(blank.period.residuals.Penalization.matrix[lambda.temp,], estimated.ATE.Penalization.matrix[lambda.temp,]))
  null.hypothesis.p.value.Penalization.vector = c(null.hypothesis.p.value.Penalization.vector, permutation.test(blank.period.residuals.Penalization.matrix[lambda.temp,], null.hypothesis.estimated.ATE.Penalization.matrix[lambda.temp,]))
}

#=====Formulation: unit level=====#
T.I.after.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.N.after.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.N.after.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.I.after.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.null.hypothesis.after.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
estimated.ATE.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATET.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATEC.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.weights.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
C.weights.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
null.hypothesis.estimated.ATE.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))

blank.period.residuals.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = (T.naught - T.prime))
T.fitted.N.before.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = T.naught)
C.fitted.N.before.UnitLevel.matrix = matrix(NA, nrow = 0, ncol = T.naught)

estimated.p.value.UnitLevel.vector = c()
null.hypothesis.p.value.UnitLevel.vector = c()
for(xi.temp in c(1:length(changing.xi)))
{
  T.weights.UnitLevel.matrix = rbind(T.weights.UnitLevel.matrix, my.results.UnitLevel.changing.xi[[xi.temp]]$Treatment.weights)
  C.weights.UnitLevel.matrix = rbind(C.weights.UnitLevel.matrix, my.results.UnitLevel.changing.xi[[xi.temp]]$Control.weights)
  
  T.I.after.UnitLevel.matrix = rbind(T.I.after.UnitLevel.matrix, t(my.results.UnitLevel.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  C.N.after.UnitLevel.matrix = rbind(C.N.after.UnitLevel.matrix, t(my.results.UnitLevel.changing.xi[[xi.temp]]$Control.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  T.N.after.UnitLevel.matrix = rbind(T.N.after.UnitLevel.matrix, t(my.results.UnitLevel.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  C.I.after.UnitLevel.matrix = rbind(C.I.after.UnitLevel.matrix, t(my.results.UnitLevel.changing.xi[[xi.temp]]$Control.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  T.null.hypothesis.after.UnitLevel.matrix = rbind(T.null.hypothesis.after.UnitLevel.matrix, t(my.results.UnitLevel.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.null.hypothesis.matrix[,(T.naught+1):T.total])
  
  estimated.ATE.UnitLevel.matrix = rbind(estimated.ATE.UnitLevel.matrix, T.I.after.UnitLevel.matrix[xi.temp,] - C.N.after.UnitLevel.matrix[xi.temp,])
  ATET.UnitLevel.matrix = rbind(ATET.UnitLevel.matrix, T.I.after.UnitLevel.matrix[xi.temp,] - T.N.after.UnitLevel.matrix[xi.temp,])
  ATEC.UnitLevel.matrix = rbind(ATEC.UnitLevel.matrix, C.I.after.UnitLevel.matrix[xi.temp,] - C.N.after.UnitLevel.matrix[xi.temp,])
  null.hypothesis.estimated.ATE.UnitLevel.matrix = rbind(null.hypothesis.estimated.ATE.UnitLevel.matrix, T.null.hypothesis.after.UnitLevel.matrix[xi.temp,] - C.N.after.UnitLevel.matrix[xi.temp,])
  
  T.fitted.N.before.UnitLevel.matrix = rbind(T.fitted.N.before.UnitLevel.matrix, t(my.results.UnitLevel.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.N.matrix[,1:T.naught])
  C.fitted.N.before.UnitLevel.matrix = rbind(C.fitted.N.before.UnitLevel.matrix, t(my.results.UnitLevel.changing.xi[[xi.temp]]$Control.weights) %*% Y.N.matrix[,1:T.naught])
  blank.period.residuals.UnitLevel.matrix = rbind(blank.period.residuals.UnitLevel.matrix, T.fitted.N.before.UnitLevel.matrix[xi.temp, (T.prime+1):T.naught] - C.fitted.N.before.UnitLevel.matrix[xi.temp, (T.prime+1):T.naught])
  estimated.p.value.UnitLevel.vector = c(estimated.p.value.UnitLevel.vector, permutation.test(blank.period.residuals.UnitLevel.matrix[xi.temp,], estimated.ATE.UnitLevel.matrix[xi.temp,]))
  null.hypothesis.p.value.UnitLevel.vector = c(null.hypothesis.p.value.UnitLevel.vector, permutation.test(blank.period.residuals.UnitLevel.matrix[xi.temp,], null.hypothesis.estimated.ATE.UnitLevel.matrix[xi.temp,]))
}

#=====Formulation: weakly targeted=====#
T.I.after.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.N.after.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.N.after.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.I.after.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.null.hypothesis.after.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
estimated.ATE.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATET.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATEC.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.weights.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
C.weights.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
null.hypothesis.estimated.ATE.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))

blank.period.residuals.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = (T.naught - T.prime))
T.fitted.N.before.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = T.naught)
C.fitted.N.before.WeaklyTargeted.matrix = matrix(NA, nrow = 0, ncol = T.naught)

estimated.p.value.WeaklyTargeted.vector = c()
null.hypothesis.p.value.WeaklyTargeted.vector = c()
for(xi.temp in c(1:length(changing.xi)))
{
  T.weights.WeaklyTargeted.matrix = rbind(T.weights.WeaklyTargeted.matrix, my.results.WeaklyTargeted.changing.xi[[xi.temp]]$Treatment.weights)
  C.weights.WeaklyTargeted.matrix = rbind(C.weights.WeaklyTargeted.matrix, my.results.WeaklyTargeted.changing.xi[[xi.temp]]$Control.weights)
  
  T.I.after.WeaklyTargeted.matrix = rbind(T.I.after.WeaklyTargeted.matrix, t(my.results.WeaklyTargeted.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  C.N.after.WeaklyTargeted.matrix = rbind(C.N.after.WeaklyTargeted.matrix, t(my.results.WeaklyTargeted.changing.xi[[xi.temp]]$Control.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  T.N.after.WeaklyTargeted.matrix = rbind(T.N.after.WeaklyTargeted.matrix, t(my.results.WeaklyTargeted.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  C.I.after.WeaklyTargeted.matrix = rbind(C.I.after.WeaklyTargeted.matrix, t(my.results.WeaklyTargeted.changing.xi[[xi.temp]]$Control.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  T.null.hypothesis.after.WeaklyTargeted.matrix = rbind(T.null.hypothesis.after.WeaklyTargeted.matrix, t(my.results.WeaklyTargeted.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.null.hypothesis.matrix[,(T.naught+1):T.total])
  
  estimated.ATE.WeaklyTargeted.matrix = rbind(estimated.ATE.WeaklyTargeted.matrix, T.I.after.WeaklyTargeted.matrix[xi.temp,] - C.N.after.WeaklyTargeted.matrix[xi.temp,])
  ATET.WeaklyTargeted.matrix = rbind(ATET.WeaklyTargeted.matrix, T.I.after.WeaklyTargeted.matrix[xi.temp,] - T.N.after.WeaklyTargeted.matrix[xi.temp,])
  ATEC.WeaklyTargeted.matrix = rbind(ATEC.WeaklyTargeted.matrix, C.I.after.WeaklyTargeted.matrix[xi.temp,] - C.N.after.WeaklyTargeted.matrix[xi.temp,])
  null.hypothesis.estimated.ATE.WeaklyTargeted.matrix = rbind(null.hypothesis.estimated.ATE.WeaklyTargeted.matrix, T.null.hypothesis.after.WeaklyTargeted.matrix[xi.temp,] - C.N.after.WeaklyTargeted.matrix[xi.temp,])
  
  T.fitted.N.before.WeaklyTargeted.matrix = rbind(T.fitted.N.before.WeaklyTargeted.matrix, t(my.results.WeaklyTargeted.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.N.matrix[,1:T.naught])
  C.fitted.N.before.WeaklyTargeted.matrix = rbind(C.fitted.N.before.WeaklyTargeted.matrix, t(my.results.WeaklyTargeted.changing.xi[[xi.temp]]$Control.weights) %*% Y.N.matrix[,1:T.naught])
  blank.period.residuals.WeaklyTargeted.matrix = rbind(blank.period.residuals.WeaklyTargeted.matrix, T.fitted.N.before.WeaklyTargeted.matrix[xi.temp, (T.prime+1):T.naught] - C.fitted.N.before.WeaklyTargeted.matrix[xi.temp, (T.prime+1):T.naught])
  estimated.p.value.WeaklyTargeted.vector = c(estimated.p.value.WeaklyTargeted.vector, permutation.test(blank.period.residuals.WeaklyTargeted.matrix[xi.temp,], estimated.ATE.WeaklyTargeted.matrix[xi.temp,]))
  null.hypothesis.p.value.WeaklyTargeted.vector = c(null.hypothesis.p.value.WeaklyTargeted.vector, permutation.test(blank.period.residuals.WeaklyTargeted.matrix[xi.temp,], null.hypothesis.estimated.ATE.WeaklyTargeted.matrix[xi.temp,]))
}



#======================================#
#=====Print and save to .txt files=====#
#======================================#

p.value.rejection.criteria = 0.05

#Average treatment effect
ATE.report = rbind(true.ATE,
                   estimated.ATE,
                   estimated.ATE.Cardinality.matrix,
                   estimated.ATE.WeaklyTargeted.matrix,
                   estimated.ATE.UnitLevel.matrix,
                   estimated.ATE.Penalization.matrix)
Loss.MAE.vec = c(0)
Loss.MSE.vec = c(0)
for(loss.temp in 2:nrow(ATE.report))
{
  Loss.MAE.vec = c(Loss.MAE.vec, sum(abs(ATE.report[loss.temp,] - ATE.report[1,])) / (T.total - T.naught))
}
for(loss.temp in 2:nrow(ATE.report))
{
  Loss.MSE.vec = c(Loss.MSE.vec, sum((ATE.report[loss.temp,] - ATE.report[1,]) * (ATE.report[loss.temp,] - ATE.report[1,])) / (T.total - T.naught))
}
estimated.p.values.report = c(1,
                              estimated.p.value,
                              estimated.p.value.Cardinality.vector,
                              estimated.p.value.WeaklyTargeted.vector,
                              estimated.p.value.UnitLevel.vector,
                              estimated.p.value.Penalization.vector)
rejection.report = as.numeric(estimated.p.values.report < p.value.rejection.criteria)
ATE.report.print = cbind(ATE.report, Loss.MAE.vec, Loss.MSE.vec)
my.file.name.copy = as.character(paste0("output_nonlinear/ ", repetition.RANDOM.SEED, "Different_Formulations_ATE.txt"))
write.table(cbind(round(ATE.report.print, digits = 3), estimated.p.values.report, rejection.report), file = my.file.name.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Average treatment effect if the outcome models are generated under the sharp null hypothesis
null.hypothesis.ATE.report = rbind(null.hypothesis.true.ATE,
                                   null.hypothesis.estimated.ATE,
                                   null.hypothesis.estimated.ATE.Cardinality.matrix,
                                   null.hypothesis.estimated.ATE.WeaklyTargeted.matrix,
                                   null.hypothesis.estimated.ATE.UnitLevel.matrix,
                                   null.hypothesis.estimated.ATE.Penalization.matrix)
null.hypothesis.Loss.MAE.vec = c(0)
null.hypothesis.Loss.MSE.vec = c(0)
for(loss.temp in 2:nrow(null.hypothesis.ATE.report))
{
  null.hypothesis.Loss.MAE.vec = c(null.hypothesis.Loss.MAE.vec, sum(abs(null.hypothesis.ATE.report[loss.temp,] - null.hypothesis.ATE.report[1,])) / (T.total - T.naught))
}
for(loss.temp in 2:nrow(null.hypothesis.ATE.report))
{
  null.hypothesis.Loss.MSE.vec = c(null.hypothesis.Loss.MSE.vec, sum((null.hypothesis.ATE.report[loss.temp,] - null.hypothesis.ATE.report[1,]) * (null.hypothesis.ATE.report[loss.temp,] - null.hypothesis.ATE.report[1,])) / (T.total - T.naught))
}
null.hypothesis.p.values.report = c(1,
                                    null.hypothesis.p.value,
                                    null.hypothesis.p.value.Cardinality.vector,
                                    null.hypothesis.p.value.WeaklyTargeted.vector,
                                    null.hypothesis.p.value.UnitLevel.vector,
                                    null.hypothesis.p.value.Penalization.vector)
null.hypothesis.rejection.report = as.numeric(null.hypothesis.p.values.report < p.value.rejection.criteria)
null.hypothesis.ATE.report = cbind(null.hypothesis.ATE.report, null.hypothesis.Loss.MAE.vec, null.hypothesis.Loss.MSE.vec)
my.file.name.null.hypothesis.copy = as.character(paste0("output_nonlinear/ ", repetition.RANDOM.SEED, "Different_Formulations_ATE_null_hypothesis.txt"))
write.table(cbind(round(null.hypothesis.ATE.report, digits = 3), null.hypothesis.p.values.report, null.hypothesis.rejection.report), file = my.file.name.null.hypothesis.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Average effect of treatment on the treated
ATET.report = rbind(ATET.T.weighted,
                    ATET.Cardinality.matrix,
                    ATET.WeaklyTargeted.matrix,
                    ATET.UnitLevel.matrix,
                    ATET.Penalization.matrix)
ATET.Loss.MAE.vec = c()
ATET.Loss.MSE.vec = c()
for(loss.temp in 1:nrow(ATET.report))
{
  ATET.Loss.MAE.vec = c(ATET.Loss.MAE.vec, sum(abs(ATE.report[loss.temp+1,] - ATET.report[loss.temp,])) / (T.total - T.naught))
}
for(loss.temp in 1:nrow(ATET.report))
{
  ATET.Loss.MSE.vec = c(ATET.Loss.MSE.vec, sum((ATE.report[loss.temp+1,] - ATET.report[loss.temp,]) * (ATE.report[loss.temp+1,] - ATET.report[loss.temp,])) / (T.total - T.naught))
}
ATET.report.print = cbind(ATET.report, ATE.report[-1,], ATET.Loss.MAE.vec, ATET.Loss.MSE.vec)
my.file.name.ATET.copy = as.character(paste0("output_nonlinear/ ", repetition.RANDOM.SEED, "Different_Formulations_ATET.txt"))
write.table(round(ATET.report.print, digits = 3), file = my.file.name.ATET.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Average effect of treatment on the control
ATEC.report = rbind(ATEC.C.weighted,
                    ATEC.Cardinality.matrix,
                    ATEC.WeaklyTargeted.matrix,
                    ATEC.UnitLevel.matrix,
                    ATEC.Penalization.matrix)
ATEC.Loss.MAE.vec = c()
ATEC.Loss.MSE.vec = c()
for(loss.temp in 1:nrow(ATEC.report))
{
  ATEC.Loss.MAE.vec = c(ATEC.Loss.MAE.vec, sum(abs(ATE.report[loss.temp+1,] - ATEC.report[loss.temp,])) / (T.total - T.naught))
}
for(loss.temp in 1:nrow(ATEC.report))
{
  ATEC.Loss.MSE.vec = c(ATEC.Loss.MSE.vec, sum((ATE.report[loss.temp+1,] - ATEC.report[loss.temp,]) * (ATE.report[loss.temp+1,] - ATEC.report[loss.temp,])) / (T.total - T.naught))
}
ATEC.report.print = cbind(ATEC.report, ATE.report[-1,], ATEC.Loss.MAE.vec, ATEC.Loss.MSE.vec)
my.file.name.ATEC.copy = as.character(paste0("output_nonlinear/ ", repetition.RANDOM.SEED, "Different_Formulations_ATEC.txt"))
write.table(round(ATEC.report.print, digits = 3), file = my.file.name.ATEC.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Treatment weights
Treatment.weights.report = rbind(T.weights,
                                 T.weights.Cardinality.matrix,
                                 T.weights.WeaklyTargeted.matrix,
                                 T.weights.UnitLevel.matrix,
                                 T.weights.Penalization.matrix)
my.file.name.TreatmentWeights.copy = as.character(paste0("output_nonlinear/ ", repetition.RANDOM.SEED, "Different_Formulations_TreatmentWeights.txt"))
write.table(round(Treatment.weights.report, digits = Weight.digits), file = my.file.name.TreatmentWeights.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Control weights
Control.weights.report = rbind(C.weights,
                               C.weights.Cardinality.matrix,
                               C.weights.WeaklyTargeted.matrix,
                               C.weights.UnitLevel.matrix,
                               C.weights.Penalization.matrix)
my.file.name.ControlWeights.copy = as.character(paste0("output_nonlinear/ ", repetition.RANDOM.SEED, "Different_Formulations_ControlWeights.txt"))
write.table(round(Control.weights.report, digits = Weight.digits), file = my.file.name.ControlWeights.copy, append = FALSE, sep = "\t", col.names=FALSE)

end_time <- Sys.time()
end_time - start_time








#=======================================================#
#===== TASK 2: Compare with randomized experiments =====#
#===== random assignments + different estimators =======#
#=======================================================#

#--- Last Change: 20250831 10AM added nonlinear data generating process to run this simulation
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

library(stats)
library(GA)

#================================================#
#===== Randomization 1: naive randomization =====#
#================================================#

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
  
  #--- the number of treated units is less than or equal to K.cardinality
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
  
  #--- fix the number of treated units to be equal to K.cardinality
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


#=====================================================#
#===== Randomization 2: stratified randomization =====#
#=====================================================#

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
    fitness = objective_function_GA_alternative, # Objective function
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
    my.blocks = my_MinMaxDiameter_blocking(random.seed_ = random.seed_, N.Regions_ = N.Regions_, feature.matrix_ = feature.matrix_, K.cardinality_ = K.cardinality)
    # my.blocks = my_k_means(random.seed_ = random.seed_, N.Regions_ = N.Regions_, feature.matrix_ = feature.matrix_, K.cardinality_ = K.cardinality)
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



#=============================================================#
#===== Generate random assignments and fetch the outputs =====#
#=============================================================#

#---------------------------------------------------------------#
#--- Two options here: whether we use pre-treatment outcomes ---#
#--- in stratified randomization and matching or not         ---#
#---------------------------------------------------------------#

my.feature.matrix_StA_Matching = cbind(t(Z.ob.covariates.matrix), Y.N.matrix[,1:T.naught])
# my.feature.matrix_StA_Matching = t(Z.ob.covariates.matrix)
my.feature.matrix_StA_Matching = apply(my.feature.matrix_StA_Matching, 2, function(x) (x - mean(x)) / sd(x))

#---------------------------------------#
#--- Now generate random assignments ---#
#---------------------------------------#

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

#-----------------------------#
#--- Now fetch the outputs ---#
#-----------------------------#

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
#--- regression adjustment can't use pre-treatment outcomes (dimension too large) ---#
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
my.file.name.copy = as.character(paste0("output_nonlinear/ ", repetition.RANDOM.SEED, "RandomAssignment_RADiM.txt"))
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
my.file.name.copy = as.character(paste0("output_nonlinear/ ", repetition.RANDOM.SEED, "StratifiedAssignment_StADiM.txt"))
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
my.file.name.copy = as.character(paste0("output_nonlinear/ ", repetition.RANDOM.SEED, "RandomAssignment_RARegAdj.txt"))
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
my.file.name.copy = as.character(paste0("output_nonlinear/ ", repetition.RANDOM.SEED, "RandomAssignment_RA1NN.txt"))
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
my.file.name.copy = as.character(paste0("output_nonlinear/ ", repetition.RANDOM.SEED, "RandomAssignment_RA5NN.txt"))
write.table(round(ATE.RA5NN.report.print, digits = 3), file = my.file.name.copy, append = FALSE, sep = "\t", col.names=FALSE)


end_time <- Sys.time()
end_time - start_time


