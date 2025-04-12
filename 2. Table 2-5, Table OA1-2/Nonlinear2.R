#===============================================================#
#========== Several different methods to find optimum ==========#
#===============================================================#
#========== In this file, many formulations used open ==========#
#========== source solvers since Gurobi solver cannot ==========#
#========== be called too often in cluster computing  ==========#
#===============================================================#

#--- Gurobi parameter: Threads = 1

#--- Last Change: 20221129 1PM --- a different nonlinear setup using exponential function, instead of gaussian pdf
#--- Last Change: 20221130 10PM --- Line 67 changed range.coefficients.max from 10 to 3 to reduce the magnitude of errors
#--- Last Change: 20221201 10PM --- added MSE (4 occasions in the outputs), put the new formulation earlier
#--- Last Change: 20221212 10AM --- added controllable digits and trimming 
#--- Last Change: 20250316 2PM: --- set the weights in the estimand (beta.vector) to be uniform weights 1/N.Regions
#---------------------------------- for simplicity, we do not use print_model.primitives to print the model primitives


# setwd("/project/jinglonggrp/Synthetic Experiment/")
# repetition.RANDOM.SEED = 415
repetition.RANDOM.SEED = as.numeric(Sys.getenv("SGE_TASK_ID"))


library(slam)
library(dplyr)
library(gurobi)
library(Matrix)
library(limSolve)
library(quadprog)

#=====Define variables=====#
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
range.coefficients.max = 3 #ChangeThis #20221130 10PM changed this from 10 to 3
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

Generate_Nonlinear2_Primitives <- function(range.intercept.max_ = range.intercept.max,
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

print_model.primitives <- function(model.primitives_, random.seed_)
{
  my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "ModelPrimitives_beta.txt"))
  write.table(t(c(random.seed_, model.primitives_$beta.vector)), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "ModelPrimitives_delta.txt"))
  write.table(t(c(random.seed_, model.primitives_$delta.N.vector)), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "ModelPrimitives_upsilon.txt"))
  write.table(t(c(random.seed_, model.primitives_$upsilon.I.vector)), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "ModelPrimitives_Z.txt"))
  write.table(cbind(rep(random.seed_, nrow(Z.ob.covariates.matrix)), model.primitives_$Z.ob.covariates.matrix), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "ModelPrimitives_mu.txt"))
  write.table(cbind(rep(random.seed_, nrow(mu.unob.covariates.matrix)), model.primitives_$mu.unob.covariates.matrix), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "ModelPrimitives_theta.txt"))
  write.table(cbind(rep(random.seed_, nrow(theta.ob.N.matrix)), model.primitives_$theta.ob.N.matrix), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "ModelPrimitives_gamma.txt"))
  write.table(cbind(rep(random.seed_, nrow(gamma.ob.I.matrix)), model.primitives_$gamma.ob.I.matrix), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "ModelPrimitives_lambda.txt"))
  write.table(cbind(rep(random.seed_, nrow(lambda.unob.N.matrix)), model.primitives_$lambda.unob.N.matrix), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "ModelPrimitives_eta.txt"))
  write.table(cbind(rep(random.seed_, nrow(eta.unob.I.matrix)), model.primitives_$eta.unob.I.matrix), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "ModelPrimitives_epsilon.txt"))
  write.table(cbind(rep(random.seed_, nrow(epsilon.N.matrix)), model.primitives_$epsilon.N.matrix), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "ModelPrimitives_xi.txt"))
  write.table(cbind(rep(random.seed_, nrow(xi.I.matrix)), model.primitives_$xi.I.matrix), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "ModelPrimitives_YN.txt"))
  write.table(cbind(rep(random.seed_, nrow(Y.N.matrix)), model.primitives_$Y.N.matrix), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
  
  my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "ModelPrimitives_YI.txt"))
  write.table(cbind(rep(random.seed_, nrow(Y.I.matrix)), model.primitives_$Y.I.matrix), file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
}

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

#======================================#
#=====The New Formulation for ATET=====#
#======================================#

Synthetic_Experiment_For_ATET <- function(T.prime_, r.ob.covariates.dim_, N.Regions_, Y.N.matrix, Z.ob.covariates.matrix, beta.vector_, xi_ = 1, digits_ = Weight.digits)
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

#----------------------------#
#-----Easier Formulation-----#
#-----No linear terms--------#
#----------------------------#

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

#Formulation (11) in the paper
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

#Formulation (13)
Synthetic_Experiment_Asymmetric <- function(T.prime_, r.ob.covariates.dim_, N.Regions_, Y.N.matrix_, Z.ob.covariates.matrix_, beta.vector_, xi_ = 1, K.cardinality_ = -1, digits_ = Weight.digits)
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


start_time <- Sys.time()

#====================================#
#===== Repeating multiple times =====#
#====================================#

repetition.model.primitives.SAVED = list()
repetition.output.Original.Formulation.SAVED = list()
repetition.output.Cardinality.SAVED = list()
repetition.output.Penalization.SAVED = list()
repetition.output.Asymmetric.SAVED = list()
repetition.output.ForATET.SAVED = list()
repetition.reported.Original.Formulation.SAVED = matrix(NA, nrow = 0, ncol = (4*(T.total-T.naught) + (T.naught - T.prime) + 4))


#=====================================#
#===== Generate Model Primitives =====#
#=====================================#
model.primitives = Generate_Nonlinear2_Primitives(random.seed_ = repetition.RANDOM.SEED)
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
#-----Save and print-----#
# print_model.primitives(model.primitives, random.seed_ = repetition.RANDOM.SEED)
repetition.model.primitives.SAVED[[repetition.RANDOM.SEED]] = model.primitives


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
changing.lambda = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)
for(lambda.temp in 1:length(changing.lambda))
{
  lambda.this.round = changing.lambda[lambda.temp]
  my.results.Penalization.changing.lambda[[lambda.temp]] = Synthetic_Experiment_Penalization(T.prime, r.ob.covariates.dim, N.Regions, Y.N.matrix, Z.ob.covariates.matrix, beta.vector, lambda_ = lambda.this.round)
}
#-----Save and print-----#
repetition.output.Penalization.SAVED[[repetition.RANDOM.SEED]] = my.results.Penalization.changing.lambda

#=====Formulation: asymmetric=====#
my.results.Asymmetric.changing.xi = list()
changing.xi = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)
for(xi.temp in 1:length(changing.xi))
{
  xi.this.round = changing.xi[xi.temp]
  my.results.Asymmetric.changing.xi[[xi.temp]] = Synthetic_Experiment_Asymmetric(T.prime, r.ob.covariates.dim, N.Regions, Y.N.matrix, Z.ob.covariates.matrix, beta.vector, xi_ = xi.this.round, K.cardinality_ = -1)
}
#-----Save and print-----#
repetition.output.Asymmetric.SAVED[[repetition.RANDOM.SEED]] = my.results.Asymmetric.changing.xi

#=====Formulation: ForATET=====#
my.results.ForATET.changing.xi = list()
changing.xi = c(0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000)
for(xi.temp in 1:length(changing.xi))
{
  xi.this.round = changing.xi[xi.temp]
  my.results.ForATET.changing.xi[[xi.temp]] = Synthetic_Experiment_For_ATET(T.prime, r.ob.covariates.dim, N.Regions, Y.N.matrix, Z.ob.covariates.matrix, beta.vector, xi_ = xi.this.round)
}
#-----Save and print-----#
repetition.output.ForATET.SAVED[[repetition.RANDOM.SEED]] = my.results.ForATET.changing.xi




end_time <- Sys.time()
end_time - start_time

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

#-----Save and print-----#
repetition.reported.Original.Formulation.SAVED = rbind(repetition.reported.Original.Formulation.SAVED, c(repetition.RANDOM.SEED, true.ATE, estimated.ATE, blank.period.residuals, MeanAbsoluteError, estimated.p.value, null.hypothesis.true.ATE, null.hypothesis.estimated.ATE, null.hypothesis.p.value))
my.file.name = as.character(paste0("output11_nonlinear/ ", Sys.Date(), "SyntheticExperiment.txt"))
write.table(t(c(repetition.RANDOM.SEED, true.ATE, estimated.ATE, blank.period.residuals, MeanAbsoluteError, estimated.p.value, null.hypothesis.true.ATE, null.hypothesis.estimated.ATE, null.hypothesis.p.value)), file = my.file.name, append = TRUE, sep = "\t", col.names=FALSE)

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

#=====Formulation: asymmetric=====#
T.I.after.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.N.after.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.N.after.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.I.after.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.null.hypothesis.after.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
estimated.ATE.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATET.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATEC.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.weights.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
C.weights.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
null.hypothesis.estimated.ATE.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))

blank.period.residuals.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = (T.naught - T.prime))
T.fitted.N.before.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = T.naught)
C.fitted.N.before.Asymmetric.matrix = matrix(NA, nrow = 0, ncol = T.naught)

estimated.p.value.Asymmetric.vector = c()
null.hypothesis.p.value.Asymmetric.vector = c()
for(xi.temp in c(1:length(changing.xi)))
{
  T.weights.Asymmetric.matrix = rbind(T.weights.Asymmetric.matrix, my.results.Asymmetric.changing.xi[[xi.temp]]$Treatment.weights)
  C.weights.Asymmetric.matrix = rbind(C.weights.Asymmetric.matrix, my.results.Asymmetric.changing.xi[[xi.temp]]$Control.weights)
  
  T.I.after.Asymmetric.matrix = rbind(T.I.after.Asymmetric.matrix, t(my.results.Asymmetric.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  C.N.after.Asymmetric.matrix = rbind(C.N.after.Asymmetric.matrix, t(my.results.Asymmetric.changing.xi[[xi.temp]]$Control.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  T.N.after.Asymmetric.matrix = rbind(T.N.after.Asymmetric.matrix, t(my.results.Asymmetric.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  C.I.after.Asymmetric.matrix = rbind(C.I.after.Asymmetric.matrix, t(my.results.Asymmetric.changing.xi[[xi.temp]]$Control.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  T.null.hypothesis.after.Asymmetric.matrix = rbind(T.null.hypothesis.after.Asymmetric.matrix, t(my.results.Asymmetric.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.null.hypothesis.matrix[,(T.naught+1):T.total])
  
  estimated.ATE.Asymmetric.matrix = rbind(estimated.ATE.Asymmetric.matrix, T.I.after.Asymmetric.matrix[xi.temp,] - C.N.after.Asymmetric.matrix[xi.temp,])
  ATET.Asymmetric.matrix = rbind(ATET.Asymmetric.matrix, T.I.after.Asymmetric.matrix[xi.temp,] - T.N.after.Asymmetric.matrix[xi.temp,])
  ATEC.Asymmetric.matrix = rbind(ATEC.Asymmetric.matrix, C.I.after.Asymmetric.matrix[xi.temp,] - C.N.after.Asymmetric.matrix[xi.temp,])
  null.hypothesis.estimated.ATE.Asymmetric.matrix = rbind(null.hypothesis.estimated.ATE.Asymmetric.matrix, T.null.hypothesis.after.Asymmetric.matrix[xi.temp,] - C.N.after.Asymmetric.matrix[xi.temp,])
  
  T.fitted.N.before.Asymmetric.matrix = rbind(T.fitted.N.before.Asymmetric.matrix, t(my.results.Asymmetric.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.N.matrix[,1:T.naught])
  C.fitted.N.before.Asymmetric.matrix = rbind(C.fitted.N.before.Asymmetric.matrix, t(my.results.Asymmetric.changing.xi[[xi.temp]]$Control.weights) %*% Y.N.matrix[,1:T.naught])
  blank.period.residuals.Asymmetric.matrix = rbind(blank.period.residuals.Asymmetric.matrix, T.fitted.N.before.Asymmetric.matrix[xi.temp, (T.prime+1):T.naught] - C.fitted.N.before.Asymmetric.matrix[xi.temp, (T.prime+1):T.naught])
  estimated.p.value.Asymmetric.vector = c(estimated.p.value.Asymmetric.vector, permutation.test(blank.period.residuals.Asymmetric.matrix[xi.temp,], estimated.ATE.Asymmetric.matrix[xi.temp,]))
  null.hypothesis.p.value.Asymmetric.vector = c(null.hypothesis.p.value.Asymmetric.vector, permutation.test(blank.period.residuals.Asymmetric.matrix[xi.temp,], null.hypothesis.estimated.ATE.Asymmetric.matrix[xi.temp,]))
}

#=====Formulation: ForATET=====#
T.I.after.ForATET.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.N.after.ForATET.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.N.after.ForATET.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.I.after.ForATET.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.null.hypothesis.after.ForATET.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
estimated.ATE.ForATET.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATET.ForATET.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
ATEC.ForATET.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.weights.ForATET.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
C.weights.ForATET.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
null.hypothesis.estimated.ATE.ForATET.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))

blank.period.residuals.ForATET.matrix = matrix(NA, nrow = 0, ncol = (T.naught - T.prime))
T.fitted.N.before.ForATET.matrix = matrix(NA, nrow = 0, ncol = T.naught)
C.fitted.N.before.ForATET.matrix = matrix(NA, nrow = 0, ncol = T.naught)

estimated.p.value.ForATET.vector = c()
null.hypothesis.p.value.ForATET.vector = c()
for(xi.temp in c(1:length(changing.xi)))
{
  T.weights.ForATET.matrix = rbind(T.weights.ForATET.matrix, my.results.ForATET.changing.xi[[xi.temp]]$Treatment.weights)
  C.weights.ForATET.matrix = rbind(C.weights.ForATET.matrix, my.results.ForATET.changing.xi[[xi.temp]]$Control.weights)
  
  T.I.after.ForATET.matrix = rbind(T.I.after.ForATET.matrix, t(my.results.ForATET.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  C.N.after.ForATET.matrix = rbind(C.N.after.ForATET.matrix, t(my.results.ForATET.changing.xi[[xi.temp]]$Control.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  T.N.after.ForATET.matrix = rbind(T.N.after.ForATET.matrix, t(my.results.ForATET.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  C.I.after.ForATET.matrix = rbind(C.I.after.ForATET.matrix, t(my.results.ForATET.changing.xi[[xi.temp]]$Control.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  T.null.hypothesis.after.ForATET.matrix = rbind(T.null.hypothesis.after.ForATET.matrix, t(my.results.ForATET.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.null.hypothesis.matrix[,(T.naught+1):T.total])
  
  estimated.ATE.ForATET.matrix = rbind(estimated.ATE.ForATET.matrix, T.I.after.ForATET.matrix[xi.temp,] - C.N.after.ForATET.matrix[xi.temp,])
  ATET.ForATET.matrix = rbind(ATET.ForATET.matrix, T.I.after.ForATET.matrix[xi.temp,] - T.N.after.ForATET.matrix[xi.temp,])
  ATEC.ForATET.matrix = rbind(ATEC.ForATET.matrix, C.I.after.ForATET.matrix[xi.temp,] - C.N.after.ForATET.matrix[xi.temp,])
  null.hypothesis.estimated.ATE.ForATET.matrix = rbind(null.hypothesis.estimated.ATE.ForATET.matrix, T.null.hypothesis.after.ForATET.matrix[xi.temp,] - C.N.after.ForATET.matrix[xi.temp,])
  
  T.fitted.N.before.ForATET.matrix = rbind(T.fitted.N.before.ForATET.matrix, t(my.results.ForATET.changing.xi[[xi.temp]]$Treatment.weights) %*% Y.N.matrix[,1:T.naught])
  C.fitted.N.before.ForATET.matrix = rbind(C.fitted.N.before.ForATET.matrix, t(my.results.ForATET.changing.xi[[xi.temp]]$Control.weights) %*% Y.N.matrix[,1:T.naught])
  blank.period.residuals.ForATET.matrix = rbind(blank.period.residuals.ForATET.matrix, T.fitted.N.before.ForATET.matrix[xi.temp, (T.prime+1):T.naught] - C.fitted.N.before.ForATET.matrix[xi.temp, (T.prime+1):T.naught])
  estimated.p.value.ForATET.vector = c(estimated.p.value.ForATET.vector, permutation.test(blank.period.residuals.ForATET.matrix[xi.temp,], estimated.ATE.ForATET.matrix[xi.temp,]))
  null.hypothesis.p.value.ForATET.vector = c(null.hypothesis.p.value.ForATET.vector, permutation.test(blank.period.residuals.ForATET.matrix[xi.temp,], null.hypothesis.estimated.ATE.ForATET.matrix[xi.temp,]))
}



#======================================#
#=====Print and save to .txt files=====#
#======================================#

p.value.rejection.criteria = 0.05

#Average treatment effect
ATE.report = rbind(true.ATE,
                   estimated.ATE,
                   estimated.ATE.Cardinality.matrix,
                   estimated.ATE.ForATET.matrix,
                   estimated.ATE.Asymmetric.matrix,
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
                              estimated.p.value.ForATET.vector,
                              estimated.p.value.Asymmetric.vector,
                              estimated.p.value.Penalization.vector)
rejection.report = as.numeric(estimated.p.values.report < p.value.rejection.criteria)
ATE.report.print = cbind(ATE.report, Loss.MAE.vec, Loss.MSE.vec)
my.file.name.copy = as.character(paste0("output11_nonlinear/ ", repetition.RANDOM.SEED, "Different_Optimization_Methods_ATE.txt"))
write.table(cbind(round(ATE.report.print, digits = 3), estimated.p.values.report, rejection.report), file = my.file.name.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Average treatment effect if the outcome models are generated under the sharp null hypothesis
null.hypothesis.ATE.report = rbind(null.hypothesis.true.ATE,
                                   null.hypothesis.estimated.ATE,
                                   null.hypothesis.estimated.ATE.Cardinality.matrix,
                                   null.hypothesis.estimated.ATE.ForATET.matrix,
                                   null.hypothesis.estimated.ATE.Asymmetric.matrix,
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
                                    null.hypothesis.p.value.ForATET.vector,
                                    null.hypothesis.p.value.Asymmetric.vector,
                                    null.hypothesis.p.value.Penalization.vector)
null.hypothesis.rejection.report = as.numeric(null.hypothesis.p.values.report < p.value.rejection.criteria)
null.hypothesis.ATE.report = cbind(null.hypothesis.ATE.report, null.hypothesis.Loss.MAE.vec, null.hypothesis.Loss.MSE.vec)
my.file.name.null.hypothesis.copy = as.character(paste0("output11_nonlinear/ ", repetition.RANDOM.SEED, "Different_Optimization_Methods_ATE_null_hypothesis.txt"))
write.table(cbind(round(null.hypothesis.ATE.report, digits = 3), null.hypothesis.p.values.report, null.hypothesis.rejection.report), file = my.file.name.null.hypothesis.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Average effect of treatment on the treated
ATET.report = rbind(ATET.T.weighted,
                    ATET.Cardinality.matrix,
                    ATET.ForATET.matrix,
                    ATET.Asymmetric.matrix,
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
my.file.name.ATET.copy = as.character(paste0("output11_nonlinear/ ", repetition.RANDOM.SEED, "Different_Optimization_Methods_ATET.txt"))
write.table(round(ATET.report.print, digits = 3), file = my.file.name.ATET.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Average effect of treatment on the control
ATEC.report = rbind(ATEC.C.weighted,
                    ATEC.Cardinality.matrix,
                    ATEC.ForATET.matrix,
                    ATEC.Asymmetric.matrix,
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
my.file.name.ATEC.copy = as.character(paste0("output11_nonlinear/ ", repetition.RANDOM.SEED, "Different_Optimization_Methods_ATEC.txt"))
write.table(round(ATEC.report.print, digits = 3), file = my.file.name.ATEC.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Treatment weights
Treatment.weights.report = rbind(T.weights,
                                 T.weights.Cardinality.matrix,
                                 T.weights.ForATET.matrix,
                                 T.weights.Asymmetric.matrix,
                                 T.weights.Penalization.matrix)
my.file.name.TreatmentWeights.copy = as.character(paste0("output11_nonlinear/ ", repetition.RANDOM.SEED, "Different_Optimization_Methods_TreatmentWeights.txt"))
write.table(round(Treatment.weights.report, digits = Weight.digits), file = my.file.name.TreatmentWeights.copy, append = FALSE, sep = "\t", col.names=FALSE)

#Control weights
Control.weights.report = rbind(C.weights,
                               C.weights.Cardinality.matrix,
                               C.weights.ForATET.matrix,
                               C.weights.Asymmetric.matrix,
                               C.weights.Penalization.matrix)
my.file.name.ControlWeights.copy = as.character(paste0("output11_nonlinear/ ", repetition.RANDOM.SEED, "Different_Optimization_Methods_ControlWeights.txt"))
write.table(round(Control.weights.report, digits = Weight.digits), file = my.file.name.ControlWeights.copy, append = FALSE, sep = "\t", col.names=FALSE)

end_time <- Sys.time()
end_time - start_time

#=====Debugging use=====#
# temp_time <- Sys.time()
# print("*********************************************")
# print("*********************************************")
# temp_time - start_time
# print("*********************************************")
# print("*********************************************")