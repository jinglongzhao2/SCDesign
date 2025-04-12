#==========================================================#
#========== Read and process Walmart Kaggle data ==========#
#==========================================================#

#--------------------------------------------#
#--- Use Lines 722 - 729 to run the codes ---#
#--------------------------------------------#

repetition.SEED = 1

library(tidyr)
library(dplyr)
library(slam)
library(gurobi)
library(Matrix)
library(limSolve)
library(matrixcalc)
library(quadprog)
library(gamlr)
library(EnvStats)

set.seed(123456)

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

#for confidence intervals
alpha.CI = 0.05

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



#=======================================================================================#
#========== Below we apply synthetic control methods to design the experiment ==========#
#=======================================================================================#

#==================================#
#=====The Original Formulation=====#
#==================================#

#Formulation Unconstrained in the paper
Synthetic_Experiment <- function(T.prime_, r.ob.covariates.dim_, N.Regions_, Y.N.matrix, Z.ob.covariates.matrix, beta.vector_, digits_ = Weight.digits, TimeLimit_ = 3600)
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
  params = list(Threads = 1, TimeLimit = TimeLimit_, NonConvex=2) #2 stands for non-convex optimization
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

Synthetic_Control <- function(target.vector_, X.matrix_)
{
  X.effective = X.matrix_ - matrix(data = target.vector_, nrow = nrow(X.matrix_), ncol = ncol(X.matrix_), byrow = FALSE)
  
  # #=====================================================================#
  # #=============================== NOTICE ==============================#
  # #=====Gurobi needs to check the license every time it is called. =====#
  # #=====And Gurobi does not like to be called too frequently. ==========#
  # #=====Cannot use Gurobi to solve many simple optimization. ===========#
  # #=====================================================================#
  # 
  # #=====Use Gurobi to solve the weights=====#
  # #===Preparations
  # model = list()
  # 
  # #=====Write the objective
  # #===Quadratic terms
  # Q.obj = t(X.effective) %*% X.effective
  # Q.obj.pd = as.matrix(nearPD(Q.obj)$mat)
  # #-------- Engineering-wise: make sure the Q matrix is p.s.d. ---------#
  # #---- 20231103: However, this prevents Gurobi to solve to optimal ----#
  # #-------------  Don't know why this happens, forfeit this approach ---#
  # # Q.obj.pd = round(Q.obj.pd, digits = 0)
  # # isSymmetric(Q.obj.pd)
  # # is.positive.definite(Q.obj.pd)
  # # sprintf("%.10f",Q.obj)
  # # sprintf("%.10f",Q.obj.pd)
  # 
  # avg.Q = mean(Q.obj.pd)
  # Q.obj.SizeReduced = Q.obj.pd / avg.Q
  # 
  # #=====Write the constraints
  # #===Linear constraints
  # A.matrix = matrix(1, nrow = 1, ncol = ncol(X.effective))
  # 
  # #=====Pass to Gurobi solver
  # #Objective
  # model$Q       <- Q.obj.SizeReduced
  # #Variables
  # model$vtype   <- 'C'
  # #Linear constraints
  # model$A       <- A.matrix
  # model$rhs     <- c(1)
  # model$sense   <- c('=')
  # 
  # #=====Solve the model
  # params = list(Threads = 1, TimeLimit=360, OutputFlag=0)
  # result = gurobi(model, params = params)
  # 
  # if(sum(result$x) > 1.0001)
  # {
  #   print("simplex constraint is violated")
  #   print("simplex constraint is violated")
  #   print("simplex constraint is violated")
  # }
  # 
  # returned = list(Weights = result$x, Objval = (t(result$x) %*% Q.obj %*% result$x), solver.output = result)
  
  # #============================================================================#
  # #=====Below we use the lsei package to solve constrained least squares. =====#
  # #=====https://search.r-project.org/CRAN/refmans/limSolve/html/lsei.html =====#
  # #=====The input formulation is as follows: ==================================#
  # #============================= min  ||Ax-b||^2 ==============================#
  # #============================= s.t. Ex  = f    ==============================#
  # #=============================      Gx >= h    ==============================#
  # #============================================================================#
  # 
  # AA = X.effective
  # BB = rep(0, nrow(AA))
  # EE = rep(1, ncol(AA))
  # FF = 1
  # GG = diag(1, nrow = ncol(AA))
  # HH = rep(0, ncol(AA))
  # result = lsei(E = EE, F = FF, A = AA, B = BB, G = GG, H = HH)
  # 
  # returned = list(Weights = result$X, Objval = result$solutionNorm, solver.output = result)
  
  #================================================================================#
  #=====Below we use the quadprog package to solve constrained least squares. =====#
  #=====https://cran.r-project.org/web/packages/quadprog/quadprog.pdf =============#
  #=====The input formulation is as follows: ======================================#
  #===== min  (-d'x + 1/2 x' D x)            ======================================#
  #===== s.t.                 A'x >= b       ======================================#
  #===== first 'meq' rows of A are equality constraints ===========================#
  #================================================================================#
  
  Dmat = 2 * (t(X.effective) %*% X.effective)
  pd_Dmat = as.matrix(nearPD(Dmat)$mat)
  avg.Dmat = mean(pd_Dmat)
  Dmat.SizeReduced = pd_Dmat / avg.Dmat
  
  dvec = rep(0, ncol(X.effective))
  A.eq = matrix(1, nrow = 1, ncol = ncol(X.effective))
  b.eq = c(1)
  A.nonneg = matrix(NA, nrow = ncol(X.effective), ncol = ncol(X.effective))
  b.nonneg = c()
  for(A.nonneg.temp in 1:ncol(X.effective))
  {
    A.nonneg[A.nonneg.temp, A.nonneg.temp] = 1
    b.nonneg = c(b.nonneg, 0)
  }
  A.nonneg[is.na(A.nonneg)] = 0
  
  result = solve.QP(Dmat = Dmat.SizeReduced,
                    dvec = dvec,
                    Amat = t(rbind(A.eq, A.nonneg)),
                    bvec = c(b.eq, b.nonneg),
                    meq = nrow(A.eq)) # this argument says the first "meq" rows of Amat are equalities
  
  if(sum(result$solution) > 1.0001)
  {
    print("simplex constraint is violated")
    print("simplex constraint is violated")
    print("simplex constraint is violated")
  }
  
  returned = list(Weights = result$solution, Objval = (0.5 * t(result$solution) %*% Dmat %*% result$solution), solver.output = result)
  
  return(returned)
}

#Formulation Constrained in the paper
Synthetic_Experiment_Cardinality_Constraint <- function(T.prime_, r.ob.covariates.dim_, N.Regions_, Y.N.matrix_, Z.ob.covariates.matrix_, f.vector_, K.cardinality_ = -1)
{
  if(K.cardinality_ == -1)
  {
    K.cardinality = floor(N.Regions_/2)
  }
  if(K.cardinality_ >= 1)
  {
    K.cardinality = K.cardinality_
  }
  
  f.vec = f.vector_
  
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
  
  center.vector = X.matrix %*% f.vec
  
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
    
    Synth.treated.weights = Synth.treated$Weights
    Synth.control.weights = Synth.control$Weights
    
    square.loss.treated = (X.tentative.treated %*% Synth.treated.weights - center.vector)^2
    square.loss.control = (X.tentative.control %*% Synth.control.weights - center.vector)^2
    
    loss = c(loss, (mean(square.loss.treated) + mean(square.loss.control)))
    # loss = c(loss, (Synth.treated$Objval + Synth.control$Objval))
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
  
  # #-- debugging use --#
  # plot(1:T.total, found.treated.weights %*% Y.N.matrix, type = 'l')
  # lines(1:T.total, found.control.weights %*% Y.N.matrix, type = 'l', col = 1, lty = 2)
  # lines(1:T.total, f.vec %*% Y.N.matrix, type = 'l', col = 2, lty = 3)
  
  returned = list(Treatment.weights = round(found.treated.weights, digits = 4), 
                  Control.weights = round(found.control.weights, digits = 4))
  return(returned)
}

permutation.test <- function(pre.intervention.residuals, post.intervention.residuals, permutation.SAMPLES_ = 100000, seed_ = 123456, type_ = 0) # default: 0 - all permutations; 1 - sliding window
{
  set.seed(seed_)
  
  permutation.test.vec = c(abs(pre.intervention.residuals), abs(post.intervention.residuals)) # two sided
  test.statistic = sum(tail(permutation.test.vec, length(post.intervention.residuals)))
  
  if(type_ == 0)
  {
    if(choose(length(permutation.test.vec), length(post.intervention.residuals)) <= permutation.SAMPLES_) # find all permutations
    {
      combinations.matrix = combn(1:length(permutation.test.vec), length(post.intervention.residuals)) # each column is a combination
      permutation.statistics = c()
      for(permu.temp in 1:ncol(combinations.matrix))
      {
        permutation.statistics = c(permutation.statistics, sum(permutation.test.vec[combinations.matrix[,permu.temp]]))
      }
      returned = sum(permutation.statistics >= test.statistic) / ncol(combinations.matrix)
    }
    else #randomly sample permutations
    {
      permutation.statistics = c()
      for(permu.temp in 1:permutation.SAMPLES_)
      {
        samples = sample(1:length(permutation.test.vec), length(post.intervention.residuals), replace = FALSE)
        permutation.statistics = c(permutation.statistics, sum(permutation.test.vec[samples]))
      }
      returned = sum(permutation.statistics >= test.statistic) / permutation.SAMPLES_
    }
  }
  
  if(type_ == 1)
  {
    num.permus = length(permutation.test.vec)
    num.post.intervention = length(post.intervention.residuals)
    permutation.statistics = c()
    for(temp.start.permus in 1:num.permus)
    {
      temp.indices = seq(temp.start.permus, (temp.start.permus+num.post.intervention-1), 1)
      temp.indices.truncated = temp.indices - (temp.indices > num.permus) * num.permus
      permutation.statistics = c(permutation.statistics, sum(permutation.test.vec[temp.indices.truncated]))
    }
    returned = sum(permutation.statistics >= test.statistic) / num.permus
  }
  
  return(returned)
}

quantile_blank <- function(blank.period.residuals_, phi_ = 0.95)
{
  abs.values = abs(blank.period.residuals_)
  length.values = length(blank.period.residuals_)
  target.value.number = ceiling(length.values * phi_)
  returned = sort(abs.values)[target.value.number]
  return(returned)
}

#========================================================#
#=====Warning: The following step easily takes hours=====#
#========================================================#

# Call Synthetic_Experiment function
start_time <- Sys.time()
my.results.Cardinality.changing.K = list()
changing.K = c(1:5)
for(K.temp in changing.K)
{
  my.results.Cardinality.changing.K[[K.temp]] = Synthetic_Experiment_Cardinality_Constraint(T.prime, r.ob.covariates.dim, N.Regions, Y.N.matrix, Z.ob.covariates.matrix, f.vector, K.cardinality_ = K.temp)
}
end_time <- Sys.time()
end_time - start_time


#---------------------------------------------------#
#--- We assume that interventions have no effect ---#
#---------------------------------------------------#
Y.I.matrix = Y.N.matrix

population.N.before = t(f.vector) %*% Y.N.matrix[,1:T.naught]

population.N.after = t(f.vector) %*% Y.N.matrix[,(T.naught+1):T.total]
population.I.after = t(f.vector) %*% Y.I.matrix[,(T.naught+1):T.total]

true.ATE = population.I.after - population.N.after


#===================================#
#===== Constrained formulation =====#
#===================================#

T.I.after.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.N.after.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.N.after.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
C.I.after.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
estimated.ATE.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.total - T.naught))
T.weights.Cardinality.matrix = matrix(NA, nrow = 0, ncol = N.Regions)
C.weights.Cardinality.matrix = matrix(NA, nrow = 0, ncol = N.Regions)

blank.period.residuals.Cardinality.matrix = matrix(NA, nrow = 0, ncol = (T.naught - T.prime))
T.fitted.N.before.Cardinality.matrix = matrix(NA, nrow = 0, ncol = T.naught)
C.fitted.N.before.Cardinality.matrix = matrix(NA, nrow = 0, ncol = T.naught)
permutation.periods.residuals.Cardinality.matrix = matrix(NA, nrow = 0, ncol = T.total-T.prime)
residuals.Cardinality.matrix = matrix(NA, nrow = 0, ncol = T.total)
lower.CI.Cardinality.matrix = matrix(NA, nrow = 0, ncol = T.total-T.naught)
upper.CI.Cardinality.matrix = matrix(NA, nrow = 0, ncol = T.total-T.naught)

estimated.p.value.Cardinality.vector = c()
q.hat.Cardinality.vector = c()

for(K.temp in changing.K)
{
  tentative.T.weights = round(my.results.Cardinality.changing.K[[K.temp]]$Treatment.weights, digits = Weight.digits)
  tentative.C.weights = round(my.results.Cardinality.changing.K[[K.temp]]$Control.weights, digits = Weight.digits)
  
  use.T.weights = tentative.T.weights
  use.C.weights = tentative.C.weights
  
  T.weights.Cardinality.matrix = rbind(T.weights.Cardinality.matrix, use.T.weights)
  C.weights.Cardinality.matrix = rbind(C.weights.Cardinality.matrix, use.C.weights)
  
  T.I.after.Cardinality.matrix = rbind(T.I.after.Cardinality.matrix, t(use.T.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  C.N.after.Cardinality.matrix = rbind(C.N.after.Cardinality.matrix, t(use.C.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  T.N.after.Cardinality.matrix = rbind(T.N.after.Cardinality.matrix, t(use.T.weights) %*% Y.N.matrix[,(T.naught+1):T.total])
  C.I.after.Cardinality.matrix = rbind(C.I.after.Cardinality.matrix, t(use.C.weights) %*% Y.I.matrix[,(T.naught+1):T.total])
  
  estimated.ATE.Cardinality.matrix = rbind(estimated.ATE.Cardinality.matrix, T.I.after.Cardinality.matrix[K.temp,] - C.N.after.Cardinality.matrix[K.temp,])
  
  T.fitted.N.before.Cardinality.matrix = rbind(T.fitted.N.before.Cardinality.matrix, t(use.T.weights) %*% Y.N.matrix[,1:T.naught])
  C.fitted.N.before.Cardinality.matrix = rbind(C.fitted.N.before.Cardinality.matrix, t(use.C.weights) %*% Y.N.matrix[,1:T.naught])
  blank.period.residuals.Cardinality.matrix = rbind(blank.period.residuals.Cardinality.matrix, T.fitted.N.before.Cardinality.matrix[K.temp, (T.prime+1):T.naught] - C.fitted.N.before.Cardinality.matrix[K.temp, (T.prime+1):T.naught])
  residuals.Cardinality.matrix = rbind(residuals.Cardinality.matrix, c(T.fitted.N.before.Cardinality.matrix[K.temp,] - C.fitted.N.before.Cardinality.matrix[K.temp,], estimated.ATE.Cardinality.matrix[K.temp,]))
  
  #================================#
  #========== Estimation ==========#
  #================================#
  
  #===Visualization===#
  plot.l.lim = min(population.N.before, population.I.after, population.N.after, C.fitted.N.before.Cardinality.matrix[K.temp,], C.N.after.Cardinality.matrix[K.temp,], T.fitted.N.before.Cardinality.matrix[K.temp,], T.I.after.Cardinality.matrix[K.temp,])*0.99
  plot.u.lim = max(population.N.before, population.I.after, population.N.after, C.fitted.N.before.Cardinality.matrix[K.temp,], C.N.after.Cardinality.matrix[K.temp,], T.fitted.N.before.Cardinality.matrix[K.temp,], T.I.after.Cardinality.matrix[K.temp,])*1.01
  
  my.png.name = as.character(paste0("Output/Fitted_", config.name, "Constrained_K=", K.temp, "_", format(Sys.time(), "%Y%m%d"), ".png"))
  png(filename = my.png.name, width = 1536, height = 768)
  par(mar=c(9.5,4,0.5,0.5))
  plot(1:T.total, c(C.fitted.N.before.Cardinality.matrix[K.temp,], C.N.after.Cardinality.matrix[K.temp,]), ylim = c(plot.l.lim, plot.u.lim), xlab = "", xaxt = 'n', ylab = "", type = 'l', lty = 2, lwd = 3, main = "", cex.lab = 2, cex.axis = 2)
  axis(1, at = (1:T.total)[seq(1, T.total, by = 4)], labels = unique.date.names.merged[seq(1, T.total, by = 4)], las = 2, cex.axis = 1.5)
  # #---Plot the control units to verify no unit perfectly matches the treated---#
  # Y.N.control.group.matrix = Y.N.matrix[C.weights!= 0,1:T.naught]
  # for(control.temp in 1:nrow(Y.N.control.group.matrix))
  # {
  #   lines(1:T.naught, Y.N.control.group.matrix[control.temp,], type = 'l', col = rgb(0.4, 0.7, 1, alpha = 0.6), lty = 3, lwd = 3)
  # }
  #---Plot all units to verify no unit perfectly matches the treated---#
  for(control.temp in 1:nrow(Y.N.matrix))
  {
    lines(1:T.naught, Y.N.matrix[control.temp,1:T.naught], type = 'l', col = rgb(0.4, 0.7, 1, alpha = 0.6), lty = 3, lwd = 3)
  }
  lines(1:T.total, c(C.fitted.N.before.Cardinality.matrix[K.temp,], C.N.after.Cardinality.matrix[K.temp,]), type = 'l', lty = 2, lwd = 3)
  lines(1:T.total, c(T.fitted.N.before.Cardinality.matrix[K.temp,], T.I.after.Cardinality.matrix[K.temp,]), type = 'l', lty = 1, lwd = 3)
  # lines(1:T.total, colMeans(Sales.panel), type = 'l', col = 2, lty = 1, lwd = 3)
  #---Finish plotting---#
  abline(v=T.prime, col = 1, lty = 3, lwd = 4)
  abline(v=T.naught, col = 1, lty = 3, lwd = 4)
  text((T.prime/2), plot.u.lim, "fitting periods", col = rgb(0.1, 0.1, 0.1), adj = c(0, 0), cex = 2)
  text(((T.prime+T.naught)/2), plot.u.lim, "blank periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.5, 0), cex = 2)
  text(((T.naught+T.total)/2), plot.u.lim, "expt periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.25, 0), cex = 2)
  legend("topleft", legend = c("synthetic treated", "synthetic control", "individual units"), lty = c(1, 2, 3), col = c(1,1,rgb(0.4, 0.7, 1, alpha = 0.6)), lwd = 4, cex = 2)
  dev.off()
  
  #===============================#
  #========== Inference ==========#
  #===============================#
  
  #------------------------#
  #--- Permutation test ---#
  #------------------------#
  
  estimated.p.value.Cardinality.vector = c(estimated.p.value.Cardinality.vector, permutation.test(blank.period.residuals.Cardinality.matrix[K.temp,], estimated.ATE.Cardinality.matrix[K.temp,], type_ = 0))
  
  #===Visualization===#
  plot.l.lim = min(T.fitted.N.before.Cardinality.matrix[K.temp,] - C.fitted.N.before.Cardinality.matrix[K.temp,], T.I.after.Cardinality.matrix[K.temp,] - C.N.after.Cardinality.matrix[K.temp,])
  plot.u.lim = max(T.fitted.N.before.Cardinality.matrix[K.temp,] - C.fitted.N.before.Cardinality.matrix[K.temp,], T.I.after.Cardinality.matrix[K.temp,] - C.N.after.Cardinality.matrix[K.temp,])
  
  my.png.name = as.character(paste0("Output/Residuals_", config.name, "Constrained_K=", K.temp, "_", format(Sys.time(), "%Y%m%d"), ".png"))
  png(filename = my.png.name, width = 1536, height = 768)
  par(mar=c(9.5,4,0.5,0.5))
  plot(1:T.total, c(T.fitted.N.before.Cardinality.matrix[K.temp,] - C.fitted.N.before.Cardinality.matrix[K.temp,], T.I.after.Cardinality.matrix[K.temp,] - C.N.after.Cardinality.matrix[K.temp,]), ylim = c(plot.l.lim, plot.u.lim), xlab = "", xaxt = 'n', ylab = "", type = 'l', lwd = 3, main = "", cex = 3, cex.lab = 2, cex.axis = 2)
  axis(1, at = (1:T.total)[seq(1, T.total, by = 4)], labels = unique.date.names.merged[seq(1, T.total, by = 4)], las = 2, cex.axis = 1.5)
  abline(h=0, col = 8, lty = 3, lwd = 3)
  abline(v=T.prime, col = 8, lty = 3, lwd = 4)
  abline(v=T.naught, col = 8, lty = 3, lwd = 4)
  text((T.prime/2), plot.u.lim, "fitting periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.5, 0), cex = 2)
  text(((T.prime+T.naught)/2), plot.u.lim, "blank periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.5, 0), cex = 2)
  text(((T.naught+T.total)/2), plot.u.lim, "expt periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.25, 0), cex = 2)
  text(((T.naught+T.total)/2), plot.u.lim, paste0("p value = ", round(estimated.p.value.Cardinality.vector[K.temp], digits = 3)), col = rgb(0.1, 0.1, 0.1), adj = c(0.25, 2), cex = 2)
  legend("topleft", legend = c("treatment effect estimate"), col = c(1), lty = c(1), lwd = 4, cex = 2)
  dev.off()
  
  #--------------------------#
  #--- Conformal Interval ---#
  #--------------------------#
  
  q.hat.Cardinality.vector = c(q.hat.Cardinality.vector, quantile_blank(blank.period.residuals.Cardinality.matrix[K.temp,], 1-alpha.CI))
  lower.CI.Cardinality.matrix = rbind(lower.CI.Cardinality.matrix, T.I.after.Cardinality.matrix[K.temp,] - C.N.after.Cardinality.matrix[K.temp,] - q.hat.Cardinality.vector[K.temp])
  upper.CI.Cardinality.matrix = rbind(upper.CI.Cardinality.matrix, T.I.after.Cardinality.matrix[K.temp,] - C.N.after.Cardinality.matrix[K.temp,] + q.hat.Cardinality.vector[K.temp])
  experimental.periods.for.drawing.only = c((T.naught+1):T.total)
  
  #===Visualization===#
  plot.l.lim = min(T.fitted.N.before.Cardinality.matrix[K.temp,] - C.fitted.N.before.Cardinality.matrix[K.temp,], T.I.after.Cardinality.matrix[K.temp,] - C.N.after.Cardinality.matrix[K.temp,], lower.CI.Cardinality.matrix[K.temp,])
  plot.u.lim = max(T.fitted.N.before.Cardinality.matrix[K.temp,] - C.fitted.N.before.Cardinality.matrix[K.temp,], T.I.after.Cardinality.matrix[K.temp,] - C.N.after.Cardinality.matrix[K.temp,], upper.CI.Cardinality.matrix[K.temp,])
  
  my.png.name = as.character(paste0("Output/ResidualsCI_", config.name, "Constrained_K=", K.temp, "_", format(Sys.time(), "%Y%m%d"), ".png"))
  png(filename = my.png.name, width = 1536, height = 768)
  par(mar=c(9.5,4,0.5,0.5))
  plot(1:T.total, c(T.fitted.N.before.Cardinality.matrix[K.temp,] - C.fitted.N.before.Cardinality.matrix[K.temp,], T.I.after.Cardinality.matrix[K.temp,] - C.N.after.Cardinality.matrix[K.temp,]), ylim = c(plot.l.lim, plot.u.lim), xlab = "", xaxt = 'n', ylab = "", type = 'l', lwd = 3, main = "", cex = 3, cex.lab = 2, cex.axis = 2)
  axis(1, at = (1:T.total)[seq(1, T.total, by = 4)], labels = unique.date.names.merged[seq(1, T.total, by = 4)], las = 2, cex.axis = 1.5)
  lines(experimental.periods.for.drawing.only, lower.CI.Cardinality.matrix[K.temp,], type = 'l', lwd = 3, main = "", cex = 3, cex.lab = 2, cex.axis = 2, col = rgb(0.75, 0.72, 0.70, alpha = 0.6))
  lines(experimental.periods.for.drawing.only, upper.CI.Cardinality.matrix[K.temp,], type = 'l', lwd = 3, main = "", cex = 3, cex.lab = 2, cex.axis = 2, col = rgb(0.75, 0.72, 0.70, alpha = 0.6))
  polygon(c(experimental.periods.for.drawing.only, rev(experimental.periods.for.drawing.only)), c(lower.CI.Cardinality.matrix[K.temp,], rev(upper.CI.Cardinality.matrix[K.temp,])), col = rgb(0.75, 0.72, 0.70, alpha = 0.4), border = NA)
  abline(h=0, col = 8, lty = 3, lwd = 3)
  abline(v=T.prime, col = 8, lty = 3, lwd = 4)
  abline(v=T.naught, col = 8, lty = 3, lwd = 4)
  text((T.prime/2), plot.u.lim, "fitting periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.5, 0), cex = 2)
  text(((T.prime+T.naught)/2), plot.u.lim, "blank periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.5, 0), cex = 2)
  text(((T.naught+T.total)/2), plot.u.lim, "expt periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.25, 0), cex = 2)
  text(((T.naught+T.total)/2), plot.u.lim, paste0("p value = ", round(estimated.p.value.Cardinality.vector[K.temp], digits = 3)), col = rgb(0.1, 0.1, 0.1), adj = c(0.25, 2), cex = 2)
  legend("topleft", legend = c("treatment effect estimate"), col = c(1), lty = c(1), lwd = 4, cex = 2)
  dev.off()
}

# #=======================================#
# #===== Sanity check: debugging use =====#
# #=======================================#
# 
# MSE.vec = c()
# for(K.temp in changing.K)
# {
#   MSE.vec = c(MSE.vec, sqrt(mean(c(T.fitted.N.before.Cardinality.matrix[K.temp,1:T.prime] - C.fitted.N.before.Cardinality.matrix[K.temp,1:T.prime])^2)))
# }
# MSE.vec
# 
# MSE.T.vec = c()
# MSE.C.vec = c()
# for(K.temp in changing.K)
# {
#   MSE.T.vec = c(MSE.T.vec, sqrt(mean(c(T.fitted.N.before.Cardinality.matrix[K.temp,1:T.prime] - population.N.before[1:T.prime])^2)))
#   MSE.C.vec = c(MSE.C.vec, sqrt(mean(c(C.fitted.N.before.Cardinality.matrix[K.temp,1:T.prime] - population.N.before[1:T.prime])^2)))
# }
# MSE.T.vec
# MSE.C.vec
# MSE.T.vec + MSE.C.vec

#=====================================================#
#===== Synthetic Control results: print and save =====#
#=====================================================#

ATE.report = rbind(true.ATE,
                   estimated.ATE.Cardinality.matrix)
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
ATE.report.print = cbind(ATE.report, Loss.MAE.vec, Loss.MSE.vec)
my.file.name.copy = as.character(paste0(Sys.Date(), "SC_cardinality_matrix.txt"))
write.table(round(ATE.report.print, digits = 3), file = my.file.name.copy, append = FALSE, sep = "\t", col.names=FALSE)
