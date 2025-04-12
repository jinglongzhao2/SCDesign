
#=====================================================#
#===== One basic run of synthetic control design =====#
#===== Change "noise.variance = 1" in Line 131   =====#
#=====================================================#


library(Matrix)
library(slam)
library(dplyr)
library(gurobi)
library(Matrix)

Weight.digits = 2

set.seed(123456)

#=====Define variables=====#
#===Basics
#Basics --- Scalars
N.Regions = 15
T.naught = 25
T.prime = 20
T.total = 30
r.ob.covariates.dim = 7
F.unob.covariates.dim = 11
#Population weights --- Vectors
f.vector = c()

#===Model Primitives
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

#=====Generate Random Values=====#
#===Basics
#Population weights --- Vectors
# f.temp = runif(N.Regions)
# f.vector = f.temp / sum(f.temp)
f.vector = rep(1/N.Regions, N.Regions)
#===Model Primitives
#Constants --- Vectors
#Range --- [0,10]
range.intercept.max = 20 
delta.N.vector = sort(c(range.intercept.max * runif(T.naught), range.intercept.max * runif(T.total - T.naught)))
upsilon.I.vector = c(rep(NA, T.naught), sort(range.intercept.max * runif(T.total - T.naught)))
#Covariates --- Matrices
#Range --- [0,1]
range.covariates.max = 1 
for(j in 1:N.Regions)
{
  z.temp = range.covariates.max * runif(r.ob.covariates.dim)
  Z.ob.covariates.matrix[,j] = z.temp
}
colnames(Z.ob.covariates.matrix) = c(1:N.Regions)
for(j in 1:N.Regions)
{
  mu.temp = range.covariates.max * runif(F.unob.covariates.dim)
  mu.unob.covariates.matrix[,j] = mu.temp
}
colnames(mu.unob.covariates.matrix) = c(1:N.Regions)
#Coefficients --- Matrices
#Range --- [0,10]
range.coefficients.max = 10 
for(t in 1:T.naught)
{
  theta.temp = range.coefficients.max * runif(r.ob.covariates.dim)
  theta.ob.N.matrix[,t] = theta.temp
}
for(t in (T.naught+1):T.total)
{
  theta.temp = range.coefficients.max * runif(r.ob.covariates.dim)
  theta.ob.N.matrix[,t] = theta.temp
}
colnames(theta.ob.N.matrix) = c(1:T.total)

for(t in 1:T.naught)
{
  gamma.temp = rep(NA, r.ob.covariates.dim)
  gamma.ob.I.matrix[,t] = gamma.temp
}
for(t in (T.naught+1):T.total)
{
  gamma.temp = range.coefficients.max * runif(r.ob.covariates.dim)
  gamma.ob.I.matrix[,t] = gamma.temp
}
colnames(gamma.ob.I.matrix) = c(1:T.total)

for(t in 1:T.naught)
{
  lambda.temp = range.coefficients.max * runif(F.unob.covariates.dim)
  lambda.unob.N.matrix[,t] = lambda.temp
}
for(t in (T.naught+1):T.total)
{
  lambda.temp = range.coefficients.max * runif(F.unob.covariates.dim)
  lambda.unob.N.matrix[,t] = lambda.temp
}
colnames(lambda.unob.N.matrix) = c(1:T.total)

for(t in 1:T.naught)
{
  eta.temp = rep(NA, F.unob.covariates.dim)
  eta.unob.I.matrix[,t] = eta.temp
}
for(t in (T.naught+1):T.total)
{
  eta.temp = range.coefficients.max * runif(F.unob.covariates.dim)
  eta.unob.I.matrix[,t] = eta.temp
}
colnames(eta.unob.I.matrix) = c(1:T.total)

#Noises --- Matrices
#Normal Distribution ~ N(0,1)
noise.variance = 1 #ChangeThis
for(t in 1:T.naught)
{
  epsilon.temp = rnorm(N.Regions, 0, noise.variance)
  epsilon.N.matrix[,t] = epsilon.temp
}
for(t in (T.naught+1):T.total)
{
  epsilon.temp = rnorm(N.Regions, 0, noise.variance)
  epsilon.N.matrix[,t] = epsilon.temp
}
colnames(epsilon.N.matrix) = c(1:T.total)

for(t in 1:T.naught)
{
  xi.temp = rep(NA, N.Regions)
  xi.I.matrix[,t] = xi.temp
}
for(t in (T.naught+1):T.total)
{
  xi.temp = rnorm(N.Regions, 0, noise.variance)
  xi.I.matrix[,t] = xi.temp
}
colnames(xi.I.matrix) = c(1:T.total)

#=====Generate Potential Outcomes=====#
#Potential Outcomes in hindsight --- Matrices
#A Factor Model
for(j in 1:N.Regions)
{
  for(t in 1:T.total)
  {
    Y.N.matrix[j,t] = delta.N.vector[t] + theta.ob.N.matrix[,t] %*% Z.ob.covariates.matrix[,j] + lambda.unob.N.matrix[,t] %*% mu.unob.covariates.matrix[,j] + epsilon.N.matrix[j,t]
  }
}

for(j in 1:N.Regions)
{
  for(t in 1:T.total)
  {
    Y.I.matrix[j,t] = upsilon.I.vector[t] + gamma.ob.I.matrix[,t] %*% Z.ob.covariates.matrix[,j] + eta.unob.I.matrix[,t] %*% mu.unob.covariates.matrix[,j] + xi.I.matrix[j,t]
  }
}

plot(1:T.total, Y.N.matrix[1,])
lines(1:T.total, Y.N.matrix[2,], col = 2, type = 'p')



#=====================================#
#===== Now run synthetic control =====#
#=====================================#

#Formulation "Unconstrained" version in the paper
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





#=================================================#
#===== Call Synthetic_Experiment function ========#
#===== Warning: This step easily takes hours =====#
#=================================================#

start_time <- Sys.time()

my.result = Synthetic_Experiment(T.prime, r.ob.covariates.dim, N.Regions, Y.N.matrix, Z.ob.covariates.matrix, f.vector)

end_time <- Sys.time()
end_time - start_time


#=====Fetch the outputs=====#
#===Rounding is somewhat needed. ===========================#
#===Since the QCQP is solved numerically, the unselected ===#
#===weights are not exactly zero ~10^(-8) error ============#
T.weights = round(my.result$Treatment.weights, digits = Weight.digits)
C.weights = round(my.result$Control.weights, digits = Weight.digits)


population.N.before = t(f.vector) %*% Y.N.matrix[,1:T.naught]
T.fitted.N.before = t(T.weights) %*% Y.N.matrix[,1:T.naught]
C.fitted.N.before = t(C.weights) %*% Y.N.matrix[,1:T.naught]

population.N.after = t(f.vector) %*% Y.N.matrix[,(T.naught+1):T.total]
population.I.after = t(f.vector) %*% Y.I.matrix[,(T.naught+1):T.total]
T.I.after = t(T.weights) %*% Y.I.matrix[,(T.naught+1):T.total]
C.N.after = t(C.weights) %*% Y.N.matrix[,(T.naught+1):T.total]


true.ATE = population.I.after - population.N.after
estimated.ATE = T.I.after - C.N.after
ATET.T.weighted = t(T.weights) %*% Y.I.matrix[,(T.naught+1):T.total] - t(T.weights) %*% Y.N.matrix[,(T.naught+1):T.total]
ATEC.C.weighted = t(C.weights) %*% Y.I.matrix[,(T.naught+1):T.total] - t(C.weights) %*% Y.N.matrix[,(T.naught+1):T.total]

blank.period.residuals = T.fitted.N.before[(T.prime+1):T.naught] - C.fitted.N.before[(T.prime+1):T.naught]

permutation.test.vec = c(abs(blank.period.residuals), abs(estimated.ATE))
test.statistic = sum(tail(permutation.test.vec, T.total - T.naught))
combinations.matrix = combn(1:(T.total-T.prime), T.total - T.naught) #each column is a combination
permutation.statistics = c()
for(permu.temp in 1:ncol(combinations.matrix))
{
  permutation.statistics = c(permutation.statistics, sum(permutation.test.vec[combinations.matrix[,permu.temp]]))
}
estimated.p.value = sum(permutation.statistics >= test.statistic) / ncol(combinations.matrix)

# #===Visualization: all===#
# plot.l.lim = min(population.N.before, population.I.after, population.N.after)*0.95
# plot.u.lim = max(population.N.before, population.I.after, population.N.after)*1.05
# 
# my.png.name = as.character(paste0("Noise Variance =", noise.variance,".png"))
# png(filename = my.png.name, width = 1536, height = 768)
# par(mar=c(4.5,5,0,0.5))
# plot(1:T.total, c(population.N.before, rep(NA, (T.total - T.naught))), ylim = c(plot.l.lim, plot.u.lim), xlab = "Horizon", ylab = "Value of Interests", type = 'l', lwd = 3, main = "", cex.lab = 2, cex.axis = 2) #main = paste0("N(0,", noise.variance,") Noise")
# lines(1:T.total, c(rep(NA, (T.naught)), population.I.after), type = 'l', col = rgb(0.3, 0.1, 0.1), lwd = 3, lty = 2)
# lines(1:T.total, c(rep(NA, (T.naught-1)), population.N.before[T.naught], population.N.after), type = 'l', col = rgb(0.1, 0.3, 0.1), lwd = 3, lty = 4)
# lines(1:T.total, c(T.fitted.N.before, rep(NA, (T.total - T.naught))), type = 'l', col = rgb(0.8, 0, 0.4), lwd = 3, lty = 5)
# lines(1:T.total, c(C.fitted.N.before, rep(NA, (T.total - T.naught))), type = 'l', col = rgb(0.4, 0.8, 0.4), lwd = 3, lty = 6)
# lines(1:T.total, c(rep(NA, (T.naught)), T.I.after), type = 'l', col = rgb(1, 0, 0.5), lwd = 3, lty = 2)
# lines(1:T.total, c(rep(NA, (T.naught-1)), C.fitted.N.before[T.naught], C.N.after), type = 'l', col = rgb(0.5, 1, 0.5), lwd = 3, lty = 4)
# abline(v=20, col = 8, lty = 3, lwd = 2)
# abline(v=25, col = 8, lty = 3, lwd = 2)
# text((T.prime/2), plot.u.lim, "Fitting periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.5, 0), cex = 2)
# text(((T.prime+T.naught)/2), plot.u.lim, "Blank periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.5, 0), cex = 2)
# text(((T.naught+T.total)/2), plot.u.lim, "Experimental periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.4, 0), cex = 2)
# # text((T.naught/2), plot.u.lim, "Pre-intervention", col = rgb(0.1, 0.1, 0.1), adj = c(0.5, 0), cex = 1)
# # text(((T.naught+T.total)/2), plot.u.lim, "Post-intervention", col = rgb(0.1, 0.1, 0.1), adj = c(0.4, 0), cex = 1)
# legend("bottomright", legend = c("(Pre-Intervention) Observed All-in-Control Outcomes", "(Pre-Intervention) Counterfactual Treatment Outcomes", "(Pre-Intervention) Counterfactual Control Outcomes", "(Post-Intervention) Counterfactual All-in-Treatment Outcomes", "(Post-Intervention) Counterfactual All-in-Control Outcomes", "(Post-Intervention) Observed Treatment Outcomes", "(Post-Intervention) Observed Control Outcomes"),
#        col = c(rgb(0, 0, 0), rgb(0.8, 0, 0.4), rgb(0.4, 0.8, 0.4), rgb(0.3, 0.1, 0.1), rgb(0.1, 0.3, 0.1), rgb(1, 0, 0.5), rgb(0.5, 1, 0.5)), lty = c(1, 5, 6, 2, 4, 2, 4), lwd = 4, cex = 1.6)
# text(T.total, plot.u.lim, paste0("p value =", round(estimated.p.value, digits = 3)), col = rgb(0, 0, 0), adj = c(1, 5), cex = 2)
# dev.off()

#===Visualization:lean===#
#===Estimation===#
plot.l.lim = (min(population.N.before, population.I.after, population.N.after)-3)*0.85
plot.u.lim = (max(population.N.before, population.I.after, population.N.after)+3)*1.15

my.png.name = as.character(paste0("ObservedData_NoiseVariance=", noise.variance,".png"))
png(filename = my.png.name, width = 1536, height = 768)
par(mar=c(4.5,5,0.5,0.5))
plot(1:T.total, c(C.fitted.N.before, C.N.after), ylim = c(plot.l.lim, plot.u.lim), xlab = "time", ylab = "", type = 'l', lty = 2, lwd = 3, main = "", cex.lab = 2, cex.axis = 2) #main = paste0("N(0,", noise.variance,") Noise")
lines(1:T.total, c(T.fitted.N.before, T.I.after), type = 'l', lty = 1, lwd = 3)
# #---Plot the control units to verify no unit perfectly match the treated---#
# Y.N.control.group.matrix = Y.N.matrix[C.weights!= 0,1:T.naught]
# for(control.temp in 1:nrow(Y.N.control.group.matrix))
# {
#   lines(1:T.naught, Y.N.control.group.matrix[control.temp,], type = 'l', col = rgb(0.4, 0.7, 1, alpha = 0.6), lty = 3, lwd = 3)
# }
#---Plot the control units to verify no unit perfectly match the treated---#
for(control.temp in 1:nrow(Y.N.matrix))
{
  lines(1:T.naught, Y.N.matrix[control.temp,1:T.naught], type = 'l', col = rgb(0.4, 0.7, 1, alpha = 0.6), lty = 3, lwd = 3)
}
#---Finish plotting---#
abline(v=T.prime, col = 8, lty = 3, lwd = 2)
abline(v=T.naught, col = 8, lty = 3, lwd = 2)
text((T.prime/2), plot.u.lim, "fitting periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.5, 0), cex = 2)
text(((T.prime+T.naught)/2), plot.u.lim, "blank periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.5, 0), cex = 2)
text(((T.naught+T.total)/2), plot.u.lim, "experimental periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.4, 0), cex = 2)
legend("bottomleft", legend = c("synthetic treated", "synthetic control", "individual units"), lty = c(1, 2, 3), col = c(1,1,rgb(0.4, 0.7, 1, alpha = 0.6)), lwd = 4, cex = 2)
dev.off()

#===Inference===#
plot.l.lim = min(T.fitted.N.before - C.fitted.N.before, T.I.after - C.N.after)
plot.u.lim = max(T.fitted.N.before - C.fitted.N.before, T.I.after - C.N.after)+3

my.png.name = as.character(paste0("Residuals_NoiseVariance=", noise.variance,".png"))
png(filename = my.png.name, width = 1536, height = 768)
par(mar=c(4.5,5,0.5,0.5))
plot(1:T.total, c(T.fitted.N.before - C.fitted.N.before, T.I.after - C.N.after), ylim = c(plot.l.lim, plot.u.lim), xlab = "time", ylab = "", type = 'l', lwd = 3, main = "", cex = 3, cex.lab = 2, cex.axis = 2)
abline(h=0, col = 8, lty = 3, lwd = 3)
abline(v=T.prime, col = 8, lty = 3, lwd = 2)
abline(v=T.naught, col = 8, lty = 3, lwd = 2)
text((T.prime/2), plot.u.lim, "fitting periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.5, 0), cex = 2)
text(((T.prime+T.naught)/2), plot.u.lim, "blank periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.5, 0), cex = 2)
text(((T.naught+T.total)/2), plot.u.lim, "experimental periods", col = rgb(0.1, 0.1, 0.1), adj = c(0.4, 0), cex = 2)
text(((T.naught+T.total)/2), plot.u.lim, paste0("p value =", round(estimated.p.value, digits = 3)), col = rgb(0.1, 0.1, 0.1), adj = c(0.4, 2), cex = 2)
legend("bottomleft", legend = c("treatment effect estimate"), col = c(1), lty = c(1), lwd = 4, cex = 2)
dev.off()

#===Save to .txt file
write.table(rbind(true.ATE,estimated.ATE), file = paste0("Synthetic_Experiments_N(0,", noise.variance,")_Noise.txt"), append = TRUE, sep = "\t", col.names=FALSE)

