#=========================================#
#===== Read in all the files printed =====#
#=========================================#

#--- Last Change: 20250912 9PM --- combined the codes with those analyzing randomization results
#--- Last Change: 20250823 8AM --- added the number of treated units
#--- Last Change: 20250316 2PM --- changed the MSE into RMSE (take the square root) to make it comparable with MAE


library(tidyverse)

start_time <- Sys.time()

#----------------------------------#
#----- different formulations -----#
#----------------------------------#

list_of_files_ATE = list.files(path = "./output/", recursive = TRUE, pattern = "\\_ATE.txt$", full.names = TRUE)
datalist_ATE = lapply(list_of_files_ATE, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_ATE = lapply(datalist_ATE, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_ATE = Reduce("+", datalist_ATE)/ length(datalist_ATE)

list_of_files_ATE_null_hypothesis = list.files(path = "./output/", recursive = TRUE, pattern = "\\_ATE_null_hypothesis.txt$", full.names = TRUE)
datalist_ATE_null_hypothesis = lapply(list_of_files_ATE_null_hypothesis, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_ATE_null_hypothesis = lapply(datalist_ATE_null_hypothesis, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_ATE_null_hypothesis = Reduce("+", datalist_ATE_null_hypothesis)/ length(datalist_ATE_null_hypothesis)

list_of_files_ATEC = list.files(path = "./output/", recursive = TRUE, pattern = "\\_ATEC.txt$", full.names = TRUE)
datalist_ATEC = lapply(list_of_files_ATEC, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_ATEC = lapply(datalist_ATEC, function(mat) {
  mat[, 13] = sqrt(mat[, 13])  # Apply square root to 13th column
  return(mat)
})
average_outputs_ATEC = Reduce("+", datalist_ATEC)/ length(datalist_ATEC)

list_of_files_ATET = list.files(path = "./output/", recursive = TRUE, pattern = "\\_ATET.txt$", full.names = TRUE)
datalist_ATET = lapply(list_of_files_ATET, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_ATET = lapply(datalist_ATET, function(mat) {
  mat[, 13] = sqrt(mat[, 13])  # Apply square root to 13th column
  return(mat)
})
average_outputs_ATET = Reduce("+", datalist_ATET)/ length(datalist_ATET)

# ATET_example = sapply(datalist_ATET, function(x) x[14, 2])
# ATE_example = sapply(datalist_ATE, function(x) x[1, 2])
# 
# x_limits = range(c(ATET_example, ATE_example))
# 
# hist(ATET_example,
#      col = rgb(0, 0, 1, 0.5), # Blue with 50% transparency
#      main = "Histogram of ATET_example and ATE_example",
#      xlab = "Value",
#      ylab = "Frequency",
#      xlim = x_limits,
#      breaks = 30)
# 
# hist(ATE_example,
#      col = rgb(1, 0, 0, 0.5), # Red with 50% transparency
#      add = TRUE,
#      breaks = 30)

list_of_files_TreatmentWeights = list.files(path = "./output/", recursive = TRUE, pattern = "\\_TreatmentWeights.txt$", full.names = TRUE)
datalist_TreatmentWeights = lapply(list_of_files_TreatmentWeights, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_NumTreatmentWeights = lapply(datalist_TreatmentWeights, function(mat)
{
  # Count the number of non-zero weights of each row
  returned = apply(mat[, -1], 1, function(row)
  {
    sum(row != 0)
  })
  return(returned)
})
average_outputs_NumTreatmentWeights = Reduce("+", datalist_NumTreatmentWeights) / length(datalist_NumTreatmentWeights)

list_of_files_ControlWeights = list.files(path = "./output/", recursive = TRUE, pattern = "\\_ControlWeights.txt$", full.names = TRUE)
datalist_ControlWeights = lapply(list_of_files_ControlWeights, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_NumControlWeights = lapply(datalist_ControlWeights, function(mat)
{
  # Count the number of non-zero weights of each row
  returned = apply(mat[, -1], 1, function(row)
  {
    sum(row != 0)
  })
  return(returned)
})
average_outputs_NumControlWeights = Reduce("+", datalist_NumControlWeights) / length(datalist_NumControlWeights)

#--------------------------------------------#
#----- randomized treatment assignments -----#
#--------------------------------------------#

#--- random assignment + difference-in-means ---#
list_of_files_RADiM = list.files(path = "./output/", recursive = TRUE, pattern = "\\_RADiM.txt$", full.names = TRUE)
datalist_RADiM = lapply(list_of_files_RADiM, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_RADiM = lapply(datalist_RADiM, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_RADiM = Reduce("+", datalist_RADiM)/ length(datalist_RADiM)

#--- stratified assignment + difference-in-means ---#
list_of_files_StADiM = list.files(path = "./output/", recursive = TRUE, pattern = "\\_StADiM.txt$", full.names = TRUE)
datalist_StADiM = lapply(list_of_files_StADiM, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_StADiM = lapply(datalist_StADiM, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_StADiM = Reduce("+", datalist_StADiM)/ length(datalist_StADiM)

#--- random assignment + regression adjustment ---#
list_of_files_RARegAdj = list.files(path = "./output/", recursive = TRUE, pattern = "\\_RARegAdj.txt$", full.names = TRUE)
datalist_RARegAdj = lapply(list_of_files_RARegAdj, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_RARegAdj = lapply(datalist_RARegAdj, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_RARegAdj = Reduce("+", datalist_RARegAdj)/ length(datalist_RARegAdj)

#--- random assignment + 1 nearest neighbor matching ---#
list_of_files_RA1NN = list.files(path = "./output/", recursive = TRUE, pattern = "\\_RA1NN.txt$", full.names = TRUE)
datalist_RA1NN = lapply(list_of_files_RA1NN, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_RA1NN = lapply(datalist_RA1NN, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_RA1NN = Reduce("+", datalist_RA1NN)/ length(datalist_RA1NN)

#--- random assignment + 5 nearest neighbor matching ---#
list_of_files_RA5NN = list.files(path = "./output/", recursive = TRUE, pattern = "\\_RA5NN.txt$", full.names = TRUE)
datalist_RA5NN = lapply(list_of_files_RA5NN, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_RA5NN = lapply(datalist_RA5NN, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_RA5NN = Reduce("+", datalist_RA5NN)/ length(datalist_RA5NN)



#================================================#
#===== check which ones are missing because =====#
#===== some servers are occasionally down   =====#
#================================================#

list_of_files_ATE_NoSpaces = gsub(" ", "", list_of_files_ATE)
list_of_files_ATE_Trimmed = sapply(list_of_files_ATE_NoSpaces, substring, 10)
list_of_files_ATE_seeds = as.numeric(gsub(".*?([0-9]+).*", "\\1", list_of_files_ATE_Trimmed))
list_of_files_ATE_seeds = sort(list_of_files_ATE_seeds)
failed.attempts = setdiff(c(1:1000), list_of_files_ATE_seeds)
print(failed.attempts)

list_of_files_RADiM_NoSpaces = gsub(" ", "", list_of_files_RADiM)
list_of_files_RADiM_Trimmed = sapply(list_of_files_RADiM_NoSpaces, substring, 10)
list_of_files_RADiM_seeds = as.numeric(gsub(".*?([0-9]+).*", "\\1", list_of_files_RADiM_Trimmed))
list_of_files_RADiM_seeds = sort(list_of_files_RADiM_seeds)
failed.attempts = setdiff(c(1:1000), list_of_files_RADiM_seeds)
print(failed.attempts)



#==============================================================================#
#===== Check if there are non-zero weights for both treatment and control =====#
#==============================================================================#

list_of_files_TreatmentWeights = list.files(path = "./output/", recursive = TRUE, pattern = "\\_TreatmentWeights.txt$", full.names = TRUE)
datalist_TreatmentWeights = lapply(list_of_files_TreatmentWeights, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))

list_of_files_ControlWeights = list.files(path = "./output/", recursive = TRUE, pattern = "\\_ControlWeights.txt$", full.names = TRUE)
datalist_ControlWeights = lapply(list_of_files_ControlWeights, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))

for(weights.temp in 1:1000)
{
  TreatmentWeights.temp = datalist_TreatmentWeights[[weights.temp]][,-1]
  ControlWeights.temp = datalist_ControlWeights[[weights.temp]][,-1]

  TreatmentNonZero = (TreatmentWeights.temp > 0)
  ControlNonZero = (ControlWeights.temp > 0)

  sum.NonZero = TreatmentNonZero + ControlNonZero
  if(max(sum.NonZero) > 1)
  {
    print(weights.temp)
  }
}



#=========================================#
#===== now print out all the results =====#
#=========================================#

#----------------------------------#
#----- different formulations -----#
#----------------------------------#

average_outputs_ATE_with_Num = cbind(average_outputs_ATE, c(0,average_outputs_NumTreatmentWeights))
average_outputs_ATE_null_hypothesis_with_Num = cbind(average_outputs_ATE_null_hypothesis, c(0,average_outputs_NumTreatmentWeights))
average_outputs_ATEC_with_Num = cbind(average_outputs_ATEC, average_outputs_NumControlWeights)
average_outputs_ATET_with_Num = cbind(average_outputs_ATET, average_outputs_NumTreatmentWeights)

average_outputs_ATE_with_Num_ConsUncons = average_outputs_ATE_with_Num[c(1:9), ]
average_outputs_ATE_with_Num_Others = average_outputs_ATE_with_Num[-c(2:9), ]

my.file.name = as.character(paste0(Sys.Date(), "Average_Different_Optimization_Methods_ATE_ConsUncons.txt"))
write.table(average_outputs_ATE_with_Num_ConsUncons, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)

my.file.name = as.character(paste0(Sys.Date(), "Average_Different_Optimization_Methods_ATE_Others.txt"))
write.table(average_outputs_ATE_with_Num_Others, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)

my.file.name = as.character(paste0(Sys.Date(), "Average_Different_Optimization_Methods_ATE_null_hypothesis.txt"))
write.table(average_outputs_ATE_null_hypothesis_with_Num, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)

# my.file.name = as.character(paste0(Sys.Date(), "Average_Different_Optimization_Methods_ATEC.txt"))
# write.table(average_outputs_ATEC_with_Num, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)

my.file.name = as.character(paste0(Sys.Date(), "Average_Different_Optimization_Methods_ATET.txt"))
write.table(average_outputs_ATET_with_Num, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)

#--------------------------------------------#
#----- randomized treatment assignments -----#
#--------------------------------------------#

# #--- write to files ---#
# my.file.name = as.character(paste0(Sys.Date(), "Summary_RADiM.txt"))
# write.table(average_outputs_RADiM, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
# 
# my.file.name = as.character(paste0(Sys.Date(), "Summary_StADiM.txt"))
# write.table(average_outputs_StADiM, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
# 
# my.file.name = as.character(paste0(Sys.Date(), "Summary_RARegAdj.txt"))
# write.table(average_outputs_RARegAdj, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
# 
# my.file.name = as.character(paste0(Sys.Date(), "Summary_RA1NN.txt"))
# write.table(average_outputs_RA1NN, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)
# 
# my.file.name = as.character(paste0(Sys.Date(), "Summary_RA5NN.txt"))
# write.table(average_outputs_RA5NN, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)

#--- compare with randomization: report RMSE only ----#
average_outputs_compare_with_randomization = cbind(average_outputs_ATE[3:9,8],
                                                   average_outputs_RADiM[2:8,8], 
                                                   average_outputs_StADiM[2:8,8],
                                                   average_outputs_RARegAdj[2:8,8],
                                                   average_outputs_RA1NN[2:8,8],
                                                   average_outputs_RA5NN[2:8,8])
my.file.name = as.character(paste0(Sys.Date(), "Compare_With_Randomization.txt"))
write.table(average_outputs_compare_with_randomization, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)



end_time <- Sys.time()
end_time - start_time
