#=========================================#
#===== Read in all the files printed =====#
#=========================================#

#--- Last Change: 20250313 changed the MSE into RMSE (take the square root)
#--- Last Change: 20250218 included random assignment + regressin adjustment

library(tidyverse)

#--- random assignment + difference-in-means ---#
list_of_files_RADiM = list.files(path = "./output11_Randomization/", recursive = TRUE, pattern = "\\_RADiM.txt$", full.names = TRUE)
datalist_RADiM = lapply(list_of_files_RADiM, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_RADiM = lapply(datalist_RADiM, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_RADiM = Reduce("+", datalist_RADiM)/ length(datalist_RADiM)

#--- stratified assignment + difference-in-means ---#
list_of_files_StADiM = list.files(path = "./output11_Randomization/", recursive = TRUE, pattern = "\\_StADiM.txt$", full.names = TRUE)
datalist_StADiM = lapply(list_of_files_StADiM, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_StADiM = lapply(datalist_StADiM, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_StADiM = Reduce("+", datalist_StADiM)/ length(datalist_StADiM)

#--- random assignment + regression adjustment ---#
list_of_files_RARegAdj = list.files(path = "./output11_Randomization/", recursive = TRUE, pattern = "\\_RARegAdj.txt$", full.names = TRUE)
datalist_RARegAdj = lapply(list_of_files_RARegAdj, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_RARegAdj = lapply(datalist_RARegAdj, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_RARegAdj = Reduce("+", datalist_RARegAdj)/ length(datalist_RARegAdj)

#--- random assignment + 1 nearest neighbor matching ---#
list_of_files_RA1NN = list.files(path = "./output11_Randomization/", recursive = TRUE, pattern = "\\_RA1NN.txt$", full.names = TRUE)
datalist_RA1NN = lapply(list_of_files_RA1NN, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_RA1NN = lapply(datalist_RA1NN, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_RA1NN = Reduce("+", datalist_RA1NN)/ length(datalist_RA1NN)

#--- random assignment + 5 nearest neighbor matching ---#
list_of_files_RA5NN = list.files(path = "./output11_Randomization/", recursive = TRUE, pattern = "\\_RA5NN.txt$", full.names = TRUE)
datalist_RA5NN = lapply(list_of_files_RA5NN, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_RA5NN = lapply(datalist_RA5NN, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_RA5NN = Reduce("+", datalist_RA5NN)/ length(datalist_RA5NN)

#--------------------------------------------#
#--- check which ones are missing because ---#
#--- some servers are occasionally down   ---#
#--------------------------------------------#

list_of_files_RADiM_NoSpaces = gsub(" ", "", list_of_files_RADiM)
list_of_files_RADiM_Trimmed = sapply(list_of_files_RADiM_NoSpaces, substring, 12)
list_of_files_RADiM_seeds = as.numeric(gsub(".*?([0-9]+).*", "\\1", list_of_files_RADiM_Trimmed))
list_of_files_RADiM_seeds = sort(list_of_files_RADiM_seeds)
failed.attempts = setdiff(c(1:1000), list_of_files_RADiM_seeds)
print(failed.attempts)


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
average_outputs_compare_with_randomization = cbind(average_outputs_RADiM[,8], 
                                                   average_outputs_StADiM[,8],
                                                   average_outputs_RARegAdj[,8],
                                                   average_outputs_RA1NN[,8],
                                                   average_outputs_RA5NN[,8])
my.file.name = as.character(paste0(Sys.Date(), "Summary_compare_with_randomization.txt"))
write.table(average_outputs_compare_with_randomization, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)

