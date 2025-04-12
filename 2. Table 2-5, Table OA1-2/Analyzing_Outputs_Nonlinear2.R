#=========================================#
#===== Read in all the files printed =====#
#=========================================#

#--- Last Change: 20250316 2PM --- changed the MSE into RMSE to make it comparable with MAE (take the square root)

library(tidyverse)

start_time <- Sys.time()

list_of_files_ATE = list.files(path = "./output11_nonlinear/", recursive = TRUE, pattern = "\\_ATE.txt$", full.names = TRUE)
datalist_ATE = lapply(list_of_files_ATE, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_ATE = lapply(datalist_ATE, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_ATE = Reduce("+", datalist_ATE)/ length(datalist_ATE)

list_of_files_ATE_null_hypothesis = list.files(path = "./output11_nonlinear/", recursive = TRUE, pattern = "\\_ATE_null_hypothesis.txt$", full.names = TRUE)
datalist_ATE_null_hypothesis = lapply(list_of_files_ATE_null_hypothesis, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_ATE_null_hypothesis = lapply(datalist_ATE_null_hypothesis, function(mat) {
  mat[, 8] = sqrt(mat[, 8])  # Apply square root to 8th column
  return(mat)
})
average_outputs_ATE_null_hypothesis = Reduce("+", datalist_ATE_null_hypothesis)/ length(datalist_ATE_null_hypothesis)

list_of_files_ATEC = list.files(path = "./output11_nonlinear/", recursive = TRUE, pattern = "\\_ATEC.txt$", full.names = TRUE)
datalist_ATEC = lapply(list_of_files_ATEC, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_ATEC = lapply(datalist_ATEC, function(mat) {
  mat[, 13] = sqrt(mat[, 13])  # Apply square root to 13th column
  return(mat)
})
average_outputs_ATEC = Reduce("+", datalist_ATEC)/ length(datalist_ATEC)

list_of_files_ATET = list.files(path = "./output11_nonlinear/", recursive = TRUE, pattern = "\\_ATET.txt$", full.names = TRUE)
datalist_ATET = lapply(list_of_files_ATET, function(x) read.table(x, header = FALSE, sep = "\t", dec = "."))
datalist_ATET = lapply(datalist_ATET, function(mat) {
  mat[, 13] = sqrt(mat[, 13])  # Apply square root to 13th column
  return(mat)
})
average_outputs_ATET = Reduce("+", datalist_ATET)/ length(datalist_ATET)


list_of_files_ATE_NoSpaces = gsub(" ", "", list_of_files_ATE)
list_of_files_ATE_Trimmed = sapply(list_of_files_ATE_NoSpaces, substring, 22)
list_of_files_ATE_seeds = as.numeric(gsub(".*?([0-9]+).*", "\\1", list_of_files_ATE_Trimmed))
list_of_files_ATE_seeds = sort(list_of_files_ATE_seeds)
failed.attempts = setdiff(c(1:1000), list_of_files_ATE_seeds)
print(failed.attempts)



my.file.name = as.character(paste0(Sys.Date(), "Nonlinear2_Average_Different_Optimization_Methods_ATE.txt"))
write.table(average_outputs_ATE, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)

my.file.name = as.character(paste0(Sys.Date(), "Nonlinear2_Average_Different_Optimization_Methods_ATE_null_hypothesis.txt"))
write.table(average_outputs_ATE_null_hypothesis, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)

my.file.name = as.character(paste0(Sys.Date(), "Nonlinear2_Average_Different_Optimization_Methods_ATEC.txt"))
write.table(average_outputs_ATEC, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)

my.file.name = as.character(paste0(Sys.Date(), "Nonlinear2_Average_Different_Optimization_Methods_ATET.txt"))
write.table(average_outputs_ATET, file = my.file.name, append = FALSE, sep = "\t", col.names=FALSE)



end_time <- Sys.time()
end_time - start_time
