# SCDesign
Synthetic Control Design

This is a list of all provided materials' names and usages. All the running times are reported using a personal laptop with the following specifications:

	Operating System: Windows 11 Enterprise version 23H2
	Processor: Intel(R) Core(TM) Ultra 5 135H   3.60 GHz
	RAM: 16.0 GB (15.6 GB usable)
	System type: 64-bit operating system, x64-based processor 


There are three folders, each containing files to generate a part of the Figures or Tables in the paper.

Folder names, file names, usages, and running times:

1. Folder: Figure 4-5, Figure OA5-OA8
1.1. File: SCdesign_LazyRun.R
---- Generates one simulation environment as described in Section 5 of our paper.
---- In generating the simulation environment, need to specify "noise.variance = 1", "noise.variance = 5", or "noise.variance = 10" in Line 131.
---- Produces Figures 4-5 as described in Section 5.1 of our paper under "noise.variance = 1"; produces Figures OA5-OA8 in Section OA7 in the Online Appendix under "noise.variance = 5" and "noise.variance = 10".
---- Running time: for each parameter of "noise.variance = 1", "noise.variance = 5", or "noise.variance = 10", running time is less than 1.5 hours.
---- Output files: 
------ (1) ObservedData_NoiseVariance=1.png: Figure 4
------ (2) Residuals_NoiseVariance=1.png: Figure 5
------ (3) ObservedData_NoiseVariance=5.png: Figure OA5
------ (4) Residuals_NoiseVariance=5.png: Figure OA6
------ (5) ObservedData_NoiseVariance=10.png: Figure OA7
------ (6) Residuals_NoiseVariance=10.png: Figure OA8





2. Folder: Table 2-5, Table OA1-2
2.1. File: Different_optimization_methods.R
---- Generates one simulation environment as described in Section 5 of our paper.
---- This file is run through online cluster computing, which sends 1000 files each with a different seed to the server using Line 29 "repetition.RANDOM.SEED = as.numeric(Sys.getenv("SGE_TASK_ID"))". 
---- For a test run, set a fixed random seed by commenting out Line 29 "repetition.RANDOM.SEED = as.numeric(Sys.getenv("SGE_TASK_ID"))" and setting "repetition.RANDOM.SEED = 1".
---- For each random seed (e.g., repetition.RANDOM.SEED = 1), produces 6 intermediate files into the "output11" folder.
---- Running time: for each parameter (e.g., repetition.RANDOM.SEED = 1), running time is less than 12 hours.
---- Intermediate output files:
------ (1) 1Different_Optimization_Methods_ATE.txt
------ (2) 1Different_Optimization_Methods_ATE_null_hypothesis.txt
------ (3) 1Different_Optimization_Methods_ATEC.txt: not useful, for sanity check only
------ (4) 1Different_Optimization_Methods_ATET.txt
------ (5) 1Different_Optimization_Methods_ControlWeights.txt: not useful, for sanity check only
------ (6) 1Different_Optimization_Methods_TreatmentWeights.txt: not useful, for sanity check only

2.2. File: Analyzing_Outputs.R
---- After running Different_optimization_methods.R 1000 times using 1000 random seeds "repetition.RANDOM.SEED = 1" ~ "repetition.RANDOM.SEED = 1000", there are 6000 files in the "output11" folder. Run this code to combine all the intermediate files.
---- The output files contain more rows than in the paper. We only present five parameters beta = 0.01, 0.1, 1, 10, 100, five parameters xi = 0.01, 0.1, 1, 10, 100, and five parameters lambda = 0.01, 0.1, 1, 10, 100 in the paper. But the output files contain eight parameters beta = 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, eight parameters xi = 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, and eight parameters lambda = 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000.
---- Running time: less than 1 minute.
---- Output files:
------ (1) "DATE"Average_Different_Optimization_Methods_ATE.txt: Table 2
------ (2) "DATE"Average_Different_Optimization_Methods_ATE_null_hypothesis.txt: Table 4
------ (3) "DATE"Average_Different_Optimization_Methods_ATEC.txt: not useful, for sanity check only
------ (4) "DATE"Average_Different_Optimization_Methods_ATET.txt: Table OA1

2.3. File: Nonlinear2.R
---- Generates one simulation environment as described in Section 5.2.2 of our paper.
---- This file is run through online cluster computing, which sends 1000 files each with a different seed to the server using Line 21 "repetition.RANDOM.SEED = as.numeric(Sys.getenv("SGE_TASK_ID"))". 
---- For a test run, set a fixed random seed by commenting out Line 21 "repetition.RANDOM.SEED = as.numeric(Sys.getenv("SGE_TASK_ID"))" and setting "repetition.RANDOM.SEED = 1".
---- For each random seed (e.g., repetition.RANDOM.SEED = 1), produces 6 intermediate files into the "output11_nonlinear" folder.
---- Running time: for each parameter (e.g., repetition.RANDOM.SEED = 1), running time is less than 18 hours.
---- Intermediate output files:
------ (1) 1Different_Optimization_Methods_ATE.txt
------ (2) 1Different_Optimization_Methods_ATE_null_hypothesis.txt: not useful, for sanity check only
------ (3) 1Different_Optimization_Methods_ATEC.txt: not useful, for sanity check only
------ (4) 1Different_Optimization_Methods_ATET.txt
------ (5) 1Different_Optimization_Methods_ControlWeights.txt: not useful, for sanity check only
------ (6) 1Different_Optimization_Methods_TreatmentWeights.txt: not useful, for sanity check only

2.4. File: Analyzing_Outputs_Nonlinear2.R
---- After running Nonlinear2.R 1000 times using 1000 random seeds "repetition.RANDOM.SEED = 1" ~ "repetition.RANDOM.SEED = 1000", there are 6000 files in the "output11_nonlinear" folder. Run this code to combine all the intermediate files.
---- The output files contain more rows than in the paper. We only present five parameters beta = 0.01, 0.1, 1, 10, 100, five parameters xi = 0.01, 0.1, 1, 10, 100, and five parameters lambda = 0.01, 0.1, 1, 10, 100 in the paper. But the output files contain eight parameters beta = 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, eight parameters xi = 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, and eight parameters lambda = 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000.
---- Running time: less than 1 minute.
---- Output files:
------ (1) "DATE"Nonlinear2_Average_Different_Optimization_Methods_ATE.txt: Table 3
------ (2) "DATE"Nonlinear2_Average_Different_Optimization_Methods_ATE_null_hypothesis.txt: not useful, for sanity check only
------ (3) "DATE"Nonlinear2_Average_Different_Optimization_Methods_ATEC.txt: not useful, for sanity check only
------ (4) "DATE"Nonlinear2_Average_Different_Optimization_Methods_ATET.txt: Table OA2

2.5. File: rerand.py
---- Created by Max Cytrynbaum.
---- Algorithm for stratification as described in Section 5.2.4 of our paper.
---- Do not need to run this file. This file will be used by RandomAssignment.R to do stratified randomization.

2.6. File: RandomAssignment.R 
---- Generates one simulation environment as described in Section 5.2.4 of our paper.
---- This file is run through online cluster computing, which sends 1000 files each with a different seed to the server using Line 35 "repetition.RANDOM.SEED = as.numeric(Sys.getenv("SGE_TASK_ID"))". 
---- For a test run, set a fixed random seed by commenting out Line 35 "repetition.RANDOM.SEED = as.numeric(Sys.getenv("SGE_TASK_ID"))" and setting "repetition.RANDOM.SEED = 1".
---- For each random seed (e.g., repetition.RANDOM.SEED = 1), produces 5 intermediate files into the "output11_Randomization" folder.
---- Running time: for each parameter (e.g., repetition.RANDOM.SEED = 1), running time is less than 2 hours.
---- Intermediate output files:
------ (1) 1RandomAssignment_RA1NN.txt
------ (2) 1RandomAssignment_RA5NN.txt
------ (3) 1RandomAssignment_RADiM.txt
------ (4) 1RandomAssignment_RARegAdj.txt
------ (5) 1RandomAssignment_StADiM.txt

2.7. File: Analyzing_Outputs_RandomAssignment.R 
---- After running RandomAssignment.R 1000 times using 1000 random seeds "repetition.RANDOM.SEED = 1" ~ "repetition.RANDOM.SEED = 1000", there are 5000 files in the "output11_Randomization" folder. Run this code to combine all the intermediate files.
---- Running time: less than 1 minute.
---- Output file:
------ (1) "DATE"Summary_compare_with_randomization.txt: Table 5





3. Folder: Figure 2-3, OA1-4, Table 1
3.1. File: Walmart.csv
---- Data file from Kaggle

3.2. File: Walmart_LazyRun.R
---- Reads Walmart.csv data and conducts synthetic control design as described in Section 4 of our paper.
---- Produces 2 figures of the constrained formulation for each parameter of m-bar = 1, 2, 3, 4, 5, into the "output" folder. How these figures are selected are described in Section 4 of our paper.
---- Running time: for all 5 parameters, totaling less than 6 hours.
---- Produces some intermediate outputs that will be used in Table 1 of our paper.
---- Intermediate output file:
------ (1) "DATE"SC_cardinality_matrix.txt: reports RMSE of synthetic control design.
---- Output files (e.g., m-bar = 1):
------ (1) Fitted_Y_original_no_covariates_weekly_Constrained_K=1_"DATE".png: Figure OA1
------ (2) ResidualsCI_Y_original_no_covariates_weekly_Constrained_K=1_"DATE".png: Figure OA2
------ (3) Fitted_Y_original_no_covariates_weekly_Constrained_K=2_"DATE".png: Figure 1
------ (4) ResidualsCI_Y_original_no_covariates_weekly_Constrained_K=2_"DATE".png: Figure 2
------ (5) Fitted_Y_original_no_covariates_weekly_Constrained_K=3_"DATE".png: Figure OA3
------ (6) ResidualsCI_Y_original_no_covariates_weekly_Constrained_K=3_"DATE".png: Figure OA4
------ (7) Fitted_Y_original_no_covariates_weekly_Constrained_K=4_"DATE".png: not useful
------ (8) ResidualsCI_Y_original_no_covariates_weekly_Constrained_K=4_"DATE".png: not useful
------ (9) Fitted_Y_original_no_covariates_weekly_Constrained_K=5_"DATE".png: not useful
------ (10) ResidualsCI_Y_original_no_covariates_weekly_Constrained_K=6_"DATE".png: not useful

3.3. File: rerand.py
---- Created by Max Cytrynbaum.
---- Algorithm for stratification as described in Section 5.2.4 of our paper.
---- Do not need to run this file. This file will be used by Walmart_Randomization.R to do stratified randomization.

3.4. File: Walmart_Randomization.R
---- Generates one simulation as described in Section 4 of our paper.
---- This file is run through online cluster computing, which sends 1000 files each with a different seed to the server using Line 14 "repetition.RANDOM.SEED = as.numeric(Sys.getenv("SGE_TASK_ID"))". 
---- For a test run, set a fixed random seed by commenting out Line 14 "repetition.RANDOM.SEED = as.numeric(Sys.getenv("SGE_TASK_ID"))" and setting "repetition.RANDOM.SEED = 1".
---- For each random seed (e.g., repetition.RANDOM.SEED = 1), produces 5 intermediate files into the "output_Randomization" folder.
---- Running time: for each parameter (e.g., repetition.RANDOM.SEED = 1), running time is less than 2 hours.
---- Intermediate output files:
------ (1) 1RandomAssignment_RA1NN.txt
------ (2) 1RandomAssignment_RA5NN.txt
------ (3) 1RandomAssignment_RADiM.txt
------ (4) 1RandomAssignment_RARegAdj.txt: Regression adjustment is the same as difference in means estimation because there is no observed covariates.
------ (5) 1RandomAssignment_StADiM.txt

3.5. File: Analyzing_Outputs_RandomAssignment.R 
---- After running Walmart_Randomization.R 1000 times using 1000 random seeds "repetition.RANDOM.SEED = 1" ~ "repetition.RANDOM.SEED = 1000", there are 5000 files in the "output_Randomization" folder. Run this code to combine all the intermediate files.
---- Run this code to also combine the intermediate file that Walmart_LazyRun.R generated for synthetic control design.
---- Running time: less than 1 minute.
---- Output file:
------ (1) "DATE"Summary_compare_with_randomization.txt: Table 1
