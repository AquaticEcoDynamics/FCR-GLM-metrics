# FCR-GLM-metrics
FCR-GLMv3.3 GLM-AED simulation calibrated with PEST.

## Calibrated_models
This repo includes the setup of all models based on the best parameter set attained through calibration. The sims only run on Linux. The time period of the simulation is from 01-12-2016 to 31-12-2019 (calibration period). The models are run by changing directory to the folder where the model files are located (cd …/FCR-GLM-metrics/Calibrated_models/Deepm*) and running the following command: ./glm+_latest/glm+.

- Deepm2_naive: The naive model with deep mixing configuration set to 2. 
- Deepm2_exm_weight1: The system-inspired model with weighting scheme 1 and deep mixing configuration 2. 
- Deepm2_exm_weight2: The system-inspired model with weighting scheme 2 and deep mixing configuration 2. 
- Deepm2_exm_weight3: The system-inspired model with weighting scheme 3 and deep mixing configuration 2. 
- Deepm1_naive: The naive model with deep mixing configuration set to 1. 
- Deepm1_exm_weight1: The system-inspired model with weighting scheme 1 and deep mixing configuration 1. 
- Deepm1_exm_weight2: The system-inspired model with weighting scheme 2 and deep mixing configuration 1. 
- Deepm1_exm_weight3: The system-inspired model with weighting scheme 3 and deep mixing configuration 1. 

## Validation
This repo includes the setup of all models based on the best parameter set attained through calibration. The sims only run on Linux. The time period of the simulation is from 19-07-2015 to 02-12-2016 (validation period). The models are run by changing directory to the folder where the model files are located (cd …/FCR-GLM-metrics/Validation/valid_*) and running the following command: ./glm+_latest/glm+.

- valid_naive_deepm2: The naive model with deep mixing configuration set to 2. 
- valid_exm_deepm2_w1: The system-inspired model with weighting scheme 1 and deep mixing configuration 2.
- valid_exm_deepm2_w2: The system-inspired model with weighting scheme 2 and deep mixing configuration 2.
- valid_exm_deepm2_w3: The system-inspired model with weighting scheme 3 and deep mixing configuration 2.
- valid_naive_deepm1: The naive model with deep mixing configuration set to 1. 
- valid_exm_deepm1_w1: The system-inspired model with weighting scheme 1 and deep mixing configuration 1.
- valid_exm_deepm1_w2: The system-inspired model with weighting scheme 2 and deep mixing configuration 1.
- valid_exm_deepm1_w3: The system-inspired model with weighting scheme 3 and deep mixing configuration 1.

## Observations
- bathymetry.csv - Bathymetry dataset of the Falling Creek Reservoir, which was used for the calculation of Schmidt stability.
- CleanedObsTemp.csv - Observed temperature profiles data for the Falling Creek Reservoir retrieved from the Environmental Data Initiative Repository (Carey et al., 2022).
- CleanedObsOxy.csv - Observed dissolved oxygen profiles data for the Falling Creek Reservoir retrieved from the Environmental Data Initiative Repository (Carey et al., 2022).
- Ice_Data_2013_2022.csv - Ice cover data for the Falling Creek Reservoir retrieved from the Environmental Data Initiative Repository (Carey & Breef-Pilz, 2022).
- Obs_SS.csv - Observed Schmidt stability data calculated from the observed temperature profile and bathymetry data by ‘Rscripts/Observed_EXM.R’.
- obs_td.csv - Observed thermocline depth data calculated from the observed temperature profile data by ‘Rscripts/Observed_EXM.R’.
- mom_observed.csv - Observed metalimnetic oxygen minima data calculated from the observed dissolved oxygen profile data by ‘Rscripts/Observed_EXM.R’.
- anoxia_observed.csv - Number of anoxic layers calculated from the observed dissolved oxygen profile data by ‘Rscripts/Observed_EXM.R’.
- error_stats.csv - MEF of all models predicting temperature, dissolved oxygen, thermocline depth, Schmidt stability, metalimnetic oxygen minima and the number of anoxic layers. 
- obs_oxy_interpolated.csv - Interpolation of the observed dissolved oxygen data. Spatial interpolation by 0.1m and temporal interpolation by day. 
- obs_temp_inerpolated.csv - Interpolation of the observed temperature data. Spatial interpolation by 0.1m and temporal interpolation by day. 

## PEST_calibration_runs
This folder includes the setup of the PEST calibration for all models. The calibration only runs on Linux. The following commands start the calibration process:
1. Changing directory to the folder where the model files are located (cd .../FCR-GLM-metrics/PEST_calibration_runs/PEST_*).
2. ./bin/pestpp-glm glm3.pst
Note: The working directory and the path for saving the csv files in the Rscript PestFCR.R has to be set to the working directory in step 1. The correct format of the date column in the input files is  ‘yyyy-mm-dd’ for inflow and outflow files, and ‘yyyy-mm-dd HH:MM:SS’ for the met file. All the four deep mixing 2 models include an additional file in their calibration setup folders called glm3.csv. This is an output file of the PEST calibration which records the value of phi in each iteration. The glm3.csv files were used to create ‘Results/Figure9’ by ‘Rscripts/Figure9_Stacked_plot_Phi.R’.

- PEST_routine_deepm2: This folder includes the setup of the PEST calibration for the deep mixing 2 naive model. 
- PEST_EXM_weight1_deepm2: This folder includes the setup of the PEST calibration for the deep mixing 2 with weighting scheme 1 based on the system-inspired approach.  
- PEST_EXM_weight2_deepm2: This folder includes the setup of the PEST calibration for the deep mixing 2 with weighting scheme 2 based on the system-inspired approach.  
- PEST_EXM_weight3_deepm2: This folder includes the setup of the PEST calibration for the deep mixing 2 with weighting scheme 3 based on the system-inspired approach. 
- PEST_routine_deepm1 folder: This folder includes the setup of the PEST calibration for the deep mixing 1 naive model. 
- PEST_EXM_weight1_deepm1: This folder includes the setup of the PEST calibration for the deep mixing 1 with weighting scheme 1 based on the system-inspired approach.  
- PEST_EXM_weight2_deepm1: This folder includes the setup of the PEST calibration for the deep mixing 1 with weighting scheme 2 based on the system-inspired approach.  
- PEST_EXM_weight3_deepm1: This folder includes the setup of the PEST calibration for the deep mixing 1 with weighting scheme 3 based on the system-inspired approach. 

