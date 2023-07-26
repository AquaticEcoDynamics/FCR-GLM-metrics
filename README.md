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
- error_stats.csv - MEF of all models predicting temperature, dissolved oxygen, thermocline depth, Schmidt stability, metalimnetic oxygen minima and the number of anoxic layers during the calibration and validation periods. The MEF values in column 'Calibration.deepm2' were used to create Table 2.
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

Table 1 includes the weighting of the extra metrics observation groups for each deep mixing 2 model. These weights can be found for each model in the corresponding control files (glm3.pst) in folders 'PEST_EXM_weight1_deepm2', 'PEST_EXM_weight2_deepm2', 'PEST_EXM_weight3_deepm2' in the third column of the *observation data section. 

## Uncertainty_analysis
This folder includes the setup of the uncertainty analysis using pestpp-ies, a model-independent iterative smoother. The following commands start the uncertainty analysis:
1. Changing directory to the 'Uncertainty analysis folder' (cd .../FCR-GLM_metrics/Uncertainty_analysis).
2. python3 workflow.py (this initiates the pestpp-ies process with 10 workers on the local computer).

The prior and posterior model output files used for creating Figure 8 are also included in the folder ('glm3_reweight_ies.0.obs.csv' = prior, 'glm3_reweight_ies.3.obs.csv' = posterior).

## R-scripts
This repo includes all the R scripts for the calculation of observed extra metrics, for creating figures and calculating error metrics. 

- Figure3_Countour_plots.R -  This script creates the contour plots of modelled, observed, and the difference between modelled and observed temperature and dissolved oxygen profiles based on the naive calibration model with deep mixing 2 (‘Results/Figure 3’). 
- Figure4_Extra_metrics_naive.R - This script provides visual comparison of the observed and modelled system-metrics including thermocline depth during the stratified period, ice cover presence and absence in the winter period, Schmidt stability, sediment temperature in zone 2  with modelled and observed water temperatures at 5 m depth in zone 2, spatial and temporal extent of anoxia and the metalimnetic oxygen minimum (‘Results/Figure 4’). The modelled values are based on the naive model with deep mixing configuration set to 2. Additionally, it includes calculation of the MEF of the thermocline depth, Schmidt stability, anoxia, and the metalimnetic oxygen minimum where the modelled values are based on the naive model with deep mixing 2. The calculated MEF values are saved to ‘observations/error_stats.csv’.
- Figure5_Oxy_temp_6panel_plos.R - This script provides visual comparison of the observed and modelled temperature and dissolved oxygen in the epilimnion, metalimnion, and hypolimnion. The model predictions are based on the three system-inspired models (PEST_exm_w1, PEST_exm_w2, PEST_exm_w3) with different extra metrics weighting schemes (‘Results/Figure 5’). Additionally, it includes calculation of the MEF of the temperature and dissolved oxygen, where the modelled values are predicted by the three system-inspired models (PEST_exm_w1, PEST_exm_w2, PEST_exm_w3). The calculated MEF values are saved to ‘observations/error_stats.csv’.
- Figure6_Extra_metrics_w1_2_3.R - This script provides visual comparison of the observed system-metrics and the system-metrics predicted by the three models with different weighting schemes based on the system-inspired approach (PEST_exm_w1, PEST_exm_w2, PEST_exm_w3; ‘Results/Figure 6’). The metrics include thermocline depth, Schmidt stability, anoxia, and the metalimnetic oxygen minimum. Additionally, it includes calculation of the MEF of the thermocline depth, Schmidt stability, anoxia, and the metalimnetic oxygen minimum where the modelled values are predicted by the three models with different weighting schemes based on the system-inspired approach (PEST_exm_w1, PEST_exm_w2, PEST_exm_w3). The calculated MEF values are saved to ‘observations/error_stats.csv’. Note: The last section of the script, Anoxia, is to be run in three different runs based on the instructions in the script.
- Figure7_error_plot.R - This script provides visual comparison of the MEF between models with two different mixing configurations (deep mixing 1 and 2) during the calibration and validation period (‘Results/Figure 7’). The MEF values are read from the table ‘observations/error_stats.csv’, which was populated by different scripts. Additionally, it includes the calculation of the MEF of temperature and dissolved oxygen for both the naive deep mixing 1 and naive deep mixing 2 models.
- Figure8_Uncertainty_plot.R - This script creates fan plots to visualise the prediction uncertainty of the thermocline depth and the metalimentic oxygen minima pre and post calibration (‘Results/Figure 8’). The prediction uncertainty includes parameter uncertainty and uncertainty due to measurement noise. The uncertainty analysis is performed on the PEST_exm_w2 model with the deep mixing configuration set to 2.
- Figure9_Stacked_plot_Phi.R - This script creates stacked plots that show the convergence of the objective function in the case of the deep mixing 2 naive calibration and in the case of the deep mixing 2 system-inspired approach with different extra metrics weighting schemes during the calibration process (‘Results/Figure 9’).
- Observed_EXM.R - This script calculates the extra metrics thermocline depth, and Schmidt stability from the observed temperature depth profile data and calculates the extra metrics anoxia, and metalimnetic oxygen minima from the observed dissolved oxygen depth profile data. The csv files created (‘obs_td.csv’,  ‘obs_SS.csv’, ‘mom_observed.csv’, ‘anoxia_observed.csv’) are saved in the ‘observations’ folder.
- Deepm1_MEF.R - This script includes the calculation of the MEF for all deep mixing 1 models, specifically MEF of extra metrics (thermocline depth, Schmidt stability, anoxia, and metalimnetic oxygen minima) predicted by the deep mixing 1 naive model, additionally MEF of temperature, dissolved oxygen and extra metrics predicted by the three deep mixing 1 models based on the system-inspired approach with different weighting schemes. The calculated MEF values are saved to ‘observations/error_stats.csv’.  Note: The last section of the script, Anoxia, is to be run in three different runs based on the instructions in the script.
- Appendix_B_plots.R - This script creates all plots included in the Supporting information, which visualise the performance of the deep mixing 1 models. 
- Rolling_average_noise.R - This script includes the calculation of the measurement noise for observation groups temperature, oxygen, thermocline depth, Schmidt stability, metalimnetic oxygen minimum, and the number of anoxic layers.
- Validation_MEF.R - This script includes the calculation of the MEF of all models in the validation period. The MEF values are added to 'observations/error_stats.csv'.

## Results
This folder includes all the figures created by the R scripts. Their numbering matches the Figure captions in the manuscript including the figures in the Supporting information.


