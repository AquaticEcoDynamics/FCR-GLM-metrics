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

