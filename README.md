# Human PBPK model for PFOA in R

The purpose of the present work is to translate an existing physiologically based pharmacokinetic (PBPK) model for perfluorooctanoic acid (PFOA) to R, which will be used to study human exposure of PFOA in blood from diet and cosmetics.

### Why PFOA
PFOA is a part of a large group of compounds called per- and polyfluoroalkyl substances (PFASs), which are widely used chemically synthesized compounds that are water, stain, and oil repellent. They are highly persistent and several of the compounds, including PFOA, accumulates in the food chain. The health effects of PFOA in humans are still being studied, but exposure are associated with a range of health effects, including developmental and reproductive toxicity, immune system effects, and cancer. Increased knowledge of human exposure and internal dose in humans are important to link exposure of PFOA to the adverse health effects.

### The PBPK model
Several PBPK model for PFOA have been published, but the model that have been used the most was first published by Loccisano et al in 2011 (https://www.sciencedirect.com/science/article/abs/pii/S0273230010002242?via%3Dihub). This model has been widely used and was partly modified and used by the European Food Safety Authority (EFSA) in 2020 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7507523/). EFSA modeled PFOA using the Berkeley Madonna software in adults and breastfeeding children and the code was published in the EFSA opinion. 

We have only used the adult part of the model and translated it from Berkeley Madonna to R. The model was modified by including an dermal exposure rout, a fecal excretion of PFOA, and a few other modifications to the kidney compartment. The first work using the PBPK model is published open access in the paper "Comparison of aggregated exposure to PFOA  from diet  and personal care products with concentrations in serum using a PBPK model – results from the Norwegian biomonitoring study in EuroMix" by Husøy et al 2023, where you can find information on the use of the model and the first results. 
Challenges and weaknesses in the PBPK model is discussed in the paper.


## Files

### Files used for the exposure estimate for PFOA from the diet
The R code to run the dietary exposure assessment  
[PFOA_food_version1.R](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Code/PFOA_food_version1.R)  
The data files as input to the dietary exposure  
[3-SumPFAS_food_conc_LB.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/3-SumPFAS_food_conc_LB.csv)  
[3-SumPFAS_food_conc_MB.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/3-SumPFAS_food_conc_MB.csv)  
[3-SumPFAS_food_conc_UB.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/3-SumPFAS_food_conc_UB.csv)  
[foodintake_dummy_day1.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/foodintake_dummy_day1.csv)  
[foodintake_dummy_day2.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/foodintake_dummy_day2.csv)  
[EuroMix_dummy_sex_weight.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/EuroMix_dummy_sex_weight.csv)

### Files used for the exposure estimates from the PCPs
The R code to run the exposure assessment from PCPs  
[PFOA_PCP_version1.R](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Code/PFOA_PCP_version1.R)  
The data files as input to the exposure estimates from PCPs  
[2_SumPCPPFAS_LB_050122.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/2_SumPCPPFAS_LB_050122.csv)  
[2_SumPCPPFAS_MB_050122.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/2_SumPCPPFAS_MB_050122.csv)  
[2_SumPCPPFAS_UB_050122.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/2_SumPCPPFAS_UB_050122.csv)  
[AmountsPerApplication_female.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/AmountsPerApplication_female.csv)  
[AmountsPerApplication_male.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/AmountsPerApplication_male.csv)  
[EEF_distributions_female.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/EEF_distributions_female.csv)  
[EEF_distributions_male.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/EEF_distributions_male.csv)  
[PCP_frequency_dummy.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/PCP_frequency_dummy.csv)  
[EuroMix_dummy_sex_weight.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/EuroMix_dummy_sex_weight.csv)  

### Files used for running the PBPK model
The R code to run the PBPK model  
[PBPK_PFOA_version1_150323.R](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Code/PBPK_PFOA%20_version1_150323.R)  
The data files as input to the PBPK modeling from the exposure  
[SumPFOA_LB_food_PCP.csv](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Data/foodintake_dummy_day1.csv)

### Files used for the sensitivity analyses of the PBPK model
[PBPK_PFOA_pksensi_version1.R](https://github.com/TrineHusoy/PBPK_PFOA/blob/main/Code/PBPK_PFOA_pksensi_version1.R)  

## Exposure assesment from diet and PCPs
Establish the folders "Code", "Data" and "Results" in your directory and copy the relevant files to the appropriate folders. Open the R or Rmd code in RStudio and insert your work directory in the code. The data files will be uploaded and create a new folder under "Results" with the date of the day to store the results. Run the exposure assessment for the PCP first as the output will be merged with the results from the food later. This can be input to the PBPK model.

## Running the PBPK model
Normally I would create a R project linked to github as a start of my work, but you can also run the PBPK model just by open it in R. Upload the dummy data file from the Data folder and run the model. The model includes exposure from both diet and personal care products, and you can run one at the time by setting the exposure to zero for one of them. Run the model from ###START PBPK## to ###END PBPK###. 

## Sensitivity analyses of the PBPK model using pksensi package in R

The PBPK model in the sensitivity analysis is presented separately for now. Here only one concentration for the diet and PCP is tested in each run, and the values have to be keyed in manually. We have used R package pksensi developed by Hsieh et al 2020 (https://pubmed.ncbi.nlm.nih.gov/33426260/), which is a tool for global sensitivity analysis of PBPK models. The repository for the published version of pksensi you can find here here https://github.com/ElsevierSoftwareX/SOFTX_2020_29.
We have included the parameters from the PBPK model of PFOA that we believe is the most important ones, and the output will give you an uncertainty evaluation and a sensitivity analyse of all the included parameters.  



