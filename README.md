# Human PBPK model for PFOA in R

The purpose of the present work is to translate an existing physiologically based pharmacokinetic (PBPK) model for perfluorooctanoic acid (PFOA) to R, which will be used to study human exposure of PFOA in blood from diet and cosmetics.

### Why PFOA
PFOA is a part of a large group of compounds called per- and polyfluoroalkyl substances (PFASs), which constitute a large group of widely used chemically synthesized compounds that are water, stain, and oil repellent. They are highly persistent and several of the compounds, including PFOA, accumulates in the food chain. The health effects of PFOA in humans are still being studied, but exposure are associated with a range of health effects, including developmental and reproductive toxicity, immune system effects, and cancer. Increased knowledge of human exposure and internal dose in humans are important to link exposure of PFOA to the adverse health effects.

### The PBPK model
Several PBPK model for PFOA have been published, but the model that have been used the most was first published by Loccisano et al in 2011 (https://www.sciencedirect.com/science/article/abs/pii/S0273230010002242?via%3Dihub). This model has been widely used and was partly modified and used by the European Food Safety Authority (EFSA) in 2020 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7507523/). EFSA modeled PFOA using the Berkeley Madonna software in adults and breastfeeding children and the code was published in the EFSA opinion.

We have only used the adult part of the model and translated it from Berkeley Madonna to R. The model was modified by including an dermal exposure rout, a fecal excretion of PFOA, and a few other modifications to the kidney compartment. The first work using the PBPK model is published open access in the paper "Comparison of aggregated exposure to PFOA  from diet  and personal care products with concentrations in serum using a PBPK model – results from the Norwegian biomonitoring study in EuroMix" by Husøy et al 2023, where you can find information on the use of the model and the first results. 
Challenges and weaknesses in the PBPK model is discussed in the paper.

## Files

## Running the PBPK model

## Sensitivity analyses of the PBPK model using pksensi package in R