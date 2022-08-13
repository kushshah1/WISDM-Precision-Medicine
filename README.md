# WISDM-Precision-Medicine
Individualized treatment rule (ITR) estimation for the WISDM T1D clinical trial (Pratley et al., [2020](https://jamanetwork.com/journals/jama/article-abstract/2767159)).

## Introduction
This repository contains an end-to-end precision medicine workflow which estimates an individualized treatment rule (ITR) for older adults with type 1 diabetes (T1D). In this trial, continuous glucose monitoring (CGM) was compared with traditional blood glucose monitoring (BGM) for its effect on hypoglycemia (reduction in time spent in hypoglycemia, defined as blood glucose < 70 mg/dL).

## Features of the Workflow
- Data Processing:
  - Compilation of variables, treatment information, and clinical outcome (and omission of patients lost to follow-up)
  - Selection of final features for analysis
- Model Fitting:
  - Implementation of interpretable precision medicine algorithms: Decision List (Zhang et al., [2018](https://www.tandfonline.com/doi/pdf/10.1080/01621459.2017.1345743?casa_token=S770pwuOCIAAAAAA:Uw1DUpN2Rz5luTVEbH7SHH9kFH__RuHgcrnqIGljnqN2ddv-prxMhHzt6VrTgweZVuVIK16gEUQH3w)) and Policy Tree (Zhou et al., [2022](https://pubsonline.informs.org/doi/pdf/10.1287/opre.2022.2271?casa_token=rxm57uTuoy0AAAAA:_Q9Ymz5Ah90nrfdtdDeBqC1R0vYieFr4hzhx0UDVzJN6_E00w1F5zM5hEymU57Da5CwRRg_eMiP1))
- Model Evaluation:
  - Empirical value function approximation for resulting ITRs
  - Nested K-fold cross validation for parameter tuning, optimal model selection, and optimal model evaluation on held-out test set
- Compilation of Results:
  - Optimal decision rule, defined by patient coefficient of variation (%CV)
  - Characteristics of study participants, stratified by decision rule subgroup
  - Visualization of differences in CGM vs. BGM treatment effects, stratified by baseline coefficient of variation (%CV)
  - Training/validation set estimates of potential decision rules, along with held-out test set evaluation of optimal rule
  
## Data Availability

All data used in this analysis is publicly available through the Jaeb Center for Health Research (JCHR) and can be accessed in two ways:
1. Direct link: https://public.jaeb.org/dataset/564
2. Visit https://public.jaeb.org/datasets/diabetes and click `WISDMPublicDataset.zip' 

## Running the Workflow
1. Create the following empty folders from the root directory: `/data` and `/data/study data`
2. Download study data, unzip, and copy all text files into the `study data` folder
3. Run `C00data_processing.R`, which will output a clean dataset (`dat_clean_XAY_full.rds`) into `/data`
4. Run `C20decision_list.R` and `C21policy_tree.R`, which will output the following into `/data`
    - `decision_list_CV.rds`: Cross validation results for decision lists with different depth parameters
    - `policy_tree_CV.rds`: Cross validation results for policy trees with different depth parameters
    - `optimal_CV.rds`: Held-out test set cross validation result for optimal decision rule, CGM-only rule, and BGM-only rule
    - `dat_opt.rds`: Dataset with column of optimal predicted treatments added
5. Run `C30results.R`, which compiles all outputted files and calculates final results

## File Structure
- `scripts/C00data_processing.R`: Processes all data and outputs a clean dataset for analysis
  - Dependencies: Study data in `data/study data`
  - Output: `dat_clean_XAY_full.rds`
- `scripts/C20decision_list.R`: Runs functions for decision list fitting and value estimation
  - Dependencies: `F10CV_base.R`, `F20decision_list.R`
  - Output: `decision_list_CV.rds`
- `scripts/C21policy_tree.R`: Runs functions for policy tree fitting and value estimation
  - Dependencies: `F10CV_base.R`, `F21policy_tree.R`
  - Outputs: `policy_tree_CV.rds`, `optimal_CV.rds`, `dat_opt.rds`
  - Plots optimal decision rule (a policy tree with depth = 1)
- `scripts/C30results.R`:
  - Dependencies: Study data in `data/study data`
  - Outputs:
    - [figure] Difference in CGM vs. BGM treatment effect on outcome (by baseline %CV)
    - [table] Characteristics of study particiapnts (stratified by decision rule subgroup)
    - [table] All cross validation results for parameter tuning and optimal final model
- `scripts/F10CV_base.R`: Supporting file containing functions to create splits for inner and outer training/testing folds (nested K-fold cross validation)
- `scripts/F20decision_list.R`: Supporting file containing functions to fit decision list algorithm and calculate cross validated value estimates of resulting decision rule
- `scripts/F21policy_tree.R`: Supporting file containing functions to fit policy tree algorithm and calculate cross validated value estimates of resulting decision rule
