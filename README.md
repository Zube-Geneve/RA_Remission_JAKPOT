# Machine Learning to Predict Remission in Rheumatoid Arthritis (JAK-pot)

This repository contains R code used to develop and evaluate machine learning models predicting clinical remission between 6 and 24 months in rheumatoid arthritis, as described in:

**Salis Z. et al.**  
*Machine Learning to Predict Remission Between Six and 24 Months in Rheumatoid Arthritis: Insights from the JAK-pot Collaboration*  
**Arthritis & Rheumatology**, 2026.

---

## Data availability and governance

The JAK-pot collaboration is based on data contributions from multiple independent rheumatology registries.

- Data are hosted and analysed centrally in Geneva.
- Data remain the property of the individual contributing registries.
- Data cannot be shared through this repository.
- Any external data request must be addressed directly to each contributing registry and requires approval from all data owners.

This repository therefore provides **code only**, not data.

---

## Repository contents

- `32_ML_for_JAKPOT_CDAI_R1_v23.2_PUBLIC_full_or_top10.R`  
  Main analysis pipeline:
  - data preprocessing
  - internal and external validation
  - XGBoost model training
  - performance evaluation
  - SHAP-based interpretability (optional)

No patient-level data, derived datasets, or trained model artifacts are included.

---

## Model variants

The pipeline supports two model specifications, controlled by the environment variable `JAKPOT_MODEL_TYPE`:

### 1) Full model (default)

Uses the complete set of predictors described in the manuscript.

macOS/Linux:
```bash
export JAKPOT_MODEL_TYPE="full"
```

Windows (PowerShell):
```powershell
$env:JAKPOT_MODEL_TYPE="full"
```

If not specified, `"full"` is used by default.

### 2) Simplified 10-variable model

The reduced model uses the following 10 baseline variables:

- `PGA`
- `TJC28`
- `Prev_btsDMARD`
- `HAQ0`
- `PhGA`
- `Age`
- `CDAI0`
- `CRP`
- `Pain`
- `AntiCCPevernever`

These correspond to the top-ranked predictors identified in the primary analysis.

macOS/Linux:
```bash
export JAKPOT_MODEL_TYPE="top10"
```

Windows (PowerShell):
```powershell
$env:JAKPOT_MODEL_TYPE="top10"
```

The same preprocessing, imputation strategy, splitting, tuning framework, and validation procedures are used. Only the feature set differs.

---

## Requirements

R version **≥ 4.2** recommended.

Core packages used include:
- `tidyverse`, `data.table`
- `caret`
- `xgboost`
- `mice`
- `pROC`
- `SHAPforxgboost`
- `here`

Install missing packages manually if needed:

```r
install.packages(c(
  "tidyverse", "data.table", "caret", "xgboost",
  "mice", "pROC", "SHAPforxgboost", "here"
))
```

---

## How to run

Because data cannot be distributed, you must have authorized access to the JAK-pot dataset.

### 1) Set required environment variable

You must define the location of the main `.RData` file:

Windows (PowerShell):
```powershell
$env:JAKPOT_RDATA_FILE="D:\secure_path\clean_data.RData"
```

macOS/Linux:
```bash
export JAKPOT_RDATA_FILE="/secure_path/clean_data.RData"
```

Optional environment variables:
- `JAKPOT_HAQ_SCQM_FILE` (optional corrected HAQ file)
- `JAKPOT_OUT_DIR` (where logs/figures would be written)
- `JAKPOT_MODEL_DIR` (where model binaries would be written if enabled)
- `JAKPOT_MODEL_TYPE` (`"full"` or `"top10"`)

**Note:** If `JAKPOT_OUT_DIR` / `JAKPOT_MODEL_DIR` are not set, the script writes to temporary directories (`tempdir()`), to reduce the risk of accidentally committing outputs to Git.

### 2) Run the script

```r
source("32_ML_for_JAKPOT_CDAI_R1_v23.2_PUBLIC_full_or_top10.R")
```

---

## Expected objects in the `.RData` file

The script assumes the `.RData` file contains:

- `jakpot.data_eff`
- `jakpot.data_effbaseline`
- `jakpot.data_ae`

These correspond to the cleaned and harmonized datasets used for model development.

---

## Hyperparameter tuning

The pipeline includes optional hyperparameter tuning for the XGBoost model.

Tuning is controlled inside the script by:

```r
RUN_TUNING <- TRUE
```

- If `RUN_TUNING = TRUE`, cross-validated grid search is performed. Runtime can increase substantially depending on hardware.
- If `RUN_TUNING = FALSE`, pre-specified hyperparameters are used. Results may differ slightly from the tuned run.

For quick test runs/debugging, temporarily set:
```r
RUN_TUNING <- FALSE
```

---

## Safety mechanisms (privacy)

The public script includes safety flags to prevent accidental export of sensitive material:

```r
EXPORT_PATIENT_LEVEL   <- FALSE
EXPORT_DERIVED_OUTPUTS <- FALSE
SAVE_MODEL_ARTIFACTS   <- FALSE
```

These should remain `FALSE` for public use.

---

## Reproducibility notes

Results can depend on:
- exact data version and harmonization
- random seeds
- package versions / computing environment

For strict reproducibility, consider using `renv` to lock package versions.

---

## Contact

**Zubeyir Salis**  
University of Montpellier, INSERM U1046 (PhyMedExp)  
Email: zubeyir.salis@inserm.fr
