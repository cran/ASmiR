# ASmiR R package

R interface to the supplied ASmiR trained models for classifying plant miRNA sequences as associated or not associated with **cold, drought, heat, or salt stress**.

> **Validation note:** The package bundles the supplied ASmiR model objects and feature rankings and wraps the original PseKNC workflow. Validate predictions against suitable reference results before using them for scientific conclusions.

## Requirements

- R >= 4.1
- `Biostrings` (Bioconductor), `e1071` (CRAN)
- Python 3 available as `python` on PATH. You can set `ASMiR_PYTHON` to the full path of your Python executable.

```r
install.packages(c("remotes", "e1071"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Biostrings")
remotes::install_local("C:/path/to/ASmiR", dependencies = TRUE)
```

## Predict with the bundled example

```r
library(ASmiR)
fasta <- system.file("extdata", "examples", "gene_mirna.fasta", package = "ASmiR")
asmir_validate_sequences(fasta)
features <- asmir_features(fasta)
dim(features)  # expected 10 x 1379 for bundled example

cold <- asmir_predict(fasta, stress = "cold")
print(cold)

all_stresses <- asmir_predict(fasta, stress = "all")
write.csv(all_stresses, "ASmiR_predictions.csv", row.names = FALSE)
```

## Workflow fidelity

For each stress, the package generates five PseKNC feature blocks with k=1..5 and parameters lambda=3 and w=0.2, combines them into 1,379 columns, selects the first N indices in the corresponding SVM-RFE file, and calls the matching bundled SVM model. N is cold=380, drought=272, heat=174, salt=340. RNA U is converted to DNA T only for the original PseKNC feature-generation step. Input validation accepts A/U/G/C only.

## Output

`Name`, `Stress`, `p_stress`, `Probability`, and `Label`. `p_stress` is the raw probability for SVM class `Y`; `Probability` reproduces the original workflow's confidence-like score `max(p_stress, 1-p_stress)`. Label thresholds are 0.5 by default.


## Version 0.1.1 note

The SVM prediction step calls e1071's `predict.svm` method directly via
`getFromNamespace()` to avoid the S3 dispatch error encountered with the bundled
serialized SVM models. Restart R before reinstalling this version. Full
prediction parity with the original ASmiR workflow should still be verified
against known reference predictions.
