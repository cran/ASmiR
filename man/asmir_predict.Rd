\name{asmir_predict}
\alias{asmir_predict}
\title{Predict Abiotic Stress Association of miRNA Sequences}
\description{Predicts whether input plant miRNA sequences are associated with cold, drought, heat, or salt stress using the bundled ASmiR models. The function generates PseKNC features and applies the relevant stress-specific model.}
\usage{
asmir_predict(
  x,
  stress = c("cold", "drought", "heat", "salt", "all"),
  threshold = 0.5,
  model_file = NULL,
  feature_file = NULL
)
}
\arguments{
  \item{x}{Input supported by \code{\link{asmir_read_fasta}}.}
  \item{stress}{Stress model to use: \code{cold}, \code{drought}, \code{heat}, \code{salt}, or \code{all}.}
  \item{threshold}{Numeric probability threshold between 0 and 1 for assigning the stress-associated class.}
  \item{model_file}{Optional path to an RDS model file, used when predicting a single stress.}
  \item{feature_file}{Optional path to a feature-ranking file, used when predicting a single stress.}
}
\value{A data frame containing sequence names, stress labels, predicted stress probabilities, confidence-like probabilities, and predicted class labels.}
