\name{asmir_model_files}
\alias{asmir_model_files}
\title{Locate ASmiR Model Files}
\description{Returns paths to the bundled stress-specific support vector machine models and feature-ranking files used by ASmiR.}
\usage{
asmir_model_files(stress = c("cold", "drought", "heat", "salt", "all"))
}
\arguments{
  \item{stress}{Stress model to locate: \code{cold}, \code{drought}, \code{heat}, \code{salt}, or \code{all}.}
}
\value{A named character vector containing paths to the model and feature-ranking files for the requested stress model or models.}
