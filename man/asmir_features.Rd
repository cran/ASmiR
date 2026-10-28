\name{asmir_features}
\alias{asmir_features}
\title{Generate ASmiR PseKNC Features}
\description{Generates the PseKNC feature representation used by ASmiR for validated plant miRNA sequences. This function invokes the bundled Python feature-generation script and requires Python 3.}
\usage{
asmir_features(x)
}
\arguments{
  \item{x}{Input supported by \code{\link{asmir_read_fasta}}.}
}
\value{A numeric matrix with 1,379 feature columns and one row per input sequence.}
