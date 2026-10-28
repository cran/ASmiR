\name{asmir_validate_sequences}
\alias{asmir_validate_sequences}
\title{Validate ASmiR RNA Sequences}
\description{Checks whether input sequences are non-empty and contain only the RNA nucleotide symbols A, U, G, and C.}
\usage{
asmir_validate_sequences(x)
}
\arguments{
  \item{x}{Input supported by \code{\link{asmir_read_fasta}}.}
}
\value{A data frame containing sequence identifiers, sequence lengths, validity flags, and validation messages.}
