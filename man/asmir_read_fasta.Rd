\name{asmir_read_fasta}
\alias{asmir_read_fasta}
\title{Read ASmiR FASTA Input}
\description{Converts a FASTA file, compatible data frame, or character vector of sequences into the standardized sequence data frame used by ASmiR.}
\usage{
asmir_read_fasta(x)
}
\arguments{
  \item{x}{A FASTA file path, a data frame containing sequence identifier and sequence columns, or a character vector of RNA sequences.}
}
\value{A data frame with \code{sequence_id} and \code{sequence} columns.}
