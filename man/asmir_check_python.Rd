\name{asmir_check_python}
\alias{asmir_check_python}
\title{Check Python Runtime}
\description{Checks whether a Python executable is available for the ASmiR PseKNC feature-generation workflow. The executable can be configured with the \code{ASMiR_PYTHON} environment variable.}
\usage{
asmir_check_python()
}
\value{Returns \code{TRUE} invisibly when Python can be invoked; otherwise raises an error.}
