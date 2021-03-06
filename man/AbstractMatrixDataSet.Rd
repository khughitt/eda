% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AbstractMatrixDataSet.R
\name{AbstractMatrixDataSet}
\alias{AbstractMatrixDataSet}
\title{An R6 class representing collection of related datasets}
\description{
An R6 class representing collection of related datasets
}
\section{Usage}{
\preformatted{# TODO
}
}

\section{Arguments}{

\itemize{
\item \code{dataset}: A list of datasets (matrices, data frames, etc.), each of
which shared some column / row identifiers with the first entry in
the list.
}
}

\section{Fields}{

\itemize{
\item \code{dat}: Primary dataset
\item \code{row_data}: List of additional data keyed on row identifiers
\item \code{col_data}: List of additional data keyed on column identifiers
}
}

\section{Methods}{

\itemize{
\item \code{cross_cor(key1=1, key2=2, meas='pearson')}: Computes cross-dataset
correlation matrix between rows in two specified datasets.
\item \code{print()}: Prints an overview of the object instance.
}
}

\section{Examples}{
\preformatted{TODO
}
}

