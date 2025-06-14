% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scsajr.R
\name{cast_xy_table}
\alias{cast_xy_table}
\title{Construct a contingency table from paired x and y vectors}
\usage{
cast_xy_table(x, y, i)
}
\arguments{
\item{x}{A vector (atomic) of length n, representing row labels.}

\item{y}{A vector of length n, representing column labels.}

\item{i}{A vector of length n, where \code{i[j]} is the value to place at row = \code{x[j]}, column = \code{y[j]}.}
}
\value{
A matrix of dimension \verb{length(unique(x)) × length(unique(y))}, with row names set to sorted unique values of \code{x}
and column names set to sorted unique values of \code{y}.
For each position \verb{(r, c)}, the entry is the value from \code{i[j]} where \code{x[j] == r} and \code{y[j] == c}. If no such \code{j} exists, the entry is \code{NA}.
}
\description{
Given two vectors \code{x} and \code{y} of equal length, and a corresponding vector \code{i} of values,
this function creates a matrix whose rows correspond to each unique value in \code{x},
columns correspond to each unique value in \code{y},
and entries are taken from \code{i} at matching \verb{(x, y)} pairs.
}
\examples{
\dontrun{
x <- c("A", "B", "A", "C")
y <- c("X", "Y", "X", "Z")
i <- c(1, 2, 3, 4)
# unique(x) = A, B, C; unique(y) = X, Y, Z
# (A, X) appears twice with values 1 and 3; the latter overwrites
result <- cast_xy_table(x, y, i)
#      X  Y  Z
#  A   3 NA NA
#  B  NA  2 NA
#  C  NA NA  4
}

}
