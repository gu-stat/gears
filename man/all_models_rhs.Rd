% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/all_models_rhs.R
\name{all_models_rhs}
\alias{all_models_rhs}
\title{Generates all possible combinations of right-hand side variables.}
\usage{
all_models_rhs(
  y.name = NULL,
  y.max.lags = NULL,
  x.names = NULL,
  x.max.lags = NULL,
  x.fixed.names = NULL,
  x.fixed.lags = NULL,
  x.interaction.names = NULL,
  x.interaction.lags = NULL,
  use.intercept = c("both", "without", "with")
)
}
\arguments{
\item{y.name}{The name of the Y (left-hand side) variable. If NULL (default),
the function creates a temporary name.}

\item{y.max.lags}{A numeric value that gives the maximum number of lags of
the Y (left-hand side) variable (see Details). Can be NULL (default) if
the past values of the Y variable are not included in the model.}

\item{x.names}{List with names of the X (right-hand side) variables that
have a maximum number of lags (see Details). Can be NULL (default) if
univariate model or if your model does not have variables of this type.}

\item{x.max.lags}{List of numeric values that give the maximum number of lags
of the X (right-hand side) variables. Can be NULL (default).}

\item{x.fixed.names}{List with names of the X (right-hand side) variables
that have a fixed number of lags (see Details).
Can be NULL (default) if univariate model or if your model does not have
variables of this type.}

\item{x.fixed.lags}{List of numeric values that give the fixed number of lags
of the variables in x.fixed.names. Can be NULL (default).}

\item{x.interaction.names}{List of character vectors with names of the
variables to be included as interaction terms (see Details).
Can be NULL (default) if your model does not have interaction terms.}

\item{x.interaction.lags}{List of numeric vectors with lags of the
variables to be included as interaction terms. List and numeric vectors
should have the same length as the ones in x.interaction.names.
Can be NULL (default) if your model does not have interaction terms.}

\item{use.intercept}{If \code{use.intercept == "both"} (default), the
function returns all possible model equations with and without intercept.
If \code{use.intercept == "without"}, then the function returns all
possible right-hand side equations without intercept. If
\code{use.intercept == "with"}, then only the right-hand side equations
with intercept are returned.}
}
\value{
A character vector containing all possible combinations of the input
    variables and lags that form the right-hand side of the model equations.
}
\description{
Takes in the names and lags of the variables, and generates all possible
combinations of these to form the right-hand side of the model equations.
}
\details{
If \code{y.max.lags} equals to a number, then all lags of
    \code{y.name} up to this number (and starting at 0) will be included
    in the list of variables to enter the right-hand side of the model
    equations. For example, if \code{y.max.lags = 2}, then
    \ifelse{html}{\out{Y<sub>t</sub>, Y<sub>t-1</sub>, Y<sub>t-2</sub>
    }}{\eqn{Y_t, Y_{t-1},Y_{t-2}}} will be included in the list of variables
    to enter the right-hand side of the equation.
}
\keyword{internal}
