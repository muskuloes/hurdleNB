% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hurdleNB_ll.R
\name{ldth0}
\alias{ldth0}
\title{Log-likelihood derivatives w.r.t. \eqn{\vartheta_0}.}
\usage{
ldth0(g, y, th0, v)
}
\arguments{
\item{g}{\eqn{\gamma}, a numeric vector,}

\item{y}{\eqn{y}, a numeric vector,}

\item{th0}{\eqn{\vartheta_0}, a numeric,}

\item{v}{\code{v}, a list containing \eqn{\kappa}, \eqn{\tau},
\code{lg}, \code{ind} and \code{ii}.}
}
\value{
A list of the first and second derivatives of the
log-likelihood w.r.t. \eqn{\vartheta_0}.
}
\description{
Log-likelihood derivatives w.r.t. \eqn{\vartheta_0}.
}
