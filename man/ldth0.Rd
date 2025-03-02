% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hurdleNB_ll.R
\name{ldth0}
\alias{ldth0}
\title{Log-likelihood derivatives w.r.t. θ₀.}
\usage{
ldth0(g, y, th0, v)
}
\arguments{
\item{g}{\itemize{
\item 𝛄, a numeric vector,
}}

\item{y}{\itemize{
\item 𝐲, a numeric vector,
}}

\item{th0}{\itemize{
\item θ₀, a numeric,
}}

\item{v}{\itemize{
\item v, a list containing 𝛋, 𝛕 and lg.
}}
}
\value{
A list of the first and second derivatives of the
log-likelihood w.r.t. θ₀.
}
\description{
Log-likelihood derivatives w.r.t. θ₀.
}
