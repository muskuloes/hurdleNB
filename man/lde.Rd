% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hurdleNB_ll.R
\name{lde}
\alias{lde}
\title{Log-likelihood derivatives w.r.t. 𝛈.}
\usage{
lde(eta, level = 2)
}
\arguments{
\item{eta}{\itemize{
\item 𝛈, a numeric vector,
}}

\item{level}{\itemize{
\item == 0 - first and second derivatives,
>  0 - derivatives needed for quasi-Newton,
>  1 - derivatives need for full Newton.
}}
}
\value{
A list of derivatives of the log-likelihood w.r.t. 𝛈 (eta).
}
\description{
Log-likelihood derivatives w.r.t. 𝛈.
}
