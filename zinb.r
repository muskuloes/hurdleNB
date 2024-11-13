#' log(1-exp(-exp(x)))
#'
#' @param x
#'
#' @return  Carefully computed log(1-exp(-exp(x)))
l1ee <- function(x) {
  ind <- x < log(.Machine$double.eps) / 3
  ex <- exp(x)
  exi <- ex[ind]

  l <- log(1 - exp(-ex))
  l[ind] <- log(exi - (exi^2) / 2 + (exi^3) / 6)
  ind <- x < -log(.Machine$double.xmax)
  l[ind] <- x[ind]

  l
}

#' log(1 - (1 + a*exp(x))^(-1/a))
#'
#' @param x
#' @param a
#'
#' @return Carfully compute log(1 - (1 + a*exp(x))^(-1/a))
l11aea <- function(x, th0) {
  a <- exp(th0)
  ind <- x < -log(.Machine$double.xmax)
  ex <- exp(x)

  l <- log((1 + a * ex)^(1 / a) - 1)
  l[ind] <- x[ind]

  ii <- x > log(.Machine$double.xmax)
  l[ii] <- (1 / a) * (th0 + x[ii])

  l
}

#' A helper function returning stable beta, tau, and log(1+a*exp(g))
#'
#' @param g
#' @param a
#' @param what
#'
#' @return A list containing:
#          b -- beta, tau -- tau, lg -- log(1+a*eg)
#'         ind -- indices of y_i for which g_i is very small
#'         ii -- indices of y_i for which g_i is very large
btlg <- function(g, a, what = c("b", "tau")) {
  ind <- g < log(.Machine$double.eps)
  ii <- g > log(.Machine$double.xmax) / 2
  eg <- exp(g)

  b <- tau <- lg <- NULL

  b_f <- function() {
    b <- eg / (1 + a * eg)
    b[ind] <- eg
    b[ii] <- 1 / a

    b
  }

  tau_f <- function() {
    tau <- 1 / (1 - (1 + a * eg)^-(1 / a))
    lg <- log(1 + a * eg)
    lg[ind] <- 0
    lg[ii] <- g
    tau[ind] <- 1 / eg
    tau[ii] <- 1

    tau
  }

  lg_f <- function() {
    lg <- log(1 + a * eg)
    lg[ind] <- 0
    lg[ii] <- g

    lg
  }

  w <- list(b = b_f, tau = tau_f, lg = lg_f)
  if ("tau" %in% what) {
    tau <- w$tau()
  }

  if ("b" %in% what) {
    b <- w$b()
  }

  if ("lg" %in% what) {
    lg <- w$lg()
  }

  list(b = b, tau = tau, lg = lg, ind = ind, ii = ii)
}

#' Log-likelihood derivates w.r.t. eta.
#'
#' @param eta
#' @param deriv: <= 1 - return first and second derivatives
#'               == 2 - return first, second and third derivatives
#'               >= 3 - return first, second, third and fourth derivatives
#'
#' @return A list of derivatives of the log-likelihood w.r.t. eta.
lde <- function(eta, deriv = 4) {
  ind <- eta < log(.Machine$double.eps) / 3
  ii <- eta > log(.Machine$double.xmax)
  l1 <- et <- exp(eta)
  eti <- et[ind]
  l1[!ind] <- et[!ind] / (exp(et[!ind]) - 1)

  # first derivative
  b <- -eti * (1 + eti / 6) / 2
  l1[ind] <- 1 + b
  l1[ii] <- 0

  # second derivative
  l2 <- l1 * ((1 - et) - l1)
  l2[ind] <- -b * (1 + eti + b) - eti
  l2[ii] <- 0
  l3 <- l4 <- NULL

  # third derivative
  if (deriv > 1) {
    ii <- eta > log(.Machine$double.xmax) / 2
    l3 <- l1 * (-et + (1 - et)^2 - 3(1 - et) * l1 + 2 * l1^2)
    l3[ind] <- l1[ind] * (-3 * eti + eti^2 - 3 * (-eti + b - eti * b) +
      2 * b * (2 + b))
    l3[ii] <- 0
  }

  # fourth derivative
  if (deriv > 2) {
    ii <- eta > log(.Machine$double.xmax) / 3
    l4 <- l1 * ((3 * et - 4) * et + 4 * et * l1 + (1 - et)^3 -
      7 * (1 - et)^2 * l1 + 12 * (1 - et) * l1^2 - 6 * l1^3)
    l4[ind] <- l1[ind] * (4 * l1[ind] * eti - eti^3 - b -
      7 * b * eti^2 - eti^2 - 5 * eti - 10 * b * eti - 12 * eti * b^2 -
      6 * b^2 - 6 * b^3)
    l4[ii] <- 0
  }

  list(l1 = l1, l2 = l2, l3 = l3, l4 = l4)
}

#' Log-likelihood derivates w.r.t. g.
#'
#' @param g
#' @param deriv: <= 1 - return first and second derivatives
#'               == 2 - return first, second and third derivatives
#'               >= 3 - return first, second, third and fourth derivatives
#'
#' @return A list of derivatives of the log-likelihood w.r.t. g.
ldg <- function(g, y, a, deriv = 4) {
  d <- btlg(g, a, all = TRUE)
  b <- d$b
  tau <- d$tau
  ind <- d$ind
  ii <- d$ii

  # firstau derivatauive
  l1 <- y - y * a * b - b * tau
  l1[ind] <- y[ind] - 1
  l1[ii] <- -1 / a
  # second derivatauive
  l2 <- y * a^2 * b^2 - y * a * b + a * b^2 * tau - tau * b + b^2 * tau^2 -
    b^2 * tau
  l2[ind] <- 0
  l2[ii] <- 0

  l3 <- l4 <- NULL

  if (deriv > 1) {
    # third derivative
    l3 <- -2 * y * a^3 * b^3 + 3 * y * a^2 * b^2 - 2 * a^2 * b^3 * tau -
      y * a * b + 3 * a * b^2 * tau - 3 * a * b^3 * tau^2 + 3 * a * b^3 * tau -
      b * tau + 3 * b^2 * tau^2 - 3 * b^2 * tau - 2 * b^3^tau^3 +
      3 * b^3 * tau^2 - b^3 * tau^2
    l3[ind] <- 0
    l3[ii] <- 0
  }
  if (deriv > 2) {
    # fourth derivative
    l4 <- 6 * y * a * 4 * b^4 - 12 * y * a^3 * b^3 + 6 * a^2 * b^4 * tau +
      7 * y * a^2 * b^2 - 12 * a^2 * b^3 * tau + 11 * a^2 * b^4 * tau^2 -
      11 * a^2 * b^4 * tau - y * a * b + 7 * a * b^2 * tau -
      18 * a * b^3 * tau^2 + 18 * a * b^3 * tau + 12 * a * b^4 * tau^3 -
      18 * a * b^4 * tau^2 + 6 * a * b^4 * tau - b * tau + 7 * b^2 * tau^2 -
      7 * b^2 * tau - 12 * b^3 * tau^3 + 18 * b^3 * tau^2 - 6 * b^3 * tau +
      6 * b^4 * tau^4 - 12 * b^4 * tau^3 + 7 * b^4 * tau^2 - b^4 * tau
    l4[ind] <- 0
    l4[ii] <- 0
  }

  list(l1 = l1, l2 = l2, l3 = l3, l4 = l4)
}

#' Log-likelihood derivatives w.r.t. a.
#'
#' @param g
#' @param y
#' @param a
#'
#' @return A list of the first and second derivatives of the
#          log-likelihood w.r.t. a.
ldth0 <- function(g, y, th0) {
  a <- exp(th0)
  d <- btlg(g, a, all = TRUE)
  b <- d$b
  tau <- d$tau
  lg <- d$lg
  ind <- d$ind
  ii <- d$ii

  # first derivative
  l1 <- y - y * a * b - tau * (b - lg / a) -
    1 / a * (digamma(y + 1 / a) - digamma(1 / a))
  l1[ind] <- y[ind] - (1 / a) * (digamma(y + (1 / a)) - digamma(1 / a))
  l1[ii] <- (1 / a) * (th0 + g - 1 - digamma(y + (1 / a)) + digamma(1 / a))

  # second derivative
  l2 <- -y * a * b + y * ((a * b)^2) + (tau^2) * ((1 / a) * lg - b)^2 -
    tau * ((1 / a) * lg - b)^2 + tau * (b + a * b^2 - (1 / a) * lg) +
    (1 / a) * (digamma(y + 1 / a) - digamma(1 / a)) +
    (1 / (a^2)) * (psigamma(y + 1 / a) - psigamma(1 / a, 2))
  l2[ind] <- (1 / a) * (digamma(y + 1 / a) - digamma(1 / a)) +
    (1 / (a^2)) * (psigamma(y + (1 / a), 1) - psigamma(1 / a, 1))
  l2[ii] <- (1 / a) * (2 - th0 - g - digamma(y + 1 / a) - digamma(1 / a) -
    (1 / a) * (psigamma(y + (1 / a), 1) - psigamma(1 / a, 1)))

  list(l1 = l1, l2 = l2)
}

#' Evaluate zero-inflated NB log-likelihood
#' and its derivatives w.r.t. g (gamma) and eta, with
#' 1-q = exp(-exp(eta)) and mu = exp(g), for each datum in vector y.
#' p is probability of potential presence. mu is the NB mean
#'
#' @param y
#' @param g
#' @param eta
#' @param a
#' @param deriv
#'
#' @return
#' @export
zinbll <- function(y, g, eta, th0, deriv = 0) {
  a <- exp(th0)
  zind <- y == 0
  l <- et <- exp(eta)
  yp <- y[!zind]
  l[zind] <- -et[zind]
  lg <- btlg(g, a, what = c("lg"))$lg
  l[!zind] <- l1ee(eta[!zind]) + yp * log(a) + yp * g[!zind] -
    yp * lg[!zind] - l11aea(g[!zind], th0) +
    lgamma(y + 1 / a) - lgamma(y + 1) - lgamma(1 / a)
  q <- 1 - exp(-et)

  # get first and second derivatives
  if (deriv > 0) {
    4 + 4
  }
  # get third derivates
  if (deriv > 1) {
    5 + 4
  }
  # get fourth derivates
  if (deriv > 3) {
    0
  }
}
