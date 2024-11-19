#' log(1-exp(-exp(x))).
#'
#' @param x - A numeric vector.
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

#' log(1 - (1 + αeˣ)^(-1/ α)).
#'
#' @param x   - A numeric vector.
#' @param th0 - θ₀, a numeric.
#'
#' @return Carfully compute log(1 - (1 + exp(θ₀)exp(x))^(-1/(exp(θ₀)))).
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

#' A helper function returning stable 𝛃, 𝞃, and log(1+αeᵞ).
#'
#' @param g    - 𝛄, a numeric vector.
#' @param a    - α, a numeric.
#' @param what - A character vector specifying what to return.
#'
#' @return A list containing:
#          b -- 𝛃, tau -- 𝛕, lg -- log(1+αeᵞ).
#'         ind -- indices of yᵢ for which γᵢ is very small.
#'         ii -- indices of yᵢ for which γᵢ is very large.
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

#' Log-likelihood derivates w.r.t. 𝛈.
#'
#' @param eta   - 𝛈, a numeric vector.
#' @param deriv - <= 1 - first and second derivatives.
#'                == 2 - first, second and third derivatives.
#'                >= 3 - first, second, third and fourth derivatives.
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

#' Log-likelihood derivates w.r.t. 𝛄.
#'
#' @param g     - 𝛄, a numeric vector.
#' @param y     - 𝐲, a numeric vector.
#' @param a     - α, a numeric.
#' @param deriv - <= 1 - first and second derivatives.
#'                == 2 - first, second and third derivatives.
#'                >= 3 - first, second, third and fourth derivatives.
#'
#' @return A list of derivatives of the log-likelihood w.r.t. g.
ldg <- function(g, y, a, deriv = 4) {
  d <- btlg(g, a, c("b", "tau"))
  b <- d$b
  tau <- d$tau
  ind <- d$ind
  ii <- d$ii

  # first derivative
  l1 <- y - y * a * b - b * tau
  l1[ind] <- y[ind] - 1
  l1[ii] <- -1 / a

  # second derivatauive
  l2 <- y * (a^2) * (b^2) - y * a * b + a * (b^2) * tau - tau * b +
    (b^2) * (tau^2) - (b^2) * tau
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

#' Log-likelihood derivatives w.r.t. θ₀.
#'
#' @param g   - 𝛄, a numeric vector
#' @param y   - 𝐲, a numeric vector
#' @param th0 - θ₀, a numeric
#'
#' @return A list of the first and second derivatives of the
#          log-likelihood w.r.t. θ₀.
ldth0 <- function(g, y, th0) {
  a <- exp(th0)
  d <- btlg(g, a, c("b", "tau"))
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
#' p is probability of potential presence. mu is the NB mean.
#'
#' @param y     - 𝐲, a numeric vector
#' @param g     - 𝛄, a numeric vector
#' @param eta   - 𝛈, a numeric vector
#' @param th0   - θ₀, a numeric
#' @param deriv - 0 - eval.
#'                1 - gradient and Hessian.
#'                2 - third derivatives.
#'                4 - fourth derivatives.
#'
#' @return ZINB log-likelihood and its derivatives.
#' @export
zinbll <- function(y, g, eta, th0, deriv = 0) {
  a <- exp(th0)
  zind <- y == 0
  l <- et <- exp(eta)
  yp <- y[!zind]
  l[zind] <- -et[zind]
  d <- btlg(g, a, what = c("b", "tau", "lg"))
  b <- d$b
  tau <- d$tau
  lg <- d$lg
  l[!zind] <- l1ee(eta[!zind]) + yp * log(a) + yp * g[!zind] -
    yp * lg[!zind] - l11aea(g[!zind], th0) +
    lgamma(y + 1 / a) - lgamma(y + 1) - lgamma(1 / a)
  q <- 1 - exp(-et)

  n <- length(y)

  l1 <- El2 <- l2 <- l3 <- l4 <- NULL

  # get first and second derivatives
  if (deriv > 0) {
    # order l_g, l_e
    l1 <- matrix(0, n, 2)
    le <- lde(eta, deriv)
    lg <- ldg(g, y, a, deriv)
    l1[!zind, 1] <- lg$l1[!zind] # l_gamma, y>0
    l1[zind, 2] <- l[zind] # l_eta, y==0
    l1[!zind, 2] <- le$l1[!zind] # l_eta, y>0

    El2 <- l2 <- matrix(0, n, 3)

    # order l_gg, l_eg, l_ee
    l2[!zind, 1] <- lg$l2[!zind] # l_gg, y>0
    l2[zind, 3] <- l[zind] # l_ee, y==0
    l2[!zind, 3] <- le$l2[!zind] # l_ee, y>0
    El2[, 1] <- q * (q * tau * exp(g) * ((a^2) * b^2 - a * b) +
      a * (b^2) * tau - tau * b + (b^2) * (tau^2) - (b^2)(tau)) # E[l_gg]
    El2[, 3] <- -(1 - q) * et + q * le$l2 # E[l_ee]
  }
  # get third derivates
  if (deriv > 1) {
    # order l_ggg, l_gge, l_gee, l_eee
    l3 <- matrix(0, n, 4)
    l3[!zind, 1] <- lg$l3[!zind]
    l3[!zind, 4] <- le$l3[!zind]
    l3[zind, 4] <- l[zind]
  }
  # get fourth derivates
  if (deriv > 3) {
    # order l_gggg, l_ggge, l_ggee, l_geee, l_eeee
    l4 <- matrix(0, n, 5)
    l4[!zind, 1] <- lg$l4[!zind]
    l4[!zind, 5] <- le$l4[!zind]
    l4[zind, 5] <- l[zind]
  }

  list(l = l, l1 = l1, l2 = l2, l3 = l3, l4 = l4, El2 = El2)
}
