# zero-inflated negative binomial extended family for mgcv

# family - name of family character string
# link - name of link character string
# linkfun - the link function
# linkinv - the inverse link function
# mu.eta - d mu/d eta function (derivative of inverse link wrt eta)
# note: for standard links this information is supplemented using
#       function fix.family.link.extended.family with functions
#       gkg where k is 2,3 or 4 giving the kth derivative of the
#       link over the first derivative of the link to the power k.
#       for non standard links these functions must be supplied.
# dev.resids - function computing deviance residuals.
# Dd - function returning derivatives of deviance residuals w.r.t. mu and theta.
# aic - function computing twice -ve log likelihood for 2df to be added to.
# initialize - expression to be evaluated in gam.fit4 and initial.spg
#              to initialize mu or eta.
# preinitialize - optional function of y and family, returning a list with
#                 optional elements Theta - intitial Theta and y - modified
#                 y for use in fitting (see e.g. ocat and betar)
# postproc - function with arguments family, y, prior.weights, fitted,
#            linear.predictors, offset, intercept to call after fit to
#            compute (optionally) the label for the family, deviance and
#            null deviance. See ocat for simple example and betar or ziP
#            for complicated. Called in estimate.gam.
# ls - function to evaluated log saturated likelihood and derivatives w.r.t.
#      phi and theta for use in RE/ML optimization. If deviance used is just
#      -2 log-lik. can just return zeroes.
# validmu, valideta - functions used to test whether mu/eta are valid.
# n.theta - number of theta parameters.
# no.r.sq - optional TRUE/FALSE indicating whether r^2 can be computed
#           for family
# ini.theta - function for initializing theta.
# putTheta, getTheta - functions for storing and retriving theta values in
#                      function environment.
# rd - optional function for simulating response data from fitted model.
# residuals - optional function for computing residuals.
# predict - optional function for predicting from model, called by predict.gam.
# scale - < 0 to estimate. ignored if NULL.

#' Title
#'
#' @param y
#' @param g
#' @param eta
#' @param a
#' @param deriv
#'
#' @return
#' @export
zinbll <- function(y, g, eta, a, deriv = 0) {
  zind <- y == 0
  yp <- y[!zind]
  l <- et <- exp(eta)
  l[zind] <- -et[zind]
  l[!zind] <- logterm1(eta[!zind]) + yp * log(a) + yp * g[!zind] -
    yp * log(1 + exp(g[!zind])) - logterm2(g[!zind], a) + lgamma(y +
      1 / a) - lgamma(y + 1) - lgamma(1 / a)
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

  }
}

#' A helper function
#'
#' @param g
#' @param a
#' @param all
#'
#' @return A list containing:
#          b -- beta, t -- tau, lg -- log(1+a*eg)
#'         ind -- indices of y_i for which g_i is very small
#'         ii -- indices of y_i for which g_i is very large
btlg <- function(g, a, all = FALSE) {
  ind <- g < log(.Machine$double.eps)
  ii <- g > log(.Machine$double.xmax) / 2
  eg <- exp(g)

  t <- 1 / (1 - (1 + a * eg)^-(1 / a))
  lg <- log(1 + a * eg)
  lg[ind] <- 0
  lg[ii] <- g
  t[ind] <- 1 / eg
  t[ii] <- 1

  b <- eg / (1 + a * eg)
  b[ind] <- eg
  b[ii] <- 1 / a

  lg <- NULL
  if (all) {
    lg <- log(1 + a * eg)
    lg[ind] <- 0
    lg[ii] <- g
  }

  list(b = b, t = t, lg = lg, ind = ind, ii = ii)
}

#' Title
#'
#' @param x
#'
#' @return
#' @export
loglterm1 <- function(x) {
  ind <- x < log(.Machine$double.eps) / 3
  ex <- exp(x)
  exi <- ex[ind]

  l <- log(1 - exp(-ex))
  l[ind] <- log(exi - (exi^2) / 2 + (exi^3) / 6)
  ind <- x < -log(.Machine$double.xmax)
  l[ind] <- x[ind]

  l
}

#' Title
#'
#' @param x
#' @param a
#'
#' @return
#' @export
loglterm2 <- function(x, a) {
  ind <- x < -log(.Machine$double.xmax)
  ex <- exp(x)
  exi <- ex[ind]

  l <- log((1 + a * ex)^(1 / a) - 1)
  l[ind] <- x[ind]
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
  t <- d$t
  ind <- d$ind
  ii <- d$ii

  # first derivative
  l1 <- y - y * a * b - b * t
  l1[ind] <- y[ind] - 1
  l1[ii] <- -1 / a
  # second derivative
  l2 <- y * a^2 * b^2 - y * a * b + a * b^2 * t - t * b + b^2 * t^2 - b^2 * t
  l2[ind] <- 0
  l2[ii] <- 0

  l3 <- l4 <- NULL

  if (deriv > 1) {
    # third derivative
    l3 <- -2 * y * a^3 * b^3 + 3 * y * a^2 * b^2 - 2 * a^2 * b^3 * t -
      y * a * b + 3 * a * b^2 * t - 3 * a * b^3 * t^2 + 3 * a * b^3 * t -
      b * t + 3 * b^2 * t^2 - 3 * b^2 * t - 2 * b^3^t^3 + 3 * b^3 * t^2 -
      b^3 * t^2
    l3[ind] <- 0
    l3[ii] <- 0
  }
  if (deriv > 2) {
    # fourth derivative
    l4 <- 6 * y * a * 4 * b^4 - 12 * y * a^3 * b^3 + 6 * a^2 * b^4 * t +
      7 * y * a^2 * b^2 - 12 * a^2 * b^3 * t + 11 * a^2 * b^4 * t^2 -
      11 * a^2 * b^4 * t - y * a * b + 7 * a * b^2 * t - 18 * a * b^3 * t^2 +
      18 * a * b^3 * t + 12 * a * b^4 * t^3 - 18 * a * b^4 * t^2 +
      6 * a * b^4 * t - b * t + 7 * b^2 * t^2 - 7 * b^2 * t - 12 * b^3 * t^3 +
      18 * b^3 * t^2 - 6 * b^3 * t + 6 * b^4 * t^4 - 12 * b^4 * t^3 +
      7 * b^4 * t^2 - b^4 * t
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
lda <- function(g, y, a) {
  d <- btlg(g, a, all = TRUE)
  b <- d$b
  t <- d$t
  lg <- d$lg
  ind <- d$ind
  ii <- d$ii

  # TODO: revise this, interaction of t and lg as g -> -\infty is not zero
  # first derivative
  l1 <- y / a - y * b - b * t / a + t * lg / a^2 -
    1 / a^2 * (digamma(y + 1 / a) - digamma(1 / a))
  l1[ind] <- y[ind] / a - 1 / a -
    1 / a^2 * (digamma(y + 1 / a) - digamma(1 / a))
  l1[ii] <- -1 / a^2 + g / a^2 -
    1 / a^2 * (digamma(y + 1 / a) - digamma(1 / a))

  # second derivative
  l2 <- y * b^2 + b^2 * t / a + 2 * b * t / 2 - 2 * t * lg / a^2 -
    y / a^2 + t^2 / a^2 * (b - lg / a)^2 - t / a^2 * (b - lg / a)^2 -
    2 / a^3 * (digamma(1 / a) - digamma(y + 1 / a)) -
    1 / a^4 * (psigamma(1 / a, 2) - psigamma(y + 1 / a))
  l2[ind] <- 2 / a^2 - y / a^2 -
    2 / a^3 * (digamma(1 / a) - digamma(y + 1 / a)) -
    1 / a^4 * (psigamma(1 / a, 2) - psigamma(y + 1 / a))
  l2[ii] <- (3 - 2 * g) / a^3 -
    2 / a^3 * (digamma(1 / a) - digamma(y + 1 / a)) -
    1 / a^4 * (psigamma(1 / a, 2) - psigamma(y + 1 / a))

  list(l1 = l1, l2 = l2)
}
