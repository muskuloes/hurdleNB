#' Evaluate hurdle negative binomial log-likelihood
#' and its derivatives.
#' @description
#' Evaluate hurdle negative binomial (NB) log-likelihood
#' and its derivatives w.r.t. \eqn{\gamma} (g) and \eqn{\eta} (eta), with
#' \eqn{1 - q = e^{-e^{\eta}}} and \eqn{\mu = e^{\gamma}}, for each datum in vector \eqn{y}.
#' \eqn{q} is the probability of potential presence. \eqn{\mu} is the NB mean.
#'
#' @param y     - \eqn{y}, a numeric vector,
#' @param g     - \eqn{\gamma}, a numeric vector,
#' @param eta   - \eqn{\eta}, a numeric vector,
#' @param th0   - \eqn{\vartheta_0}, a numeric,
#' @param level   \itemize{
#'                 \item \eqn{== 0} - eval,
#'                 \item \eqn{> 0} - derivatives for estimating \eqn{\beta} and \eqn{\rho} using quasi-Newton,
#'                 \item \eqn{> 1} - derivatives for estimating \eqn{\beta} and \eqn{\rho} using full Newton.}
#'
#' @return hurdleNB log-likelihood and its derivatives.
#' @export
hurdleNB_ll <- function(y, g, eta, th0, level = 0) {
  a <- exp(th0)
  v <- ktlg(g, a, what = c("k", "lg", "tau"))
  k <- v$k
  tau <- v$tau
  lg <- v$lg
  zind <- y == 0
  l <- et <- exp(eta)
  yp <- y[!zind]
  l[zind] <- -et[zind]
  l[!zind] <- yp * g[!zind] + yp * log(a) -
    yp * lg[!zind] + l1ee(eta[!zind]) - l11aea(g[!zind], th0) -
    lgamma(yp + 1) - lgamma(1 / a) + lgamma(yp + 1 / a)
  q <- 1 - exp(-et)

  l1 <- El2 <- l2 <- l3 <- l4 <- NULL
  if (level >= 1) {
    n <- length(y)

    l_e <- lde(eta, level)
    l_g <- ldg(g, y, a, v, level)

    # order ∂ℓ/∂𝛄, ∂ℓ/∂𝛈, ∂ℓ/∂ϑ₀.
    l1 <- matrix(0, n, 3)
    l1[!zind, 1] <- l_g$l1[!zind]
    l1[zind, 2] <- l[zind]
    l1[!zind, 2] <- l_e$l1[!zind]
    l1[, 3] <- NaN

    # order ∂²ℓ/∂𝛄², ∂²ℓ/∂𝛈∂𝛄, ∂²ℓ/∂𝛈², ∂²ℓ/∂𝛄∂ϑ₀, ∂²ℓ/∂ϑ₀².
    l2 <- matrix(0, n, 5)
    l2[!zind, 1] <- l_g$l2[!zind]
    l2[zind, 3] <- l[zind]
    l2[!zind, 3] <- l_e$l2[!zind]
    l2[, 4] <- NaN
    l2[, 5] <- NaN

    # order 𝔼[∂²ℓ/∂𝛄²], 𝔼[∂²ℓ/∂𝛈∂𝛄], 𝔼[∂²ℓ/∂𝛈²].
    El2 <- matrix(0, n, 3)
    El2[, 1] <- q * (tau * exp(g) * ((a^2) * k^2 - a * k) + a * (k^2) * tau -
      tau * k + k^2 * tau^2 - k^2 * tau)
    El2[, 3] <- -(1 - q) * et + q * l_e$l2

    l_dgth0 <- ldgth0(g, y, th0, v, level)

    l1[zind, 3] <- 0
    l1[!zind, 3] <- l_dgth0$l1[!zind]
    l2[zind, 4] <- 0
    l2[!zind, 4] <- l_dgth0$l_gth0[!zind]

    # order ∂³ℓ/∂𝛄³, ∂³ℓ/∂𝛄²∂𝛈, ∂³ℓ/∂𝛄∂𝛈², ∂³ℓ/∂𝛈³, ∂³ℓ/∂𝛄²∂ϑ₀, ∂³ℓ/∂𝛄∂ϑ₀².
    l3 <- matrix(0, n, 6)
    l3[!zind, 1] <- l_g$l3[!zind]
    l3[!zind, 4] <- l_e$l3[!zind]
    l3[zind, 4] <- l[zind]
    l3[!zind, 5] <- l_dgth0$l_ggth0[!zind]
    l3[, 6] <- NaN
  }

  if (level >= 2) {
    l2[zind, 5] <- 0
    l2[!zind, 5] <- l_dgth0$l2[!zind]
    l3[zind, 6] <- 0
    l3[!zind, 6] <- l_dgth0$l_gth0th0[!zind]

    # order ∂⁴ℓ/∂𝛄⁴, ∂⁴ℓ/∂𝛄³∂𝛈, ∂⁴ℓ/∂𝛄²∂𝛈², ∂⁴ℓ/∂𝛄∂𝛈³, ∂⁴ℓ/∂𝛈⁴, ∂⁴ℓ/∂𝛄³∂ϑ₀,
    # ∂⁴ℓ/∂𝛄²∂ϑ₀².
    l4 <- matrix(0, n, 7)
    l4[!zind, 1] <- l_g$l4[!zind]
    l4[!zind, 5] <- l_e$l4[!zind]
    l4[zind, 5] <- l[zind]
    l4[!zind, 6] <- l_dgth0$l_gggth0[!zind]
    l4[!zind, 7] <- l_dgth0$l_ggth0th0[!zind]
  }

  list(l = l, l1 = l1, l2 = l2, l3 = l3, l4 = l4, El2 = El2)
}

#' \eqn{\log(1 - e^{-e^{\eta}})}.
#'
#' @param eta - \eqn{\eta}, a numeric vector.
#'
#' @return  Carefully computed \eqn{\log(1 - e^{-e^{\eta}})}.
l1ee <- function(eta) {
  ind <- eta < log(.Machine$double.eps) / 3
  ex <- exp(eta)
  exi <- ex[ind]

  l <- log(1 - exp(-ex))
  l[ind] <- log(exi - (exi^2) / 2 + (exi^3) / 6)
  ind <- eta < -log(.Machine$double.xmax)
  l[ind] <- eta[ind]

  l
}

#'
#' \eqn{\log\left(\left(1 + \alpha e^{\gamma}\right)^{\frac{1}{\alpha}} - 1\right)}
#'
#' @param g   - \eqn{\gamma}, a numeric vector,
#' @param th0 - \eqn{\vartheta_0}, a numeric.
#'
#' @return Carefully computed
#' \eqn{\log\left(\left(1 + e^{\vartheta_0} e^{\gamma}\right)^{\frac{1}{e^{\vartheta_0}}} - 1\right)}.
l11aea <- function(g, th0) {
  a <- exp(th0)
  ind <- g < -log(.Machine$double.xmax)
  eg <- exp(g)

  l <- log((1 + a * eg)^(1 / a) - 1)
  l[ind] <- g[ind]

  ii <- g > log(.Machine$double.xmax)
  l[ii] <- (1 / a) * (th0 + g[ii])

  l
}

#' A helper function returning stable
#' \eqn{\kappa}, \eqn{\tau}, and \eqn{\log(1 + \alpha e^{\gamma})}.
#'
#' @param g    - \eqn{\gamma}, a numeric vector,
#' @param a    - \eqn{\alpha}, a numeric,
#' @param what - A character vector specifying what to return.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{k} – \eqn{\kappa},
#'   \item \code{tau} – \eqn{\tau},
#'   \item \code{lg} – \eqn{\log(1 + \alpha e^{\gamma})},
#'   \item \code{ind} – indices of \eqn{y_i} for which \eqn{\gamma_i} is very small,
#'   \item \code{ii} – indices of \eqn{y_i} for which \eqn{\gamma_i} is very large.
#' }
ktlg <- function(g, a, what = c("k", "tau")) {
  ind <- g < log(.Machine$double.eps)
  ii <- g > log(.Machine$double.xmax) / 2
  eg <- exp(g)

  k <- tau <- lg <- NULL

  k_f <- function() {
    k <- eg / (1 + a * eg)
    k[ind] <- eg[ind]
    k[ii] <- 1 / a

    k
  }

  tau_f <- function() {
    tau <- 1 / (1 - (1 + a * eg)^(-(1 / a)))
    tau[ind] <- 1 / eg[ind]
    tau[ii] <- 1

    tau
  }

  lg_f <- function() {
    lg <- log(1 + a * eg)
    lg[ind] <- 0
    lg[ii] <- g[ii]

    lg
  }

  w <- list(k = k_f, tau = tau_f, lg = lg_f)
  if ("tau" %in% what) {
    tau <- w$tau()
  }

  if ("k" %in% what) {
    k <- w$k()
  }

  if ("lg" %in% what) {
    lg <- w$lg()
  }

  list(k = k, tau = tau, lg = lg, ind = ind, ii = ii)
}

#' Log-likelihood derivatives w.r.t. \eqn{\eta}.
#'
#' @param eta   - \eqn{\eta}, a numeric vector,
#' @param level
#' \itemize{
#'   \item \eqn{==0} – first and second derivatives,
#'   \item \eqn{> 0} – derivatives needed for quasi-Newton,
#'   \item \eqn{> 1} – derivatives needed for full Newton.
#' }
#'
#' @return A list of derivatives of the log-likelihood w.r.t. \eqn{\eta} (eta).
lde <- function(eta, level = 2) {
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

  if (level > 0) {
    # third derivative
    ii <- eta > log(.Machine$double.xmax) / 2
    l3 <- l1 * (-et + (1 - et)^2 - 3 * (1 - et) * l1 + 2 * l1^2)
    l3[ind] <- l1[ind] * (-3 * eti + eti^2 - 3 * (-eti + b - eti * b) +
      2 * b * (2 + b))
    l3[ii] <- 0
  }

  if (level > 1) {
    # fourth derivative
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

#' Log-likelihood derivatives w.r.t. \eqn{\gamma}.
#'
#' @param g     - \eqn{\gamma}, a numeric vector,
#' @param y     - \eqn{y}, a numeric vector,
#' @param a     - \eqn{\alpha}, a numeric,
#' @param v     - \code{v}, a list containing \eqn{\kappa}, \eqn{\tau}, \code{ind} and \code{ii},
#' @param level
#' \itemize{
#'   \item \eqn{0} – first and second derivatives.
#'   \item \eqn{> 0} – derivatives needed for quasi-Newton.
#'   \item \eqn{> 1} – derivatives needed for full Newton.
#' }
#'
#' @return A list of derivatives of the log-likelihood w.r.t. \eqn{\gamma} (g).
ldg <- function(g, y, a, v, level = 2) {
  k <- v$k
  tau <- v$tau
  ind <- v$ind
  ii <- v$ii

  # first derivative
  l1 <- -a * k * y - k * tau + y
  l1[ind] <- y[ind] - 1
  l1[ii] <- -1 / a

  # second derivative
  l2 <- a^2 * k^2 * y + a * k^2 * tau - a * k * y + k^2 * tau^2 - k^2 * tau -
    k * tau
  l2[ind] <- 0
  l2[ii] <- 0

  l3 <- l4 <- NULL

  if (level > 0) {
    # third derivative
    l3 <- -2 * a^3 * k^3 * y - 2 * a^2 * k^3 * tau + 3 * a^2 * k^2 * y -
      3 * a * k^3 * tau^2 + 3 * a * k^3 * tau + 3 * a * k^2 * tau - a * k * y -
      2 * k^3 * tau^3 + 3 * k^3 * tau^2 - k^3 * tau + 3 * k^2 * tau^2 -
      3 * k^2 * tau - k * tau
    l3[ind] <- 0
    l3[ii] <- 0
  }
  if (level > 1) {
    # fourth derivative
    l4 <- 6 * a^4 * k^4 * y + 6 * a^3 * k^4 * tau - 12 * a^3 * k^3 * y +
      11 * a^2 * k^4 * tau^2 - 11 * a^2 * k^4 * tau - 12 * a^2 * k^3 * tau +
      7 * a^2 * k^2 * y + 12 * a * k^4 * tau^3 - 18 * a * k^4 * tau^2 +
      6 * a * k^4 * tau - 18 * a * k^3 * tau^2 + 18 * a * k^3 * tau +
      7 * a * k^2 * tau - a * k * y + 6 * k^4 * tau^4 - 12 * k^4 * tau^3 +
      7 * k^4 * tau^2 - k^4 * tau - 12 * k^3 * tau^3 + 18 * k^3 * tau^2 -
      6 * k^3 * tau + 7 * k^2 * tau^2 - 7 * k^2 * tau - k * tau
    l4[ind] <- 0
    l4[ii] <- 0
  }

  list(l1 = l1, l2 = l2, l3 = l3, l4 = l4)
}

#' Log-likelihood derivatives w.r.t. \eqn{\vartheta_0}.
#'
#' @param g     - \eqn{\gamma}, a numeric vector,
#' @param y     - \eqn{y}, a numeric vector,
#' @param th0   - \eqn{\vartheta_0}, a numeric,
#' @param v     - \code{v}, a list containing \eqn{\kappa}, \eqn{\tau}, \code{lg}, \code{ind} and \code{ii}.
#'
#' @return A list of the first and second derivatives of the
#'          log-likelihood w.r.t. \eqn{\vartheta_0}.
ldth0 <- function(g, y, th0, v) {
  a <- exp(th0)
  k <- v$k
  tau <- v$tau
  w <- (v$lg / a) - k
  ind <- v$ind
  ii <- v$ii

  l1 <- l2 <- NULL

  # first derivative
  l1 <- -a * k * y + tau * w + y -
    1 / a * (digamma(y + 1 / a) - digamma(1 / a))
  l1[ind] <- y[ind] - (1 / a) * (digamma(y[ind] + (1 / a)) - digamma(1 / a))
  l1[ii] <- (1 / a) * (th0 + g[ii] - 1 - digamma(y[ii] + (1 / a)) +
    digamma(1 / a))

  # second derivative
  l2 <- a^2 * k^2 * y + a * k^2 * tau - a * k * y + tau^2 * w^2 - tau * w^2 -
    tau * w + (1 / a) * (digamma(y + 1 / a) - digamma(1 / a)) +
    (1 / (a^2)) * (psigamma(y + 1 / a, 1) - psigamma(1 / a, 1))
  l2[ind] <- (1 / a) * (digamma(y[ind] + 1 / a) - digamma(1 / a)) +
    (1 / (a^2)) * (psigamma(y[ind] + (1 / a), 1) - psigamma(1 / a, 1))
  l2[ii] <- (1 / a) * (2 - th0 - g[ii] - digamma(y[ii] + 1 / a) -
    digamma(1 / a) - (1 / a) * (psigamma(y[ii] + (1 / a), 1) -
      psigamma(1 / a, 1)))

  list(l1 = l1, l2 = l2)
}

#' Log-likelihood (\eqn{\ell}) derivatives w.r.t. \eqn{\vartheta_0},
#' and mixed derivatives of \eqn{\ell} w.r.t. \eqn{\gamma} and \eqn{\vartheta_0}.
#'
#' @param g     - \eqn{\gamma}, a numeric vector,
#' @param y     - \eqn{y}, a numeric vector,
#' @param th0   - \eqn{\vartheta_0}, a numeric,
#' @param v     - \code{v}, a list containing \eqn{\kappa}, \eqn{\tau}, \code{lg}, \code{ind} and \code{ii},
#' @param level
#' \itemize{
#'   \item \eqn{==0} – list of NULLs (not needed for estimating parameters),
#'   \item \eqn{> 0} – derivatives needed for quasi-Newton,
#'   \item \eqn{> 1} – derivatives needed for full Newton.
#' }
#'
#' @return A list of the first, second and mixed derivatives of the
#'          log-likelihood w.r.t. \eqn{\vartheta_0} and \eqn{\gamma}.
ldgth0 <- function(g, y, th0, v, level = 2) {
  a <- exp(th0)
  k <- v$k
  tau <- v$tau
  w <- (v$lg / a) - k
  ind <- v$ind
  ii <- v$ii

  l1 <- l2 <- NULL
  l_gth0 <- l_ggth0 <- l_gth0th0 <- l_gggth0 <- l_ggth0th0 <- NULL

  if (level > 0) {
    # ∂ℓ/ϑ₀
    l1 <- -a * k * y + tau * w + y -
      1 / a * (digamma(y + 1 / a) - digamma(1 / a))
    l1[ind] <- y[ind] - (1 / a) * (digamma(y[ind] + (1 / a)) - digamma(1 / a))
    l1[ii] <- (1 / a) * (th0 + g[ii] - 1 - digamma(y[ii] + (1 / a)) +
      digamma(1 / a))

    # ∂²ℓ/∂𝛄∂ϑ₀
    l_gth0 <- a^2 * k^2 * y + a * k^2 * tau - a * k * y - k * tau^2 * w +
      k * tau * w
    l_gth0[ind] <- 0
    l_gth0[ii] <- 0

    # ∂³ℓ/∂𝛄²∂ϑ₀
    l_ggth0 <- -2 * a^3 * k^3 * y - 2 * a^2 * k^3 * tau + 3 * a^2 * k^2 * y -
      2 * a * k^3 * tau^2 + 2 * a * k^3 * tau + a * k^2 * tau^2 * w -
      a * k^2 * tau * w + 2 * a * k^2 * tau - a * k * y + 2 * k^2 * tau^3 * w -
      3 * k^2 * tau^2 * w + k^2 * tau * w - k * tau^2 * w + k * tau * w
    l_ggth0[ind] <- 0
    l_ggth0[ii] <- 1 / a
  }
  if (level > 1) {
    # ∂²ℓ/∂ϑ₀²
    l2 <- a^2 * k^2 * y + a * k^2 * tau - a * k * y + tau^2 * w^2 - tau * w^2 -
      tau * w + (1 / a) * (digamma(y + 1 / a) - digamma(1 / a)) +
      (1 / (a^2)) * (psigamma(y + 1 / a, 1) - psigamma(1 / a, 1))
    l2[ind] <- (1 / a) * (digamma(y[ind] + 1 / a) - digamma(1 / a)) +
      (1 / (a^2)) * (psigamma(y[ind] + (1 / a), 1) - psigamma(1 / a, 1))
    l2[ii] <- (1 / a) * (2 - th0 - g[ii] - digamma(y[ii] + 1 / a) -
      digamma(1 / a) - (1 / a) * (psigamma(y[ii] + (1 / a), 1) -
        psigamma(1 / a, 1)))

    # ∂³ℓ/∂𝛄∂ϑ₀²
    l_gth0th0 <- -2 * a^3 * k^3 * y - 2 * a^2 * k^3 * tau + 3 * a^2 * k^2 * y -
      a * k^3 * tau^2 + a * k^3 * tau + 2 * a * k^2 * tau^2 * w -
      2 * a * k^2 * tau * w + a * k^2 * tau - a * k * y - 2 * k * tau^3 * w^2 +
      3 * k * tau^2 * w^2 + k * tau^2 * w - k * tau * w^2 - k * tau * w
    l_gth0th0[ind] <- 0
    l_gth0th0[ii] <- -1 / a

    # ∂⁴ℓ/∂𝛄³∂ϑ₀
    l_gggth0 <- 6 * a^4 * k^4 * y + 6 * a^3 * k^4 * tau - 12 * a^3 * k^3 * y +
      9 * a^2 * k^4 * tau^2 - 9 * a^2 * k^4 * tau - 2 * a^2 * k^3 * tau^2 * w +
      2 * a^2 * k^3 * tau * w - 10 * a^2 * k^3 * tau + 7 * a^2 * k^2 * y +
      6 * a * k^4 * tau^3 - 9 * a * k^4 * tau^2 + 3 * a * k^4 * tau -
      6 * a * k^3 * tau^3 * w + 9 * a * k^3 * tau^2 * w - 9 * a * k^3 * tau^2 -
      3 * a * k^3 * tau * w + 9 * a * k^3 * tau + 3 * a * k^2 * tau^2 * w -
      3 * a * k^2 * tau * w + 4 * a * k^2 * tau - a * k * y -
      6 * k^3 * tau^4 * w + 12 * k^3 * tau^3 * w - 7 * k^3 * tau^2 * w +
      k^3 * tau * w + 6 * k^2 * tau^3 * w - 9 * k^2 * tau^2 * w +
      3 * k^2 * tau * w - k * tau^2 * w + k * tau * w
    l_gggth0[ind] <- 0
    l_gggth0[ii] <- 0

    # ∂⁴ℓ/∂𝛄²∂ϑ₀²
    l_ggth0th0 <- 6 * a^4 * k^4 * y + 6 * a^3 * k^4 * tau - 12 * a^3 * k^3 * y +
      7 * a^2 * k^4 * tau^2 - 7 * a^2 * k^4 * tau - 4 * a^2 * k^3 * tau^2 * w +
      4 * a^2 * k^3 * tau * w - 8 * a^2 * k^3 * tau + 7 * a^2 * k^2 * y +
      2 * a * k^4 * tau^3 - 3 * a * k^4 * tau^2 + a * k^4 * tau -
      8 * a * k^3 * tau^3 * w + 12 * a * k^3 * tau^2 * w - 3 * a * k^3 * tau^2 -
      4 * a * k^3 * tau * w + 3 * a * k^3 * tau + 2 * a * k^2 * tau^3 * w^2 -
      3 * a * k^2 * tau^2 * w^2 + 3 * a * k^2 * tau^2 * w +
      a * k^2 * tau * w^2 - 3 * a * k^2 * tau * w + 2 * a * k^2 * tau -
      a * k * y + 6 * k^2 * tau^4 * w^2 - 12 * k^2 * tau^3 * w^2 -
      2 * k^2 * tau^3 * w + 7 * k^2 * tau^2 * w^2 + 3 * k^2 * tau^2 * w -
      k^2 * tau * w^2 - k^2 * tau * w - 2 * k * tau^3 * w^2 +
      3 * k * tau^2 * w^2 + k * tau^2 * w - k * tau * w^2 - k * tau * w
    l_ggth0th0[ind] <- 0
    l_ggth0th0[ii] <- 0
  }

  list(
    l1 = l1, l2 = l2, l_gth0 = l_gth0, l_ggth0 = l_ggth0, l_gth0th0 = l_gth0th0,
    l_gggth0 = l_gggth0, l_ggth0th0 = l_ggth0th0
  )
}
