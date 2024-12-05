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

#' log(1 - (1 + αeˣ)^(-1/α)).
#'
#' @param x   - A numeric vector,
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
#' @param g    - 𝛄, a numeric vector,
#' @param a    - α, a numeric,
#' @param what - A character vector specifying what to return.
#'
#' @return A list containing:
#          k -- 𝛋, tau -- 𝛕, lg -- log(1+αeᵞ),
#'         ind -- indices of yᵢ for which γᵢ is very small,
#'         ii -- indices of yᵢ for which γᵢ is very large.
btlg <- function(g, a, what = c("k", "tau")) {
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
    tau <- 1 / (1 - (1 + a * eg)^-(1 / a))
    lg <- log(1 + a * eg)
    lg[ind] <- 0
    lg[ii] <- g[ii]
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

#' Log-likelihood derivates w.r.t. 𝛈.
#'
#' @param eta   - 𝛈, a numeric vector,
#' @param deriv - <= 1 - first and second derivatives,
#'                == 2 - first, second and third derivatives,
#'                >= 3 - first, second, third and fourth derivatives.
#'
#' @return A list of derivatives of the log-likelihood w.r.t. 𝛈 (eta).
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
    l3 <- l1 * (-et + (1 - et)^2 - 3 * (1 - et) * l1 + 2 * l1^2)
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
#' @param g     - 𝛄, a numeric vector,
#' @param y     - 𝐲, a numeric vector,
#' @param a     - α, a numeric,
#' @param deriv - <= 1 - first and second derivatives,
#'                == 2 - first, second and third derivatives,
#'                >= 3 - first, second, third and fourth derivatives.
#'
#' @return A list of derivatives of the log-likelihood w.r.t. 𝛄 (g).
ldg <- function(g, y, a, deriv = 4) {
  d <- btlg(g, a, c("k", "tau"))
  k <- d$k
  tau <- d$tau
  ind <- d$ind
  ii <- d$ii

  # first derivative
  l1 <- -a * k * y - k * tau + y
  l1[ind] <- y[ind] - 1
  l1[ii] <- -1 / a

  # second derivatauive
  l2 <- a^2 * k^2 * y + a * k^2 * tau - a * k * y + k^2 * tau^2 - k^2 * tau -
    k * tau
  l2[ind] <- 0
  l2[ii] <- 0

  l3 <- l4 <- NULL

  if (deriv > 1) {
    # third derivative
    l3 <- -2 * a^3 * k^3 * y - 2 * a^2 * k^3 * tau + 3 * a^2 * k^2 * y -
      3 * a * k^3 * tau^2 + 3 * a * k^3 * tau + 3 * a * k^2 * tau - a * k * y -
      2 * k^3 * tau^3 + 3 * k^3 * tau^2 - k^3 * tau + 3 * k^2 * tau^2 -
      3 * k^2 * tau - k * tau
    l3[ind] <- 0
    l3[ii] <- 0
  }
  if (deriv > 2) {
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

#' Log-likelihood derivatives w.r.t. θ₀.
#'
#' @param g     - 𝛄, a numeric vector
#' @param y     - 𝐲, a numeric vector
#' @param th0   - θ₀, a numeric
#'
#' @return A list of the first and second derivatives of the
#'          log-likelihood w.r.t. θ₀.
ldth0 <- function(g, y, th0) {
  a <- exp(th0)
  d <- btlg(g, a, c("k", "lg", "tau"))
  k <- d$k
  tau <- d$tau
  w <- d$lg / a - k
  ind <- d$ind
  ii <- d$ii

  l1 <- l2 <- NULL

  # first derivative
  l1 <- -a * k * y + tau * w + y -
    1 / a * (digamma(y + 1 / a) - digamma(1 / a))
  l1[ind] <- y[ind] - (1 / a) * (digamma(y + (1 / a)) - digamma(1 / a))
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

#' Mixed derivatives of l w.r.t. 𝛄 and θ₀.
#'
#' @param g     - 𝛄, a numeric vector.
#' @param y     - 𝐲, a numeric vector.
#' @param th0   - θ₀, a numeric.
#' @param deriv - <= 1 - second mixed derivatives,
#'                <= 3 - second and third mixed derivatives,
#'                 > 3 - second, third and fourth mixed derivatives.
#' @return A list of mixed derivatives of l w.r.t. 𝛄 and θ₀.
ldgth0 <- function(g, y, th0, deriv = 4) {
  a <- exp(th0)
  d <- btlg(g, a, c("k", "lg", "tau"))
  k <- d$k
  tau <- d$tau
  lg <- d$lg
  w <- (lg / a) - k
  ind <- d$ind
  ii <- d$ii

  l_gth0 <- l_ggth0 <- l_gth0th0 <- l_gggth0 <- l_ggth0th0 <- NULL

  # ∂²ℓ/∂𝛄∂θ₀
  l_gth0 <- a^2 * k^2 * y + a * k^2 * tau - a * k * y - k * tau^2 * w +
    k * tau * w
  l_gth0[ind] <- 0
  l_gth0[ii] <- 0

  if (deriv > 1) {
    # ∂³ℓ/∂𝛄²∂θ₀
    l_ggth0 <- -2 * a^3 * k^3 * y - 2 * a^2 * k^3 * tau + 3 * a^2 * k^2 * y -
      2 * a * k^3 * tau^2 + 2 * a * k^3 * tau + a * k^2 * tau^2 * w -
      a * k^2 * tau * w + 2 * a * k^2 * tau - a * k * y + 2 * k^2 * tau^3 * w -
      3 * k^2 * tau^2 * w + k^2 * tau * w - k * tau^2 * w + k * tau * w
    l_ggth0[ind] <- 0
    l_ggth0[ii] <- 1 / a

    # ∂³ℓ/∂𝛄∂θ₀²
    l_gth0th0 <- -2 * a^3 * k^3 * y - 2 * a^2 * k^3 * tau + 3 * a^2 * k^2 * y -
      a * k^3 * tau^2 + a * k^3 * tau + 2 * a * k^2 * tau^2 * w -
      2 * a * k^2 * tau * w + a * k^2 * tau - a * k * y - 2 * k * tau^3 * w^2 +
      3 * k * tau^2 * w^2 + k * tau^2 * w - k * tau * w^2 - k * tau * w
    l_gth0th0[ind] <- 0
    l_gth0th0[ii] <- 0
  }
  if (deriv > 3) {
    # ∂⁴ℓ/∂𝛄³∂θ₀
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

    # ∂⁴ℓ/∂𝛄²∂θ₀²
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
    l_gth0 = l_gth0, l_ggth0 = l_ggth0, l_gth0th0 = l_gth0th0,
    l_gggth0 = l_gggth0, l_ggth0th0 = l_ggth0th0
  )
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
  d <- btlg(g, a, what = c("k", "lg", "tau"))
  k <- d$k
  tau <- d$tau
  lg <- d$lg
  l[!zind] <- l1ee(eta[!zind]) + yp * log(a) + yp * g[!zind] -
    yp * lg[!zind] - l11aea(g[!zind], th0) +
    lgamma(yp + 1 / a) - lgamma(yp + 1) - lgamma(1 / a)
  q <- 1 - exp(-et)

  n <- length(y)

  l1 <- El2 <- l2 <- l3 <- l4 <- NULL

  # first and second derivatives.
  if (deriv > 0) {
    # order ∂ℓ/∂𝛄, ∂ℓ/∂𝛈, ∂ℓ/∂θ₀.
    l1 <- matrix(0, n, 3)
    l_e <- lde(eta, deriv)
    l_g <- ldg(g, y, a, deriv)
    l_dth0 <- ldth0(g, y, th0)
    l_dgth0 <- ldgth0(g, y, th0, deriv)

    l1[!zind, 1] <- l_g$l1[!zind] # ∂ℓ/∂𝛄, y>0
    l1[zind, 2] <- l[zind] # ∂ℓ/∂𝛈, y==0
    l1[!zind, 2] <- l_e$l1[!zind] # ∂ℓ/∂𝛈, y>0
    l1[!zind, 3] <- l_dth0$l1[!zind] # ∂ℓ/∂θ₀, y>0

    # order ∂²ℓ/∂𝛄², ∂²ℓ/∂𝛈∂𝛄, ∂²ℓ/∂𝛈², ∂²ℓ/∂𝛄∂θ₀, ∂²ℓ/∂θ₀².
    l2 <- matrix(0, n, 5)
    # order E[∂²ℓ/∂𝛄²], E[∂²ℓ/∂𝛈∂𝛄], E[∂²ℓ/∂𝛈²].
    El2 <- matrix(0, n, 3)

    l2[!zind, 1] <- l_g$l2[!zind] # ∂²ℓ/∂𝛄², y>0
    l2[zind, 3] <- l[zind] # ∂²ℓ/𝛈², y==0
    l2[!zind, 3] <- l_e$l2[!zind] # ∂²ℓ/∂𝛈², y>0
    l2[!zind, 4] <- l_dgth0$l_gth0[!zind] # ∂²ℓ/∂𝛄θ₀, y>0
    l2[!zind, 5] <- l_dth0$l2[!zind] # ∂²ℓ/∂θ₀², y>0
    El2[, 1] <- q * (q * tau * exp(g) * ((a^2) * k^2 - a * k) +
      a * (k^2) * tau - tau * k + (k^2) * (tau^2) - (k^2) * (tau)) # E[∂²ℓ/∂𝛄²]
    El2[, 3] <- -(1 - q) * et + q * l_e$l2 # E[∂²ℓ/∂𝛈²]
  }

  # third derivates.
  if (deriv > 1) {
    # order ∂³ℓ/∂𝛄³, ∂³ℓ/∂𝛄²∂𝛈, ∂³ℓ/∂𝛄∂𝛈², ∂³ℓ/∂𝛈³, ∂³ℓ/∂𝛄²∂θ₀, ∂³ℓ/∂𝛄∂θ₀².
    l3 <- matrix(0, n, 6)
    l3[!zind, 1] <- l_g$l3[!zind]
    l3[!zind, 4] <- l_e$l3[!zind]
    l3[zind, 4] <- l[zind]
    l3[!zind, 5] <- l_dgth0$l_ggth0[!zind]
    l3[!zind, 6] <- l_dgth0$l_gth0th0[!zind]
  }

  # fourth derivatives
  if (deriv > 3) {
    # order ∂⁴ℓ/∂𝛄⁴, ∂⁴ℓ/∂𝛄³∂𝛈, ∂⁴ℓ/∂𝛄²∂𝛈², ∂⁴ℓ/∂𝛄∂𝛈³, ∂⁴ℓ/∂𝛈⁴, ∂⁴ℓ/∂𝛄³∂θ₀,
    # ∂⁴ℓ/∂𝛄²∂θ₀².
    l4 <- matrix(0, n, 7)
    l4[!zind, 1] <- l_g$l4[!zind]
    l4[!zind, 5] <- l_e$l4[!zind]
    l4[zind, 5] <- l[zind]
    l4[!zind, 6] <- l_dgth0$l_gggth0[!zind]
    l4[!zind, 7] <- l_dgth0$l_ggth0th0[!zind]
  }

  list(l = l, l1 = l1, l2 = l2, l3 = l3, l4 = l4, El2 = El2)
}
