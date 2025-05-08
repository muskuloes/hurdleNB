#  An S3 object of class `extended.family` consists of:
#            family          - name of family character string.
#            link            - name of link character string.
#            linkfun         - the link function.
#            linkinv         - the inverse link function.
#            dev.resids      - function computing deviance residuals.
#            Dd              - function returning derivatives of deviance
#                              residuals w.r.t. 𝝻, 𝛄 and 𝝷.
#            rd              - optional function simulating response data from
#                              fitted model.
#            residuals       - optional function for computing residuals.
#            aic             - function computing twice -ve log-likelihood for
#                              2df to be added to.
#            mu.eta          - d𝝻/d𝛈 function (derivative of inverse link with
#                              respect to 𝛈).
#            gkg             - supplements `mu.eta` using function
#                              `fix.family.link.extended.family` with functions
#                              gkg, where k = 2, 3, or 4 giving the kth deriva-
#                              tive of the link over the first derivative of
#                              the link to the power k. For non standard links
#                              these functions must be supplied.
#            preinitialize   - optional function of y and family, returning a
#                              list with optional elements:
#                                Theta - initial 𝛉,
#                                y     - modified 𝐲 for use in fitting
#            initialize      - expression to be evaluated in `gam.fit4` and
#                              `initial.spg` (see `mgcv`) to initialize 𝝲 or 𝛈.
#            postproc        - function with arguments `family`, `y`,
#                              `prior.weights`, `fitted`, `linear.predictors`,
#                              `offset`, `intercept` to call after fit to
#                              compute (optionally) the label for the family,
#                              deviance and null deviance.
#            ls              - function to evaluate log saturated log-likeli-
#                              hood and derivatives w.r.t. ϕ and 𝛉 for use in
#                              RE/ML optimization analytically. If deviance
#                              used is just -2*log-likelihood can just return
#                              zeroes.
#            no.r.sq         - optional TRUE/FALSE indicating whether r² can
#                              computed for family.
#            validmu         - function used to test whether 𝛄 (𝝻) is valid.
#            valideta        - function used to test whether 𝛈 is valid.
#            n.theta         - number of 𝛉 parameters.
#            predict         - optional function for predicting from model,
#                              call by `predict.gam`.
#            ini.theta       - function for initializing 𝛉.
#            putTheta        - function for storing 𝛉 values in function
#                              environment.
#            getTheta        - function for retrieving theta values in
#                              function environment.
#            saturated.ll    - optional function to compute saturated
#                              log-likelihood by Newton method when no
#                              analytic solution exists.
#            scale           - < 0 to estimate. Ignored if NULL.


#' Hurdle negative binomial extended family for use with mgcv's `gam` or
#' `bam`.
#'
#' @details
#' \eqn{(\vartheta_1, b + e^{\vartheta_2}, \vartheta_0)} are the `theta`
#' parameters reported in the model summary.
#'
#' @param theta \eqn{\vartheta}, a numeric vector containing the 3 parameters
#'              of the model, \eqn{\vartheta_0}, \eqn{\vartheta_1},
#'               \eqn{\vartheta_2},
#' @param link  Link function name, a character string or function name.
#'              Only the identity -- \eqn{\log}-\eqn{\log} link and the
#'              \eqn{\log} link for \eqn{q} and \eqn{\mu}, respectively,
#'              are currently supported,
#' @param b     A non-negative numeric parameter, specifying the minimum
#'              dependence of the zero-inflation rate on the linear predictor.
#'
#' @return An object of class `extended.family`.
#'
#' @export
#'
#' @examples
#' library(mgcv)
#' library(hurdleNB)
#'
#' set.seed(1)
#' n <- 400
#' dat <- gamSim(1, n = n)
#' dat$y <- rhurdleNB(dat$f / 4 - 1)
#'
#' m <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3),
#'   data = dat, family = hurdleNB()
#' )
#'
#' m$outer.info # check convergence!
#' plot(m, page = 1)
#' plot(m, pages = 1, unconditional = TRUE)
hurdleNB <- function(theta = NULL, link = "identity", b = 0) {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  }
  if (linktemp %in% c("identity")) {
    stats <- make.link(linktemp)
  } else {
    stop(linktemp, " link not available for hurdleNB")
  }

  n_theta <- 3
  iniTheta <- c(0, 0, 1)
  if (!is.null(theta)) {
    # fixed theta supplied
    iniTheta <- c(theta[1], theta[2], theta[3])
    n_theta <- 0 # no thetas to estimate
  }

  env <- new.env(parent = environment(hurdleNB))

  if (b < 0) b <- 0

  assign(".b", b, envir = env)
  assign(".Theta", iniTheta, envir = env)

  getTheta <- function(trans = FALSE) {
    th <- get(".Theta")
    if (trans) {
      th[2] <- get(".b") + exp(th[2])
    }

    th
  }

  putTheta <- function(theta) {
    assign(".Theta", theta, envir = environment(sys.function()))
  }

  validmu <- function(mu) {
    all(is.finite(mu))
  }

  dev.resids <- function(y, g, wt, theta = NULL) {
    if (is.null(theta)) {
      theta <- get(".Theta")
    }
    b <- get(".b")
    eta <- theta[1] + (b + exp(theta[2])) * g

    -2 * hurdleNB_ll(y, g, eta, theta[3], level = 0)$l
  }

  Dd <- function(y, g, theta, wt = NULL, level = 0) {
    if (is.null(theta)) {
      theta <- get(".Theta")
    }

    b <- get(".b")
    lin <- lind(g, theta, level, b)

    deriv <- 1
    if (level == 1) deriv <- 2 else if (level > 1) deriv <- 4
    z <- hurdleNB_ll(y, g, lin$eta, theta[3], level = deriv)
    n <- length(y)
    if (is.null(wt)) wt <- rep(1, n)

    oo <- list()
    oo$Dmu <- -2 * wt * (z$l1[, 1] + z$l1[, 2] * lin$eta_g)
    oo$Dmu2 <- -2 * wt * (z$l2[, 1] + z$l2[, 3] * (lin$eta_g^2) +
      z$l1[, 2] * lin$eta_gg)
    # print("z$l2[, 1][20291]")
    # print(z$l2[, 1][20291])
    # print("z$l2[, 3][20291]")
    # print(z$l2[, 3][20291])
    # print("lin$eta_g^2")
    # print((lin$eta_g^2))
    # print("(z$l2[, 3] * (lin$eta_g^2))[20291]")
    # print((z$l2[, 3] * (lin$eta_g^2))[20291])
    # print("z$l1[, 2][20291]")
    # print(z$l1[, 2][20291])
    # print("(lin$eta_gg)")
    # print((lin$eta_gg))
    # print("(z$l1[, 2] * lin$eta_gg)[20291]")
    # print((z$l1[, 2] * lin$eta_gg)[20291])
    oo$EDmu2 <- -2 * wt * (z$El2[, 1] + z$El2[, 3] * (lin$eta_g^2))

    # quasi-Newton
    if (level > 0) {
      oo$Dth <- oo$Dmuth <- oo$Dmu2th <- matrix(0, n, 3)

      oo$Dth[, 1:2] <- -2 * wt * z$l1[, 2] * lin$eta_th
      # dϑ₀.
      oo$Dth[, 3] <- -2 * wt * z$l1[, 3]

      oo$Dmuth[, 1:2] <- -2 * wt * (z$l2[, 3] * lin$eta_th * lin$eta_g +
        z$l1[, 2] * lin$eta_gth)
      # dϑ₀.
      oo$Dmuth[, 3] <- -2 * wt * z$l2[, 4]

      oo$Dmu2th[, 1:2] <- -2 * wt * (z$l3[, 4] * lin$eta_th * (lin$eta_g^2) +
        z$l2[, 3] * (2 * lin$eta_gth * lin$eta_g + lin$eta_gg * lin$eta_th) +
        z$l1[, 2] * lin$eta_ggth)
      # dϑ₀.
      oo$Dmu2th[, 3] <- -2 * wt * z$l3[, 5]

      oo$Dmu3 <- -2 * wt * (z$l3[, 1] + z$l3[, 4] * (lin$eta_g^3) +
        3 * z$l2[, 3] * lin$eta_g * lin$eta_gg + z$l1[, 2] * lin$eta_ggg)
    }

    # full Newton
    if (level > 1) {
      eta_thth <- matrix(0, n, 3)
      eta_thth[, 1] <- lin$eta_th[, 1]^2
      eta_thth[, 2] <- lin$eta_th[, 1] * lin$eta_th[, 2]
      eta_thth[, 3] <- lin$eta_th[, 2]^2

      oo$Dmu2th2 <- oo$Dmuth2 <- oo$Dth2 <- matrix(0, n, 6)

      # contains dϑ₁dϑ₁, dϑ₁dϑ₂, dϑ₂dϑ₂.
      Dth2 <- -2 * wt * (z$l2[, 3] * eta_thth +
        z$l1[, 2] * lin$eta_th2)

      oo$Dth2[, 1:2] <- Dth2[, 1:2]
      # dϑ₁dϑ₀.
      oo$Dth2[, 3] <- 0
      # dϑ₂dϑ₂.
      oo$Dth2[, 4] <- Dth2[, 3]
      # dϑ₂dϑ₀.
      oo$Dth2[, 5] <- 0
      # dϑ₀dϑ₀.
      oo$Dth2[, 6] <- -2 * wt * z$l2[, 5]

      eta_gthth <- matrix(0, n, 3)
      eta_gthth[, 1] <- 2 * lin$eta_gth[, 1] * lin$eta_th[, 1]
      eta_gthth[, 2] <- lin$eta_gth[, 1] * lin$eta_th[, 2] +
        lin$eta_gth[, 2] * lin$eta_th[, 1]
      eta_gthth[, 3] <- 2 * lin$eta_gth[, 2] * lin$eta_th[, 2]

      # contains dϑ₁dϑ₁, dϑ₁dϑ₂, dϑ₂dϑ₂.
      Dmuth2 <- -2 * wt * (z$l3[, 4] * eta_thth * lin$eta_g +
        z$l2[, 3] * (lin$eta_th2 * lin$eta_g + eta_gthth) +
        z$l1[, 2] * lin$eta_gth2)
      oo$Dmuth2[, 1:2] <- Dmuth2[, 1:2]
      # dϑ₁dϑ₀.
      oo$Dmuth2[, 3] <- 0
      # dϑ₂dϑ₂.
      oo$Dmuth2[, 4] <- Dmuth2[, 3]
      # dϑ₂dϑ₀.
      oo$Dmuth2[, 5] <- 0
      # dϑ₀dϑ₀.
      oo$Dmuth2[, 6] <- -2 * wt * z$l3[, 6]

      oo$Dmu3th <- matrix(0, n, 3)
      oo$Dmu3th[, 1:2] <- -2 * wt * (z$l4[, 5] * lin$eta_th * (lin$eta_g^3) +
        3 * z$l3[, 4] * (lin$eta_g^2 * lin$eta_gth +
          lin$eta_th * lin$eta_g * lin$eta_gg) +
        z$l2[, 3] * (3 * lin$eta_gth * lin$eta_gg +
          3 * lin$eta_g * lin$eta_ggth + lin$eta_th * lin$eta_ggg) +
        z$l1[, 2] * lin$eta_gggth)
      # dϑ₀.
      oo$Dmu3th[, 3] <- -2 * wt * z$l4[, 6]

      eta_gthgth <- matrix(0, n, 3)
      eta_gthgth[, 1] <- 2 * (lin$eta_gth[, 1]^2)
      eta_gthgth[, 2] <- 2 * (lin$eta_gth[, 1] * lin$eta_gth[, 2])
      eta_gthgth[, 3] <- 2 * (lin$eta_gth[, 2]^2)
      eta_ggthth <- matrix(0, n, 3)
      eta_ggthth[, 1] <- 2 * lin$eta_th[, 1] * lin$eta_ggth[, 1]
      eta_ggthth[, 2] <- lin$eta_th[, 1] * lin$eta_ggth[, 2] +
        lin$eta_th[, 2] * lin$eta_ggth[, 1]
      eta_ggthth[, 3] <- 2 * lin$eta_th[, 2] * lin$eta_ggth[, 2]

      # contains dϑ₁dϑ₁, dϑ₁dϑ₂, dϑ₂dϑ₂.
      Dmu2th2 <- -2 * wt * (z$l4[, 5] * eta_thth * (lin$eta_g^2) +
        z$l3[, 4] * (lin$eta_th2 * (lin$eta_g^2) +
          2 * eta_gthth * lin$eta_g + eta_thth * lin$eta_gg) +
        z$l2[, 3] * (eta_gthgth + 2 * lin$eta_g * lin$eta_gth2 + eta_ggthth +
          lin$eta_th2 * lin$eta_gg) + z$l1[, 2] * lin$eta_ggth2)
      oo$Dmu2th2[, 1:2] <- Dmu2th2[, 1:2]
      # dϑ₁dϑ₀.
      oo$Dmu2th2[, 3] <- 0
      # dϑ₂dϑ₂.
      oo$Dmu2th2[, 4] <- Dmu2th2[, 3]
      # dϑ₂dϑ₀.
      oo$Dmu2th2[, 5] <- 0
      # dϑ₀dϑ₀.
      oo$Dmu2th2[, 6] <- -2 * wt * z$l4[, 7]

      oo$Dmu4 <- -2 * wt * (z$l4[, 1] + z$l4[, 5] * (lin$eta_g^4) +
        6 * z$l3[, 4] * (lin$eta_g^2) * lin$eta_gg +
        z$l2[, 3] * (3 * (lin$eta_gg^2) + 4 * lin$eta_g * lin$eta_ggg) +
        z$l1[, 2] * lin$eta_gggg)
    }

    oo
  }

  aic <- function(y, g, theta = NULL, wt, dev) {
    if (is.null(theta)) {
      theta <- get(".Theta")
    }
    b <- get(".b")
    eta <- theta[1] + (b + exp(theta[2])) * g

    sum(-2 * wt * hurdleNB_ll(y, g, eta, theta[3], level = 0)$l)
  }

  ls <- function(y, wt, theta, scale) {
    list(
      ls = 0,
      lsth1 = c(0, 0, 0),
      LSTH1 = matrix(0, length(y), 3),
      lsth2 = matrix(0, 3, 3)
    )
  }

  initialize <- expression({
    if (any(y < 0)) {
      stop("negative values not allowed for the hurdleNB family")
    }
    # if (all.equal(y, round(y)) != TRUE) {
    #   stop("Non-integer response variables are not allowed with hurdleNB")
    # }
    if ((min(y) == 0 && max(y) == 1)) {
      stop("using hurdleNB for binary data makes no sense")
    }

    mustart <- log(y + (y == 0) / 5)
  })

  postproc <- function(family, y, prior.weights, fitted, linear.predictors,
                       offset, intercept) {
    posr <- list()
    posr$family <- paste("Hurdle negative binomial(",
      paste(round(family$getTheta(TRUE), 3), collapse = ","), ")",
      sep = ""
    )
    lf <- family$saturated.ll(y, family, prior.weights)
    l2 <- family$dev.resids(y, linear.predictors, prior.weights)
    posr$deviance <- sum(l2 - lf)

    fnull <- function(g, family, y, wt) {
      sum(family$dev.resids(y, rep(g, length(y)), wt))
    }
    meany <- mean(y)
    posr$null.deviance <- optimize(fnull,
      interval = c(meany / 5, meany * 3),
      family = family, y = y, wt = prior.weights
    )$objective - sum(lf)

    posr
  }

  rd <- function(g, wt, scale) {
    rhurdleNB(g, get(".Theta"), get(".b"))
  }

  saturated.ll <- function(y, family, wt = rep(1, length(y))) {
    pind <- y > 0
    wt <- wt[pind]
    y <- y[pind]
    g <- log(y)
    keep_on <- TRUE
    theta <- family$getTheta()
    r <- family$Dd(y, g, theta, wt)
    l <- family$dev.resids(y, g, wt, theta)
    lmax <- max(abs(l))
    ucov <- abs(r$Dmu) > lmax * 1e-7
    k <- 0
    
    # eta = theta[1] + exp(theta[2])*g
    
    # print("Initial values")
    print("y[20291]")
    print(y[20291])
    # print("g")
    # print(g)
    # print("eta")
    # print(eta)
    # print("lmax")
    # print(lmax)
    print("l[20291]")
    print(l[20291])
    # print("theta")
    # print(theta)
    # print("g[20291]")
    # print(g[20291])
    # print("r$Dmu[20291]")
    # print(r$Dmu[20291])
    # print("r$Dmu2[20291]")
    # print(r$Dmu2[20291])

    while (keep_on) {
      print(k)
      step <- -r$Dmu / (r$Dmu2+1e-8)
      step[!ucov] <- 0
      g1 <- g + step
      l1 <- family$dev.resids(y, g1, wt, theta)
      ind <- l1 > l & ucov
      kk <- 0

      print("Loop")
      print("very small y")
      print(y[which(r$Dmu2 == 0)])
      # print("y[20291]")
      # print(y[20291])
      # print("step[20291]")
      # print(step[20291])
      # print("g1[20291]")
      # print(g1[20291])
      print("l1[20291]")
      print(l1[20291])
      # print("ind[20291]")
      # print(ind[20291])
      # print("r$Dmu[20291]")
      # print(r$Dmu[20291])
      # print("r$Dmu2[20291]")
      # print(r$Dmu2[20291])
      # print("kk")
      # print(kk)
      # print("where NA")
      # print(which(is.na(ind)))
      
      while (sum(ind) > 0 && kk < 50) {
        step[ind] <- step[ind] / 2
        g1 <- g + step
        l1 <- family$dev.resids(y, g1, wt, theta)
        ind <- l1 > l & ucov
        kk <- kk + 1
        
        # print("inner")
        # eta = theta[1] + exp(theta[2])*g1
        
        # print("inner\n")
        # print("eta\n")
        # print(eta[which(is.na(ind))])
        # print("theta\n")
        # print(theta)
        # print("l1\n")
        # print(l1[which(is.na(ind))])
        # print("l\n")
        # print(l[which(is.na(ind))])
        # print("y\n")
        # print(y[which(is.na(ind))])
        # print("g1\n")
        # print(g1[which(is.na(ind))])
        # print("theta\n")
      }

      g <- g1
      l <- l1

      r <- family$Dd(y, g, theta, wt)
      ucov <- abs(r$Dmu) > lmax * 1e-7
      k <- k + 1
      if (all(!ucov) || k == 100) keep_on <- FALSE
    }
    l1 <- rep(0, length(pind))
    l1[pind] <- l

    l1
  }

  residuals <- function(object, type = c("deviance", "working", "reponse")) {
    if (type == "working") {
      res <- object$residuals
    } else if (type == "response") {
      res <- object$y - mgcv::predict.gam(object, type = "response")
    } else if (type == "deviance") {
      y <- object$y
      g <- object$linear.predictors
      wts <- object$prior.weights
      res <- object$family$dev.resids(y, g, wts)

      res <- res - object$family$saturated.ll(y, object$family, wts)
      fv <- mgcv::predict.gam(object, type = "response")
      s <- attr(res, "sign")
      if (is.null(s)) s <- sign(y - fv)
      res <- as.numeric(s * sqrt(pmax(res, 0)))
    }

    res
  }

  predict <- function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
                      beta = NULL, off = NULL, Vb = NULL) {
    theta <- family$getTheta()

    if (is.null(eta)) {
      discrete <- is.list(X)
      g <- off + if (discrete) {
        mgcv::Xbd(X$Xd, beta,
          k = X$kd, ks = X$ks, ts = X$ts, dt = X$dt,
          v = X$v, qc = X$qc, drop = X$drop
        )
      } else {
        drop(X %*% beta)
      }

      if (se) {
        se <- if (discrete) {
          sqrt(pmax(0, mgcv::diagXVXd(X$Xd, Vb,
            k = X$kd, ks = X$ks, ts = X$ts, dt = X$dt,
            v = X$v, qc = X$qc, drop = X$drop, nthreads = 1
          )))
        } else {
          sqrt(pmax(0, rowSums((X %*% Vb) * X)))
        }
      } else {
        se <- NULL
      }
    } else {
      se <- NULL
      g <- eta
    }

    b <- get(".b")
    eta <- theta[1] + (b + exp(theta[2])) * g
    et <- exp(eta)
    q <- 1 - exp(-et)
    fv <- lambda <- exp(g)
    v <- ktlg(g, exp(theta[3]), what = c("b", "tau"))
    mu <- v$tau * lambda
    # the above should handle limiting behaviour of g already,
    # but just in case we have the lines below.
    mu[v$ii] <- lambda[v$ii]
    mu[v$ind] <- 1

    fv <- list(q * mu)

    if (is.null(se)) {
      return(fv)
    } else {
      dq_dg <- exp(-et) * et * (b + exp(theta[2]))
      dmu_dg <- lambda * (v$k * v$tau - v$k * v$tau^2 + v$tau)
      dmu_dg[v$ind] <- v$k[v$ind]
      dmu_dg[v$ii] <- lambda[v$ii]

      fv[[2]] <- abs(dq_dg * mu + dmu_dg * q) * se
      names(fv) <- c("fit", "se.fit")

      return(fv)
    }
  }

  environment(aic) <- environment(Dd) <- environment(dev.resids) <-
    environment(getTheta) <- environment(predict) <- environment(putTheta) <-
    environment(rd) <- environment(saturated.ll) <- env

  structure(list(
    family = "hurdle negative binomial", link = linktemp,
    linkfun = stats$linkfun, linkinv = stats$linkinv, dev.resids = dev.resids,
    Dd = Dd, rd = rd, residuals = residuals, aic = aic, mu.eta = stats$mu.eta,
    g2g = stats$g2g, g3g = stats$g3g, g4g = stats$g4g, initialize = initialize,
    postproc = postproc, ls = ls, no.r.sq = TRUE, validmu = validmu,
    valideta = stats$valideta, n.theta = n_theta, predict = predict,
    ini.theta = iniTheta, putTheta = putTheta, getTheta = getTheta,
    saturated.ll = saturated.ll
  ), class = c("extended.family", "family"))
}

#' \eqn{\eta = \vartheta_1 + (b + e^{\vartheta_2})\gamma}.
#'
#' @param g     \eqn{\gamma}, a numeric vector,
#' @param theta \eqn{\vartheta}, a numeric vector,
#' @param level A numeric, indicating whether to return deriv
#'              w.r.t. \eqn{\vartheta_1} and \eqn{\vartheta_2},
#' @param b     A numeric.
#'
#' @return A list with \eqn{\eta = \vartheta_1 + (b + e^{\vartheta_2})\gamma}
#'         and its derivatives with respect to \eqn{\gamma}, \eqn{\vartheta_1},
#'         and \eqn{\vartheta_2}.
lind <- function(g, theta, level = 0, b = 0) {
  theta[2] <- exp(theta[2])
  r <- list(eta = theta[1] + (b + theta[2]) * g)
  r$eta_g <- b + theta[2]
  r$eta_gg <- 0

  if (level) {
    n <- length(g)
    r$eta_gggth <- r$eta_ggth <- r$eta_gth <- r$eta_th <- matrix(0, n, 2)
    r$eta_th[, 1] <- 1 # d𝛈/dϑ₁
    r$eta_th[, 2] <- theta[2] * g # d𝛈/dϑ₂
    r$eta_gth[, 2] <- theta[2] # d²𝛈/d𝛄dϑ₂
    r$eta_gggg <- r$eta_ggg <- 0 # d⁴𝛈/d𝛄⁴, d³𝛈/d𝛄³
    # order dϑ₁dϑ₁, dϑ₁dϑ₂, dϑ₂dϑ₂.
    r$eta_ggth2 <- r$eta_gth2 <- r$eta_th2 <- matrix(0, n, 3)
    r$eta_th2[, 3] <- theta[2] * g
    r$eta_gth2[, 3] <- theta[2]
  }

  r
}

#' Generate hurdleNB random variables.
#'
#' @param g     \eqn{\gamma}, a numeric vector,
#' @param theta \eqn{\vartheta}, a numeric vector,
#' @param b     A numeric.
#'
#' @return hurdle negative binomial random variables.
#' @export
rhurdleNB <- function(g, theta = c(-2, 0.3, 2), b = 0) {
  y <- g
  n <- length(y)
  lambda <- exp(g)

  eta <- theta[1] + (b + exp(theta[2])) * g
  p <- 1 - exp(-exp(eta))
  ind <- p > runif(n)
  y[!ind] <- 0

  # generate from zero-truncated NB, given y > 0
  a <- exp(theta[3])
  prob <- 1 / (1 + a * lambda)
  y[ind] <- actuar::rztnbinom(sum(ind), (1 / a), prob[ind])

  y
}
