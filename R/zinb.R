#' Zero-Inflated Negative Binomial extended family for mgcv
#'
#' @param theta - 𝛉, a numeric vector containing the 3 parameters of the model,
#'                θ₀, θ₁, θ₂.
#' @param link  - link function name, a character string or function name.
#' @param b     - a numeric parameter.
#'
#' @return An S3 object of type c("extended family", "family") consisting of:
#'            family          - name of family character string.
#'            link            - name of link character string.
#'            linkfun         - the link function.
#'            linkinv         - the inverse link function.
#'            dev.resids      - function computing deviance residuals.
#'            Dd              - function returning derivates of deviance
#'                              residuals w.r.t. 𝝻, 𝛄 and 𝝷.
#'            rd              - optional function simulating response data from
#'                              fitted model.
#'            residuals       - optional function for computing residuals.
#'            aic             - function computing twice -ve log-likelihood for
#'                              2df to be added to.
#'            mu.eta          - d𝝻/d𝛈 function (derivative of inverse link with
#'                              respect to 𝛈).
#'            gkg             - supplements `mu.eta` using function
#'                              `fix.family.link.extended.family` with functions
#'                              gkg, where k = 2, 3, or 4 giving the kth deriva-
#'                              tive of the link over the first derivative of
#'                              the link to the power k. For non standard links
#'                              these functions must be supplied.
#'            preinitialize   - optional function of y and family, returning a
#'                              list with optional elements:
#'                                Theta - initial 𝛉,
#'                                y     - modified 𝐲 for use in fitting
#'            initialize      - expression to be evaluated in `gam.fit4` and
#'                              `initial.spg` (see `mgcv`) to initialize 𝝲 or 𝛈.
#'            postproc        - function with arguments `family`, `y`,
#'                              `prior.weights`, `fitted`, `linear.predictors`,
#'                              `offset`, `intercept` to call after fit to
#'                              compute (optionally) the label for the family,
#'                              deviance and null deviance.
#'            ls              - function to evaluate log saturated log-likeli-
#'                              hood and derivates w.r.t. ϕ and 𝛉 for use in
#'                              RE/ML optimization analytically. If deviance
#'                              used is just -2*log-likelihood can just return
#'                              zeroes.
#'            no.r.sq         - optional TRUE/FALSE indicating whether r² can
#'                              computed for family.
#'            validmu         - function used to test whether 𝛄 (𝝻) is valid.
#'            valideta        - function used to test whether 𝛈 is valid.
#'            n.theta         - number of 𝛉 parameters.
#'            predict         - optional function for predicting from model,
#'                              call by `predict.gam`.
#'            ini.theta       - function for initializing 𝛉.
#'            putTheta        - function for storing 𝛉 values in function
#'                              environment.
#'            getTheta        - function for retrieving theta values in
#'                              function environment.
#'            saturated.ll    - optional function to compute saturated
#'                              log-likelihood by Newton method when no
#'                              analytic solution exists.
#'            scale           - < 0 to estimate. Ignored if NULL.
#' @export
ziNB <- function(theta = NULL, link = "identity", b = 0) {
  linktemp <- substitute(link)
  if (!is.character(linktemp)) {
    linktemp <- deparse(linktemp)
  }
  if (linktemp %in% c("identity")) {
    stats <- make.link(linktemp)
  } else {
    stop(linktemp, " link not available for zero-inflated NB")
  }

  n_theta <- 3
  ini_theta <- c(0, 0, 0)
  if (!is.null(theta)) {
    # fixed theta supplied
    ini_theta <- c(theta[1], theta[2], theta[3])
    n_theta <- 0
  }

  env <- new.env(parent = environment(ziNB))

  if (b < 0) {
    b <- 0
    assign(".b", b, envir = env)
  }
  assign(".Theta", ini_theta, envir = env)

  get_theta <- function(trans = FALSE) {
    th <- get(".Theta")
    if (trans) {
      th[2] <- get(".b") + exp(th[2])
    }
    th
  }

  put_theta <- function(theta) {
    assign(".Theta", theta, envir = environment(sys.function()))
  }

  validmu <- function(mu) {
    all(is.finite(mu))
  }

  dev_resids <- function(y, g, wt, theta = NULL) {
    # CAUTION wt parameter not used at all
    if (is.null(theta)) {
      theta <- get(".Theta")
    }
    b <- get(".b")
    eta <- theta[1] + (b + exp(theta[2])) * g

    -2 * zinbll(y, g, eta, theta[3], deriv = 0)$l
  }

  Dd <- function(y, mu, theta, wt = NULL, level = 0) {
    if (is.null(theta)) {
      theta <- get(".Theta")
    }

    deriv <- 1

    if (level == 1) {
      deriv <- 2
    } else if (level > 1) {
      deriv <- 4
    }

    b <- get(".b")
  }

  aic <- function(y, g, theta = NULL, wt, dev) {
    if (is.null(theta)) {
      theta <- get(".Theta")
    }
    b <- get(".b")
    eta <- theta[1] + (b + exp(theta[2])) * g
    sum(-2 * wt * zinbll(y, g, eta, theta[3], 0))
  }

  # CAUTION not sure about this
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
      stop("negative values not allowed for the zero-inflated NB family")
    }
    if (all.equal(y, round(y)) != TRUE) {
      stop("Non-integer response variables are not allowed with ziNB")
    }
    if ((min(y) == 0 && max(y) == 1)) {
      stop("using ziNB for binary data makes no sense")
    }

    mustart <- log(y + (y == 0) / 5)
  })

  postproc <- function(family, y, prior_weights, fitted, linear_predictors,
                       offset, intercept) {
    posr <- list()
    posr$family <- paste("Zero-Inflated Negative Binomial(",
      paste(round(family$getTheta(TRUE), 3), collapse = ","), ")",
      sep = ""
    )
    lf <- family$saturated.ll(y, family, prior_weights)
    l2 <- family$dev.resids(y, linear_predictors, prior_weights)
    posr$deviance <- sum(l2 - lf)

    fnull <- function(g, family, y, wt) {
      sum(family$dev.resids(y, rep(g, length(y)), wt))
    }
    meany <- mean(y)
    posr$null.deviance <- optimize(fnull,
      interval = c(meany / 5, meany * 3),
      family = family, y = y, wt = prior_weights
    )$objective - sum(lf)

    posr
  }

  rd <- function(g, wt, scale) {
    rzip <- function(g, theta) {
      y <- g
      n <- length(y)
      lambda <- exp(g)
      mlam <- max(c(lambda[is.finite(lambda)], .Machine$double.eps^.2))
      lambda[!is.finite(lambda)] <- mlam
      b <- get(".b")
      eta <- theta[1] + (b + exp(theta[2])) * g
      p <- 1 - exp(-exp(eta))
      ind <- p > runif(n)
      y[!ind] <- 0

      # generate from zero-truncated NB, given y>0
      a <- exp(theta[3])
      prob <- 1 / (1 + a * exp(g))
      y[ind] <- actuar::rztnbinom(sum(ind), (1 / a), prob[ind])

      y
    }

    rzip(g, get(".Theta"))
  }

  saturated_ll <- function(y, family, wt = rep(1, length(y))) {
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

    while (keep_on) {
      step <- -r$Dmu / r$Dmu2
      step[!ucov] <- 0
      mu1 <- mu + step
      l1 <- family$dev.resids(y, mu1, wt, theta)
      ind <- l1 > l & ucov
      kk <- 0

      while (step(ind) > 0 && k < 50) {
        step[ind] <- step[ind] / 2
        mu1 <- mu + step
        l1 <- family$dev.resids(y, mu1, wt, theta)
        ind <- l1 > l & ucov
        kk <- kk + 1
      }

      mu <- mu1
      l <- l1

      ucov <- abs(r$Dmu) > lmax * 1e-7
      k <- k + 1
      if (all(!ucov) || k == 100) {
        keep_on <- FALSE
      }
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
      fv <- mgcv::predict.gam(object(object, type = "response"))
      s <- attr(res, "sign")
      if (is.null(s)) {
        s <- sign(y - fv)
      }
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
    fv <- lambda <- exp(gamma)
    d <- zinb::btlg(g, theta[3], what = c("b", "tau"))
    mu <- d$tau * lambda
    # the above should handle limiting behaviour of g already,
    # but just in case we have the lines below.
    mu[d$ii] <- lambda
    mu[d$ind] <- 1

    fv <- list(q * mu)

    if (is.null(se)) {
      return(fv)
    } else {
      dq_dg <- exp(-et) * et * (b + exp(theta[2]))
      dmu_dg <- lambda * (d$b * d$tau - d$b * d$tau^2 + d$tau)
      dmu_dg[d$ind] <- d$b
      dmu_dg[d$ii] <- lambda

      fv[[2]] <- abs(dq_dg * mu + dmu_dg * q) * se
      names(fv) <- c("fit", "se.fit")
    }
  }


  environment(aic) <- environment(Dd) <- environment(dev_resids) <-
    environment(get_theta) <- environment(predict) <- environment(put_theta) <-
    environment(rd) <- environment(saturated_ll)

  structure(list(
    family = "zero-inflated negative binomial", link = linktemp,
    linkfun = stats$linkfun, linkinv = stats$linkinv, dev.resids = dev_resids,
    Dd = Dd, rd = rd, residuals = residuals, aic = aic, mu.eta = stats$mu.eta,
    g2g = stats$g2g, g3g = stats$g3g, g4g = stats$g4g, initialize = initialize,
    postproc = postproc, ls = ls, no.r.sq = TRUE, validmu = validmu,
    valideta = stats$valideta, n.theta = n_theta, predict = predict,
    ini.theta = ini_theta, putTheta = put_theta, getTheta = get_theta,
    saturated.ll = saturated_ll
  ), class = c("extended.family", "family"))
}
