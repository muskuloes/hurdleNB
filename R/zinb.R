#' Zero-Inflated Negative Binomial extended family for mgcv
#'
#' @param theta - 𝛉, a numeric vector containing the 3 parameters of the model,
#'                θ₀,θ₁,θ₂.
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

  n.theta <- 3
  iniTheta <- c(0, 0, 0)
  if (!is.null(theta)) {
    # fixed theta supplied
    iniTheta <- c(theta[1], theta[2], theta[3])
    n.theta <- 0
  }

  env <- new.env(parent = environment(ziNB))

  if (b < 0) {
    b <- 0
    assign(".b", b, envir = env)
  }
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
    # CAUTION wt parameter not used at all
    if (is.null(theta)) {
      theta <- get(".Theta")
    }
    b <- get(".b")
    eta <- theta[1] + (b + exp(theta[2])) * g

    -2 * zinbll(y, g, eta, theta[3], deriv = 0)$l
  }

  Dd <- function(y, mu, theta, wt = NULL, level = 0) {

  }
}
