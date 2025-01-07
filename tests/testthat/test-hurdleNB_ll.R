test_that("hurdleNB_ll works for scalar y at level 0", {
  level <- 0
  z <- hurdleNB_ll(y = 0, g = 10, eta = 0, th0 = log(0.5), level = level)

  # y==0: ℓ, ∂ℓ/∂𝛄, ∂ℓ/∂𝛈, ∂ℓ/∂θ⁰, ∂²ℓ/∂𝛄², ∂²ℓ/∂𝛈∂𝛄, ∂²ℓ/∂𝛈², ∂²ℓ/∂𝛄∂θ₀,
  # ∂²ℓ/∂θ₀².
  expect_equal(z$l, -1)
  expect_equal(z$l1[, 1], 0)
  expect_equal(z$l1[, 2], -1)
  expect_equal(z$l1[, 3], NaN)
  expect_equal(z$l2[, 1], 0)
  expect_equal(z$l2[, 2], 0)
  expect_equal(z$l2[, 3], -1)
  expect_equal(z$l2[, 4], NaN)
  expect_equal(z$l2[, 5], NaN)


  # ∂³ℓ, ∂⁴ℓ.
  expect_null(z$l3)
  expect_null(z$l4)

  z <- hurdleNB_ll(y = 25, g = 10, eta = 0, th0 = log(0.5), level = level)

  # y>0: ℓ, ∂ℓ/∂𝛄, ∂ℓ/∂𝛈, ∂ℓ/∂θ⁰, ∂²ℓ/∂𝛄², ∂²ℓ/∂𝛈∂𝛄, ∂²ℓ/∂𝛈², ∂²ℓ/∂𝛄∂θ₀,
  # ∂²ℓ/∂θ₀².
  expect_equal(round(z$l, 5), -15.81674)
  expect_equal(round(z$l1[, 1], 5), -1.99755)
  expect_equal(round(z$l1[, 2], 5), 0.58198)
  expect_equal(z$l1[, 3], NaN)
  expect_equal(round(z$l2[, 1], 5), -0.00245)
  expect_equal(z$l2[, 2], 0)
  expect_equal(round(z$l2[, 3], 5), -0.33870)
  expect_equal(z$l2[, 4], NaN)
  expect_equal(z$l2[, 5], NaN)

  # 𝔼[∂²ℓ/∂𝛄²], 𝔼[∂²ℓ/∂𝛈∂𝛄], 𝔼[∂²ℓ/∂𝛈²].
  expect_equal(round(z$El2[, 1], 5), -1.26413)
  expect_equal(z$El2[, 2], 0)
  expect_equal(round(z$El2[, 3], 5), -0.58198)

  # ∂³ℓ, ∂⁴ℓ.
  expect_null(z$l3)
  expect_null(z$l4)
})

test_that("hurdleNB_ll works for vector y at level 0", {
  level <- 0

  # example generated using:
  #> x0 <- runif(18)
  #> x1 <- runif(18)
  #> x2 <- runif(18)
  #> x3 <- runif(18)
  #> g <- x0 + 2 * x1 + 3 * x2 + 4 * x3
  #> y <- rhurdleNB(g, theta = c(-2, 0.3, 2), b = 0)
  # the expectations are then gotten using sympy.
  g <- c(
    1.3130418, 1.1740350, 1.7447116, 1.0350755, 0.8175248, 1.2184900,
    1.3780124, 1.1322345, 1.9907406, 0.7240296, 0.8100226, 1.1046235,
    0.9087231, 1.7945166, 0.8435220, 1.4348064, 1.0780594, 0.5603018
  )
  y <- c(15, 0, 60, 0, 7, 0, 0, 1, 78, 0, 0, 4, 0, 7, 0, 1, 0, 0)
  n <- length(y)
  eta <- -2 + exp(0.3) * g
  th0 <- 2

  z <- hurdleNB_ll(y, g, eta, th0, level = level)
  zind <- y == 0

  # y==0: ℓ, ∂ℓ/∂𝛄, ∂ℓ/∂𝛈, ∂²ℓ/∂𝛄², ∂²ℓ/∂𝛈².
  expect_equal(z$l[zind], -exp(eta[zind]))
  expect_equal(z$l1[zind, 1], rep(0, sum(zind)))
  expect_equal(z$l1[zind, 2], -exp(eta[zind]))
  expect_equal(z$l2[zind, 1], rep(0, sum(zind)))
  expect_equal(z$l2[zind, 3], -exp(eta[zind]))

  l <- c(
    -4.86116, -6.74632, -4.38402, -2.18796, -6.95198, -3.48723, -3.63684,
    -2.01474
  )
  l_g <- c(
    0.16857, 1.05483, -0.00144, -0.32947, 1.09841, -0.20166, -0.17410,
    -0.31890
  )
  l_e <- c(
    0.65408, 0.45088, 0.80983, 0.72025, 0.31547, 0.72936, 0.42404,
    0.60301
  )
  l_gg <- c(
    -0.43937, -1.29564, -0.28844, 0.03413, -1.33727, -0.08930, -0.09360, 0.03550
  )
  l_ee <- c(
    -0.29469, -0.39553, -0.17641, -0.24793, -0.41127, -0.24106, -0.40266,
    -0.32668
  )

  # y>0: ℓ, ∂ℓ/∂𝛄, ∂ℓ/∂𝛈, ∂²ℓ/∂𝛄², ∂²ℓ/∂𝛈².
  expect_equal(round(z$l[!zind], 5), l)
  expect_equal(round(z$l1[!zind, 1], 5), l_g)
  expect_equal(round(z$l1[!zind, 2], 5), l_e)
  expect_equal(round(z$l2[!zind, 1], 5), l_gg)
  expect_equal(round(z$l2[!zind, 3], 5), l_ee)

  # ∂ℓ/∂θ₀, ∂²ℓ/∂𝛈∂𝛄, ∂²ℓ/∂𝛄∂θ₀, ∂²ℓ/∂θ₀².
  expect_equal(z$l1[, 3], rep(NaN, n))
  expect_equal(z$l2[, 2], rep(0, n))
  expect_equal(z$l2[, 4], rep(NaN, n))
  expect_equal(z$l2[, 5], rep(NaN, n))

  El_gg <- c(
    -0.15194, -0.13553, -0.20149, -0.11977, -0.09709, -0.14073, -0.15971,
    -0.13070, -0.22336, -0.08825, -0.09636, -0.12755, -0.10626, -0.20648,
    -0.09965, -0.16649, -0.12455, -0.07415
  )
  El_ee <- c(
    -0.52095, -0.46607, -0.64311, -0.41111, -0.33042, -0.48378,
    -0.54558, -0.44942, -0.62722, -0.29884, -0.32781, -0.43845, -0.36315,
    -0.64689, -0.33957, -0.56607, -0.42796, -0.24875
  )

  # 𝔼[∂²ℓ/∂𝛄²], 𝔼[∂²ℓ/∂𝛈∂𝛄], 𝔼[∂²ℓ/∂𝛈²].
  expect_equal(round(z$El2[, 1], 5), El_gg)
  expect_equal(z$El2[, 2], rep(0, n))
  expect_equal(round(z$El2[, 3], 5), El_ee)

  # ∂³ℓ, ∂⁴ℓ.
  expect_null(z$l3)
  expect_null(z$l4)
})

test_that("hurdleNB_ll works for scalar y at level 1", {
  level <- 1
  z <- hurdleNB_ll(y = 0, g = 10, eta = 0, th0 = log(0.5), level = level)

  # y==0: ∂ℓ/∂θ₀, ∂²ℓ/∂𝛄∂θ₀, ∂²ℓ/∂θ₀², ∂³ℓ/∂𝛄³, ∂³ℓ/∂𝛄²∂𝛈, ∂³ℓ/∂𝛄∂𝛈², ∂³ℓ/∂𝛈³,
  # ∂³ℓ/∂𝛄²θ₀, ∂³ℓ/∂𝛄θ₀².
  expect_equal(z$l1[, 3], 0)
  expect_equal(z$l2[, 4], 0)
  expect_equal(z$l2[, 5], NaN)
  expect_equal(z$l3[, 1], 0)
  expect_equal(z$l3[, 2], 0)
  expect_equal(z$l3[, 3], 0)
  expect_equal(z$l3[, 4], -1)
  expect_equal(z$l3[, 5], 0)
  expect_equal(z$l3[, 6], NaN)

  # ∂⁴ℓ.
  expect_null(z$l4)

  z <- hurdleNB_ll(y = 25, g = 10, eta = 0, th0 = log(0.5), level = level)

  # y>0: ∂ℓ/∂θ₀, ∂²ℓ/∂𝛄∂θ₀, ∂²ℓ/∂θ₀², ∂³ℓ/∂𝛄³, ∂³ℓ/∂𝛄²∂𝛈, ∂³ℓ/∂𝛄∂𝛈², ∂³ℓ/∂𝛈³,
  # ∂³ℓ/∂𝛄²θ₀, ∂³ℓ/∂𝛄θ₀².
  expect_equal(round(z$l1[, 3], 5), 10.90750)
  expect_equal(round(z$l2[, 4], 5), 1.99737)
  expect_equal(round(z$l2[, 5], 5), NaN)
  expect_equal(round(z$l3[, 1], 5), 0.00245)
  expect_equal(z$l3[, 2], 0)
  expect_equal(z$l3[, 3], 0)
  expect_equal(round(z$l3[, 4], 5), -0.18775)
  expect_equal(round(z$l3[, 5], 5), 0.00263)
  expect_equal(round(z$l3[, 6], 5), NaN)

  # ∂⁴ℓ.
  expect_null(z$l4)
})

test_that("hurdleNB_ll works for vector y at level 1", {
  level <- 1
  g <- c(
    1.3130418, 1.1740350, 1.7447116, 1.0350755, 0.8175248, 1.2184900,
    1.3780124, 1.1322345, 1.9907406, 0.7240296, 0.8100226, 1.1046235,
    0.9087231, 1.7945166, 0.8435220, 1.4348064, 1.0780594, 0.5603018
  )
  y <- c(15, 0, 60, 0, 7, 0, 0, 1, 78, 0, 0, 4, 0, 7, 0, 1, 0, 0)
  n <- length(y)
  eta <- -2 + exp(0.3) * g
  th0 <- 2

  z <- hurdleNB_ll(y, g, eta, th0, level = level)
  zind <- y == 0

  l_eee0 <- c(
    -0.66020, -0.54728, -0.70103, -0.86946, -0.35964, -0.40390,
    -0.46146, -0.42259, -0.57997, -0.28832
  )

  # y==0: ∂ℓ/∂θ₀, ∂²ℓ/∂𝛄∂θ₀, ∂²ℓ/∂θ₀², ∂³ℓ/∂𝛄³, ∂³ℓ/∂𝛈³, ∂³ℓ/∂𝛄²θ₀, ∂³ℓ/∂𝛄θ₀².
  expect_equal(z$l1[zind, 3], rep(0, sum(zind)))
  expect_equal(z$l2[zind, 4], rep(0, sum(zind)))
  expect_equal(z$l2[zind, 5], rep(NaN, sum(zind)))
  expect_equal(z$l3[zind, 1], rep(0, sum(zind)))
  expect_equal(round(z$l3[zind, 4], 5), l_eee0)
  expect_equal(z$l3[zind, 5], rep(0, sum(zind)))
  expect_equal(z$l3[zind, 6], rep(NaN, sum(zind)))

  l_t0 <- c(
    -0.00161, 0.72745, -0.10116, -0.09930, 0.75448, -0.19927, -0.20057,
    -0.06606
  )
  l_gt0 <- c(
    -0.36428, -1.21840, -0.21635, 0.10820, -1.25893, -0.01540,
    -0.01613, 0.11123
  )
  l_ggg <- c(
    0.44515, 1.26825, 0.29758, 0.00720, 1.31895, 0.12039, 0.12148,
    0.00215
  )
  l_eee <- c(
    -0.19543, -0.11781, -0.14913, -0.18551, 0.03868, -0.18296,
    -0.09379, -0.1921
  )
  l_ggt0 <- c(
    0.45057, 1.27287, 0.30432, 0.01305, 1.32323, 0.12631, 0.12603,
    0.00732
  )

  # y>0: ∂ℓ/∂θ₀, ∂²ℓ/∂𝛄∂θ₀, ∂²ℓ/∂θ₀², ∂³ℓ/∂𝛄³, ∂³ℓ/∂𝛈³, ∂³ℓ/∂𝛄²θ₀, ∂³ℓ/∂𝛄θ₀².
  expect_equal(round(z$l1[!zind, 3], 5), l_t0)
  expect_equal(round(z$l2[!zind, 4], 5), l_gt0)
  expect_equal(round(z$l2[!zind, 5], 5), rep(NaN, sum(!zind)))
  expect_equal(round(z$l3[!zind, 1], 5), l_ggg)
  expect_equal(round(z$l3[!zind, 4], 5), l_eee)
  expect_equal(round(z$l3[!zind, 5], 5), l_ggt0)
  expect_equal(round(z$l3[!zind, 6], 5), rep(NaN, sum(!zind)))

  # ∂³ℓ/∂𝛄²∂𝛈, ∂³ℓ/∂𝛄∂𝛈².
  expect_equal(z$l3[, 2], rep(0, sum(n)))
  expect_equal(z$l3[, 3], rep(0, sum(n)))

  # ∂⁴ℓ.
  expect_null(z$l4)
})

test_that("hurdleNB_ll works for scalar y at level 2", {
  level <- 2
  z <- hurdleNB_ll(y = 0, g = 10, eta = 0, th0 = log(0.5), level = level)

  expect_false(is.null(z$l1))

  # y==0: ∂²ℓ/∂θ₀², ∂³ℓ/∂𝛄θ₀², ∂ℓ⁴/∂𝛄⁴, ∂⁴ℓ/∂𝛄³∂𝛈, ∂⁴ℓ/∂𝛄²∂𝛈², ∂⁴ℓ/∂𝛄∂𝛈³,
  # ∂ℓ⁴/∂𝛈⁴, ∂⁴ℓ/∂𝛄³∂θ₀, ∂⁴ℓ/∂𝛄²∂θ₀².
  expect_equal(z$l2[, 5], 0)
  expect_equal(z$l3[, 6], 0)
  expect_equal(z$l4[, 1], 0)
  expect_equal(z$l4[, 1], 0)
  expect_equal(z$l4[, 2], 0)
  expect_equal(z$l4[, 3], 0)
  expect_equal(z$l4[, 4], 0)
  expect_equal(z$l4[, 5], -1)
  expect_equal(z$l4[, 6], 0)
  expect_equal(z$l4[, 7], 0)

  z <- hurdleNB_ll(y = 25, g = 10, eta = 0, th0 = log(0.5), level = level)

  expect_false(is.null(z$l1))

  # y>0: ∂²ℓ/∂θ₀², ∂³ℓ/∂𝛄θ₀², ∂ℓ⁴/∂𝛄⁴, ∂⁴ℓ/∂𝛄³∂𝛈, ∂⁴ℓ/∂𝛄²∂𝛈², ∂⁴ℓ/∂𝛄∂𝛈³,
  # ∂ℓ⁴/∂𝛈⁴, ∂⁴ℓ/∂𝛄³∂θ₀, ∂⁴ℓ/∂𝛄²∂θ₀².
  expect_equal(round(z$l2[, 5], 5), -11.33667)
  expect_equal(round(z$l3[, 6], 5), -1.99701)
  expect_equal(round(z$l4[, 1], 5), -0.00245)
  expect_equal(z$l4[, 2], 0)
  expect_equal(z$l4[, 3], 0)
  expect_equal(z$l4[, 4], 0)
  expect_equal(round(z$l4[, 5], 5), 0.08452)
  expect_equal(round(z$l4[, 6], 5), -0.00263)
  expect_equal(round(z$l4[, 7], 5), -0.00299)
})

test_that("hurdleNB_ll works for vector y at level 2", {
  level <- 2
  g <- c(
    1.3130418, 1.1740350, 1.7447116, 1.0350755, 0.8175248, 1.2184900,
    1.3780124, 1.1322345, 1.9907406, 0.7240296, 0.8100226, 1.1046235,
    0.9087231, 1.7945166, 0.8435220, 1.4348064, 1.0780594, 0.5603018
  )
  y <- c(15, 0, 60, 0, 7, 0, 0, 1, 78, 0, 0, 4, 0, 7, 0, 1, 0, 0)
  n <- length(y)
  eta <- -2 + exp(0.3) * g
  th0 <- 2

  z <- hurdleNB_ll(y, g, eta, th0, level = level)
  zind <- y == 0

  expect_false(is.null(z$l1))

  l_eeee0 <- c(
    -0.66020, -0.54728, -0.70103, -0.86946, -0.35964, -0.40390,
    -0.46146, -0.42259, -0.57997, -0.28832
  )

  # y==0: ∂²ℓ/∂θ₀², ∂³ℓ/∂𝛄θ₀², ∂ℓ⁴/∂𝛄⁴, ∂ℓ⁴/∂𝛈⁴, ∂⁴ℓ/∂𝛄³∂θ₀, ∂⁴ℓ/∂𝛄²∂θ₀².
  expect_equal(z$l2[zind, 5], rep(0, sum(zind)))
  expect_equal(z$l3[zind, 6], rep(0, sum(zind)))
  expect_equal(z$l4[zind, 1], rep(0, sum(zind)))
  expect_equal(round(z$l4[zind, 5], 5), l_eeee0)
  expect_equal(z$l4[zind, 6], rep(0, sum(zind)))
  expect_equal(z$l4[zind, 7], rep(0, sum(zind)))

  l_t0t0 <- c(
    -0.15981, -0.85996, -0.07941, -0.06316, -0.88569, 0.02104,
    0.04354, -0.08402
  )
  l_gt0t0 <- c(
    0.37124, 1.18930, 0.23081, -0.06429, 1.23750, 0.04929, 0.04200,
    -0.07328
  )
  l_gggg <- c(
    -0.39412, -1.16052, -0.24334, -0.01949, -1.23220, -0.11289,
    -0.12133, -0.01404
  )
  l_eeee <- c(
    -0.00934, 0.32878, -0.09545, -0.06549, 0.60723, -0.07093, 0.38621,
    0.05373
  )
  l_gt0t0t0 <- c(
    -0.39631, -1.16207, -0.24652, -0.02201, -1.23346, -0.11547,
    -0.12281, -0.01603
  )
  l_ggt0t0 <- c(
    -0.40693, -1.17118, -0.25950, -0.03342, -1.24189, -0.12700,
    -0.13177, -0.02618
  )

  # y>0: ∂²ℓ/∂θ₀², ∂³ℓ/∂𝛄θ₀², ∂ℓ⁴/∂𝛄⁴, ∂ℓ⁴/∂𝛈⁴, ∂⁴ℓ/∂𝛄³∂θ₀, ∂⁴ℓ/∂𝛄²∂θ₀².
  expect_equal(round(z$l2[!zind, 5], 5), l_t0t0)
  expect_equal(round(z$l3[!zind, 6], 5), l_gt0t0)
  expect_equal(round(z$l4[!zind, 1], 5), l_gggg)
  expect_equal(round(z$l4[!zind, 5], 5), l_eeee)
  expect_equal(round(z$l4[!zind, 6], 5), l_gt0t0t0)
  expect_equal(round(z$l4[!zind, 7], 5), l_ggt0t0)

  # ∂⁴ℓ/∂𝛄³∂𝛈, ∂⁴ℓ/∂𝛄²∂𝛈², ∂⁴ℓ/∂𝛄∂𝛈³.
  expect_equal(z$l4[, 2], rep(0, n))
  expect_equal(z$l4[, 3], rep(0, n))
  expect_equal(z$l4[, 4], rep(0, n))
})

test_that("ktlg works", {
  g <- c(393.54581, 4.76799, -63.07129)
  a <- exp(2)

  # test 𝛋, 𝛕.
  v <- ktlg(g, a)

  ind <- c(FALSE, FALSE, TRUE)
  ii <- c(TRUE, FALSE, FALSE)
  k <- c(1 / exp(2), 0.13518, exp(g[3]))
  tau <- c(1, 1.66687, 1 / exp(g[3]))

  expect_null(v$lg)
  expect_equal(v$ind, ind)
  expect_equal(v$ii, ii)

  v$k[2] <- round(v$k[2], 5)
  v$tau[2] <- round(v$tau[2], 5)
  expect_equal(v$k, k)
  expect_equal(v$tau, tau)

  # test lg, log(1+αexp(𝛄)).
  v <- ktlg(g, a, what = c("lg"))

  lg <- c(g[1], 6.76914, 0)

  expect_null(v$k)
  expect_null(v$tau)

  v$lg[2] <- round(v$lg[2], 5)
  expect_equal(v$lg, lg)
})
