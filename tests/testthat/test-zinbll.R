test_that("zinbll works for scalar y at level=0", {
  # for y == 0
  z <- zinbll(y = 0, g = 10, eta = 0, th0 = log(.5), level = 0)
  expect_equal(z$l, -1)
  expect_equal(z$l1[, 1], 0)
  expect_equal(z$l1[, 2], -1)
  expect_equal(z$l1[, 3], NaN)
  expect_equal(z$l2[, 1], 0)
  expect_equal(z$l2[, 2], 0)
  expect_equal(z$l2[, 3], -1)
  expect_equal(z$l2[, 4], NaN)
  expect_equal(z$l2[, 5], NaN)

  # for y > 0
  z <- zinbll(y = 25, g = 10, eta = 0, th0 = log(.5), level = 0)
  expect_equal(round(z$l, 5), -15.81674)
  expect_equal(round(z$l1[, 1], 5), -1.99755)
  expect_equal(round(z$l1[, 2], 5), 0.58198)
  expect_equal(z$l1[, 3], NaN)
  expect_equal(round(z$l2[, 1], 5), -0.00245)
  expect_equal(z$l2[, 2], 0)
  expect_equal(round(z$l2[, 3], 5), -0.33870)
  expect_equal(z$l2[, 4], NaN)
  expect_equal(z$l2[, 5], NaN)
})

test_that("zinbll works for vector y at level=0", {
  # example generated using:
  #> x0 <- runif(18)
  #> x1 <- runif(18)
  #> x2 <- runif(18)
  #> x3 <- runif(18)
  #> g <- x0 + 2 * x1 + 3 * x2 + 4 * x3
  #> y <- rzinb(g, theta = c(-2, .3, 2), b = 0)
  g <- c(
    1.3130418, 1.1740350, 1.7447116, 1.0350755, 0.8175248, 1.2184900,
    1.3780124, 1.1322345, 1.9907406, 0.7240296, 0.8100226, 1.1046235,
    0.9087231, 1.7945166, 0.8435220, 1.4348064, 1.0780594, 0.5603018
  )
  y <- c(15, 0, 60, 0, 7, 0, 0, 1, 78, 0, 0, 4, 0, 7, 0, 1, 0, 0)
  eta <- -2 + exp(.3) * g
  th0 <- 2
  z <- zinbll(y, g, eta, th0, level = 0)
  zind <- y == 0
  expect_equal(z$l[zind], -exp(eta[zind]))
  l <- c(
    -4.86116, -6.74632, -4.38402, -2.18796, -6.95198, -3.48723, -3.63684,
    -2.01474
  )
  l1_g <- c(
    0.16857, 1.05483, -0.00144, -0.32947, 1.09841, -0.20166, -0.17410,
    -0.31890
  )
  l1_eta <- c(
    0.65408, 0.45088, 0.80983, 0.72025, 0.31547, 0.72936, 0.42404,
    0.60301
  )
  l2_g <- c(
    -0.43937, -1.29564, -0.28844, 0.03413, -1.33727, -0.08930, -0.09360, 0.03550
  )
  l2_eta <- c(
    -0.29469, -0.39553, -0.17641, -0.24793, -0.41127, -0.24106, -0.40266,
    -0.32668
  )
  expect_equal(round(z$l[!zind], 5), l)
  expect_equal(round(z$l1[!zind, 1], 5), l1_g)
  expect_equal(round(z$l1[!zind, 2], 5), l1_eta)
  expect_equal(z$l1[, 3], rep(NaN, 18))
  expect_equal(round(z$l2[!zind, 1], 5), l2_g)
  expect_equal(z$l2[, 2], rep(0, 18))
  expect_equal(round(z$l2[!zind, 3], 5), l2_eta)
  expect_equal(z$l2[, 4], rep(NaN, 18))
  expect_equal(z$l2[, 5], rep(NaN, 18))
})
