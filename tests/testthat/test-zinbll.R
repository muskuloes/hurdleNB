test_that("zinbll works for scalar y at level=0", {
  # for y==0: ℓ, ∂ℓ/∂𝛄, ∂ℓ/∂𝛈, ∂ℓ/∂θ⁰, ∂²ℓ/∂𝛄², ∂²ℓ/∂𝛈∂𝛄, ∂²ℓ/∂𝛈², ∂²ℓ/∂𝛄∂θ₀,
  # ∂²ℓ/∂θ₀².
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

  q <- 1 - exp(-exp(0))
  El_eta2 <- (1 - q) * (-exp(0)) + q * -0.33870

  # 𝔼[∂²ℓ/∂𝛄²], 𝔼[∂²ℓ/∂𝛈∂𝛄], 𝔼[∂²ℓ/∂𝛈²].
  expect_equal(round(z$El2[, 1], 5), -1.26413)
  expect_equal(z$El2[, 2], 0)
  expect_equal(round(z$El2[, 3], 5), round(El_eta2, 5))

  # ∂³ℓ, ∂⁴ℓ.
  expect_equal(z$l3, NULL)
  expect_equal(z$l4, NULL)

  # for y>0: ℓ, ∂ℓ/∂𝛄, ∂ℓ/∂𝛈, ∂ℓ/∂θ⁰, ∂²ℓ/∂𝛄², ∂²ℓ/∂𝛈∂𝛄, ∂²ℓ/∂𝛈², ∂²ℓ/∂𝛄∂θ₀,
  # ∂²ℓ/∂θ₀².
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

  # 𝔼[∂²ℓ/∂𝛄²], 𝔼[∂²ℓ/∂𝛈∂𝛄], 𝔼[∂²ℓ/∂𝛈²].
  expect_equal(round(z$El2[, 1], 5), -1.26413)
  expect_equal(z$El2[, 2], 0)
  expect_equal(round(z$El2[, 3], 5), round(El_eta2, 5))

  # ∂³ℓ, ∂⁴ℓ.
  expect_equal(z$l3, NULL)
  expect_equal(z$l4, NULL)
})

test_that("zinbll works for vector y at level=0", {
  # example generated using:
  #> x0 <- runif(18)
  #> x1 <- runif(18)
  #> x2 <- runif(18)
  #> x3 <- runif(18)
  #> g <- x0 + 2 * x1 + 3 * x2 + 4 * x3
  #> y <- rzinb(g, theta = c(-2, .3, 2), b = 0)
  # the expectations are then gotten using sympy.
  g <- c(
    1.3130418, 1.1740350, 1.7447116, 1.0350755, 0.8175248, 1.2184900,
    1.3780124, 1.1322345, 1.9907406, 0.7240296, 0.8100226, 1.1046235,
    0.9087231, 1.7945166, 0.8435220, 1.4348064, 1.0780594, 0.5603018
  )

  n <- 18
  y <- c(15, 0, 60, 0, 7, 0, 0, 1, 78, 0, 0, 4, 0, 7, 0, 1, 0, 0)
  eta <- -2 + exp(.3) * g
  th0 <- 2
  z <- zinbll(y, g, eta, th0, level = 0)
  zind <- y == 0

  # for y==0: ℓ, ∂ℓ/∂𝛄, ∂ℓ/∂𝛈, ∂²ℓ/∂𝛄², ∂²ℓ/∂𝛈².
  expect_equal(z$l[zind], -exp(eta[zind]))
  expect_equal(z$l1[zind, 1], rep(0, sum(zind)))
  expect_equal(z$l1[zind, 2], -exp(eta[zind]))
  expect_equal(z$l2[zind, 1], rep(0, sum(zind)))
  expect_equal(z$l2[zind, 3], -exp(eta[zind]))

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

  # for y>0: ℓ, ∂ℓ/∂𝛄, ∂ℓ/∂𝛈, ∂²ℓ/∂𝛄², ∂²ℓ/∂𝛈².
  expect_equal(round(z$l[!zind], 5), l)
  expect_equal(round(z$l1[!zind, 1], 5), l1_g)
  expect_equal(round(z$l1[!zind, 2], 5), l1_eta)
  expect_equal(round(z$l2[!zind, 1], 5), l2_g)
  expect_equal(round(z$l2[!zind, 3], 5), l2_eta)

  # ∂ℓ/∂θ₀, ∂²ℓ/∂𝛈∂𝛄, ∂²ℓ/∂𝛄∂θ₀, ∂²ℓ/∂θ₀².
  expect_equal(z$l1[, 3], rep(NaN, n))
  expect_equal(z$l2[, 2], rep(0, n))
  expect_equal(z$l2[, 4], rep(NaN, n))
  expect_equal(z$l2[, 5], rep(NaN, n))

  El_g2 <- c(
    -0.15194, -0.13553, -0.20149, -0.11977, -0.09709, -0.14073, -0.15971,
    -0.13070, -0.22336, -0.08825, -0.09636, -0.12755, -0.10626, -0.20648,
    -0.09965, -0.16649, -0.12455, -0.07415
  )
  q <- 1 - exp(-exp(eta))
  nl2_eta <- rep(0, n)
  nl2_eta[zind] <- c(
    -0.25849, -0.22421, -0.26993, -0.31183, -0.15835, -0.17491, -0.19549,
    -0.18171, -0.23455, -0.13035
  )
  nl2_eta[!zind] <- l2_eta
  El_eta2 <- (1 - q) * (-exp(eta)) + q * nl2_eta

  # 𝔼[∂²ℓ/∂𝛄²], 𝔼[∂²ℓ/∂𝛈∂𝛄], 𝔼[∂²ℓ/∂𝛈²].
  expect_equal(round(z$El2[, 1], 5), El_g2)
  expect_equal(z$El2[, 2], rep(0, n))
  expect_equal(round(z$El2[, 3], 4), round(El_eta2, 4))

  # ∂³ℓ, ∂⁴ℓ.
  expect_equal(z$l3, NULL)
  expect_equal(z$l4, NULL)
})
