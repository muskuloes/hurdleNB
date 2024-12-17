test_that("zinb is the required structure", {
  model <- zinb()

  expect_s3_class(model, "extended.family")
  expect_s3_class(model, "family")
  expect_equal(model$family, "zero-inflated negative binomial")
  expect_type(model$link, "character")
  expect_type(model$linkfun, "closure")
  expect_type(model$linkinv, "closure")
  expect_type(model$dev.resids, "closure")
  expect_type(model$Dd, "closure")
  expect_type(model$rd, "closure")
  expect_type(model$residuals, "closure")
  expect_type(model$aic, "closure")
  expect_type(model$mu.eta, "closure")
  expect_null(model$g2g)
  expect_null(model$g3g)
  expect_null(model$g4g)
  expect_type(model$initialize, "expression")
  expect_type(model$postproc, "closure")
  expect_type(model$ls, "closure")
  expect_true(model$no.r.sq)
  expect_type(model$validmu, "closure")
  expect_type(model$valideta, "closure")
  expect_type(model$n.theta, "double")
  expect_equal(model$n.theta, 3)
  expect_type(model$predict, "closure")
  expect_vector(model$ini.theta, ptype = double(), size = 3)
  expect_type(model$putTheta, "closure")
  expect_type(model$getTheta, "closure")
  expect_type(model$saturated.ll, "closure")
})

test_that("lind works", {
  # level 0
  r <- lind(5, c(1, 2, 0), level = 0, b = 0)

  expect_equal(round(r$eta, 5), 37.94528)
  expect_equal(round(r$eta_g, 5), 7.38906)
  expect_equal(r$eta_gg, 0)
  expect_null(r$eta_th)
  expect_null(r$eta_gth)
  expect_null(r$eta_ggth)
  expect_null(r$eta_gggth)
  expect_null(r$eta_ggg)
  expect_null(r$eta_gggg)
  expect_null(r$eta_th2)
  expect_null(r$eta_gth2)
  expect_null(r$eta_ggth2)

  # level 1
  r <- lind(25, c(8, 2.4, 6), level = 1, b = 2)

  expect_equal(round(r$eta, 5), 333.57941)
  expect_equal(round(r$eta_g, 5), 13.02318)
  expect_equal(r$eta_gg, 0)
  expect_equal(r$eta_th[, 1], 1)
  expect_equal(round(r$eta_th[, 2], 5), 275.57941)
  expect_equal(r$eta_gth[, 1], 0)
  expect_equal(round(r$eta_gth[, 2], 5), 11.02318)
  expect_equal(r$eta_ggth[, 1], 0)
  expect_equal(r$eta_ggth[, 2], 0)
  expect_equal(r$eta_gggth[, 1], 0)
  expect_equal(r$eta_gggth[, 2], 0)
  expect_equal(r$eta_ggg, 0)
  expect_equal(r$eta_gggg, 0)
  expect_equal(r$eta_th2[, 1], 0)
  expect_equal(r$eta_th2[, 2], 0)
  expect_equal(round(r$eta_th2[, 3], 5), 275.57941)
  expect_equal(r$eta_gth2[, 1], 0)
  expect_equal(r$eta_gth2[, 2], 0)
  expect_equal(round(r$eta_gth2[, 3], 5), 11.02318)
  expect_equal(r$eta_ggth2[, 1], 0)
  expect_equal(r$eta_ggth2[, 2], 0)
  expect_equal(r$eta_ggth2[, 3], 0)
})

test_that("dev.resids works", {
  model <- zinb(b = 0.1)
  dr <- model$dev.resids(25, 10, 1, theta = c(0, log(0.02), log(0.5)))

  expect_equal(round(dr, 5), 30.78976)
})

test_that("aic works", {
  g <- c(
    1.3130418, 1.1740350, 1.7447116, 1.0350755, 0.8175248, 1.2184900,
    1.3780124, 1.1322345, 1.9907406, 0.7240296, 0.8100226, 1.1046235,
    0.9087231, 1.7945166, 0.8435220, 1.4348064, 1.0780594, 0.5603018
  )
  y <- c(15, 0, 60, 0, 7, 0, 0, 1, 78, 0, 0, 4, 0, 7, 0, 1, 0, 0)

  model <- zinb()
  aic <- model$aic(y, g, theta = c(-2, 0.3, 2), 1)

  expect_equal(round(aic, 5), 79.12820)
})

test_that("ls works", {
  model <- zinb()
  ls <- model$ls(0, 1, c(0, 0, 0), scale = 1)

  expect_equal(ls$ls, 0)
  expect_equal(ls$lsth1, c(0, 0, 0))
  expect_equal(ls$LSTH1, matrix(0, 1, 3))
  expect_equal(ls$lsth2, matrix(0, 3, 3))
})

test_that("Dd works at level 0", {
  model <- zinb(b = 0)
  oo <- model$Dd(
    y = 25, g = 10, theta = c(0, 0, log(0.5)), wt = 1,
    level = 0
  )

  expect_equal(round(oo$Dmu, 5), 3.99510)
  expect_equal(round(oo$Dmu2, 5), 0.00490)
  expect_equal(round(oo$EDmu2, 5), 3.99964)
  expect_null(oo$Dth)
  expect_null(oo$Dmuth)
  expect_null(oo$Dmu2th)
  expect_null(oo$Dmu3)
  expect_null(oo$Dth2)
  expect_null(oo$Dmuth2)
  expect_null(oo$Dmu3th)
  expect_null(oo$Dmu2th2)
  expect_null(oo$Dmu4)
})

test_that("Dd works at level 1", {
  model <- zinb(b = 0.3)
  oo <- model$Dd(
    y = 25, g = 0.5, theta = c(8, log(8), log(0.5)), wt = 1,
    level = 1
  )

  expect_equal(round(oo$Dmu, 5), -24.82311)
  expect_equal(round(oo$Dmu2, 5), 12.79751)
  expect_equal(round(oo$EDmu2, 5), 1.58087)
  expect_equal(round(oo$Dth[, 1], 5), 0)
  expect_equal(round(oo$Dth[, 2], 5), 0)
  expect_equal(round(oo$Dth[, 3], 5), -16.84326)
  expect_equal(round(oo$Dmuth[, 1], 5), 0)
  expect_equal(round(oo$Dmuth[, 2], 5), 0)
  expect_equal(round(oo$Dmuth[, 3], 5), 11.54815)
  expect_equal(round(oo$Dmu2th[, 1], 5), 0)
  expect_equal(round(oo$Dmu2th[, 2], 5), 0)
  expect_equal(round(oo$Dmu2th[, 3], 5), 0.44347)
  expect_equal(round(oo$Dmu3, 5), 1.36434)
  expect_null(oo$Dth2)
  expect_null(oo$Dmuth2)
  expect_null(oo$Dmu3th)
  expect_null(oo$Dmu2th2)
  expect_null(oo$Dmu4)
})

test_that("Dd works at level 2", {
  model <- zinb(b = 0.7)
  oo <- model$Dd(
    y = 12, g = 0.5, theta = c(0.1, log(0.001), log(2)), wt = 1,
    level = 2
  )

  expect_equal(round(oo$Dmu, 5), -4.68089)
  expect_equal(round(oo$Dmu2, 5), 4.49828)
  expect_equal(round(oo$EDmu2, 5), 1.39044)
  expect_equal(round(oo$Dth[, 1], 5), -0.82534)
  expect_equal(round(oo$Dth[, 2], 5), -0.00041)
  expect_equal(round(oo$Dth[, 3], 5), -2.47044)
  expect_equal(round(oo$Dmuth[, 1], 5), 0.56801)
  expect_equal(round(oo$Dmuth[, 2], 5), -0.00054)
  expect_equal(round(oo$Dmuth[, 3], 5), 3.62484)
  expect_equal(round(oo$Dmu2th[, 1], 5), 0.08115)
  expect_equal(round(oo$Dmu2th[, 2], 5), 0.00118)
  expect_equal(round(oo$Dmu2th[, 3], 5), -2.42007)
  expect_equal(round(oo$Dmu3, 5), -2.20606)
  expect_equal(round(oo$Dth2[, 1], 5), 0.81029)
  expect_equal(round(oo$Dth2[, 2], 5), 0.00041)
  expect_equal(round(oo$Dth2[, 3], 5), -0.00041)
  expect_equal(round(oo$Dth2[, 4], 5), 0)
  expect_equal(round(oo$Dth2[, 5], 5), 0)
  expect_equal(round(oo$Dth2[, 6], 5), 2.02972)
  expect_equal(round(oo$Dmuth2[, 1], 5), 0.11576)
  expect_equal(round(oo$Dmuth2[, 2], 5), 0.00087)
  expect_equal(round(oo$Dmuth2[, 3], 5), -0.00054)
  expect_equal(round(oo$Dmuth2[, 4], 5), 0)
  expect_equal(round(oo$Dmuth2[, 5], 5), 0)
  expect_equal(round(oo$Dmuth2[, 6], 5), -2.01351)
  expect_equal(round(oo$Dmu3th[, 1], 5), -0.28299)
  expect_equal(round(oo$Dmu3th[, 2], 5), 0.00010)
  expect_equal(round(oo$Dmu3th[, 3], 5), -0.23094)
  expect_equal(round(oo$Dmu2th2[, 1], 5), -0.40369)
  expect_equal(round(oo$Dmu2th2[, 2], 5), 0.00003)
  expect_equal(round(oo$Dmu2th2[, 3], 5), 0.00118)
  expect_equal(round(oo$Dmu2th2[, 4], 5), 0)
  expect_equal(round(oo$Dmu2th2[, 5], 5), 0)
  expect_equal(round(oo$Dmu2th2[, 6], 5), 0.01762)
  expect_equal(round(oo$Dmu4, 5), -0.45853)
})
