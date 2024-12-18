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
  expect_equal(round(oo$Dth2[, 3], 5), 0)
  expect_equal(round(oo$Dth2[, 4], 5), -0.00041)
  expect_equal(round(oo$Dth2[, 5], 5), 0)
  expect_equal(round(oo$Dth2[, 6], 5), 2.02972)
  expect_equal(round(oo$Dmuth2[, 1], 5), 0.11576)
  expect_equal(round(oo$Dmuth2[, 2], 5), 0.00087)
  expect_equal(round(oo$Dmuth2[, 3], 5), 0)
  expect_equal(round(oo$Dmuth2[, 4], 5), -0.00054)
  expect_equal(round(oo$Dmuth2[, 5], 5), 0)
  expect_equal(round(oo$Dmuth2[, 6], 5), -2.01351)
  expect_equal(round(oo$Dmu3th[, 1], 5), -0.28299)
  expect_equal(round(oo$Dmu3th[, 2], 5), 0.00010)
  expect_equal(round(oo$Dmu3th[, 3], 5), -0.23094)
  expect_equal(round(oo$Dmu2th2[, 1], 5), -0.40369)
  expect_equal(round(oo$Dmu2th2[, 2], 5), 0.00003)
  expect_equal(round(oo$Dmu2th2[, 3], 5), 0)
  expect_equal(round(oo$Dmu2th2[, 4], 5), 0.00118)
  expect_equal(round(oo$Dmu2th2[, 5], 5), 0)
  expect_equal(round(oo$Dmu2th2[, 6], 5), 0.01762)
  expect_equal(round(oo$Dmu4, 5), -0.45853)
})

test_that("Dd works for vector y at level 2", {
  model <- zinb()
  g <- c(
    1.3130418, 1.1740350, 1.7447116, 1.0350755, 0.8175248, 1.2184900,
    1.3780124, 1.1322345, 1.9907406, 0.7240296, 0.8100226, 1.1046235,
    0.9087231, 1.7945166, 0.8435220, 1.4348064, 1.0780594, 0.5603018
  )
  y <- c(15, 0, 60, 0, 7, 0, 0, 1, 78, 0, 0, 4, 0, 7, 0, 1, 0, 0)
  n <- length(y)
  zind <- y == 0

  oo <- model$Dd(
    y = y, g = g, theta = c(-2, 0.3, 2), wt = 1,
    level = 2
  )

  D0mu <- c(
    1.78234, 1.4775, 1.89257, 2.34731, 0.97092, 1.09042, 1.24582,
    1.14086, 1.56577, 0.77839
  )
  D0mu2 <- c(
    2.40591, 1.99442, 2.5547, 3.16854, 1.3106, 1.47191, 1.68168,
    1.54, 2.11356, 1.05072
  )
  D0th1 <- c(
    1.32039, 1.09456, 1.40205, 1.73893, 0.71927, 0.8078, 0.92293,
    0.84517, 1.15995, 0.57665
  )
  D0th2 <- c(
    2.09253, 1.52933, 2.30608, 3.23462, 0.70297, 0.88327, 1.13211,
    0.96234, 1.68799, 0.43613
  )
  D0muth1 <- c(
    1.78234, 1.4775, 1.89257, 2.34731, 0.97092, 1.09042, 1.24582,
    1.14086, 1.56577, 0.77839
  )
  D0muth2 <- c(
    4.60697, 3.54188, 5.00545, 6.71359, 1.91983, 2.2827, 2.77401,
    2.43988, 3.84431, 1.36711
  )
  D0mu2th1 <- c(
    2.40591, 1.99442, 2.5547, 3.16854, 1.3106, 1.47191, 1.68168,
    1.54, 2.11356, 1.05072
  )
  D0mu2th2 <- c(
    8.62467, 6.77545, 9.31136, 12.23093, 3.9021, 4.55324,
    5.4262, 4.8335, 7.30284, 2.89613
  )
  D0mu3 <- c(
    3.24764, 2.69218, 3.44849, 4.27707, 1.76912, 1.98688, 2.27003,
    2.07878, 2.85301, 1.41832
  )
  D0th1th1 <- c(
    1.32039, 1.09456, 1.40205, 1.73893, 0.71927, 0.8078, 0.92293,
    0.84517, 1.15995, 0.57665
  )
  D0th1th2 <- c(
    2.09253, 1.52933, 2.30608, 3.23462, 0.70297, 0.88327, 1.13211,
    0.96234, 1.68799, 0.43613
  )
  D0th2th2 <- c(
    5.40874, 3.66611, 6.0991, 9.25141, 1.39001, 1.84904, 2.5208,
    2.0581, 4.1444, 0.766
  )
  D0muth1th1 <- c(
    1.78234, 1.4775, 1.89257, 2.34731, 0.97092, 1.09042, 1.24582,
    1.14086, 1.56577, 0.77839
  )
  D0muth1th2 <- c(
    4.60697, 3.54188, 5.00545, 6.71359, 1.91983, 2.2827, 2.77401,
    2.43988, 3.84431, 1.36711
  )
  D0muth2th2 <- c(
    14.73263, 10.55498, 16.35125, 23.56797, 4.74506, 5.97093,
    7.70492, 6.51705, 11.71721, 2.98982
  )
  D0mu3th1 <- c(
    3.24764, 2.69218, 3.44849, 4.27707, 1.76912, 1.98688, 2.27003,
    2.07878, 2.85301, 1.41832
  )
  D0mu3th2 <- c(
    14.88972, 11.83809, 16.01751, 20.78711, 7.03641, 8.13311,
    9.59464, 8.60333, 12.71082, 5.32769
  )
  D0mu2th1th1 <- c(
    2.40591, 1.99442, 2.5547, 3.16854, 1.3106, 1.47191, 1.68168,
    1.54, 2.11356, 1.05072
  )
  D0mu2th1th2 <- c(
    8.62467, 6.77545, 9.31136, 12.23093, 3.9021, 4.55324,
    5.4262, 4.8335, 7.30284, 2.89613
  )
  D0mu2th2th2 <- c(
    34.73038, 25.80422, 38.1399, 53.10676, 12.89876, 15.69448,
    19.57127, 16.9241, 28.3087, 8.77737
  )
  D0mu4 <- c(
    4.38386, 3.63407, 4.65498, 5.77345, 2.38807, 2.682, 3.06423,
    2.80606, 3.85116, 1.91454
  )

  expect_equal(round(oo$Dmu[zind], 5), D0mu)
  expect_equal(round(oo$Dmu2[zind], 5), D0mu2)
  expect_equal(round(oo$Dth[zind, 1], 5), D0th1)
  expect_equal(round(oo$Dth[zind, 2], 5), D0th2)
  expect_equal(round(oo$Dth[zind, 3], 5), rep(0, sum(zind)))
  expect_equal(round(oo$Dmuth[zind, 1], 5), D0muth1)
  expect_equal(round(oo$Dmuth[zind, 2], 5), D0muth2)
  expect_equal(round(oo$Dmuth[zind, 3], 5), rep(0, sum(zind)))
  expect_equal(round(oo$Dmu2th[zind, 1], 5), D0mu2th1)
  expect_equal(round(oo$Dmu2th[zind, 2], 5), D0mu2th2)
  expect_equal(round(oo$Dmu2th[zind, 3], 5), rep(0, sum(zind)))
  expect_equal(round(oo$Dmu3[zind], 5), D0mu3)
  expect_equal(round(oo$Dth2[zind, 1], 5), D0th1th1)
  expect_equal(round(oo$Dth2[zind, 2], 5), D0th1th2)
  expect_equal(round(oo$Dth2[zind, 4], 5), D0th2th2)
  expect_equal(round(oo$Dth2[zind, 6], 5), rep(0, sum(zind)))
  expect_equal(round(oo$Dmuth2[zind, 1], 5), D0muth1th1)
  expect_equal(round(oo$Dmuth2[zind, 2], 5), D0muth1th2)
  expect_equal(round(oo$Dmuth2[zind, 4], 5), D0muth2th2)
  expect_equal(round(oo$Dmuth2[zind, 6], 5), rep(0, sum(zind)))
  expect_equal(round(oo$Dmu3th[zind, 1], 5), D0mu3th1)
  expect_equal(round(oo$Dmu3th[zind, 2], 5), D0mu3th2)
  expect_equal(round(oo$Dmu3th[zind, 3], 5), rep(0, sum(zind)))
  expect_equal(round(oo$Dmu2th2[zind, 1], 5), D0mu2th1th1)
  expect_equal(round(oo$Dmu2th2[zind, 2], 5), D0mu2th1th2)
  expect_equal(round(oo$Dmu2th2[zind, 4], 5), D0mu2th2th2)
  expect_equal(round(oo$Dmu2th2[zind, 6], 5), rep(0, sum(zind)))
  expect_equal(round(oo$Dmu4[zind], 5), D0mu4)

  Dmu <- c(
    -2.10297, -3.32691, -2.18342, -1.28552, -3.04851, -1.56574,
    -0.7966, -0.99015
  )
  Dmu2 <- c(
    1.95266, 4.03268, 1.21977, 0.83525, 4.1733, 1.05709, 1.6546,
    1.1195
  )
  Dth1 <- c(
    -1.30816, -0.90177, -1.61966, -1.4405, -0.63095, -1.45872,
    -0.84809, -1.20602
  )
  Dth2 <- c(
    -2.31862, -2.12376, -1.78736, -2.20159, -1.69549, -2.17508,
    -2.05435, -2.3358
  )
  Dth0 <- c(
    0.00322, -1.4549, 0.20231, 0.19859, -1.50895, 0.39855,
    0.40115, 0.13211
  )
  Dmuth1 <- c(
    0.79558, 1.06781, 0.47627, 0.66933, 1.1103, 0.6508, 1.08707,
    0.88195
  )
  Dmuth2 <- c(
    -0.35573, 1.29756, -1.66073, -0.92149, 2.13194, -0.99868,
    1.48846, 0.08019
  )
  Dmuth0 <- c(
    0.72856, 2.4368, 0.43269, -0.2164, 2.51787, 0.03079,
    0.03226, -0.22247
  )
  Dmu2th1 <- c(
    0.71218, 0.42933, 0.54345, 0.67603, -0.14094, 0.66675,
    0.34179, 0.70005
  )
  Dmu2th2 <- c(
    3.41014, 3.89391, 1.88551, 2.84022, 2.61877, 2.75115, 3.76271,
    3.73686
  )
  Dmu2th0 <- c(
    -0.90114, -2.54574, -0.60864, -0.02609, -2.64646, -0.25262,
    -0.25205, -0.01464
  )
  Dmu3 <- c(
    0.07105, -1.95696, 0.13843, 0.89814, -2.82815, 0.65924,
    0.21841, 0.94067
  )
  Dth1th1 <- c(
    0.58938, 0.79105, 0.35283, 0.49586, 0.82253, 0.48212, 0.80532,
    0.65336
  )
  Dth1th2 <- c(
    1.04463, 1.86302, 0.38936, 0.75784, 2.21033, 0.71888, 1.95076,
    1.26542
  )
  Dth2th2 <- c(
    -0.46709, 2.26386, -1.35768, -1.04334, 4.24414, -1.10316,
    2.67107, 0.11505
  )
  Dth0th0 <- c(
    0.31961, 1.71991, 0.15881, 0.12631, 1.77139, -0.04208,
    -0.08707, 0.16805
  )
  Dmuth1th1 <- c(
    0.5276, 0.31806, 0.4026, 0.50081, -0.10441, 0.49394,
    0.2532, 0.51861
  )
  Dmuth1th2 <- c(
    1.73071, 1.81687, 0.92055, 1.43476, 0.82973, 1.38731, 1.70042,
    1.88639
  )
  Dmuth2th2 <- c(
    4.12192, 8.09131, -0.11928, 2.29431, 7.34522, 2.04031,
    8.24071, 5.44186
  )
  Dmuth0th0 <- c(
    -0.74248, -2.3786, -0.46162, 0.12859, -2.47501, -0.09858,
    -0.08401, 0.14655
  )
  Dmu3th1 <- c(
    0.04596, -1.61736, 0.46956, 0.32214, -2.98707, 0.34893,
    -1.89986, -0.26433
  )
  Dmu3th2 <- c(
    2.9655, -2.07044, 2.71891, 3.22997, -8.59766, 3.22036,
    -3.21801, 2.32297
  )
  Dmu3th0 <- c(
    0.79262, 2.32414, 0.49305, 0.04403, 2.46691, 0.23094, 0.24562,
    0.03206
  )
  Dmu2th1th1 <- c(
    0.03405, -1.19817, 0.34786, 0.23865, -2.21288, 0.2585,
    -1.40745, -0.19582
  )
  Dmu2th1th2 <- c(
    1.48471, -1.96315, 1.47077, 1.71679, -6.22836, 1.71895,
    -2.72575, 1.02084
  )
  Dmu2th2th2 <- c(
    10.7141, 4.17551, 5.99379, 9.33753, -11.87818, 9.05959,
    1.75066, 10.80673
  )
  Dmu2th0th0 <- c(
    0.81387, 2.34236, 0.51901, 0.06683, 2.48379, 0.254, 0.26355,
    0.05236
  )
  Dmu4 <- c(
    0.85028, 0.13785, 1.12053, 0.47383, -1.56772, 0.6968,
    -2.32188, -0.32873
  )

  expect_equal(round(oo$Dmu[!zind], 5), Dmu)
  expect_equal(round(oo$Dmu2[!zind], 5), Dmu2)
  expect_equal(round(oo$Dth[!zind, 1], 5), Dth1)
  expect_equal(round(oo$Dth[!zind, 2], 5), Dth2)
  expect_equal(round(oo$Dth[!zind, 3], 5), Dth0)
  expect_equal(round(oo$Dmuth[!zind, 1], 5), Dmuth1)
  expect_equal(round(oo$Dmuth[!zind, 2], 5), Dmuth2)
  expect_equal(round(oo$Dmuth[!zind, 3], 5), Dmuth0)
  expect_equal(round(oo$Dmu2th[!zind, 1], 5), Dmu2th1)
  expect_equal(round(oo$Dmu2th[!zind, 2], 5), Dmu2th2)
  expect_equal(round(oo$Dmu2th[!zind, 3], 5), Dmu2th0)
  expect_equal(round(oo$Dmu3[!zind], 5), Dmu3)
  expect_equal(round(oo$Dth2[!zind, 1], 5), Dth1th1)
  expect_equal(round(oo$Dth2[!zind, 2], 5), Dth1th2)
  expect_equal(round(oo$Dth2[!zind, 4], 5), Dth2th2)
  expect_equal(round(oo$Dth2[!zind, 6], 5), Dth0th0)
  expect_equal(round(oo$Dmuth2[!zind, 1], 5), Dmuth1th1)
  expect_equal(round(oo$Dmuth2[!zind, 2], 5), Dmuth1th2)
  expect_equal(round(oo$Dmuth2[!zind, 4], 5), Dmuth2th2)
  expect_equal(round(oo$Dmuth2[!zind, 6], 5), Dmuth0th0)
  expect_equal(round(oo$Dmu3th[!zind, 1], 5), Dmu3th1)
  expect_equal(round(oo$Dmu3th[!zind, 2], 5), Dmu3th2)
  expect_equal(round(oo$Dmu3th[!zind, 3], 5), Dmu3th0)
  expect_equal(round(oo$Dmu2th2[!zind, 1], 5), Dmu2th1th1)
  expect_equal(round(oo$Dmu2th2[!zind, 2], 5), Dmu2th1th2)
  expect_equal(round(oo$Dmu2th2[!zind, 4], 5), Dmu2th2th2)
  expect_equal(round(oo$Dmu2th2[!zind, 6], 5), Dmu2th0th0)
  expect_equal(round(oo$Dmu4[!zind], 5), Dmu4)

  EDmu2 <- c(
    2.20235, 1.96954, 2.74663, 1.73773, 1.39831, 2.04448, 2.30763,
    1.89919, 2.73244, 1.26552, 1.38734, 1.85294, 1.53592, 2.77038,
    1.43676, 2.39588, 1.80868, 1.05482
  )
  expect_equal(round(oo$EDmu2, 5), EDmu2)

  # dθ₁dθ₀, dθ₂dθ₀.
  expect_equal(oo$Dth2[, 3], rep(0, n))
  expect_equal(oo$Dth2[, 5], rep(0, n))
  expect_equal(oo$Dmuth2[, 3], rep(0, n))
  expect_equal(oo$Dmuth2[, 5], rep(0, n))
  expect_equal(oo$Dmu2th2[, 3], rep(0, n))
  expect_equal(oo$Dmu2th2[, 5], rep(0, n))
})
