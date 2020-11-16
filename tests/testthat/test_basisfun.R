library(lgpr)

context("kernel decompositions and fits using basis functions")

test_that("constant kernel matrix decompositions are correct", {
  m <- create_model(y ~ sex + age + age | sex + age | id + group, testdata_002)
  K <- const_kernels(m)
  dec <- const_kernels.decompositions(m)
  K_rec <- const_kernels.reconstruct(dec)
  expect_equal(length(K), 5)
  expect_equal(length(K_rec), 5)

  # Compute maximum absolute errors between recontructions and originals
  mae1 <- max(abs(K[[1]] - K_rec[[1]]))
  mae2 <- max(abs(K[[2]] - matrix(1, 96, 96)))
  mae3 <- max(abs(K[[3]] - K_rec[[3]]))
  mae4 <- max(abs(K[[4]] - K_rec[[4]]))
  mae5 <- max(abs(K[[5]] - K_rec[[5]]))
  expect_lt(mae1, 1e-6)
  expect_lt(mae2, 1e-6)
  expect_lt(mae3, 1e-6)
  expect_lt(mae4, 1e-6)
  expect_lt(mae5, 1e-6)
})

test_that("entire additive gp matrix low-rank decomposition works", {
  m <- create_model(y ~ sex + age + age | sex + age | id + group, testdata_002)
  alpha <- c(1, 1, 1, 1, 1)
  ell <- c(1, 1, 1)
  dec <- kernel_decomposition(m, alpha, ell, num_basisfun = 4)
  K_rec <- dec$V %*% diag(dec$D) %*% t(dec$V)
  expect_equal(dec$ranks, c(1, 1, 1, 11, 1))
  expect_equal(dim(dec$V), c(96, 15 * 4))
  expect_equal(length(dec$D), 15 * 4)

  # Reconstruct
  expect_equal(dim(K_rec), c(96, 96))
})

test_that("additive gp matrix decomposition works when there's one component", {
  m <- create_model(y ~ age, testdata_002)
  alpha <- c(1)
  ell <- c(1)

  # True kernel
  x_age <- m@stan_input$x_cont[1, ]
  K_se <- lgpr:::sim.kernel_se(x_age, x_age, 1, 1)

  # Reconstructions with different numbers of basis functions
  NUM_BF <- seq(2, 20, by = 2)
  L <- length(NUM_BF)
  MAE <- rep(0, L)
  for (j in seq_len(L)) {
    dec <- kernel_decomposition(m, alpha, ell,
      num_basisfun = NUM_BF[j],
      width_basisfun = 5
    )
    K_rec <- dec$V %*% diag(dec$D) %*% t(dec$V)
    MAE[j] <- mean(abs(K_se - K_rec))
  }
  logerr <- log(MAE)

  # Reconstruct
  expect_false(is.unsorted(rev(logerr)))
})

test_that("fitting basis fun approx works correctly (one component)", {
  sd <- simulate_data(
    N = 6, t_data = seq(1, 5, by = 0.6),
    relevances = c(0, 1),
    lengthscales = c(1, 0.5), t_jitter = 0.5
  )
  dat <- sd@data
  # plot(dat$age, dat$y)
  suppressWarnings({
    f1 <- lgp(y ~ age, dat, iter = 1000, chains = 1, refresh = 0)
    f2 <- lgp(y ~ age, dat,
      iter = 1000, chains = 1,
      options = list(num_basisfun = 30, width_basisfun = 2.5),
      refresh = 0
    )
  })
  # test that f2 is close to f1
  pm1 <- rstan::get_posterior_mean(f1@stan_fit)[1:3]
  pm2 <- rstan::get_posterior_mean(f2@stan_fit)[1:3]
  diff <- abs(pm1 - pm2)
  expect_lt(diff[1], 0.2)
  expect_lt(diff[2], 0.1)
  expect_lt(diff[3], 0.1)
})
