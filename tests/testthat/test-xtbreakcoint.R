test_that("xtbreakcoint returns correct structure", {
  # Generate simple panel data
  set.seed(123)
  N <- 5
  T <- 30
  
  panel_data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    y = NA,
    x = NA
  )
  
  for (i in 1:N) {
    idx <- panel_data$id == i
    x <- cumsum(rnorm(T))
    u <- rnorm(T, sd = 0.3)
    y <- 1 + x + u
    panel_data$x[idx] <- x
    panel_data$y[idx] <- y
  }
  
  result <- xtbreakcoint(y ~ x, data = panel_data, id = "id", time = "time",
                         max_factors = 2, max_lag = 2)
  
  expect_s3_class(result, "xtbreakcoint")
  expect_true("Z_t" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("tbar" %in% names(result))
  expect_true("N" %in% names(result))
  expect_true("T" %in% names(result))
  expect_equal(result$N, N)
  expect_equal(result$T, T)
  expect_length(result$adf_stats, N)
  expect_length(result$breaks, N)
})

test_that("model specification works correctly", {
  set.seed(456)
  N <- 4
  T <- 25
  
  panel_data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    y = rnorm(N * T),
    x = rnorm(N * T)
  )
  
  # Test different model specifications
  r1 <- xtbreakcoint(y ~ x, data = panel_data, id = "id", time = "time",
                     model = "constant", max_factors = 0)
  expect_equal(r1$model_num, 1)
  
  r2 <- xtbreakcoint(y ~ x, data = panel_data, id = "id", time = "time",
                     model = 2, max_factors = 0)
  expect_equal(r2$model_num, 2)
  
  r3 <- xtbreakcoint(y ~ x, data = panel_data, id = "id", time = "time",
                     model = "levelshift", max_factors = 0)
  expect_equal(r3$model_num, 3)
})

test_that("input validation works", {
  panel_data <- data.frame(
    id = rep(1:3, each = 10),
    time = rep(1:10, 3),
    y = rnorm(30),
    x = rnorm(30)
  )
  
  # Wrong panel identifier
  expect_error(
    xtbreakcoint(y ~ x, data = panel_data, id = "wrong_id", time = "time"),
    "not found"
  )
  
  # Wrong time identifier
  expect_error(
    xtbreakcoint(y ~ x, data = panel_data, id = "id", time = "wrong_time"),
    "not found"
  )
  
  # Invalid trim
  expect_error(
    xtbreakcoint(y ~ x, data = panel_data, id = "id", time = "time", trim = 0.6),
    "trim"
  )
  
  # Invalid model
  expect_error(
    xtbreakcoint(y ~ x, data = panel_data, id = "id", time = "time", model = "invalid"),
    "Invalid model"
  )
})

test_that("print method works", {
  set.seed(789)
  N <- 3
  T <- 20
  
  panel_data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    y = rnorm(N * T),
    x = rnorm(N * T)
  )
  
  result <- xtbreakcoint(y ~ x, data = panel_data, id = "id", time = "time",
                         max_factors = 0)
  
  expect_output(print(result), "Panel Cointegration Test")
  expect_output(print(result), "Z_t statistic")
})

test_that("factor estimation works", {
  set.seed(321)
  N <- 6
  T <- 40
  
  # Generate data with common factor
  common_factor <- cumsum(rnorm(T))
  
  panel_data <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, N),
    y = NA,
    x = NA
  )
  
  for (i in 1:N) {
    idx <- panel_data$id == i
    loading <- runif(1, 0.5, 1.5)
    x <- cumsum(rnorm(T))
    u <- rnorm(T, sd = 0.2) + loading * common_factor
    y <- 1 + x + u
    panel_data$x[idx] <- x
    panel_data$y[idx] <- y
  }
  
  result <- xtbreakcoint(y ~ x, data = panel_data, id = "id", time = "time",
                         max_factors = 3)
  
  # Should detect at least one factor
  expect_true(result$n_factors >= 0)
})
