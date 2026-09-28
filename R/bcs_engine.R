# Port of the GAUSS procedures of Banerjee and Carrion-i-Silvestre (2015),
# replication files "factcoint.src" (procedures factcoint, factcoint_iter,
# MQ_test) and "brkcoint.src" (procedure ADFRC), July 2013, available from
# https://sites.google.com/view/carrion-i-silvestre/code-and-data
#
# Models (GAUSS numbering, first differences estimation):
#   1 individual effects
#   2 individual effects and time trend
#   3 individual effects and level shift
#   4 individual effects, time trend and level shift
#   5 individual effects, time trend, level and trend shift (common break)
#
# One deliberate difference: for model 4 the GAUSS code selects the break by
# minimising b'b (the squared coefficient vector) instead of the sum of
# squared residuals used for every other model and described in the paper;
# here the break minimises the sum of squared residuals for all models.


#' Regressors of the first-differenced cointegrating regression
#' @keywords internal
#' @noRd
.bcs_design <- function(model, T, j, Dx) {
  DTb <- numeric(T)
  if (model >= 3) DTb[j + 1] <- 1
  one <- rep(1, T - 1)
  if (model == 1) {
    Dx
  } else if (model == 2) {
    cbind(one, Dx)
  } else if (model == 3) {
    cbind(DTb[2:T], Dx)
  } else if (model == 4) {
    cbind(one, DTb[2:T], Dx)
  } else {
    du <- c(rep(0, j), rep(1, T - j))
    cbind(one, DTb[2:T], du[2:T], Dx)
  }
}


#' OLS residuals (QR based, as GAUSS "/")
#' @keywords internal
#' @noRd
.bcs_resid <- function(y, X) {
  as.numeric(qr.resid(qr(X), y))
}


#' GAUSS procedure factcoint
#'
#' @param y T x N matrix of dependent variables.
#' @param x list of k T x N matrices (one per regressor).
#' @param model model number (1 to 5).
#' @param unknown TRUE if the break dates are estimated.
#' @param m_Tb integer N-vector of known break dates (ignored if unknown).
#' @param kmax maximum number of factors (0 = no factors).
#' @param trim trimming fraction (0.15 in the GAUSS code).
#' @return list(e = cumulated idiosyncratic residuals (T-1 x N),
#'   De = their first differences, Fhat = cumulated factors,
#'   csi = loadings, m_tbe = break dates, r = number of factors)
#' @keywords internal
#' @noRd
.bcs_factcoint <- function(y, x, model, unknown, m_Tb, kmax, trim = 0.15) {
  T <- nrow(y)
  N <- ncol(y)
  D_res <- matrix(0, T - 1, N)
  lo <- floor(trim * T)
  hi <- floor((1 - trim) * T)
  m_tbe <- m_Tb
  SSR <- NULL
  if (unknown && model == 5) SSR <- matrix(0, hi - lo + 1, N)

  Dx_i <- function(i) {
    do.call(cbind, lapply(x, function(xk) diff(xk[, i])))
  }

  for (i in seq_len(N)) {
    Dy <- diff(y[, i])
    Dx <- Dx_i(i)
    if (model <= 2) {
      D_res[, i] <- .bcs_resid(Dy, .bcs_design(model, T, 0, Dx))
      m_tbe <- rep(0L, N)
    } else if (!unknown) {
      D_res[, i] <- .bcs_resid(Dy, .bcs_design(model, T, m_Tb[i], Dx))
    } else if (model %in% c(3, 4)) {
      best <- Inf
      for (j in lo:hi) {
        r <- .bcs_resid(Dy, .bcs_design(model, T, j, Dx))
        s <- sum(r^2)
        if (s < best) {
          best <- s
          D_res[, i] <- r
          m_tbe[i] <- j
        }
      }
    } else {
      for (j in lo:hi) {
        r <- .bcs_resid(Dy, .bcs_design(model, T, j, Dx))
        SSR[j - lo + 1, i] <- sum(r^2)
      }
    }
  }

  if (unknown && model == 5) {
    jb <- which.min(rowSums(SSR)) + lo - 1
    m_tbe <- rep(as.integer(jb), N)
    for (i in seq_len(N)) {
      D_res[, i] <- .bcs_resid(diff(y[, i]), .bcs_design(5, T, jb, Dx_i(i)))
    }
  }

  r <- 0L
  csi <- NULL
  fhat <- NULL
  if (kmax > 0) {
    # Bai and Ng (2002) IC1, as in the GAUSS code
    ev <- eigen(crossprod(D_res), symmetric = TRUE)$vectors
    sig0 <- mean(colSums(D_res^2 / T))
    IC <- log(sig0)
    kk <- min(kmax, N - 1, T - 2)
    for (k in seq_len(kk)) {
      cs <- sqrt(N) * ev[, 1:k, drop = FALSE]
      fh <- D_res %*% cs / N
      De <- D_res - fh %*% t(cs)
      CT <- log(N * T / (N + T)) * k * (N + T) / (N * T)
      IC <- c(IC, log(mean(colSums(De^2 / T))) + CT)
    }
    r <- which.min(IC) - 1L
    if (r > 0) {
      csi <- sqrt(N) * ev[, 1:r, drop = FALSE]
      fhat <- D_res %*% csi / N
    }
  }
  De <- if (r > 0) D_res - fhat %*% t(csi) else D_res
  list(e = apply(De, 2, cumsum), De = De,
       Fhat = if (r > 0) apply(fhat, 2, cumsum) else NULL,
       fhat = fhat, csi = csi, m_tbe = as.integer(m_tbe), r = r)
}


#' GAUSS procedure factcoint_iter (version in factcoint.src)
#' @keywords internal
#' @noRd
.bcs_factcoint_iter <- function(y, x, model, kmax, tolerance, max_iter,
                                trim = 0.15) {
  T <- nrow(y)
  N <- ncol(y)
  zeros <- rep(0L, N)
  f0 <- .bcs_factcoint(y, x, model, TRUE, zeros, 0, trim)
  SSR_opt <- sum(f0$e^2)
  m_opt <- f0$m_tbe
  i <- 0L
  if (model >= 3) {
    while (i <= max_iter) {
      f1 <- .bcs_factcoint(y, x, model, FALSE, m_opt, kmax, trim)
      SSR_2 <- sum(f1$e^2)
      if (abs(SSR_2 - SSR_opt) > tolerance) {
        y_temp <- y[2:T, , drop = FALSE]
        if (f1$r > 0) y_temp <- y_temp - f1$Fhat %*% t(f1$csi)
        x_temp <- lapply(x, function(xk) xk[2:T, , drop = FALSE])
        f2 <- .bcs_factcoint(y_temp, x_temp, model, TRUE, zeros, 0, trim)
        m_opt <- f2$m_tbe + 1L
        SSR_opt <- SSR_2
      } else {
        break
      }
      i <- i + 1L
    }
  }
  fin <- .bcs_factcoint(y, x, model, FALSE, m_opt, kmax, trim)
  fin$m_tbe <- if (model >= 3) m_opt else zeros
  fin$iterations <- i
  fin$ssr <- sum(fin$e^2)
  fin
}


#' GAUSS procedure ADFRC (no deterministic terms)
#'
#' method = 1: general-to-specific selection from p_max, dropping the last
#' lag while |t| < 1.645; method = 0: p_max lags.
#' @keywords internal
#' @noRd
.bcs_adfrc <- function(res, method = 1, p_max = 4) {
  T <- length(res)
  res1 <- c(NA, res[-T])
  d_res <- res - res1
  lagn <- function(v, k) c(rep(NA, k), v[seq_len(length(v) - k)])
  fit <- function(j) {
    M <- cbind(d_res, res1)
    if (j > 0) for (l in seq_len(j)) M <- cbind(M, lagn(d_res, l))
    M <- M[-seq_len(j + 1), , drop = FALSE]
    Xm <- M[, -1, drop = FALSE]
    XtXi <- solve(crossprod(Xm))
    b <- as.numeric(XtXi %*% crossprod(Xm, M[, 1]))
    e <- M[, 1] - as.numeric(Xm %*% b)
    s2 <- sum(e^2) / (T - ncol(Xm))
    list(b = b, se = sqrt(diag(s2 * XtXi)))
  }
  if (method == 0) {
    f <- fit(p_max)
    return(list(t_adf = f$b[1] / f$se[1], p = p_max))
  }
  for (j in p_max:0) {
    f <- fit(j)
    if (j == 0 || abs(f$b[j + 1] / f$se[j + 1]) >= 1.645) {
      return(list(t_adf = f$b[1] / f$se[1], p = j))
    }
  }
}


#' Critical values (5% column used) of the MQ tests, GAUSS MQ_test
#' @keywords internal
#' @noRd
.bcs_mq_cv <- function(model, T) {
  if (model %in% c(1, 3)) {
    c(-13.730, -23.535, -32.296, -40.442, -48.617, -57.040, -67.465203,
      -76.042352, -83.823797, -92.623120, -100.26543, -108.17400)
  } else if (model %in% c(2, 4)) {
    c(-21.313, -31.356, -40.180, -48.421, -55.818, -64.393, -74.068345,
      -82.332640, -90.895894, -98.474676, -106.64284, -114.87589)
  } else if (T <= 75) {
    c(-24.828, -32.792, -39.703, -44.865, -47.472, -48.444)
  } else if (T <= 200) {
    c(-26.833, -36.464, -45.879, -53.251, -61.099, -67.183)
  } else {
    c(-25.697, -38.103, -45.066, -53.392, -62.404, -68.748)
  }
}


#' GAUSS procedure MQ_test (Bai and Ng, 2004)
#'
#' @param F T x r matrix of detrended common factors (levels).
#' @param parametric FALSE for the non-parametric (Bartlett) version, TRUE
#'   for the parametric (VAR filtered, BIC) version.
#' @return list(MQ, n_trends)
#' @keywords internal
#' @noRd
.bcs_mq <- function(F, model, N, parametric) {
  F <- as.matrix(F)
  r <- ncol(F)
  T <- nrow(F)
  cv <- .bcs_mq_cv(model, T)
  u_F <- svd(crossprod(F) / T^2)$u
  lagm <- function(M, k) rbind(matrix(NA, k, ncol(M)), M[seq_len(nrow(M) - k), , drop = FALSE])
  r_star <- r
  MQ <- NA_real_
  repeat {
    Yc <- F %*% u_F[, 1:r_star, drop = FALSE]
    m <- ncol(Yc)
    if (!parametric) {
      tmp <- cbind(Yc, lagm(Yc, 1))[-1, , drop = FALSE]
      Z <- tmp[, (m + 1):(2 * m), drop = FALSE]
      res <- tmp[, 1:m, drop = FALSE] - Z %*% qr.solve(Z, tmp[, 1:m, drop = FALSE])
      bigJ <- 4 * ceiling((min(T, N) / 100)^(1 / 4))
      sigma <- matrix(0, m, m)
      nr <- nrow(res)
      for (j in seq_len(bigJ)) {
        if (j >= nr) break
        sigma <- sigma + (1 - j / (bigJ + 1)) *
          crossprod(res[1:(nr - j), , drop = FALSE], res[(j + 1):nr, , drop = FALSE]) / T
      }
      A <- Yc[2:T, , drop = FALSE]
      B <- Yc[1:(T - 1), , drop = FALSE]
      PHI <- 0.5 * (crossprod(A, B) + crossprod(B, A) - T * (sigma + t(sigma))) %*%
        solve(crossprod(B))
    } else {
      pmax <- 4L * as.integer(floor((T / 100)^(1 / 4)))
      DY <- Yc - lagm(Yc, 1)
      p <- 0L
      if (pmax > 0) {
        bic <- numeric(pmax)
        for (i in seq_len(pmax)) {
          Dl <- do.call(cbind, lapply(seq_len(i), function(l) lagm(DY, l)))
          tmp <- cbind(DY, Dl)[-seq_len(pmax + 1), , drop = FALSE]
          Z <- tmp[, -(1:m), drop = FALSE]
          res <- tmp[, 1:m, drop = FALSE] - Z %*% qr.solve(Z, tmp[, 1:m, drop = FALSE])
          detO <- det(crossprod(res)) / (T - m)
          l <- (-T / 2) * (m * (1 + log(2 * pi)) + log(detO))
          bic[i] <- -2 * l / T + (m^2 * ncol(Z)) * log(T) / T
        }
        D0 <- DY[-1, , drop = FALSE]
        detO <- det(crossprod(D0)) / (T - 1)
        l <- (-T / 2) * (m * (1 + log(2 * pi)) + log(detO))
        p <- which.min(c(-2 * l / T, bic)) - 1L
      }
      if (p > 0) {
        Dl <- do.call(cbind, lapply(seq_len(p), function(l) lagm(DY, l)))
        tmp <- cbind(DY, Dl)[-seq_len(p + 1), , drop = FALSE]
        Z <- tmp[, -(1:m), drop = FALSE]
        beta <- qr.solve(Z, tmp[, 1:m, drop = FALSE])
        Yl <- do.call(cbind, lapply(seq_len(p), function(l) lagm(Yc, l)))
        tmp2 <- cbind(Yc, Yl)[-seq_len(p), , drop = FALSE]
        Yf <- tmp2[, 1:m, drop = FALSE] - tmp2[, -(1:m), drop = FALSE] %*% beta
      } else {
        Yf <- Yc
      }
      n <- nrow(Yf)
      A <- Yf[2:n, , drop = FALSE]
      B <- Yf[1:(n - 1), , drop = FALSE]
      PHI <- 0.5 * (crossprod(A, B) + crossprod(B, A)) %*% solve(crossprod(B))
    }
    v_c <- min(svd(PHI)$d)
    MQ <- T * (v_c - 1)
    if (r_star > 0 && r_star <= length(cv) && MQ < cv[r_star]) {
      r_star <- r_star - 1L
      if (r_star > 0) next
    }
    break
  }
  list(MQ = MQ, n_trends = r_star)
}
