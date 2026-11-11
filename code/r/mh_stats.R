# Step 15 / DT14 WP-E (T15-E): exact re-statement of the modehunt statistics as functions of a sorted sample, so that
# critical values are calibrated on THE SAME statistic that the test uses ("matched recalibration").
# Source of truth: modehunt 1.0.8 modeHunting (all intervals) and modeHuntingBlock (blocks); unknown support
# (lower = -Inf, upper = Inf), so the interval endpoints are sample points X[j] < X[k], k >= j + 2, and
# T_jk = sum_{i=j+1}^{k-1} (2 (X_i - X_j)/(X_k - X_j) - 1), standardised by sqrt((k-j-1)/3).
# DW08 additive correction as coded in modeHunting: Gamma_{jk} = sqrt(2 (1 + log((n - 1)/(k - j)))), n = sample size.

MH_LIB <- if (nzchar(Sys.getenv("DIPTEST_RLIB"))) Sys.getenv("DIPTEST_RLIB") else NULL
suppressMessages(library(modehunt, lib.loc = MH_LIB))

# All standardised local statistics Z_jk = T_jk / sigma_jk for 1 <= j, k <= n, k >= j + 2 (fully vectorised;
# sum_{i=j+1}^{k-1} X_i = cs[k-1] - cs[j]). Pair index tables are cached per n.
.PAIRS <- new.env()
pairs_of <- function(n) {
  key <- as.character(n)
  if (!exists(key, envir = .PAIRS)) {
    jk <- which(upper.tri(matrix(0, n, n)), arr.ind = TRUE)    # row < col
    jk <- jk[jk[, 2] - jk[, 1] >= 2, , drop = FALSE]
    assign(key, list(j = jk[, 1], k = jk[, 2], i = jk[, 2] - jk[, 1] - 1), envir = .PAIRS)
  }
  get(key, envir = .PAIRS)
}
local_z <- function(X) {
  n <- length(X)
  P <- pairs_of(n)
  cs <- cumsum(X)
  Tjk <- 2 * (cs[P$k - 1] - cs[P$j] - P$i * X[P$j]) / (X[P$k] - X[P$j]) - P$i
  list(j = P$j, k = P$k, Z = Tjk * sqrt(3 / P$i))
}

# Max statistics used by the global calibration of DW08: with Gamma and without.
dw_max <- function(X) {
  n <- length(X)
  L <- local_z(sort(X))
  g <- sqrt(2 * (1 + log((n - 1) / (L$k - L$j))))
  c(withadd = max(abs(L$Z) - g), noadd = max(abs(L$Z)))
}

# Global unimodality rejection from directional interval sets: some interval with a significant
# DECREASE (Z < -c) lies entirely to the left of some interval with a significant INCREASE (Z > c).
# Minimal intervals lose nothing, so we test all significant intervals directly: reject iff
#   min over increase intervals of X[j]  >=  min over decrease intervals of X[k]   is NOT the right test; we need
#   exists (dec, inc) with X[k_dec] <= X[j_inc]  <=>  min_dec X[k] <= max_inc X[j].
antimode_reject <- function(X, j, k, zpos, zneg) {
  if (!any(zpos) || !any(zneg)) return(FALSE)
  min(X[k[zneg]]) <= max(X[j[zpos]])
}

# ---------------- RW10 block procedure ----------------
# (a) "rw_default": the exact modehunt::modeHuntingBlock interval set I_app and block assignment (d0 = 2, m0 = 10,
#     fm = 2); structurally empty for n < 20 (n.blocks = 0) and unable to give two ordered intervals for n < 23.
# (b) "rw_all": prespecified small-sample adaptation: ALL intervals, blocks by point count k - j + 1 in
#     dyadic classes [3,3], [4,7], [8,15], [16,31], [32,63], [64,127], ordered LARGEST first (block 1 = largest), as in
#     modehunt::blocks. Both are calibrated by the RW10 block criterion (A = 10) with our own Monte Carlo.
rw_default_pairs <- function(n, d0 = 2, m0 = 10, fm = 2) {
  n.blocks <- floor(log(n / m0) / log(fm))
  if (n.blocks < 1) return(NULL)
  fd <- sqrt(fm)
  J <- c(); K <- c()
  for (block in 1:n.blocks) {
    d <- myRound(d0 * fd^(n.blocks - block))
    m <- myRound(m0 * fm^(n.blocks - block))
    minspan <- d * ceiling((m + 1) / d)
    maxspan <- d * floor(fm * m / d)
    for (j in seq(1, n - m + 1, by = d)) {
      if (min(n, j + minspan) > min(n, j + maxspan)) next
      k <- seq(min(n, j + minspan), min(n, (j + maxspan)), by = d)
      J <- c(J, rep(j, length(k))); K <- c(K, k)
    }
  }
  keep <- K - J >= 2
  J <- J[keep]; K <- K[keep]
  u <- !duplicated(cbind(J, K))
  J <- J[u]; K <- K[u]
  bs <- blocks(n, m0, fm)                      # rows: block 1 = largest sizes; columns l1, l2 on k - j + 1
  sz <- K - J + 1
  B <- rep(NA_integer_, length(J))
  for (b in 1:n.blocks) B[sz >= bs[b, 1] & sz <= bs[b, 2]] <- b   # modeHuntingBlock tests rows 1..n.blocks only
  ok <- !is.na(B)
  list(j = J[ok], k = K[ok], i = K[ok] - J[ok] - 1, block = B[ok], nb = n.blocks)
}
rw_all_pairs <- function(n) {
  P <- pairs_of(n)
  sz <- P$k - P$j + 1
  cls <- floor(log2(sz))          # 3 -> 1, 4..7 -> 2, 8..15 -> 3, ...
  nb <- max(cls)
  list(j = P$j, k = P$k, i = P$i, block = nb - cls + 1, nb = nb)   # block 1 = largest class
}
.RWP <- new.env()
rw_pairs <- function(n, variant) {
  key <- paste(variant, n)
  if (!exists(key, envir = .RWP)) assign(key, if (variant == "rw_default") rw_default_pairs(n) else rw_all_pairs(n), envir = .RWP)
  get(key, envir = .RWP)
}
rw_z <- function(X, P) {
  cs <- cumsum(X)
  Tjk <- 2 * (cs[P$k - 1] - cs[P$j] - P$i * X[P$j]) / (X[P$k] - X[P$j]) - P$i
  Tjk * sqrt(3 / P$i)
}
rw_blockmax <- function(Xraw, variant) {        # per-block max |Z| (length nb); NULL if no intervals
  X <- sort(Xraw)
  P <- rw_pairs(length(X), variant)
  if (is.null(P)) return(NULL)
  Z <- abs(rw_z(X, P))
  vapply(seq_len(P$nb), function(b) { v <- Z[P$block == b]; if (length(v)) max(v) else -Inf }, numeric(1))
}
rw_reject <- function(Xraw, variant, cvec) {     # cvec: per-block critical values on |Z|
  X <- sort(Xraw)
  P <- rw_pairs(length(X), variant)
  if (is.null(P)) return(FALSE)
  Z <- rw_z(X, P)
  cb <- cvec[P$block]
  antimode_reject(X, P$j, P$k, Z >= cb, -Z >= cb)        # modeHuntingBlock uses >=
}

dw_reject <- function(Xraw, cv) {      # cv = c(withadd, noadd)
  X <- sort(Xraw)
  n <- length(X)
  L <- local_z(X)
  g <- sqrt(2 * (1 + log((n - 1) / (L$k - L$j))))
  c(withadd = antimode_reject(X, L$j, L$k, L$Z - g > cv[1], -L$Z - g > cv[1]),
    noadd = antimode_reject(X, L$j, L$k, L$Z > cv[2], -L$Z > cv[2]))
}
