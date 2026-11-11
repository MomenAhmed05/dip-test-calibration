# Step 15 / DT14 WP-E (T15-E): Monte Carlo calibration of the multiscale comparators under U(0,1) (matched to the
# statistics in mh_stats.R). Usage: Rscript --vanilla mh_calib.R <root> <n> <M> <outdir>
# DW08 (with Gamma / without): TOLERANCE cutoff: c = T_(k), k = smallest index with
#   BetaCDF(1 - alpha; k, M - k + 1) <= DELTA_C, so that P_calib(size of the frozen cutoff > alpha) <= DELTA_C.
# RW10 block criterion (A = 10; RW10 §4; blocks indexed LARGEST first as in modehunt::blocks): on half 1 the largest a
#   with empirical P(exists l: max_l |Z| > q_l(a / (A + l)^2)) <= alpha; thresholds q_l(a) from half 1 (type-1 empirical
#   quantiles); certification gate on half 2: frozen ladder a_j = a * 1.1^(-j), j = 0..80, first a_j whose
#   one-sided CP upper limit (level DELTA_C) of the half-2 union exceedance is <= alpha (fixed sequence: the union
#   probability is nonincreasing in a). If none passes: the test never rejects (recorded).
# DELTA_C = 0.001 / (11 n * 3 alpha * 4 calibrated methods) per cutoff (Bonferroni over the frozen panel).
# Seeds: set.seed(root + n) with RNGkind("Mersenne-Twister", "Inversion", "Rejection") recorded in the output.
source(file.path(Sys.getenv("DIPTEST_R", "code/r"), "mh_stats.R"))
args <- commandArgs(trailingOnly = TRUE)
root <- as.integer(args[1]); n <- as.integer(args[2]); M <- as.integer(args[3]); outdir <- args[4]
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
set.seed(root + n)
ALPHAS <- c(0.10, 0.05, 0.01)
A <- 10
t0 <- Sys.time()
dw <- matrix(NA_real_, M, 2)
rwd <- NULL; rwa <- NULL
Pd <- rw_pairs(n, "rw_default"); Pa <- rw_pairs(n, "rw_all")
if (!is.null(Pd) && length(Pd$j) > 0) rwd <- matrix(NA_real_, M, Pd$nb)
rwa <- matrix(NA_real_, M, Pa$nb)
for (s in 1:M) {
  U <- runif(n)
  dw[s, ] <- dw_max(U)
  if (!is.null(rwd)) rwd[s, ] <- rw_blockmax(U, "rw_default")
  rwa[s, ] <- rw_blockmax(U, "rw_all")
}
DELTA_C <- 0.001 / (11 * 3 * 4)
tol_cut <- function(v, a) {
  v <- sort(v); M <- length(v)
  ks <- seq(ceiling(M * (1 - a)), M)
  ok <- ks[pbeta(1 - a, ks, M - ks + 1) <= DELTA_C]
  if (length(ok) == 0) Inf else v[min(ok)]
}
q_emp <- function(v, p) { v <- sort(v); k <- ceiling(length(v) * p); if (k < 1) -Inf else if (k > length(v)) Inf else v[k] }
block_crit <- function(Mx, alpha) {
  half <- floor(nrow(Mx) / 2)
  M1 <- Mx[1:half, , drop = FALSE]; M2 <- Mx[(half + 1):nrow(Mx), , drop = FALSE]
  L <- ncol(Mx)
  thr <- function(a) vapply(1:L, function(l) q_emp(M1[, l], 1 - a / (A + l)^2), numeric(1))
  p_union <- function(Mz, q) mean(apply(sweep(Mz, 2, q, ">"), 1, any))
  lo <- 0; hi <- alpha * (A + 1)^2      # v1.1 fix: a may exceed 1; only a/(A + l)^2 <= alpha <= 1 is needed
  for (it in 1:60) { mid <- (lo + hi) / 2; if (p_union(M1, thr(mid)) <= alpha) lo <- mid else hi <- mid }
  n2 <- nrow(M2)
  for (jj in 0:80) {
    aj <- lo * 1.1^(-jj)
    q <- thr(aj)
    k2 <- sum(apply(sweep(M2, 2, q, ">"), 1, any))
    cpu <- if (k2 >= n2) 1 else qbeta(1 - DELTA_C, k2 + 1, n2 - k2)
    if (cpu <= alpha) return(list(a_tilde = lo, a_cert = aj, ladder_step = jj, crit = q, union_half1 = p_union(M1, q),
                                  union_half2 = k2 / n2, cp_upper_half2 = cpu, certified = TRUE))
  }
  list(a_tilde = lo, crit = rep(Inf, L), certified = FALSE)
}
res <- list(root = root, n = n, M = M, rng = RNGkind(), A = A, R = R.version.string, DELTA_C = DELTA_C,
            modehunt = as.character(packageVersion("modehunt", lib.loc = MH_LIB)))
for (a in ALPHAS) {
  key <- as.character(a)
  res[[paste0("dw_", key)]] <- list(withadd = tol_cut(dw[, 1], a), noadd = tol_cut(dw[, 2], a))
  res[[paste0("rwa_", key)]] <- block_crit(rwa, a)
  res[[paste0("rwd_", key)]] <- if (is.null(rwd)) NULL else block_crit(rwd, a)
}
res$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
writeLines(deparse(res), file.path(outdir, sprintf("calib_n%d.Rdput", n)))
saveRDS(res, file.path(outdir, sprintf("calib_n%d.rds", n)))
cat(sprintf("n=%d M=%d %.1fs dw05=%.4f/%.4f rwa05=%s\n", n, M, res$seconds, res$dw_0.05$withadd, res$dw_0.05$noadd,
            paste(round(res$rwa_0.05$crit, 3), collapse = ",")))
