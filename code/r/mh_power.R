# Step 15 / DT14 WP-E (T15-E): rejection rates of dip-based and multiscale tests on COMMON samples.
# Usage: Rscript --vanilla mh_power.R <root> <n> <N> <calibdir> <cutoffs.csv> <outdir>
# Tests per sample x (size n), alpha in {.10, .05, .01}:
#   LF   : dip(x) >= C_HI[alpha] / sqrt(n)            (the "fix"; certified asymptotic constant, upper enclosure end)
#   R1   : dip(x) >= d_cert[n, alpha]                  (T15-D certified finite-n cutoff; absent cells = no test)
#   TAB  : diptest::dip.test(x)$p.value <= alpha       (conventional uniform-table calibration; NOT valid over U_c)
#   DWa  : DW08 all intervals with Gamma, matched MC calibration (mh_calib.R)
#   DWn  : DW08 all intervals without Gamma, matched MC calibration
#   RWd  : RW10 block procedure, modehunt default interval set (structurally unable to reject at small n)
#   RWa  : RW10 block criterion on all intervals (prespecified small-sample adaptation)
# Laws: nulls (diagnostic) and alternatives with documented modality (see LAWS). Seeds: set.seed(root + 1000 * n + law).
source(file.path(Sys.getenv("DIPTEST_R", "code/r"), "mh_stats.R"))
suppressMessages(library(diptest, lib.loc = MH_LIB))
args <- commandArgs(trailingOnly = TRUE)
root <- as.integer(args[1]); n <- as.integer(args[2]); N <- as.integer(args[3])
calibdir <- args[4]; cutfile <- args[5]; outdir <- args[6]
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
ALPHAS <- c(0.10, 0.05, 0.01)
C_HI <- c("0.1" = 0.6685, "0.05" = 0.7315, "0.01" = 0.861)

witness <- function(m, w, eps) { L <- 1 - eps; top <- runif(m) < w; ifelse(top, L + eps * runif(m), L * runif(m)) }
mixn <- function(m, p, d, s2 = 1) { a <- runif(m) < p; ifelse(a, rnorm(m, 0, 1), rnorm(m, d, s2)) }
mixu <- function(m, g) { a <- runif(m) < 0.5; ifelse(a, runif(m), 1 + g + runif(m)) }
# Local valley in a broad background: density on [0, 1] proportional to 1 - h * max(0, 1 - |x - 1/2| / w)
# (a piecewise-linear notch of half-width w and relative depth h at 1/2); NOT unimodal for h > 0. Exact rejection sampler.
valley <- function(m, w, h) {
  out <- numeric(0)
  while (length(out) < m) {
    x <- runif(2 * m)
    keep <- runif(2 * m) <= 1 - h * pmax(0, 1 - abs(x - 0.5) / w)
    out <- c(out, x[keep])
  }
  out[1:m]
}
LAWS <- list(
  N_U01 = function(m) runif(m),
  N_N01 = function(m) rnorm(m),
  N_EXP = function(m) rexp(m),
  N_T3 = function(m) rt(m, 3),
  N_BETA25 = function(m) rbeta(m, 2, 5),
  N_TRAP = function(m) runif(m, -1, 1) + runif(m, -0.3, 0.3),
  N_W01 = function(m) witness(m, 0.1, 0.005),
  A_N50_D22 = function(m) mixn(m, 0.5, 2.2),
  A_N50_D3 = function(m) mixn(m, 0.5, 3),
  A_N50_D4 = function(m) mixn(m, 0.5, 4),
  A_N50_D5 = function(m) mixn(m, 0.5, 5),
  A_N30_D35 = function(m) mixn(m, 0.3, 3.5),
  A_N30_D5 = function(m) mixn(m, 0.3, 5),
  A_NSC_D3 = function(m) mixn(m, 0.5, 3, 0.5),
  A_U_G05 = function(m) mixu(m, 0.5),
  A_U_G1 = function(m) mixu(m, 1.0),
  A_U_G3 = function(m) mixu(m, 3.0),
  A_N50_D8 = function(m) mixn(m, 0.5, 8),
  A_TIGHT = function(m) ifelse(runif(m) < 0.5, rnorm(m, 0, 0.05), rnorm(m, 1, 0.05)),
  A_VAL_W10_H6 = function(m) valley(m, 0.10, 0.6),
  A_VAL_W10_H10 = function(m) valley(m, 0.10, 1.0),
  A_VAL_W25_H6 = function(m) valley(m, 0.25, 0.6),
  A_VAL_W25_H10 = function(m) valley(m, 0.25, 1.0)
)

cal <- readRDS(file.path(calibdir, sprintf("calib_n%d.rds", n)))
cut <- read.csv(cutfile)
r1 <- setNames(rep(NA_real_, 3), as.character(ALPHAS))
for (a in ALPHAS) { v <- cut$d_cert[cut$n == n & abs(cut$alpha - a) < 1e-12]; if (length(v) == 1) r1[as.character(a)] <- v }

TESTS <- c("LF", "R1", "TAB", "DWa", "DWn", "RWd", "RWa")
decide <- function(x) {
  D <- dip(x)
  p <- dip.test(x)$p.value
  out <- matrix(FALSE, length(TESTS), 3, dimnames = list(TESTS, as.character(ALPHAS)))
  for (a in ALPHAS) {
    k <- as.character(a)
    out["LF", k] <- D >= C_HI[[k]] / sqrt(n)
    out["R1", k] <- if (is.na(r1[k])) NA else D >= r1[[k]]
    out["TAB", k] <- p <= a
    dwk <- cal[[paste0("dw_", k)]]
    dd <- dw_reject(x, c(dwk$withadd, dwk$noadd))
    out["DWa", k] <- dd[["withadd"]]
    out["DWn", k] <- dd[["noadd"]]
    rwd <- cal[[paste0("rwd_", k)]]
    out["RWd", k] <- if (is.null(rwd)) FALSE else rw_reject(x, "rw_default", rwd$crit)
    out["RWa", k] <- rw_reject(x, "rw_all", cal[[paste0("rwa_", k)]]$crit)
  }
  out
}

res <- list(root = root, n = n, N = N, rng = RNGkind(), R = R.version.string, r1 = r1,
            diptest = as.character(packageVersion("diptest", lib.loc = MH_LIB)), laws = list())
t0 <- Sys.time()
for (li in seq_along(LAWS)) {
  set.seed(root + 1000 * n + li)
  cnt <- array(0L, c(length(TESTS), 3), dimnames = list(TESTS, as.character(ALPHAS)))
  pair <- array(0L, c(length(TESTS), 3, 2, 2), dimnames = list(TESTS, as.character(ALPHAS), c("R1no", "R1yes"), c("Tno", "Tyes")))
  for (s in 1:N) {
    d <- decide(LAWS[[li]](n))
    cnt <- cnt + (d %in% TRUE)
    for (a in as.character(ALPHAS)) if (!is.na(d["R1", a])) for (t in TESTS) {
      i <- 1 + d["R1", a]; j <- 1 + isTRUE(d[t, a])
      pair[t, a, i, j] <- pair[t, a, i, j] + 1L
    }
  }
  res$laws[[names(LAWS)[li]]] <- list(counts = cnt, paired_vs_R1 = pair)
  cat(sprintf("n=%d %s %s %.0fs\n", n, names(LAWS)[li],
              paste(sprintf("%s=%.3f", TESTS, cnt[, "0.05"] / N), collapse = " "),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
res$seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
saveRDS(res, file.path(outdir, sprintf("power_n%d.rds", n)))
