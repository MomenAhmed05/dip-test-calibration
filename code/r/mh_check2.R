# DEV check (root 2026110827): rw_reject("rw_default") reproduces modehunt::modeHuntingBlock decisions. Not a result.
source(file.path(Sys.getenv("DIPTEST_R", "code/r"), "mh_global.R"))
source(file.path(Sys.getenv("DIPTEST_R", "code/r"), "mh_stats.R"))
set.seed(2026110827)
mism <- 0; tot <- 0; rej <- 0; errs <- 0
for (t in 1:300) {
  n <- sample(c(25, 33, 41, 49), 1)
  x <- if (t %% 2 == 0) c(rnorm(n %/% 2, 0, 1), rnorm(n - n %/% 2, runif(1, 2, 8), 1)) else runif(n)
  P <- rw_default_pairs(n)
  cv <- runif(P$nb, 0.5, 2.5)
  ref <- tryCatch(modeHuntingBlock(x, crit.vals = cv, min.int = TRUE), error = function(e) NULL)
  if (is.null(ref)) { errs <- errs + 1; next }
  a <- antimode_found(ref$Dp, ref$Dm)
  b <- rw_reject(x, "rw_default", cv)
  tot <- tot + 1; rej <- rej + a
  if (a != b) mism <- mism + 1
}
cat("rw_default decision mismatches vs modeHuntingBlock:", mism, "/", tot, " (ref rejections", rej, "; package errors", errs, ")\n")
for (n in c(8, 15, 20, 23, 30, 49)) { P <- rw_pairs(n, "rw_all"); Q <- rw_default_pairs(n)
  cat("n", n, " rw_all blocks", P$nb, " intervals", length(P$j), " | rw_default blocks", if (is.null(Q)) 0 else Q$nb,
      " intervals", if (is.null(Q)) 0 else length(Q$j), "\n") }
