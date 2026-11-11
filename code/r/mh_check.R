# DEV check (root 2026110826): (1) our local_z/dw_reject reproduce modehunt::modeHunting decisions on random samples;
# (2) timing of the max statistic at n = 49. Not a result.
source(file.path(Sys.getenv("DIPTEST_R", "code/r"), "mh_global.R"))
source(file.path(Sys.getenv("DIPTEST_R", "code/r"), "mh_stats.R"))
set.seed(2026110826)
mism <- 0
for (t in 1:300) {
  n <- sample(c(8, 12, 20, 35, 49), 1)
  x <- if (t %% 2 == 0) c(rnorm(n %/% 2, 0, 1), rnorm(n - n %/% 2, runif(1, 2, 6), 1)) else runif(n)
  cv <- c(runif(1, 0.3, 2), runif(1, 1.5, 3.5))
  ref <- modeHunting(x, crit.vals = cv, min.int = TRUE)
  a <- c(antimode_found(ref$Dp, ref$Dm), antimode_found(ref$Dp.noadd, ref$Dm.noadd))
  b <- dw_reject(x, cv)
  if (any(a != b)) mism <- mism + 1
}
cat("decision mismatches vs modeHunting:", mism, "/ 300\n")
t0 <- Sys.time(); for (i in 1:500) dw_max(runif(49)); cat("dw_max n=49: ", format((Sys.time() - t0) / 500), "per sample\n")
