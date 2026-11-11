# Summary of the T15-E calibration deposits (read-only). Usage: Rscript --vanilla mh_calib_summary.R <calibdir>
args <- commandArgs(trailingOnly = TRUE)
for (n in c(8, 10, 12, 15, 20, 25, 30, 35, 40, 45, 49)) {
  r <- readRDS(file.path(args[1], sprintf("calib_n%d.rds", n)))
  s <- character(0)
  for (a in c("0.1", "0.05", "0.01")) {
    x <- r[[paste0("rwa_", a)]]; d <- r[[paste0("rwd_", a)]]; w <- r[[paste0("dw_", a)]]
    s <- c(s, sprintf("a%s DW %.3f/%.3f RWa cert=%s step=%s u2=%.4f cp=%.4f RWd %s", a, w$withadd, w$noadd,
                      x$certified, x$ladder_step, x$union_half2, x$cp_upper_half2,
                      if (is.null(d)) "none" else sprintf("cert=%s step=%s", d$certified, d$ladder_step)))
  }
  cat("n", n, ":", paste(s, collapse = " ; "), "\n")
}
