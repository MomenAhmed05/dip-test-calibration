# Read-only digest 2 (paired-interval counts). Usage: Rscript --vanilla mh_power_digest2.R <summary.rds>
args <- commandArgs(trailingOnly = TRUE)
s <- readRDS(args[1])$table
z <- s[grepl("^A_", s$law) & s$test != "R1", ]
cat("NA paired intervals:", sum(is.na(z$d_lo)), "of", nrow(z), "\n")
print(head(z[is.na(z$d_lo), c("n", "law", "alpha", "test", "k", "d_R1", "d_lo", "d_hi")], 5))
for (t in c("LF", "DWa", "DWn", "RWa", "RWd")) for (a in c(0.1, 0.05, 0.01)) {
  w <- z[z$test == t & z$alpha == a & !is.na(z$d_lo), ]
  cat(sprintf("%-4s alpha %.2f: R1 better (lower99 > 0) in %3d / %3d cells; X better (upper99 < 0) in %d; max diff %.3f\n",
              t, a, sum(w$d_lo > 0), nrow(w), sum(w$d_hi < 0), max(w$d_R1)))
}
