# Read-only digest 3: cells where a comparator beats R1 (upper99 < 0). Usage: Rscript --vanilla mh_power_digest3.R <summary.rds>
args <- commandArgs(trailingOnly = TRUE)
s <- readRDS(args[1])$table
z <- s[s$test != "R1" & !is.na(s$d_hi) & s$d_hi < 0 & s$test != "TAB", ]
print(z[, c("n", "law", "alpha", "test", "k", "rate", "d_R1", "d_lo", "d_hi")])
r1 <- s[s$test == "R1" & s$law %in% z$law & s$n %in% z$n, c("n", "law", "alpha", "k", "rate")]
print(r1[r1$alpha %in% z$alpha, ])
