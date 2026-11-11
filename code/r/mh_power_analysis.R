# Step 15 / DT14 WP-E (T15-E): frozen analysis of the power deposits (DT14-E-PREREG v1.1 §5).
# Usage: Rscript --vanilla mh_power_analysis.R <powerdir> <out_prefix>
# Writes <out_prefix>.rds (all numbers) and <out_prefix>.md (tables). Rates: two-sided 99% Clopper-Pearson.
# Paired R1 - X: Newcombe (1998) method 10 (Wilson-score based, no continuity correction), 99%, approximate.
# Modality documentation of every law: numerical density on a grid (mode count = strict local maxima above 1e-9 of the max,
# valley depth = 1 - f(lowest antimode)/min(f(modes))), population dip = dip of the quantile sample F^-1((i - 1/2)/1e5)
# with F from a 2e5-point grid integral (numerical, labelled).
suppressMessages(library(diptest))
args <- commandArgs(trailingOnly = TRUE)
pdir <- args[1]; outp <- args[2]
NS <- c(8, 10, 12, 15, 20, 25, 30, 35, 40, 45, 49)
TESTS <- c("LF", "R1", "TAB", "DWa", "DWn", "RWd", "RWa")
ALPHAS <- c("0.1", "0.05", "0.01")
Z <- qnorm(1 - 0.01 / 2)
cp <- function(k, n) c(if (k == 0) 0 else qbeta(0.005, k, n - k + 1), if (k == n) 1 else qbeta(0.995, k + 1, n - k))
wilson <- function(k, n) { p <- k / n; den <- 1 + Z^2 / n; c((p + Z^2 / (2 * n) - Z * sqrt(p * (1 - p) / n + Z^2 / (4 * n^2))) / den,
                                                            (p + Z^2 / (2 * n) + Z * sqrt(p * (1 - p) / n + Z^2 / (4 * n^2))) / den) }
newcombe10 <- function(a, b, c, d) {          # a: both reject, b: R1 only, c: X only, d: neither
  a <- as.numeric(a); b <- as.numeric(b); c <- as.numeric(c); d <- as.numeric(d)   # fix: 32-bit integer overflow in e*f*g*h
  n <- a + b + c + d; p1 <- (a + b) / n; p2 <- (a + c) / n
  w1 <- wilson(a + b, n); w2 <- wilson(a + c, n)
  e <- a + b; f <- c + d; g <- a + c; h <- b + d
  phi <- if (min(e, f, g, h) == 0) 0 else max(-1, min(1, (a * d - b * c) / sqrt(e * f * g * h)))
  dl <- sqrt(max(0, (p1 - w1[1])^2 - 2 * phi * (p1 - w1[1]) * (w2[2] - p2) + (w2[2] - p2)^2))
  du <- sqrt(max(0, (w1[2] - p1)^2 - 2 * phi * (w1[2] - p1) * (p2 - w2[1]) + (p2 - w2[1])^2))
  c(diff = p1 - p2, lo = p1 - p2 - dl, hi = p1 - p2 + du)
}

# ---- modality documentation (densities match mh_power.R LAWS) ----
dmixn <- function(x, p, d, s) p * dnorm(x) + (1 - p) * dnorm(x, d, s)
dvalley <- function(x, w, h) ifelse(x >= 0 & x <= 1, (1 - h * pmax(0, 1 - abs(x - 0.5) / w)) / (1 - h * w), 0)
DENS <- list(
  A_N50_D22 = list(f = function(x) dmixn(x, .5, 2.2, 1), r = c(-6, 8)), A_N50_D3 = list(f = function(x) dmixn(x, .5, 3, 1), r = c(-6, 9)),
  A_N50_D4 = list(f = function(x) dmixn(x, .5, 4, 1), r = c(-6, 10)), A_N50_D5 = list(f = function(x) dmixn(x, .5, 5, 1), r = c(-6, 11)),
  A_N30_D35 = list(f = function(x) dmixn(x, .3, 3.5, 1), r = c(-6, 10)), A_N30_D5 = list(f = function(x) dmixn(x, .3, 5, 1), r = c(-6, 11)),
  A_NSC_D3 = list(f = function(x) dmixn(x, .5, 3, .5), r = c(-6, 6)),
  A_U_G05 = list(f = function(x) 0.5 * dunif(x) + 0.5 * dunif(x, 1.5, 2.5), r = c(-0.1, 2.6)),
  A_U_G1 = list(f = function(x) 0.5 * dunif(x) + 0.5 * dunif(x, 2, 3), r = c(-0.1, 3.1)),
  A_U_G3 = list(f = function(x) 0.5 * dunif(x) + 0.5 * dunif(x, 4, 5), r = c(-0.1, 5.1)),
  A_N50_D8 = list(f = function(x) dmixn(x, .5, 8, 1), r = c(-6, 14)),
  A_TIGHT = list(f = function(x) 0.5 * dnorm(x, 0, .05) + 0.5 * dnorm(x, 1, .05), r = c(-0.5, 1.5)),
  A_VAL_W10_H6 = list(f = function(x) dvalley(x, .1, .6), r = c(0, 1)), A_VAL_W10_H10 = list(f = function(x) dvalley(x, .1, 1), r = c(0, 1)),
  A_VAL_W25_H6 = list(f = function(x) dvalley(x, .25, .6), r = c(0, 1)), A_VAL_W25_H10 = list(f = function(x) dvalley(x, .25, 1), r = c(0, 1)))
modality <- function(dl) {
  x <- seq(dl$r[1], dl$r[2], length.out = 200001); f <- dl$f(x); fm <- max(f)
  # plateau-aware extrema: compress runs of equal values, then local maxima/minima of the run sequence
  rr <- rle(round(f / fm, 12)); v <- rr$values; L <- length(v)
  isMax <- vapply(seq_len(L), function(i) (i == 1 || v[i] > v[i - 1]) && (i == L || v[i] > v[i + 1]) && v[i] > 1e-9, logical(1))
  isMin <- vapply(seq_len(L), function(i) i > 1 && i < L && v[i] < v[i - 1] && v[i] < v[i + 1], logical(1))
  modes <- sum(isMax); depth <- if (any(isMin) && modes >= 2) 1 - min(v[isMin]) / min(v[isMax]) else 0
  Fx <- cumsum(f) / sum(f); q <- approx(Fx, x, xout = (seq_len(1e5) - 0.5) / 1e5, ties = "ordered", rule = 2)$y
  list(modes = modes, valley_depth = depth, pop_dip = dip(q))
}
mod <- lapply(DENS, modality)

res <- list(); rows <- list()
for (n in NS) {
  p <- readRDS(file.path(pdir, sprintf("power_n%d.rds", n)))
  for (law in names(p$laws)) {
    L <- p$laws[[law]]; N <- p$N
    for (a in ALPHAS) for (t in TESTS) {
      k <- L$counts[t, a]
      r1na <- is.na(p$r1[[a]]) && t == "R1"
      ci <- cp(k, N)
      pr <- L$paired_vs_R1[t, a, , ]
      nc <- if (sum(pr) > 0 && t != "R1") newcombe10(pr["R1yes", "Tyes"], pr["R1yes", "Tno"], pr["R1no", "Tyes"], pr["R1no", "Tno"]) else c(NA, NA, NA)
      structural <- t == "RWd" && n < 23
      rows[[length(rows) + 1]] <- data.frame(n = n, law = law, alpha = as.numeric(a), test = t, k = k, N = N, rate = if (r1na) NA else k / N,
                                             lo = ci[1], hi = ci[2], d_R1 = nc[1], d_lo = nc[2], d_hi = nc[3],
                                             status = if (r1na) "no test" else if (structural) "structural zero" else "estimated")
    }
  }
}
tab <- do.call(rbind, rows)
saveRDS(list(table = tab, modality = mod), paste0(outp, ".rds"))
fmt <- function(x) ifelse(is.na(x), "-", sprintf("%.3f", x))
con <- file(paste0(outp, ".md"), "w")
writeLines("# T15-E power tables (generated by code/r/mh_power_analysis.R; 99% CP for rates, Newcombe-10 99% for R1 - X)\n", con)
writeLines("## Modality of the alternatives (numerical)\n| law | modes | valley depth | population dip |\n|---|---|---|---|", con)
for (nm in names(mod)) writeLines(sprintf("| %s | %d | %.3f | %.4f |", nm, mod[[nm]]$modes, mod[[nm]]$valley_depth, mod[[nm]]$pop_dip), con)
for (a in ALPHAS) {
  writeLines(sprintf("\n## alpha = %s: rejection rates (TAB = NOT valid over unimodal laws)\n", a), con)
  writeLines(paste0("| law | n | ", paste(TESTS, collapse = " | "), " | R1-DWn [99%] | R1-RWa [99%] |"), con)
  writeLines(paste0("|", paste(rep("---", length(TESTS) + 4), collapse = "|"), "|"), con)
  for (law in unique(tab$law)) for (n in NS) {
    s <- tab[tab$law == law & tab$n == n & tab$alpha == as.numeric(a), ]
    r <- vapply(TESTS, function(t) { z <- s[s$test == t, ]; if (z$status == "structural zero") "0 (str.)" else fmt(z$rate) }, "")
    dn <- s[s$test == "DWn", ]; dr <- s[s$test == "RWa", ]
    writeLines(sprintf("| %s | %d | %s | %s [%s, %s] | %s [%s, %s] |", law, n, paste(r, collapse = " | "),
                       fmt(dn$d_R1), fmt(dn$d_lo), fmt(dn$d_hi), fmt(dr$d_R1), fmt(dr$d_lo), fmt(dr$d_hi)), con)
  }
}
close(con)
cat("wrote", paste0(outp, ".md"), "rows", nrow(tab), "\n")
