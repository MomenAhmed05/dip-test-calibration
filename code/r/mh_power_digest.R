# Read-only digest of the T15-E summary (for the result file). Usage: Rscript --vanilla mh_power_digest.R <summary.rds>
args <- commandArgs(trailingOnly = TRUE)
s <- readRDS(args[1]); tab <- s$table
cat("== modality ==\n"); for (nm in names(s$modality)) with(s$modality[[nm]], cat(sprintf("%-14s modes=%d depth=%.3f popdip=%.4f\n", nm, modes, valley_depth, pop_dip)))
for (a in c(0.05, 0.1, 0.01)) {
  cat(sprintf("\n== alpha %.2f: alternatives (rates %%) at n = 10/20/30/49: R1 | LF | DWa | DWn | RWa | RWd | TAB ==\n", a))
  for (law in unique(tab$law[grepl("^A_", tab$law)])) {
    out <- c()
    for (n in c(10, 20, 30, 49)) {
      z <- tab[tab$law == law & tab$n == n & tab$alpha == a, ]
      g <- function(t) { r <- z$rate[z$test == t]; if (is.na(r)) "-" else sprintf("%.1f", 100 * r) }
      out <- c(out, sprintf("n%d %s|%s|%s|%s|%s|%s|%s", n, g("R1"), g("LF"), g("DWa"), g("DWn"), g("RWa"), g("RWd"), g("TAB")))
    }
    cat(sprintf("%-14s %s\n", law, paste(out, collapse = "  ")))
  }
  cat(sprintf("== alpha %.2f: nulls, max over n of each test's rate (%%) ==\n", a))
  for (law in unique(tab$law[grepl("^N_", tab$law)])) {
    z <- tab[tab$law == law & tab$alpha == a, ]
    cat(sprintf("%-10s %s\n", law, paste(vapply(c("R1", "LF", "DWa", "DWn", "RWa", "RWd", "TAB"), function(t) sprintf("%s %.2f", t, 100 * max(z$rate[z$test == t], na.rm = TRUE)), ""), collapse = "  ")))
  }
  zz <- tab[tab$alpha == a & grepl("^A_", tab$law) & tab$test %in% c("DWa", "DWn", "RWa", "RWd", "LF") & !is.na(tab$d_R1), ]
  cat(sprintf("== alpha %.2f: paired R1 - X over all alternative cells: min diff %.4f, min lower99 %.4f, cells with lower99 > 0: %d / %d ==\n",
              a, min(zz$d_R1), min(zz$d_lo), sum(zz$d_lo > 0), nrow(zz)))
  for (t in c("LF", "DWa", "DWn", "RWa", "RWd")) { w <- zz[zz$test == t, ]; cat(sprintf("   vs %s: max diff %.3f (law %s n %d); min diff %.4f; lower99>0 in %d / %d\n", t, max(w$d_R1), w$law[which.max(w$d_R1)], w$n[which.max(w$d_R1)], min(w$d_R1), sum(w$d_lo > 0), nrow(w))) }
}
