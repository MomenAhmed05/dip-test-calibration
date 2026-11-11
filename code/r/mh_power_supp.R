# T15-E SUPPLEMENTARY (post-result; labelled supplementary): exact paired (McNemar-type) tests on the discordant
# counts, two-sided, Holm-adjusted over the family {5 valid comparators x 176 alternative cells} per alpha, familywise 1%.
# Also: (i) check that the defect-record summary and the corrected summary have identical rates; (ii) largest comparator advantage.
# Usage: Rscript --vanilla mh_power_supp.R <powerdir> <summary.rds> <old_summary.rds>
args <- commandArgs(trailingOnly = TRUE)
pdir <- args[1]; new <- readRDS(args[2])$table; old <- readRDS(args[3])$table
cat("rates identical old vs corrected:", isTRUE(all.equal(old$rate, new$rate)), " (k identical:", identical(old$k, new$k), ")\n")
NS <- c(8, 10, 12, 15, 20, 25, 30, 35, 40, 45, 49)
COMP <- c("LF", "DWa", "DWn", "RWa", "RWd")
for (a in c("0.1", "0.05", "0.01")) {
  rows <- list()
  for (n in NS) {
    p <- readRDS(file.path(pdir, sprintf("power_n%d.rds", n)))
    for (law in names(p$laws)[grepl("^A_", names(p$laws))]) for (t in COMP) {
      pr <- p$laws[[law]]$paired_vs_R1[t, a, , ]
      b <- as.numeric(pr["R1yes", "Tno"]); c <- as.numeric(pr["R1no", "Tyes"])
      pv <- if (b + c == 0) 1 else binom.test(b, b + c, 0.5)$p.value
      rows[[length(rows) + 1]] <- data.frame(n = n, law = law, test = t, b = b, c = c, p = pv)
    }
  }
  d <- do.call(rbind, rows); d$p_holm <- p.adjust(d$p, "holm")
  cat(sprintf("alpha %s: family %d tests; Holm-significant at FW 1%%: R1 better %d, comparator better %d\n", a, nrow(d),
              sum(d$p_holm <= 0.01 & d$b > d$c), sum(d$p_holm <= 0.01 & d$c > d$b)))
  for (t in COMP) cat(sprintf("   %s: R1 better %d / 176; comparator better %d\n", t, sum(d$test == t & d$p_holm <= 0.01 & d$b > d$c),
                              sum(d$test == t & d$p_holm <= 0.01 & d$c > d$b)))
  z <- new[new$alpha == as.numeric(a) & grepl("^A_", new$law) & new$test %in% COMP, ]
  i <- which.min(z$d_R1)
  cat(sprintf("   largest comparator advantage (any significance): %s %s n=%d by %.4f\n", z$test[i], z$law[i], z$n[i], -z$d_R1[i]))
}
