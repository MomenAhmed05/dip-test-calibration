# Step 15 / DT14 WP-E (T15-E): global unimodality tests derived from the modehunt multiscale procedures
# (Duembgen & Walther 2008, "DW08", all intervals; Rufibach & Walther 2010, "RW10", block procedure).
#
# Rule (S derivation, to be replayed by P): reject unimodality iff some interval with a significant DECREASE lies
# entirely to the left of some interval with a significant INCREASE (an antimode must lie between them).
# Under a unimodal density with mode m: a strict decrease on I_dec forces a point > m inside I_dec, and a strict increase
# on I_inc forces a point < m inside I_inc; with sup I_dec <= inf I_inc this is impossible. Hence the rejection event is
# contained in the event that the simultaneous confidence statement fails, which has probability <= alpha (DW08 / RW10
# guarantee, critical values from the uniform). Signs: T_jk = sum(2 X_local - 1) > 0 <=> mass to the right <=> increase.
#
# Usage from the shell: Rscript --vanilla mh_global.R <mode: dev> ...  (functions are sourced by the drivers)

MH_LIB <- if (nzchar(Sys.getenv("DIPTEST_RLIB"))) Sys.getenv("DIPTEST_RLIB") else NULL
suppressMessages(library(modehunt, lib.loc = MH_LIB))

antimode_found <- function(Dp, Dm) {
  # Dp / Dm: matrices (or NULL) with columns = interval endpoints (lower, upper) as returned by modehunt
  if (is.null(Dp) || is.null(Dm) || length(Dp) == 0 || length(Dm) == 0) return(FALSE)
  Dp <- matrix(Dp, ncol = ncol(as.matrix(Dp)))
  Dm <- matrix(Dm, ncol = ncol(as.matrix(Dm)))
  lo_inc <- Dp[, 1]
  hi_dec <- Dm[, 2]
  any(outer(hi_dec, lo_inc, "<="))
}

mh_reject <- function(x, cv, method = c("all", "block"), with_add = TRUE) {
  method <- match.arg(method)
  res <- if (method == "all") modeHunting(x, crit.vals = cv, min.int = TRUE)
         else modeHuntingBlock(x, crit.vals = cv, min.int = TRUE)
  if (with_add) antimode_found(res$Dp, res$Dm) else antimode_found(res$Dp.noadd, res$Dm.noadd)
}
