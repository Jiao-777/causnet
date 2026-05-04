#' Causnet with Stability Selection
#'
#' A variant of \code{\link{causnet}} that replaces the single p-value
#' threshold in parent identification with stability selection
#' (Meinshausen & Buhlmann, 2010), providing a finite-sample bound on
#' the per-family error rate (PFER).
#'
#' @param data A data matrix or data frame, with rows as observations
#'   and columns as variables (nodes of the network).
#' @param max_parents Maximum number of parents per node. Default 2.
#' @param score_fun Scoring function. If \code{NULL}, defaults to
#'   \code{score_bic_lm}.
#' @param tau Target PFER, i.e. the expected number of falsely selected
#'   parents to tolerate. Must be positive. Default 1.
#' @param pi_thr Stability threshold. Must lie in (0.5, 1). A candidate
#'   parent is retained if its empirical selection frequency across
#'   subsamples is at least \code{pi_thr}. Default 0.75.
#' @param B Number of subsamples. Default 100.
#' @param seed Optional random seed for reproducibility.
#'
#' @return An object of class \code{causnet}, identical in structure to
#'   the output of \code{\link{causnet}}.
#'
#' @details
#' This function uses \code{find_possible_parents_stable} for the parent
#' identification step, and inherits all downstream steps from the
#' original \code{causnet} pipeline. The user-facing parameters
#' \code{tau}, \code{pi_thr}, and \code{B} replace the original
#' \code{alpha}; the per-subsample selection size \code{q} is derived
#' internally from Theorem 1 of Meinshausen & Buhlmann (2010):
#' \deqn{q = \lfloor \sqrt{\tau \cdot p \cdot (2 \pi_{thr} - 1)} \rfloor}
#'
#' @references
#' Meinshausen, N. and Buhlmann, P. (2010). Stability selection.
#' Journal of the Royal Statistical Society: Series B, 72(4), 417-473.
#'
#' @export
causnet_stable <- function(data,
                           max_parents = 2,
                           score_fun = NULL,
                           tau = 1,
                           pi_thr = 0.75,
                           B = 100,
                           seed = NULL) {

  # ---- Argument validation ----
  if (pi_thr <= 0.5 || pi_thr >= 1) {
    stop("`pi_thr` must be in (0.5, 1).", call. = FALSE)
  }
  if (!is.numeric(tau) || tau <= 0) {
    stop("`tau` must be a positive numeric value.", call. = FALSE)
  }
  if (!is.numeric(B) || B < 1) {
    stop("`B` must be a positive integer.", call. = FALSE)
  }

  if (is.null(score_fun)) {
    score_fun <- score_bic_lm
  }

  # ---- Step 1: no-parent score ----
  ms <- no_parent_score(data, score_fun)

  # ---- Step 2: parent identification via stability selection ----
  possible_parents <- find_possible_parents_stable(
    data, tau = tau, pi_thr = pi_thr, B = B, seed = seed
  )

  # ---- Step 3 onwards: identical to original causnet pipeline ----
  possible_offspring <- find_possible_offspring(possible_parents)
  possible_parent_sets <- find_possible_parent_sets(
    possible_parents, max_parents
  )
  possible_parent_set_scores <- score_possible_parent_sets(
    data, possible_parent_sets, score_fun
  )
  bps1 <- best_possible_parent_sets(
    ms, possible_parent_sets, possible_parent_set_scores, max_parents
  )
  best_sinks <- find_best_sinks(
    possible_parents, ms, possible_offspring,
    possible_parent_sets, bps1, max_parents
  )
  bnsets <- find_best_network(best_sinks)

  out <- sink2net(
    bnsets, possible_parents, possible_parent_sets, bps1, max_parents
  )
  out[["mutliple_networks"]] <- list(
    bnsets = bnsets,
    possible_parents = possible_parents,
    possible_parent_sets = possible_parent_sets,
    bps1 = bps1,
    max_parents = max_parents
  )
  names(out[[1]])[1:2] <- c("from", "to")
  class(out) <- "causnet"
  out
}
