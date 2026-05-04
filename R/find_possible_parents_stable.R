#' Stability Selection for Parent Identification
#'
#' A stability-selection variant of \code{find_possible_parents} that
#' replaces the single p-value threshold with subsampling-based selection
#' frequencies and a per-family error rate (PFER) bound, following
#' Meinshausen & Buhlmann (2010).
#'
#' @param data A data matrix or data frame (n x p).
#' @param tau Target PFER (expected number of falsely selected parents).
#'   Default 1.
#' @param pi_thr Stability threshold. Must be in (0.5, 1). Default 0.75.
#' @param B Number of subsamples. Default 100.
#' @param seed Optional random seed.
#'
#' @return A list of length \code{ncol(data)}; element \code{i} contains
#'   the indices of stable parent candidates for node \code{i}.
#'
#' @details
#' For each node \eqn{V_i}, the function:
#' \enumerate{
#'   \item Computes \code{q} from Theorem 1 of Meinshausen & Buhlmann (2010):
#'     \deqn{q = \lfloor \sqrt{\tau \cdot (p-1) \cdot (2 \pi_{thr} - 1)} \rfloor}
#'   \item For \code{B} subsamples of size \eqn{\lfloor n/2 \rfloor}, drawn
#'     without replacement, ranks all other nodes by absolute Pearson
#'     correlation with \eqn{V_i} and retains the top \code{q}.
#'   \item Computes the empirical selection frequency \eqn{\hat\Pi_j^{(i)}}
#'     for each candidate \eqn{V_j}.
#'   \item Returns \eqn{\{j : \hat\Pi_j^{(i)} \geq \pi_{thr}\}}.
#' }
#'
#' @references
#' Meinshausen, N. and Buhlmann, P. (2010). Stability selection.
#' Journal of the Royal Statistical Society: Series B, 72(4), 417-473.
#'
#' @export
find_possible_parents_stable <- function(data,
                                         tau = 1,
                                         pi_thr = 0.75,
                                         B = 100,
                                         seed = NULL) {

  if (!is.null(seed)) set.seed(seed)
  if (pi_thr <= 0.5 || pi_thr >= 1) {
    stop("`pi_thr` must be in (0.5, 1).", call. = FALSE)
  }

  n     <- nrow(data)
  n_var <- ncol(data)
  n_sub <- floor(n / 2)

  # ---- Step 0: derive q from Theorem 1 ----
  p_candidate <- n_var - 1
  q <- floor(sqrt(tau * p_candidate * (2 * pi_thr - 1)))
  q <- max(q, 1)

  # ---- Step 1: B subsamples, accumulate selection counts ----
  selection_count <- matrix(0, nrow = n_var, ncol = n_var)

  for (b in seq_len(B)) {
    I_b <- sample.int(n, size = n_sub, replace = FALSE)
    data_sub <- data[I_b, , drop = FALSE]

    for (vertex in seq_len(n_var)) {
      scores <- numeric(n_var)
      for (parent in seq_len(n_var)) {
        if (parent != vertex) {
          scores[parent] <- abs(stats::cor(
            data_sub[, vertex, drop = TRUE],
            data_sub[, parent, drop = TRUE]
          ))
        } else {
          scores[parent] <- -Inf
        }
      }
      top_q <- order(scores, decreasing = TRUE)[seq_len(q)]
      selection_count[vertex, top_q] <-
        selection_count[vertex, top_q] + 1
    }
  }

  # ---- Step 2: empirical selection frequencies ----
  Pi <- selection_count / B

  # ---- Step 3: apply stability threshold ----
  out <- vector("list", n_var)
  for (vertex in seq_len(n_var)) {
    out[[vertex]] <- which(Pi[vertex, ] >= pi_thr)
    if (length(out[[vertex]]) == 0) out[[vertex]] <- NULL
  }

  message(sprintf(
    "[stability selection] n=%d, p=%d, q=%d, pi_thr=%.2f, tau=%.2f, B=%d",
    n, n_var, q, pi_thr, tau, B
  ))

  out
}
