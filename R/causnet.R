#' Find globally optimal causal network
#'
#' @details
#' Candidate parents are identified with a Pearson correlation test
#' (\code{cor.test()} p-value below \code{alpha}).
#'
#' Before the dynamic programme, nodes that have no candidate parents and are
#' not a candidate parent of any node are removed (the feasible set of the
#' published CausNet code). Such nodes cannot be part of any edge; the search
#' on the remaining nodes finds the same best sinks with the same scores.
#' \code{score_fun} is called on the remaining columns of \code{data}, and the
#' node indices it receives refer to those columns.
#'
#' @param data A data.frame of numeric variables.
#' @param alpha Numerical, p-value threshold to determine cutoff for possible
#'   parents.
#' @param max_parents Numeric, maximal number of parents. Defaults to 2.
#' @param score_fun A function used for scoring parent sets. For more details
#' see score_fun vignette.
#'
#' @return An object of class \code{causnet}, a list with 3 elements:
#' \code{network}, a data.frame with 3 columns \code{from}, \code{to} and
#' \code{component}; \code{n_best_parents}, the number of parents in the best
#' parent set of each sink; and \code{mutliple_networks}, internal objects
#' used by \code{\link{return_all_best_networks}} and
#' \code{\link{return_random_best_networks}}. Node indices in \code{from} and
#' \code{to} are column indices of \code{data}.
#' @export
#'
#' @examples
#' set.seed(1)
#' new_data <- simdat(n_var = 5)
#' out <- causnet(data = new_data, 0.5)
causnet <- function(data, alpha = 0.05, max_parents = 2, score_fun = NULL) {

  if (length(alpha) != 1 || !is.numeric(alpha) || alpha < 0 || alpha > 1) {
    stop("`alpha` must be a single numeric value in [0, 1].", call. = FALSE)
  }

  if (is.null(score_fun)) {
    score_fun <- score_bic_lm
  }

  possible_parents <- find_possible_parents(data, alpha)

  causnet_search(data, possible_parents, max_parents, score_fun)
}
