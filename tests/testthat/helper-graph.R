# Helpers shared by the test files.

has_two_cycle <- function(net) {
  any(paste(net$from, net$to) %in% paste(net$to, net$from))
}

is_dag <- function(net) {
  e <- unique(net[, c("from", "to")])
  while (nrow(e) > 0) {
    roots <- setdiff(e$from, e$to)
    if (length(roots) == 0) return(FALSE)
    e <- e[!e$from %in% roots, , drop = FALSE]
  }
  TRUE
}

# Three pairwise dependent variables (V1 -> V2, V1 -> V3, V2 -> V3); the best
# network is a complete DAG in which one node has two parents.
triangle_data <- function() {
  set.seed(3)
  n <- 300
  d <- data.frame(V1 = rnorm(n))
  d$V2 <- 0.5 * d$V1 + rnorm(n)
  d$V3 <- 0.5 * d$V1 + 0.5 * d$V2 + rnorm(n)
  d
}

# Linear Gaussian data from a random DAG in which every edge runs from a lower
# to a higher column index. About half of the nodes get no parents, so many
# columns are independent of all others.
sparse_dag_data <- function(p, n = 300, seed = 1, coef = 0.5) {
  set.seed(seed)
  x <- matrix(stats::rnorm(n * p), n, p)
  for (i in 2:p) {
    k <- sample(0:2, 1, prob = c(0.5, 0.25, 0.25))
    for (pa in sample.int(i - 1, min(k, i - 1))) {
      x[, i] <- x[, i] + coef * x[, pa]
    }
  }
  colnames(x) <- paste0("V", seq_len(p))
  as.data.frame(x)
}
