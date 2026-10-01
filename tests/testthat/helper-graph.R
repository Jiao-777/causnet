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
