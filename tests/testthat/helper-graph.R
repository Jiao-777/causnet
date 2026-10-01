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
