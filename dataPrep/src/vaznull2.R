## Vazquez-style null: retain species, link count and total interactions.
## Link probabilities are proportional to observed marginal interaction totals.
## Marginal totals themselves are not fixed.
vaznull.fast <- function(web) {
  probability <- outer(rowSums(web), colSums(web)) / sum(web)^2
  finalmat <- matrix(0, nrow(web), ncol(web), dimnames=dimnames(web))
  n.links <- sum(web > 0)

  ## First give every species at least one link. Do not choose a link
  ## that makes it impossible to cover all species within the link budget.
  while (any(rowSums(finalmat) == 0) || any(colSums(finalmat) == 0)) {
    missing.rows <- rowSums(finalmat) == 0
    missing.cols <- colSums(finalmat) == 0
    remaining <- n.links - sum(finalmat > 0) - 1
    feasible <- outer(missing.rows, missing.cols, function(r, c) {
      (r | c) & pmax(sum(missing.rows) - r, sum(missing.cols) - c) <= remaining
    })
    cells <- which(feasible & finalmat == 0)
    cell <- cells[sample.int(length(cells), 1, prob=probability[cells])]
    finalmat[cell] <- 1
  }

  ## Add links to reach observed connectance.
  n.add <- n.links - sum(finalmat > 0)
  if (n.add > 0) {
    cells <- which(finalmat == 0)
    add <- cells[sample.int(length(cells), n.add, prob=probability[cells])]
    finalmat[add] <- 1
  }

  ## Allocate the remaining visits only among realized links.
  n.add <- sum(web) - n.links
  if (n.add > 0) {
    cells <- which(finalmat > 0)
    add <- sample.int(length(cells), n.add, replace=TRUE, prob=probability[cells])
    finalmat[cells] <- finalmat[cells] + tabulate(add, nbins=length(cells))
  }
  finalmat
}
