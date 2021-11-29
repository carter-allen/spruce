#' Sum neighboring psis in spot i
#' @param Psi an n x 1 vector of component k psis
#' @param ai the ith row (or column) of the adjacency matrix
psi_sums <- function(ai,Psi)
{
  Psi = as.vector(Psi)
  return(sum(Psi[ai == 1]))
}

#' Sum all neighboring psis
#' @param Psi an n x 1 vector of component k psis
#' @param A an n x n adjacency matrix
get_psi_sums <- function(Psi, A)
{
  return(apply(A, 1, psi_sums, Psi = Psi))
}