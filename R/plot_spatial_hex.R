#' Plot hexagonal spatial expression data
#'
#' Allows for visualization of gene expression patterns across discrete bins
#'
#' @param coords An n x 2 data frame where n = # of cells
#' @param feature A vector of length n containing the feature to display
#' @param n_hex_bins The number of hexagonal bins to use for each of row/col
#'
#' @return a ggplot
#' @export
#' @import ggplot2
#' @examples 
#' n <- 100
#' x <- runif(n)
#' y <- runif(n)
#' z <- rnorm(n)
#' coords_df <- data.frame(X = x, Y = y)
#' hex_plt <- plot_spatial_hex(coords_df,z)
#' plot(hex_plt)
plot_spatial_hex <- function(coords,feature, n_hex_bins = 20)
{
  coords_df <- as.data.frame(coords)
  colnames(coords_df) <- c("X","Y")
  coords_df$Z <- feature
  
  hex_plt <- ggplot(data = coords_df, aes(x = .data$X, y = .data$Y,z = .data$Z)) + 
    stat_summary_hex(bins = n_hex_bins) + 
    theme_void() + 
    theme(legend.position = "none") 
}
