#' Plot delta parameters from multinomial regression model
#'
#' Allows for visualization of multinomial regression models from spatial or non-spatial models
#'
#' @param fit A model fit returned by one of the fit_*_PG model functions
#'
#' @return a ggplot
#' @export
#' @import patchwork
#' @importFrom tidyr pivot_longer separate
#' @importFrom dplyr filter
#' @importFrom tidyselect everything
#' @import ggplot2
#' @examples 
#' # plot_deltas(fit)
plot_deltas <- function(fit)
{
  DELTA <- fit$DELTA
  v <- ncol(fit$W)
  K <- fit$K
  # name columns
  c_count = 1
  c_names <- rep("C",ncol(DELTA))
  for(k in 2:K)
  {
    for(vs in 1:v)
    {
      c_names[c_count] <- paste0("k",k,"_p",vs)
      c_count = c_count + 1
    }
  }
  colnames(DELTA) <- c_names
  DELTA <- as.data.frame(DELTA) 
  D_df_long <- tidyr::pivot_longer(DELTA,cols = tidyselect::everything(),names_to = "delta",values_to = "draw") 
  D_df_long <- tidyr::separate(data = D_df_long,col = "delta",into = c("k","p"),remove = FALSE)
  
  D_df_long_p1 <- dplyr::filter(D_df_long,p == "p1")
  D1_add <- data.frame(delta = rep("k1_p1",nrow(DELTA)),
                       k = rep("k1",nrow(DELTA)),
                       p = rep("p1",nrow(DELTA)),
                       draw = rep(0,nrow(DELTA)))
  D_df_long_p1 <- rbind(D_df_long_p1,D1_add)
  d_plot_p <- ggplot(data = D_df_long_p1, aes(x = k, y = draw, fill = delta)) +
    geom_boxplot() + 
    geom_hline(yintercept = 0,linetype = "dashed") + 
    theme_classic() + 
    ylab("p = 1") + 
    xlab(NULL) + 
    annotate("text", x = 1,y = 0.1,label = "Reference") + 
    theme(legend.position = "none")
  d_plots <- d_plot_p
  for(vs in 2:v)
  {
    D_df_long_p <- dplyr::filter(D_df_long,p == paste0("p",vs))
    D_add <- data.frame(delta = rep("k1_p1",nrow(DELTA)),
                         k = rep("k1",nrow(DELTA)),
                         p = rep(paste0("p",vs),nrow(DELTA)),
                         draw = rep(0,nrow(DELTA)))
    D_df_long_p <- rbind(D_df_long_p,D_add)
    d_plot_p <- ggplot(data = D_df_long_p, aes(x = k, y = draw, fill = delta)) +
      geom_boxplot() + 
      geom_hline(yintercept = 0,linetype = "dashed")+ 
      theme_classic() + 
      ylab(paste0("p = ",vs)) + 
      xlab(NULL)+ 
      theme(legend.position = "none")
    if(vs == v) d_plot_p <- d_plot_p + xlab("Cluster")
    d_plots <- d_plots + d_plot_p
  }
  d_plots <- d_plots + patchwork::plot_layout(nrow = v)
  return(d_plots)
}