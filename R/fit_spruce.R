#' Fit spruce Bayesian spatial mixture model
#'
#' This function allows you to detect sub-populations single-sample spatial transcriptomics experiments.
#' @param seurat_obj An integrated Seurat object 
#' @param K The number of sub-populations to infer. Each should be present in each sample.
#' @param emb Either one of "PCs", "HVGs", or "SVGs" OR a matrix with custom embeddings. If the latter, rows should be sorted as in meta data of Seurat object.
#' @param n_dim The number of dimensions to use if emb is specified as one of "PCs", "HVGs", or "SVGs". Ignored if emb is a matrix of custom embeddings.
#' @param r Spatial smoothing parameter. Should be greater than 0 with larger values enforcing stronger prior spatial association.
#' @param MCAR Logical. Include multivariate CAR random intercepts in gene expression model?
#' @param CAR Logical. Include univariate CAR random intercepts in multinomial gene expression model?
#' @param smooth Logical. Use manual spatial smoothing controled by r parameter?
#' @param nsim Number of total MCMC iterations to conduct. 
#' @param burn Number of initial MCMC iterations to discard as burn in. The number of saved iterations is nsim-burn
#' @param z_init Initialized cluster allocation vector to aid in MCMC convergence. If NULL z_init will be set using hierarchical clustering. 
#'
#' @keywords spatial transcriptomics Bayesian
#' @import Seurat
#' @export
#' @return A list of MCMC samples, including the MAP estimate of cluster indicators (z)
#' 
fit_spruce <- function(seurat_obj,
                       K,
                       emb = "PCs",
                       n_dim = 8,
                       r = 3,
                       MCAR = TRUE,
                       CAR = TRUE,
                       smooth = TRUE,
                       nsim = 2000,
                       burn = 1000,
                       z_init = NULL)
{
  if(emb == "PCs")
  {
    # check PCs are present
    if(!is.null(seurat_obj@reductions$pca))
    {
      Y <- seurat_obj@reductions$pca@cell.embeddings[,1:n_dim]
      rownames(Y) <- rownames(seurat_obj@reductions$pca@cell.embeddings)
    }
    else
    {
      message("Error: No PCA reductions found in the supplied Seurat object. Please add PCA embeddings to seurat_obj$reductions$pca@cell.embeddings.")
      return(NULL)
    }
  }
  else if(emb == "HVGs")
  {
    hvgs <- VariableFeatures(seurat_obj)[1:n_dim]
    Y <- t(seurat_obj@assays$SCT@scale.data[hvgs,])
  }
  else if(emb == "SVGs")
  {
    svgs <- SpatiallyVariableFeatures(seurat_obj)[1:n_dim]
    Y <- t(seurat_obj@assays$SCT@scale.data[svgs,])
  }
  else if(is.matrix(emb))
  {
    Y <- emb
  }
  else
  {
    message("Error: emb should be one of (PCs, HVGs, SVGs) or a matrix of custom embeddings with rows in same order as in seurat_obj@meta.data")
  }
  
  # compile coordinates
  L = length(seurat_obj@images)
  coords = NULL
  for(l in 1:L)
  {
    coords_x_l <- seurat_obj@images[[l]]@coordinates$col
    coords_y_l <- seurat_obj@images[[l]]@coordinates$row
    if(l > 1)
    {
      coords_x_l <- coords_x_l + max(seurat_obj@images[[l-1]]@coordinates$col) + 50
    }
    coords_l <- data.frame(x = coords_x_l,
                           y = coords_y_l)
    rownames(coords_l) <- rownames(seurat_obj@images[[l]]@coordinates)
    coords <- rbind(coords,coords_l)
  }
  
  # make W matrix
  W <- matrix(1,
              nrow = nrow(Y),
              ncol = 1)
  
  # check dimensions match
  if(nrow(Y) != nrow(coords))
  {
    message("Number of rows (cells) of embedding must match that of spatial coordinates")
    return(NULL)
  }
  else
  {
    # order coords according to rownames of Y
    coords <- coords[rownames(Y),]
    N <- nrow(Y)
    meta <- seurat_obj@meta.data[rownames(Y),]
  }
  
  # dispatch to spruce functions
  print("Dispatching to appropriate model fit function")
  if(emb %in% c("HVGs","SVGs"))
  {
    print("Fitting MSN model to account for skewness of HVGs/SVGs")
    fit <- spruce::fit_msn_PG_smooth(Y = Y,
                                     W = W,
                                     coords_df = coords,
                                     K = K,
                                     r = r,
                                     nsim = nsim, 
                                     burn = burn,
                                     z_init = z_init)
  }
  else
  {
    if((MCAR == TRUE) & (CAR == TRUE) & (smooth == TRUE))
    {
      fit <- spruce::fit_mvn_PG_CAR_MCAR_smooth(Y = Y,
                                                W = W,
                                                coords_df = coords,
                                                K = K,
                                                r = r,
                                                nsim = nsim, 
                                                burn = burn,
                                                z_init = z_init)
    }
    else if((MCAR == TRUE) & (CAR == TRUE) & (smooth == FALSE))
    {
      fit <- spruce::fit_mvn_PG_CAR_MCAR(Y = Y,
                                         W = W,
                                         coords_df = coords,
                                         K = K,
                                         nsim = nsim, 
                                         burn = burn,
                                         z_init = z_init)
    }
    else if((MCAR == TRUE) & (CAR == FALSE) & (smooth == TRUE))
    {
      fit <- spruce::fit_mvn_PG_MCAR_smooth(Y = Y,
                                            W = W,
                                            coords_df = coords,
                                            K = K,
                                            r = r,
                                            nsim = nsim, 
                                            burn = burn,
                                            z_init = z_init)
    }
    else if((MCAR == TRUE) & (CAR == FALSE) & (smooth == FALSE))
    {
      fit <- spruce::fit_mvn_PG_MCAR(Y = Y,
                                     W = W,
                                     coords_df = coords,
                                     K = K,
                                     nsim = nsim, 
                                     burn = burn,
                                     z_init = z_init)
    }
    else if((MCAR == FALSE) & (CAR == TRUE) & (smooth == TRUE))
    {
      fit <- spruce::fit_mvn_PG_CAR_smooth(Y = Y,
                                           W = W,
                                           coords_df = coords,
                                           K = K,
                                           r = r,
                                           nsim = nsim, 
                                           burn = burn,
                                           z_init = z_init)
    }
    else if((MCAR == FALSE) & (CAR == TRUE) & (smooth == FALSE))
    {
      fit <- spruce::fit_mvn_PG_CAR(Y = Y,
                                    W = W,
                                    coords_df = coords,
                                    K = K,
                                    nsim = nsim, 
                                    burn = burn,
                                    z_init = z_init)
    }
    else if((MCAR == FALSE) & (CAR == FALSE) & (smooth == TRUE))
    {
      fit <- spruce::fit_mvn_PG_smooth(Y = Y,
                                       W = W,
                                       coords_df = coords,
                                       K = K,
                                       r = r,
                                       nsim = nsim, 
                                       burn = burn,
                                       z_init = z_init)
    }
    else
    {
      fit <- spruce::fit_mvn_PG(Y = Y,
                                W = W,
                                K = K,
                                nsim = nsim, 
                                burn = burn,
                                z_init = z_init)
    }
  }
  return(fit)
}

