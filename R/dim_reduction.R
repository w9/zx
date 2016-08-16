#' @import Rtsne
#' @import dplyr
#' @export
dim_reduction <-
  function(x,
           what='all',
					 zp=T,
           transpose=T,
           id_col_name='id') {
  title <- sprintf('DR of %s (%s)', deparse(substitute(x)), what)

  if (transpose) { x <- t(x) }

  output <- data_frame(row_num=1:nrow(x))
  if (!is.null(rownames(x))) {
    output[[id_col_name]] <- rownames(x)
  }

  if (any(what == 'all')) {
    what <- c('pca', 'pca_scale', 'mds_cor', 'tsne', 'tsne_cor', 'tsne_abs_cor')
  }

  df_list <- list()
  coord_list <- list()
  for (method in what) {
    message(sprintf('* doing %s ... ', method))

    if (method == 'pca') {
      pr <- prcomp(x)$x
      df <- data_frame(pc1=pr[,1], pc2=pr[,2], pc3=pr[,3])
    } else if (method == 'pca_scale') {
      non_zeros <- apply(x, 2, function(x)var(x)!=0)
      if (sum(non_zeros) < ncol(x)) message(sprintf('%d of %d features have variance > 0.', sum(non_zeros), ncol(x)))
      pr <- prcomp(x[,non_zeros], scale.=T)$x
      df <- data_frame(pc_scale1=pr[,1], pc_scale2=pr[,2], pc_scale3=pr[,3])
    } else if (method == 'mds_cor') {
      mds <- cmdscale((1-cor(t(x)))^3, k=3)
      df <- data_frame(mds1=mds[,1], mds2=mds[,2], mds3=mds[,3])
    } else if (method == 'tsne') {
      ret <- Rtsne(x, dims=3, perplexity=min(30, floor((nrow(x)-1)/3)))$Y
      df <- data_frame(tsne1=ret[,1], tsne2=ret[,2], tsne3=ret[,3])
    } else if (method == 'tsne_cor') {
      ret <- Rtsne((1-cor(t(x)))^3, dims=3, is_distance=T, perplexity=min(30, floor((nrow(x)-1)/3)))$Y
      df <- data_frame(tsne_cor1=ret[,1], tsne_cor2=ret[,2], tsne_cor3=ret[,3])
    } else if (method == 'tsne_abs_cor') {
      ret <- Rtsne((abs(cor(t(x))))^3, dims=3, is_distance=T, perplexity=min(30, floor((nrow(x)-1)/3)))$Y
      df <- data_frame(tsne_cor1=ret[,1], tsne_cor2=ret[,2], tsne_cor3=ret[,3])
    }

    df_list[[method]] <- df
    coord_list[[method]] <- colnames(df)
  }

  out_df <- bind_cols(df_list)

  if (zp) {
    zp(out_df) %>%
      zp_coords_(coord_list)
  } else {
    out_df
  }
}
