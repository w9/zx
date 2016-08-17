#' @import Rtsne
#' @import dplyr
#' @import RDRToolbox
#' @export
dim_reduction <-
  function(x,
           what='all',
					 only_df=F,
           transpose=T,
           title=NULL,
           id_col_name='id') {
  if (is.null(title)) {
    title <- deparse(substitute(x))
  }

  if (transpose) { x <- t(x) }

  info <- data_frame(row_num=1:nrow(x))
  if (!is.null(rownames(x))) {
    info[[id_col_name]] <- rownames(x)
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
      df <- data_frame(pr[,1], pr[,2], pr[,3])
    } else if (method == 'pca_scale') {
      non_zeros <- apply(x, 2, function(x)var(x)!=0)
      if (sum(non_zeros) < ncol(x)) message(sprintf('%d of %d features have variance > 0.', sum(non_zeros), ncol(x)))
      pr <- prcomp(x[,non_zeros], scale.=T)$x
      df <- data_frame(pr[,1], pr[,2], pr[,3])
    } else if (method == 'mds_cor') {
      mds <- cmdscale((1-cor(t(x)))^3, k=3)
      df <- data_frame(mds[,1], mds[,2], mds[,3])
    } else if (method == 'tsne') {
      ret <- Rtsne(x, dims=3, perplexity=min(30, floor((nrow(x)-1)/3)))$Y
      df <- data_frame(ret[,1], ret[,2], ret[,3])
    } else if (method == 'tsne_cor') {
      ret <- Rtsne((1-cor(t(x)))^3, dims=3, is_distance=T, perplexity=min(30, floor((nrow(x)-1)/3)))$Y
      df <- data_frame(ret[,1], ret[,2], ret[,3])
    } else if (method == 'tsne_abs_cor') {
      ret <- Rtsne((1-abs(cor(t(x))))^3, dims=3, is_distance=T, perplexity=min(30, floor((nrow(x)-1)/3)))$Y
      df <- data_frame(ret[,1], ret[,2], ret[,3])
    } else if (method == 'isomap') {
      ret <- Isomap(t(x), 3)[[1]]
      df <- data_frame(ret[,1], ret[,2], ret[,3])
    }

    colnames(df) <- paste0(method, 1:3)
    df_list[[method]] <- df
    coord_list[[method]] <- colnames(df)
  }

  out_df <- bind_cols(info, df_list)

  if (only_df) {
    out_df
  } else {
    zp(out_df) %>%
      zp_options(title = sprintf('%s - (%s)', title, what %>% str_join(', '))) %>%
      zp_coords_(coord_list)
  }
}
