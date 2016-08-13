#' @import Rtsne
#' @import dplyr
#' @export
dim_reduction <-
  function(x,
           what='all',
					 zp=T,
           transpose=T,
           id_col_name='id',
           additional=NULL) {
  if (transpose) { x <- t(x) }

  output <- data_frame(row_num=1:nrow(x))
  if (!is.null(rownames(x))) {
    output[[id_col_name]] <- rownames(x)
  }
  

  if (what == 'all' || 'pca' %in% what) {
    message('* doing pca ... ')
    pr <- prcomp(x)$x
    output <- output %>% mutate(pc1=pr[,1], pc2=pr[,2], pc3=pr[,3])
  }

  if (what == 'all' || 'pca_scale' %in% what) {
    message('* doing pca_scale ... ')
    non_zero_rows <- apply(x, 1, function(x)var(x)!=0)
    if (sum(non_zero_rows) < nrow(x)) message(sprintf('%d of %d rows have variance > 0.', sum(non_zero_rows), nrow(x)))
    pr <- prcomp(x[non_zero_rows,], scale.=T)$x
    output <- output %>% mutate(pc1=pr[,1], pc2=pr[,2], pc3=pr[,3])
  }

  if (what == 'all' || 'mds_cor' %in% what) {
    message('* doing mds_cor ... ')
    mds <- cmdscale((1-cor(t(x)))^3, k=3)
    output <- output %>% mutate(mds1=mds[,1], mds2=mds[,2], mds3=mds[,3])
  }

  if (what == 'all' || 'tsne' %in% what) {
    message('* doing tsne ... ')
    ret <- Rtsne(x, dims=3, perplexity=min(30, floor((nrow(x)-1)/3)))$Y
    output <- output %>% mutate(tsne1=ret[,1], tsne2=ret[,2], tsne3=ret[,3])
  }

  if (what == 'all' || 'tsne_cor' %in% what) {
    message('* doing tsne_cor ... ')
    ret <- Rtsne((1-cor(t(x)))^3, dims=3, is_distance=T, perplexity=min(30, floor((nrow(x)-1)/3)))$Y
    output <- output %>% mutate(tsne_cor1=ret[,1], tsne_cor2=ret[,2], tsne_cor3=ret[,3])
  }

  if (!is.null(additional)) {
    output <- output %>% bind_cols(additional)
  }

  if (zp) {
    zp_output <- zp(output)

		if ('pca' %in% what) {
      zp_output <- zp_output %>% zp_coord(pc1, pc2, pc3)
		}

		if ('mds_cor' %in% what) {
      zp_output <- zp_output %>% zp_coord(mds1, mds2, mds3)
		}

		if ('tsne' %in% what) {
      zp_output <- zp_output %>% zp_coord(tsne1, tsne2, tsne3)
		}

		if ('tsne_cor' %in% what) {
      zp_output <- zp_output %>% zp_coord(tsne_cor1, tsne_cor2, tsne_cor3)
		}

    zp_output
  } else{
    output
  }
}
