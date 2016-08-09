#' Generate sample order from consensus matrix
#'
#' @import dplyr
#' @export
get_order_from_consensus <- function(consensus_mat) {
  consensus_mat %>% dist %>% hclust %>% as.dendrogram %>% order.dendrogram
}

#' Generate Consensus Heatmap
#'
#' @import ggplot2
#' @import scales
#' @import reshape2
#' @import magrittr
#' @import dplyr
#' @export
consensus_heatmap <-
  function(consensus_mat, phe=NULL, sample_names=NULL) {
    if (!is.null(sample_names)) {
      rownames(consensus_mat) <- sample_names
      colnames(consensus_mat) <- sample_names
    } else if (is.null(rownames(sample_names)) || is.null(colnames(sample_names))) {
      rownames(consensus_mat) <- 1:nrow(consensus_mat) %>% paste0('x', .)
      colnames(consensus_mat) <- 1:ncol(consensus_mat) %>% paste0('x', .)
    }

    dendro_order <- get_order_from_consensus(consensus_mat)

    samples_ordered <- rownames(consensus_mat)[dendro_order]

    ggdat <- consensus_mat %>%
      melt(c('s1', 's2'), value.name='consensus') %>%
      tbl_df %>%
      print %>%
      mutate(s1=factor(s1, samples_ordered)) %>%
      mutate(s2=factor(s2, samples_ordered))

    if (!is.null(phe)) {
      phe_ordered <- factor(phe[dendro_order])
      ggdat_colorbar <- data_frame(sample=samples_ordered, phe=phe)
    } else {
      ggdat_colorbar <- data_frame(sample=character(), phe=character())
    }

    ggplot() +
      geom_raster(aes(x=s1, y=s2, fill=consensus), ggdat) +
      geom_rug(aes(x=sample, color=phe), ggdat_colorbar) +
      labs(x='sample', y='sample') +
      ggtitle('Clustering consensus') +
      scale_fill_continuous(name=wrap_format(10)('percentage of runs that cluster the pair together'), labels=percent_format()) +
      scale_color_discrete(name='group') +
      theme(axis.text.x=element_text(angle=90))
	}
