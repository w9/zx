#' Generate Consensus Heatmap
#'
#' @param consensus_mat_ The consensus matrix. The row and column names of this matrix will be used in the plot.
#'   Example: \code{NMF::consensus(...)}.
#' @param phe_ Grouping vector.
#' @import ggplot2
#' @import scales
#' @import reshape2
#' @import dplyr
#' @export
consensus_heatmap <-
  function(consensus_mat_, phe_=NULL) {
    output <- list()

    dendro_order <- consensus_mat_ %>% dist %>% hclust %>% as.dendrogram %>% order.dendrogram
    samples_ordered <- rownames(consensus_mat_)[dendro_order]

    output$samples_ordered <- samples_ordered

    output$ggdat <- consensus_mat_[dendro_order, dendro_order] %>%
      melt(c('s1', 's2'), value.name='consensus') %>%
      tbl_df %>%
      mutate(s1=factor(s1, samples_ordered)) %>%
      mutate(s2=factor(s2, samples_ordered))

    if (!is.null(phe_)) {
      phe_ordered <- factor(phe_, phe_[dendro_order])
      output$ggdat_colorbar <- data_frame(sample=samples_ordered, phe=phe_)
    } else {
      output$ggdat_colorbar <- data_frame(sample=character(), phe=character())
    }

		output$plot <-
			function(title_='Clustering consensus', color_bar_scale_name_='group', axis_label_=F) {
				ggtheme <- 
					if (axis_label_) {
						theme(axis.text.x=element_text(angle=90))
					} else {
						theme_void() +
							theme(axis.title=element_text(),
										legend.title=element_text())
					}
				ggplot() +
					geom_raster(aes(x=s1, y=s2, fill=consensus), output$ggdat) +
					geom_rug(aes(x=sample, color=phe), output$ggdat_colorbar) +
					labs(x='sample', y='sample') +
					ggtitle(title_) +
					scale_fill_continuous(name=wrap_format(10)('percentage of runs that cluster the pair together'), labels=percent_format()) +
					scale_color_discrete(name=color_bar_scale_name_) +
					ggtheme
			}

		output
	}
