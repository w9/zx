#' @export
sink_reset <- function(){
    for(i in seq_len(sink.number())){
        sink(NULL)
    }
}

#' @import rstudioapi
#' @export
browse <- function(x, use_viewer=F, effect=F...) {
	temp_f <- sprintf('%s.txt', tempfile())
	sink(temp_f)

	if (effect) {
		x
	} else {
		print(x, ...)
	}
	sink_reset()

	if (use_viewer) {
		viewer(temp_f)
	} else {
		browseURL(temp_f)
	}
}

#' @export
corner <- function(x, n=5, m=10) {
	print(x[1:n, 1:m])
}


#' @import dplyr
#' @export
sort_row <- function(x, decreasing=T) {
  x[apply(x, 1, mean) %>% sort(decreasing=decreasing) %>% names,]
}

#' @import dplyr
#' @export
sort_col <- function(x, decreasing=T) {
  x[, apply(x, 2, mean) %>% sort(decreasing=decreasing) %>% names]
}

#' @export
apply_if <- function(x, p, f) {
  if (p) {
    . %>% f
  } else {
    identity
  }
}
    

#' @import Rtsne
#' @import dplyr
#' @export
dim_reduction <-
  function(x,
           what='pca',
					 zp=F,
           transpose=T,
           id_col_name='id',
           additional=NULL) {
  if (transpose) { x <- t(x) }

  output <- data_frame(row_num=1:nrow(x))
  if (!is.null(rownames(x))) {
    output[[id_col_name]] <- rownames(x)
  }
  
  if (what == 'all') {
    what <- c('pca', 'mds_cor', 'tsne', 'tsne_cor')
  }

  if ('pca' %in% what) {
    message('* doing pca ... ')
    pr <- prcomp(x)$x
    output <- output %>% mutate(pc1=pr[,1], pc2=pr[,2], pc3=pr[,3])
  }

  if ('mds_cor' %in% what) {
    message('* doing mds_cor ... ')
    mds <- cmdscale((1-cor(t(x)))^3, k=3)
    output <- output %>% mutate(mds1=mds[,1], mds2=mds[,2], mds3=mds[,3])
  }

  if ('tsne' %in% what) {
    message('* doing tsne ... ')
    ret <- Rtsne(x, dims=3, perplexity=min(30, floor((ncol(rwl)-1)/3)))$Y
    output <- output %>% mutate(tsne1=ret[,1], tsne2=ret[,2], tsne3=ret[,3])
  }

  if ('tsne_cor' %in% what) {
    message('* doing tsne_cor ... ')
    ret <- Rtsne((1-cor(t(x)))^3, dims=3, is_distance=T, perplexity=min(30, floor((ncol(rwl)-1)/3)))$Y
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

#' @export
hash_vec <- function(hash_table, vec) {
  out_vec <- vector(length = length(vec))
  for (i in 1:length(vec)) {
    out_vec[i] <- hash_table[[vec[i]]]
  }
  out_vec
}


#' @export
bind_tbls <- function(list_of_dfs, id_col_name='list_id') {
  col_names <- lapply(list_of_dfs, colnames)
  common_col_names <- Reduce(function(x, y)intersect(x, y), col_names)
  
  result_df <- data.frame(row.names=1:sum(sapply(list_of_dfs, nrow)))
  for (col_name in common_col_names) {
    message(col_name)
    
    new_col <- NULL
    for (df in list_of_dfs) {
      if (class(df[[col_name]])=='factor') {
        message(col_name, '  is factor, convert to character')
        new_col <- c(new_col, as.character(df[[col_name]]))
      } else {
        new_col <- c(new_col, df[[col_name]])
      }
    }
    result_df[[col_name]] <- new_col 
  }
  
  list_names <- names(list_of_dfs)
  if (!is.null(list_names)) {
    message('sadfasf')
    new_col <- NULL
    for (list_name in list_names) {
      new_col <- c(new_col, rep(list_name, nrow(list_of_dfs[[list_name]])))
    }
    result_df[[id_col_name]] <- new_col
  }
  
  result_df
}



#' @export
scale_ <- function(v, min_v=0, max_v=1) {
  (v-min(v))/(max(v)-min(v)) * (max_v-min_v) + min_v
}

#' @export
get_slope <- function(x, y) lm(y ~ x, data.frame(x, y))$coefficients['x']

#' @export
rand_measure <- function(a, b) {
  mean(apply(combn(1:length(b),2), 2, function(x)(b[x[1]]==b[x[2]])==(a[x[1]]==a[x[2]])))
}



#' @export
pca <- function(pr=NULL, rwl=NULL, rw=NULL, phe=NULL, labels=F) {
  library(ggplot2)
  
  if (is.null(pr)) {
    if (is.null(rwl)) {
      if (is.null(rw)) stop('At least one of rw, rwl, or pr needs to be presented.')
      rwl <- log_trans(rw)
    }
    rwl <- rwl[apply(rwl, 1, function(x)any(x>0)),]
    pr <- prcomp(t(rwl))
  }
  
  pcLabel <- function(i)paste0('PC',i,' (variance explained = ', sprintf('%.2f', pr$sdev[i]/sum(pr$sdev)*100), '%)')
  
  ggdat <- data.frame(n=1:length(pr$x[,1]), x=pr$x[,1], y=pr$x[,2])
  
  p <- ggplot(ggdat)
    
  if (is.null(phe)) {
    p <- p + geom_point(aes(x=x, y=y))
  } else {
    p <- p + geom_point(aes(x=x, y=y, color=phe)) + scale_color_brewer(type='qual', palette=6)
  }
  
  if (labels) {
    library(FField)
    jittered <- FFieldPtRep(cbind(pr$x[,1], pr$x[,2]))
    p <- p + geom_text(aes(x=jittered$x, y=jittered$y, label=n), alpha=0.2)
  }
  
  p <- p + labs(x=pcLabel(1),y=pcLabel(2))
  
  p
}

#' @export
vj <- function(rwl=NULL, rw=NULL, phe=NULL, point_alpha=0.1, violin_alpha=0.4, jitter_width=0.4) {
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  
  if (is.null(rwl)) {
    if (is.null(rw)) stop('At least one of rw or rwl needs to be presented.')
    rwl <- log_trans(rw)
  }
  
  if (is.null(phe)) {
    ggdat <- melt(data.matrix(rwl), varnames=c('gene', 'sample')) %>%
      mutate(sample=factor(sample), gene=factor(gene)) %>%
      filter(value>=.Machine$double.eps)
    p <- ggplot(ggdat) + geom_point(aes(x=sample, y=value), position=position_jitter(width=jitter_width), alpha=point_alpha)
  } else {
    names(phe) <- colnames(rwl)
    ggdat <- melt(data.matrix(rwl), varnames=c('gene', 'sample')) %>%
      mutate(sample=factor(sample), gene=factor(gene)) %>%
      filter(value>=.Machine$double.eps) %>%
      mutate(phe=phe[sample])
    p <- ggplot(ggdat) + geom_point(aes(x=sample, y=value, color=phe), position=position_jitter(width=jitter_width), alpha=point_alpha)
  }
  
  p <- p + geom_violin(aes(x=sample, y=value), alpha=violin_alpha) +
    labs(x='Sample', 'Log Expression')
  
  p
}


#' @export
smoothDensity <- function(rwl=NULL, rw=NULL) {
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  
  if (is.null(rwl)) {
    if (is.null(rw)) stop('At least one of rw or rwl needs to be presented.')
    rwl <- log_trans(rw)
  }
  
  ggdat <- melt(data.matrix(rwl), varnames=c('gene', 'sample')) %>%
    filter(value>=.Machine$double.eps)
  
  ggplot(ggdat) +
    geom_density(aes(color=sample, x=value))
}


#' @export
zero_percentage_distribution <- function(rw=NULL, jitter_width=0.05) {
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  
  ggdat <- melt(data.matrix(rw), varnames=c('gene', 'sample')) %>%
    group_by(sample) %>%
    summarize(pct_zero=mean(value==0))
  
  ggplot(ggdat) +
    geom_point(aes(x=1, y=pct_zero), position=position_jitter(width=jitter_width), alpha=0.5) +
    geom_violin(aes(x=1, y=pct_zero), alpha=0.5)
}

########## PACKAGES ##########



#' @export
nmf_ <- function(rwl=NULL, rw=NULL, rank=2:5, nrun=30, method="brunet", .options="p32v3", seed=12345) {
  library(NMF)

  if (is.null(rwl)) {
    if (is.null(rw)) stop('At least one of rw or rwl needs to be presented.')
    rwl <- log_trans(rw)
  }

  nmf(rwl, rank=rank, nrun=nrun, method=method, .options=.options, seed=seed)
}
  

#' @export
edgeR_ <- function(rw, phe) {
  library(edgeR)

  cds <- DGEList(rw, group=phe)
  cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
  cds <- calcNormFactors( cds )
  cds <- estimateCommonDisp( cds )
  ret <- exactTest(cds)$table
  ret <- ret[order(ret$PValue),]
}
