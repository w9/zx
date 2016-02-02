hash_vec <- function(hash_table, vec) {
  out_vec <- vector(length = length(vec))
  for (i in 1:length(vec)) {
    out_vec[i] <- hash_table[[vec[i]]]
  }
  out_vec
}

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



scale_ <- function(v, min_v=0, max_v=1) {
  (v-min(v))/(max(v)-min(v)) * (max_v-min_v) + min_v
}

get_slope <- function(x, y) lm(y ~ x, data.frame(x, y))$coefficients['x']

rand_measure <- function(a, b) {
  mean(apply(combn(1:length(b),2), 2, function(x)(b[x[1]]==b[x[2]])==(a[x[1]]==a[x[2]])))
}


log_trans <- function (rw) 
{
  rw[rw == 0] <- NaN
  rw <- log(rw)
}


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


nmf_ <- function(rwl=NULL, rw=NULL, rank=2:5, nrun=30, method="brunet", .options="p32v3", seed=12345) {
  library(NMF)

  if (is.null(rwl)) {
    if (is.null(rw)) stop('At least one of rw or rwl needs to be presented.')
    rwl <- log_trans(rw)
  }

  nmf(rwl, rank=rank, nrun=nrun, method=method, .options=.options, seed=seed)
}
  
edgeR_ <- function(rw, phe) {
  library(edgeR)

  cds <- DGEList(rw, group=phe)
  cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
  cds <- calcNormFactors( cds )
  cds <- estimateCommonDisp( cds )
  ret <- exactTest(cds)$table
  ret <- ret[order(ret$PValue),]
}
