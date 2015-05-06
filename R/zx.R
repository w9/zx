pca_ <- function(pr=NULL, rwl=NULL, rw=NULL, phe=NULL) {
  library(ggplot2)
  library(FField)
  
  if (is.null(pr)) {
    if (is.null(rwl)) {
      if (is.null(rw)) stop('At least one of rw, rwl, or pr needs to be presented.')
      rwl <- log(rw + 1)
    }
    pr <- prcomp(t(rwl))
  }
  
  pcLabel <- function(i)paste0('PC',i,' (variance explained = ', sprintf('%.2f', pr$sdev[i]/sum(pr$sdev)*100), '%)')
  
  jittered <- FFieldPtRep(cbind(pr$x[,1], pr$x[,2]))
  
  ggdat <- data.frame(n=1:length(pr$x[,1]), x=pr$x[,1], y=pr$x[,2], x.t=jittered$x, y.t=jittered$y)
  ggplot(ggdat) +
    geom_point(aes(x=x, y=y, color=phe)) +
    geom_text(aes(x=x.t, y=y.t, label=n), alpha=0.2) +
    scale_color_brewer(type='qual', palette=6) +
    labs(x=pcLabel(1),y=pcLabel(2))
}

vj_ <- function(rwl=NULL, rw=NULL, phe=NULL, point_alpha=0.1, violin_alpha=0,4, jitter_width=0.4) {
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  
  if (is.null(rwl)) {
    if (is.null(rw)) stop('At least one of rw or rwl needs to be presented.')
    rwl <- log(rw + 1)
  }
  
  ggdat <- melt(data.matrix(rwl), varnames=c('gene', 'sample')) %>%
    filter(value>.Machine$double.eps)
  
  ggplot(ggdat) +
    geom_point(aes(x=sample, y=value, color=phe), position=position_jitter(width=jitter_width), alpha=point_alpha) +
    geom_violin(aes(x=sample, y=value), alpha=violin_alpha)
}

########## PACKAGES ##########


nmf_ <- function(rwl=NULL, rw=NULL, rank=2:5, nrun=30, method="brunet", .options="p32v3", seed=12345) {
  library(NMF)

  if (is.null(rwl)) {
    if (is.null(rw)) stop('At least one of rw or rwl needs to be presented.')
    rwl <- log(rw + 1)
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
