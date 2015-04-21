pca <- function(pr=NULL, rwl=NULL, rw=NULL, phe=NULL) {
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
  
  if (is.null(phe)) {
    ggdat <- data.frame(n=1:length(pr$x[,1]), x=pr$x[,1], y=pr$x[,2], x.t=jittered$x, y.t=jittered$y)
    ggplot(ggdat) +
      geom_point(aes(x=x, y=y)) +
      geom_text(aes(x=x.t, y=y.t, label=n), alpha=0.2) +
      labs(x=pcLabel(1),y=pcLabel(2))
  } else {
    ggdat <- data.frame(n=1:length(pr$x[,1]), x=pr$x[,1], y=pr$x[,2], x.t=jittered$x, y.t=jittered$y, phe=phe)
    ggplot(ggdat) +
      geom_point(aes(x=x, y=y, color=phe)) +
      geom_text(aes(x=x.t, y=y.t, label=n), alpha=0.2) +
      scale_color_brewer(type='qual', palette=6) +
      labs(x=pcLabel(1),y=pcLabel(2))
  }
}


