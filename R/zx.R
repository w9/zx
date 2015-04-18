pca <- function(mat, phe, logged=F) {
  library(ggplot2)
  library(FField)
  
  if (!logged) {
    rwl <- log(mat+1)
  } else {
    rwl <- mat
  }
  
  pr <- prcomp(t(rwl))
  
  pcLabel <- function(i)paste0('PC',i,' (variance explained = ', sprintf('%.2f', pr$sdev[i]/sum(pr$sdev)*100), '%)')
  
  jittered <- FFieldPtRep(cbind(pr$x[,1], pr$x[,2]))
  
  ggdat <- data.frame(n=1:length(pr$x[,1]), x=pr$x[,1], y=pr$x[,2], x.t=jittered$x, y.t=jittered$y)
  ggplot(ggdat) +
    geom_point(aes(x=x, y=y)) +
    geom_text(aes(x=x.t, y=y.t, label=1:nrow(ggdat)), alpha=0.2) +
    labs(x=pcLabel(1),y=pcLabel(2))
}
