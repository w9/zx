pca <- function(mat, phe, logged=F) {
  library(ggplot2)
  
  if (!logged) mat <- log(mat+1)
  
  pr <- prcomp(t(mat))
  out$pr <- pr
  ggdat <- data.frame(x=pr$x[,1], y=pr$x[,2], phe=phe)
  out$plot <- ggplot() +
    geom_point(aes(x=x, y=y, color=phe)) +
    scale_color_brewer(type='qual', palette=6)
  print(out$plot)
  
  out
}