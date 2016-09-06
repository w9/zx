#' @import stringr
#' @import readr
#' @import dplyr
#' @export
read_mat <- function(filename, format='auto', row_names=T, col_names=T, ...) {
  if (format == 'auto') {
    if (filename %>% str_detect('\\.csv$')) {
      format <- 'csv'
    } else {
      format <- 'tsv'
    }
  }
  
  if (format == 'csv') {
    tb <- read_csv(filename, col_names=col_names, ...)
  } else {
    tb <- read_tsv(filename, col_names=col_names, ...)
  }
  
  tb <- tb %>% as.data.frame
  
  if (row_names) {
    rownames(tb) <- tb[[1]]
    tb[[1]] <- NULL
  }

  rw <- tb %>% data.matrix

  rw
}

#' @export
corner <- function(x, n=5, m=10) {
	x[1:n, 1:m]
}


#' @import purrr
#' @import dplyr
#' @export
row_sort <- function(x, p=mean, decreasing=T) {
  x[apply(x, 1, as_function(p)) %>% order(decreasing=decreasing),]
}

#' @import purrr
#' @import dplyr
#' @export
col_sort <- function(x, p=mean, decreasing=T) {
  x[, apply(x, 2, as_function(p)) %>% order(decreasing=decreasing)]
}

#' @import purrr
#' @import dplyr
#' @export
row_filter <- function(x, p) {
  x[apply(x, 1, as_function(p)),]
}

#' @import purrr
#' @import dplyr
#' @export
col_filter <- function(x, p) {
  x[, apply(x, 2, as_function(p))]
}
