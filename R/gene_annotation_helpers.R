#' Anntate genes using NCBI gene summary
#'
#' @export
#' @import xml2
#' @import httr
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import AnnotationDbi
#' @import dplyr
gene_annotation <-
  function(genes_, organism_, md_filename_=NULL, html_filename_=NULL, verbose_level_=1, show_in_viewer_=T) {
    gene_ranking <- data_frame(symbol=genes_, rank=1:length(genes_))
    symbol2eg <- switch(organism_,
                        mouse = org.Mm.eg.db::org.Mm.egSYMBOL2EG,
                        human = org.Hs.eg.db::org.Hs.egSYMBOL2EG,
                        stop(sprintf('Error: Unrecognized organism %s.', organism_)))

    mapped_genes <- intersect(genes_, mappedkeys(symbol2eg))
    genes_entrez <- symbol2eg[mapped_genes] %>% as.list %>% unlist

    if (verbose_level_ >= 2) message(sprintf('genes_entrez = (length %d)', length(genes_entrez)))
    if (verbose_level_ >= 2) print(head(genes_entrez, 50))

    gene_summary_xml <- paste0(genes_entrez, collapse=',') %>%
      sprintf('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s', .) %>% GET

    if (verbose_level_ >= 3) message('gene_summary_xml = ')
    if (verbose_level_ >= 3) print(gene_summary_xml)

    gene_annotation <- data_frame(symbol=character(), full_name=character(), summary=character())
    for (record in gene_summary_xml %>% content %>% xml_find_all('./*/DocumentSummary')) {
      gene_symbol      <- record %>% xml_find_all('./Name')        %>% xml_contents %>% as.character %>% ifelse(length(.) > 0, ., '')
      gene_description <- record %>% xml_find_all('./Description') %>% xml_contents %>% as.character %>% ifelse(length(.) > 0, ., '')
      gene_summary     <- record %>% xml_find_all('./Summary')     %>% xml_contents %>% as.character %>% ifelse(length(.) > 0, ., '')

      gene_annotation <- gene_annotation %>% add_row(symbol=gene_symbol, full_name=gene_description, summary=gene_summary)
    }
    gene_annotation <-  gene_annotation %>% left_join(gene_ranking, by='symbol') %>% select(rank, symbol, full_name, summary)

    temp_prefix <- tempfile()
    md_filename   <- ifelse(is.null(md_filename_  ), sprintf('%s.md'  , temp_prefix), md_filename_  )
    html_filename <- ifelse(is.null(html_filename_), sprintf('%s.html', temp_prefix), html_filename_)

    if (verbose_level_ >= 2) message(sprintf('md_filename = %s', md_filename))
    if (verbose_level_ >= 2) message(sprintf('html_filename = %s', html_filename))

    with(gene_annotation, paste0('# **', rank, '** ', symbol, '\n\n***', full_name, '***\n\n', summary, '\n\n')) %>% write(md_filename)

    rmarkdown::render(md_filename, output_file=basename(html_filename), output_dir=dirname(html_filename), quiet=ifelse(verbose_level_ >= 3, F, T))
    if (show_in_viewer_) rstudioapi::viewer(html_filename)

    gene_annotation
  }
