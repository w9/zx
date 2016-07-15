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
  function(genes_, organism_, format_='markdown', genes_with_summary_only_=F) {
    gene_ranking <- data_frame(symbol=genes_, rank=1:length(genes_))
    symbol2eg <- switch(organism_,
                        mouse = org.Mm.eg.db::org.Mm.egSYMBOL2EG,
                        human = org.Hs.eg.db::org.Hs.egSYMBOL2EG,
                        stop(sprintf('Error: Unrecognized organism %s.', organism_)))

    mapped_genes <- intersect(genes_, mappedkeys(symbol2eg))
    genes_entrez <- symbol2eg[mapped_genes] %>% as.list %>% unlist

    gene_summary_xml <- paste0(genes_entrez, collapse=',') %>%
      sprintf('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s', .) %>% GET

    gene_annotation <- data_frame(symbol=character(), summary=character())
    for (record in gene_summary_xml %>% content %>% xml_child %>% xml_children) {
      gene_symbol      <- record %>% xml_find_all('./Name')        %>% xml_contents %>% as.character
      gene_description <- record %>% xml_find_all('./Description') %>% xml_contents %>% as.character %>% sprintf('*%s*', .)
      gene_summary     <- record %>% xml_find_all('./Summary')     %>% xml_contents %>% as.character
      if (length(gene_description) > 0) {
        if (length(gene_summary) > 0) {
          gene_annotation <- gene_annotation %>% add_row(symbol=gene_symbol, summary=paste(gene_description, gene_summary))
        } else if (!genes_with_summary_only_) {
          gene_annotation <- gene_annotation %>% add_row(symbol=gene_symbol, summary=gene_description)
        }
      }
    }
    gene_annotation <-  gene_annotation %>% left_join(gene_ranking, by='symbol') %>% select(rank, symbol, summary)

    if (format_ == 'data_frame') {
      gene_annotation
    } else if (format_ == 'markdown') {
      temp_prefix <- tempfile()
      md_file <- sprintf('%s.md', temp_prefix)
      html_file <- sprintf('%s.html', temp_prefix)
      paste0('# **', gene_annotation$rank, '** ', gene_annotation$symbol, '\n\n', gene_annotation$summary, '\n\n') %>% write(md_file)
      rmarkdown::render(md_file)
      rstudioapi::viewer(html_file)
    }
  }
