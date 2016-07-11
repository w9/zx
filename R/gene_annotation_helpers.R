#' Anntate genes using NCBI gene summary
#'
#' @import dplyr
#' @import xml2
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
gene_annotation <-
  function(genes_, organism_='mouse', format_='markdown') {
    symbol2eg <- switch(organism_,
                        mouse = org.Mm.eg.db,
                        human = org.Hs.eg.db,
                        stop(sprintf('Error: Unrecognized organism %s.', organism_)))

    mapped_genes <- intersect(genes_, mappedkeys(symbol2eg))
    genes_entrez <- symbol2eg[mapped_genes] %>% as.list %>% unlist

    gene_summary_xml <- paste0(genes_entrez, collapse=',') %>%
      sprintf('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=%s', .) %>% GET

    gene_annotation <- data_frame(symbol=character(), summary=character())
    for (record in gene_summary_xml %>% content %>% xml_child %>% xml_children) {
      gene_symbol <- record %>% xml_find_all('./Name') %>% xml_contents %>% as.character
      gene_summary <- record %>% xml_find_all('./Summary') %>% xml_contents %>% as.character
      if (length(gene_summary) > 0) {
        gene_annotation <- gene_annotation %>% add_row(symbol=gene_symbol, summary=gene_summary)
      }
    }

    if (format_ == 'data_frame') {
      gene_annotation
    } else if (format_ == 'markdown') {
      temp_prefix <- tempfile()
      md_file <- sprintf('%s.md', temp_prefix)
      html_file <- sprintf('%s.html', temp_prefix)
      paste0('# **', 1:nrow(gene_annotation), '** ', gene_annotation$symbol, '\n\n', gene_annotation$summary, '\n\n') %>% write(md_file)
      rmarkdown::render(md_file)
      rstudioapi::viewer(html_file)
    }
  }