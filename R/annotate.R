

#' annotate the `snps` with genes
#' @param snps a `tibble` with at least chrom and pos columns
#'   the chrom column is expected in `1,2,...X,Y` format
#' @param all_genes a `GRanges` object with the genes, named with `ENTREZID`
#' @return the snps `tibble` extended with annotated genes
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom stringr str_c
#' @importFrom dplyr pull inner_join left_join select rename_all
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr nest
#' @importFrom IRanges IRanges nearest
#' @importFrom S4Vectors queryHits subjectHits
#' @import AnnotationDbi
#' @importFrom snakecase to_snake_case
#' @importFrom magrittr set_names
annotate_genes <- function(snps, all_genes) {

  gr <- with(snps,
    GenomicRanges::GRanges(
      seqnames = stringr::str_c("chr", chrom),
      ranges = IRanges::IRanges(start = pos, width = 1))) %>%
    magrittr::set_names(dplyr::pull(snps, id))

  nms <- names(gr)
  nearest_genes <- IRanges::nearest(
    gr, all_genes, select = "all", ignore.strand = TRUE)

  snp_ids <- S4Vectors::queryHits(nearest_genes)
  gene_ids <- S4Vectors::subjectHits(nearest_genes)
  entrezids <- names(all_genes[gene_ids])

  gene_search <- AnnotationDbi::select(human,
    columns = c("SYMBOL"),
    keytype = "ENTREZID", keys = unique(entrezids)) %>%
    tibble::as_tibble() %>%
    dplyr::rename_all(list(snakecase::to_snake_case))

  out <- tibble(id = nms[snp_ids], entrezid = entrezids) %>%
    dplyr::left_join(gene_search, by = "entrezid") %>%
    dplyr::select(-entrezid)

  id <- NULL
  human <- NULL
  entrezid <- NULL
  symbol <- NULL

  return(
    snps %>%
      dplyr::inner_join(out, by = "id") %>%
      tidyr::nest(genes = c(symbol)))
}
