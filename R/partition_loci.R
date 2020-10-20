
#' Separates a `tibble` with the summary statistics of a GWAS experiment with 
#' columns id, pos, and p (snp id, snp position, and p.value)
#' @param data A `tibble` with at least id, pos and p columns
#' @param snp_dist an integer to find all SNPs within that distance from the 
#'   significant snp. By default 250 kbps.
#' @param pval_thr p.value significant threshold. By default 1e-8.
#' @return a nested `tibble` with the separated loci
#' @export
#' @importFrom dplyr filter inner_join distinct
#' @importFrom IRanges IRanges resize findOverlaps
#' @importFrom tibble tibble
#' @importFrom tidyr nest
#' @importFrom S4Vectors queryHits subjectHits
partition_by_snp <- function(data, snp_dist = 250e3, pval_thr = 1e-8) {

  sign_snps <- data %>%
    dplyr::filter(p <= pval_thr)

  if (nrow(sign_snps) > 1) {
    sign_snps <- with(sign_snps,
      IRanges::IRanges(start = pos, width = 1, names = id))
    all_snps <- with(data,
      IRanges::IRanges(start = pos, width = 1, names = id))

    sign_snps <- IRanges::resize(sign_snps, width = snp_dist, fix = "center")
    overlaps <- IRanges::findOverlaps(sign_snps, all_snps)

    outcome <- tibble::tibble(
      ref_snp = names(sign_snps)[S4Vectors::queryHits(overlaps)],
      id = names(all_snps)[S4Vectors::subjectHits(overlaps)]) %>%
      dplyr::inner_join(data, by = "id") %>%
      tidyr::nest(locus = c(pos, id, p)) %>%
      dplyr::distinct(locus)
  } else {
    warning("There are not significant SNPs, increase significant threshold")
    outcome <- NULL
  }

  pos <- id <- p <- locus <- NULL

  return(outcome)
}
