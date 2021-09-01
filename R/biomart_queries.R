#' obtain the rs ids in a range
#' @param chrom chromosome in numeric format
#' @param start starting position of the range
#' @param end ending position of the range
#' @param snp_mart a snp biomart obtained with `biomaRt::useEnsembl`
#' @param attributes a vector with the attributes to obtain. By default we get
#'   "refsnp_id", "allele", "chrom_start", and "chrom_strand"
#' @return a data.frame with the results of the query
#' @export
#' @importFrom biomaRt getBM
get_rsid <- function(chrom, start, end, snp_mart,
  attributes = c("refsnp_id", "allele", "chrom_start", "chrom_strand")) {
  biomaRt::getBM(
    attributes = attributes,
      filters = c("chr_name", "start", "end"),
      values = list(chrom, start, end),
      mart = snp_mart)
}
