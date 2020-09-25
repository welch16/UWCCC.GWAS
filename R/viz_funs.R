#' Manhattan plot for the GWAS data
#' 
#' @param gwas_data a `data.frame` with at least the column chrom, pos, and
#'  p.value
#' @param sig a numeric value with the significance level
#' @param sig_color color for the significance threshold line
#' @return a `ggplot` object with the Manhattan plot
#' @param chr_colors colors to use per chromosome, if there are two colors is 
#'  going to color odd even chromosomes
#' 
#' @export
#' @importFrom magrittr `%>%` `%<>%`
#' @importFrom dplyr arrange select group_by mutate
#' @importFrom dplyr ungroup count distinct group_split
#' @importFrom ggplot2 ggplot aes geom_point geom_hline scale_size_continuous
#' @importFrom ggplot2 scale_x_continuous labs scale_color_manual 
#' @importFrom forcats fct_collapse
manhattan_plot <- function(gwas_data, sig = 5e-8, sig_color = "red",
  chr_colors = c("grey40", "darkgrey")) {

  `%<>%` <- magrittr::`%<>%`

  gwas_data %<>%
    dplyr::arrange(chrom, pos)

  gwas_data %<>%
    dplyr::group_by(chrom) %>%
    dplyr::mutate(
      bp_pos = (pos - min(pos)) / (max(pos) - min(pos)),
      logp = - log10(p),
      chrom_col = factor(as.character(chrom))) %>%
    dplyr::ungroup()

  if (length(chr_colors) == 2) {

    chroms <- dplyr::distinct(gwas_data, chrom) %>%
      dplyr::arrange(chrom) %>%
      dplyr::mutate(
        mod = dplyr::row_number() %% 2) %>%
      dplyr::group_split(mod)

    gwas_data %<>%
      dplyr::mutate(
        chrom_col = forcats::fct_collapse(
          chrom_col,
          odds = as.character(chroms[[1]]$chrom),
          even = as.character(chroms[[2]]$chrom)))
  }

  axis_set <- gwas_data %>%
    dplyr::count(chrom) %>%
    dplyr::mutate(
      w = n / sum(n),
      w_start = dplyr::lag(w, n = 1, 0),
      w_end = cumsum(w),
      w_start = cumsum(w_start))

  gwas_data %<>%
    dplyr::left_join(
      dplyr::select(axis_set, chrom, w_start, w_end), by = "chrom")

  gwas_data %<>%
    dplyr::mutate(
      bp_pos = (w_end - w_start) * bp_pos + w_start) %>%
    dplyr::select(-w_start, -w_end)

  axis_set %<>%
    dplyr::mutate(center = w_start + w / 2)

  yl <- abs(floor(max(gwas_data$logp)) + 10)

  mplot <- gwas_data %>%
    ggplot2::ggplot(
      ggplot2::aes(x = bp_pos, y = logp, color = chrom_col,
      size = logp)) +
  ggplot2::geom_point(alpha = 0.75) +
  ggplot2::geom_hline(yintercept = -log10(sig), color = sig_color,
    linetype = "dashed") +
  ggplot2::scale_x_continuous(
      label = axis_set$chrom, breaks = axis_set$center) +
  ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, yl)) +
  ggplot2::scale_color_manual(values = chr_colors) +
  ggplot2::scale_size_continuous(range = c(0.5, 3)) +
  ggplot2::labs(
    x = NULL,
    color = "chr",
    y = expression(-log[10](p)))

  chrom <- NULL
  pos <- NULL
  p <- NULL
  n <- NULL
  w <- NULL
  w_start <- NULL
  w_end <- NULL
  bp_pos <- NULL
  logp <- NULL
  mod <- NULL
  chrom_col <- NULL


  return(mplot)


}

#' generate a qqplot for GWAS studies
#' @param gwas_data a `data.frame` with at least the column chrom, pos, and
#'  p.value
#' @export
#' @importFrom dplyr select arrange
#' @importFrom stats ppoints
#' @importFrom ggplot2 labs
qqplot_pvalue <- function(gwas_data) {

  `%<>%` <- magrittr::`%<>%`

  plot_data <- gwas_data %>%
    dplyr::select(id, p) %>%
    dplyr::arrange(p)

  nobs <- nrow(plot_data)

  plot_data %<>%
    dplyr::mutate(
      obs = -log10(p),
      exp = -log10(ppoints(nobs)))

  log10Pe <- expression(paste("expected -log"[10], plain(p.value)))
  log10Po <- expression(paste("observed -log"[10], plain(p.value)))

  qqplot <- plot_data %>%
    ggplot2::ggplot(ggplot2::aes(x = exp, y = obs, size = obs)) +
    ggplot2::geom_point(alpha = .75, shape = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    ggplot2::scale_size_continuous(range = c(0.5, 1)) +
    ggplot2::labs(x = log10Pe, y = log10Po)

  id <- NULL
  p <- NULL
  obs <- NULL

  return(qqplot)

}