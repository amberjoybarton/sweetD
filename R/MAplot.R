utils::globalVariables(c("Sample1", "Sample2", "D", "Dataset"))

#' MA plot
#'
#' This function allows you to take a log-transformed expression matrix, and specify one (in which case it will be compared to the median) or two samples, outputing an MA plot, where A = mean expression count or intensity and M = log fold difference
#' @param expression  A log-transformed expression matrix, with columns corresponding to samples and rows to genes or other features.
#' @param sample1 Character string specifiying the name of a column in your expression matrix, corresponding to a sample
#' @param sample2 Character string specifiying the name of a column in your expression matrix, corresponding to a sample. Defaults to the median if not specified.
#' @export
#' @examples
#' data(expr.raw, expr.batchcorrected, expr.normalised)
#' MAplot(expr.raw, "S8")
#' MAplot(expr.raw, "S8", "S15")
#' MAplot(expr.normalised, "S8", "S15")


MAplot = function(expression, sample1, sample2=NULL){
  Median = apply(as.matrix(expression), 1, stats::median)
  x= (as.matrix(expression)[,sample1])
  if(is.null(sample2)){y= (Median)
  name = "Median"
  } else{y= as.matrix(expression)[,sample2]
  name = sample2}
  M = x-y
  A = (x+y)/2
  df=data.frame(A,M)
  h = Hmisc::hoeffd(as.matrix(df))$D[1,2]
  alpha = 1/nrow(expression)
  ggplot2::ggplot(df, ggplot2::aes(x = A, y= M, fill = abs(M)),  environment = environment()) +ggplot2::geom_hline(yintercept = 0, colour = "grey") + ggplot2::geom_point(colour="black",pch=21, size=5) + ggplot2::theme_classic() + ggplot2::labs(x = "A", y = paste("M (Sample", sample1, "-", name, ")")) + ggplot2::scale_fill_gradient(high = "deeppink2", low = "white") + ggplot2::guides(fill =F)  + ggplot2::geom_smooth(method = stats::lm, formula = y ~ splines::bs(x, 3), se = FALSE, color = "black", size = 3) + ggplot2::facet_wrap(~paste("Hoeffding's D =", round(h, 2)))
}
