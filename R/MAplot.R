#' MA plot
#'
#' This function allows you to take a log-transformed expression matrix, and specify one (in which case it will be compared to the median) or two samples, outputing an MA plot, where A = mean expression count or intensity and M = log fold difference
#' @param expression  A log-transformed expression matrix, with columns corresponding to samples and rows to genes or other features. 
#' @param sample1 Character string specifiying the name of a column in your expression matrix, corresponding to a sample
#' @param sample2 Character string specifiying the name of a column in your expression matrix, corresponding to a sample. Defaults to the median if not specified. 
#' @export
#' @examples 
#' MAplot(expr, "S8")
#' MAplot(expr, "S8", "S15")
#' MAplot(expr.normalised, "S8", "S15)


MAplot = function(expression, sample1, sample2=NULL){
  Median = apply(as.matrix(expression), 1, median)
  x= (as.matrix(expression)[,sample1])
  if(is.null(sample2)){y= (Median)
  name = "Median"
  } else{y= as.matrix(expression)[,sample2]
  name = sample2}
  M = x-y
  A = (x+y)/2
  df=data.frame(A,M)
  h = hoeffd(as.matrix(df))$D[1,2]
  alpha = 1/nrow(expression)
  ggplot(df, aes(x = A, y= M, fill = abs(M)),  environment = environment()) +geom_hline(yintercept = 0, colour = "grey") + geom_point(colour="black",pch=21, size=5) + theme_classic() + labs(x = "A", y = paste("M (Sample", sample1, "-", name, ")")) + scale_fill_gradient(high = "deeppink2", low = "white") + guides(fill =F)  + geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE, color = "black", size = 3) + facet_wrap(~paste("Hoeffding's D =", round(h, 2)))
}