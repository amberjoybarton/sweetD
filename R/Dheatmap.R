#' Hoeffding's D statistic against every other sample
#'
#' This function allows you to take one or more log-transformed expression matrices, and calculates Hoeffding's D-statistic for each sample in the matrix against every other sample. This outputs a dataframe with the values off Hoeffding's D for each sample compared to every other sample in each expression matrix.
#' @param ...  One or more log transformed matrices, with columns corresponding to samples and rows to genes or other features.
#' @export
#' @examples
#' data(expr.raw, expr.batchcorrected, expr.normalised)
#' sweetDall(expr.raw)
#' Dstat_results_all = sweetDall(expr.raw, expr.batchcorrected, expr.normalised)


sweetDall = function(...){
  all = list(...)
  names(all) <- as.character(substitute(c(...))[-1])
  print(names(all))
  D_statistics = as.data.frame(matrix(ncol = 4, nrow = 0, data = NA ))
  colnames(D_statistics) = c("Sample1", "Sample2", "D", "Dataset")
  for(i in names(all)){
    expression = all[[i]]
    print(paste("Calculating for dataset", i))
    Median = apply(as.matrix(expression), 1, stats::median)
    for(j in 1:length(colnames(expression))){
      for(k in 1:length(colnames(expression))){
        print(paste0(round(j/length(colnames(expression))*100), "%"))
        x= (as.matrix(expression)[,j])
        y= (as.matrix(expression)[,k])
        M = x-y
        A = (x+y)/2
        df=data.frame(A,M)
        h = Hmisc::hoeffd(as.matrix(df))$D[1,2]
        newline = t(as.data.frame(c(colnames(expression)[j], colnames(expression)[k], h,i)))
        colnames(newline) = c("Sample1", "Sample2", "D", "Dataset")
        rownames(newline) = make.names(paste(colnames(expression)[j], colnames(expression)[k], i, sep = "."))
        D_statistics = rbind(D_statistics, newline)
      }
    }
  }
  D_statistics = D_statistics[,c(4,1,2, 3)]
  D_statistics$D = as.numeric(as.character(D_statistics$D))
  D_statistics = D_statistics[order(D_statistics$D, decreasing =T),]
  return(D_statistics)
}


#' Plot Hoeffding's D statistic against every other sample
#'
#' This function takes the output of sweetDall, with the values off Hoeffding's D for each sample compared to every other sample in each expression matrix, and gives a heatmap of D statistics.
#' @param Dall_output  Output from sweetDall() performed on one or more log transformed matrices, with columns corresponding to samples and rows to genes or other features.
#' @export
#' @examples
#' data(expr.raw, expr.batchcorrected, expr.normalised)
#' Dstat_results_all = sweetDall(expr.raw, expr.batchcorrected, expr.normalised)
#' sweetDheatmap(Dstat_results_all)

sweetDheatmap = function(Dall_output){
  ggplot2::ggplot(Dall_output, ggplot2::aes(x = Sample1, y= Sample2),  environment = environment())  + ggplot2::facet_wrap(~Dataset, scale = "free_x") + ggplot2::geom_tile(ggplot2::aes(fill = D), colour = "white") + ggplot2::scale_fill_gradient(low = "white", high = "deeppink2") + ggplot2::theme_classic()
}


#' sweetD example dataset
#'
#' @name expr.raw
#' @docType data
#' @keywords data
NULL


#' sweetD example raw dataset
#'
#' @name expr.raw
#' @docType data
#' @keywords data
NULL

#' sweetD example normalised dataset
#'
#' @name expr.normalised
#' @docType data
#' @keywords data
NULL

#' sweetD example batch corrected dataset
#'
#' @name expr.batchcorrected
#' @docType data
#' @keywords data
NULL
