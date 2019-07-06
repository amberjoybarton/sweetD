#' Hoeffding's D statistic against the median for an expression matrix
#'
#' This function allows you to take one or more log-transformed expression matrices, and calculates Hoeffding's D-statistic for each sample against the median for that expression matrix. This outputs a dataframe with the values off Hoeffding's D for each sample in each expression matrix.
#' @param ...  One or more log transformed matrices, with columns corresponding to samples and rows to genes or other features. 
#' @export
#' @examples 
#' Dmedian(expr)
#' DstatResults = Dmedian(expr, expr.batchcorrected, expr.normalised)


Dmedian = function(...){
  all = list(...)
  names(all) <- as.character(substitute(c(...))[-1])
  print(names(all))
  D_statistics = as.data.frame(matrix(ncol = 3, nrow = 0, data = NA ))
  colnames(D_statistics) = c("Sample", "D", "Dataset")
  for(i in names(all)){
    expression = all[[i]]
    print(paste("Calculating for dataset", i))
    Median = apply(as.matrix(expression), 1, median)
    for(j in 1:length(colnames(expression)))
    { print(paste0(round(j/length(colnames(expression))*100), "%"))
      x= (as.matrix(expression)[,j])
      y= (Median)
      M = x-y
      A = (x+y)/2
      df=data.frame(A,M)
      h = hoeffd(as.matrix(df))$D[1,2]
      newline = t(as.data.frame(c(colnames(expression)[j], h,i)))
      colnames(newline) = c("Sample", "D", "Dataset")
      D_statistics = rbind(D_statistics, newline)  
    }
    
  }
  rownames(D_statistics) = make.names(paste(D_statistics$Dataset, D_statistics$Sample))
  D_statistics = D_statistics[,c(3,1,2)]
  D_statistics$D = as.numeric(as.character(D_statistics$D))
  D_statistics = D_statistics[order(D_statistics$D, decreasing =T),]
  return(D_statistics)
}


#' Plot Hoeffding's D statistic against the median for an expression matrix
#'
#' Generates a violin plot from the output of Dmedian to quickly visualise the distribution of D statistics.
#' @param   Dmedian_output The output after having performed Dmedian() on one or more log transformed matrices, with columns corresponding to samples and rows to genes or other features. 
#' @export
#' @examples 
#' DstatResults = Dmedian(expr, expr.batchcorrected, expr.normalised)
#' Dplot(DstatResults)

Dplot = function(Dmedian_output){
  ggplot(Dmedian_output, aes(x = Dataset, y= D, fill = D),  environment = environment()) + geom_violin(fill = "gold2", colour = NA) + geom_jitter(colour="black",pch=21, size=5) + theme_classic() + theme(axis.text.x=element_blank())+ labs(x = "", y = "Hoeffding's D statistic") + scale_fill_gradient(high = "deeppink2", low = "white") + guides(fill =F) + facet_wrap(~Dataset, scale = "free_x")
}
