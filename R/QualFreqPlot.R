QualFreqPlot <- function(QualFreqBefore, QualFreqAfter, before_col = "blue", after_col = "red",
   inc_legend = TRUE, xlabel = "Quality Score", ylabel = "Relative Frequency", 
   main_title = "Frequency Distributions of Quality Scores"){

   if(round(sum(QualFreqBefore),5) != 1) {
      warning("WARNING! Before recalibration frequencies do not sum to 1.")
   }
   if(round(sum(QualFreqAfter),5) != 1) {
      warning("WARNING! After recalibration frequencies do not sum to 1.")
   }
   
   lim <- max(QualFreqBefore, QualFreqAfter)
   maxqc <- (length(QualFreqBefore) - 1)
   plot(c(0:maxqc), QualFreqBefore, 'l', col = before_col,
      xlab = xlabel, ylab = ylabel, lwd = 2, ylim = c(-0.05*lim,1.05*lim))
   title(main = list(main_title, cex=1))     
   lines(c(0:maxqc), QualFreqAfter, col = after_col, lwd = 2)
   if(inc_legend){
      legend(0, 0.95*lim, c("Before", "After"), col = c(before_col, after_col),
         lty = c(1, 1),lwd = c(2,2))
	}
		   
}