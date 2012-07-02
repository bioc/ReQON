ReadPosErrorPlot <- function(ReadPosErrors, startpos = 1, error_col = "blue", thresh = 1.5, thresh_col = "cyan",
   xlabel = "Read Position", ylabel = "# Errors", main_title = "Distribution of Errors by Read Position"){
   
    readlen <- length(ReadPosErrors)
    plot(c(startpos:(readlen-1+startpos)), ReadPosErrors, 'l', col = error_col, xlab = xlabel, ylab = ylabel) 
	title(main = list(main_title, cex=1))   
    lines(c(1:readlen), matrix(1,1,readlen)*thresh*(sum(ReadPosErrors)/readlen), col = thresh_col,
	   lwd = 2, lty = 2)
	   
}