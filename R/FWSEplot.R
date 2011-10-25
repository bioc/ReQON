FWSEplot <- function(ErrRates, QualFreq, FWSE_out = TRUE, col = "blue", min_freq = 0.001, 
   lim = c(0, length(QualFreq) - 1), xlabel = "Reported Quality", ylabel = "Empirical Quality", 
   main_title = "Reported vs. Empirical Quality"){

    if(length(ErrRates) != length(QualFreq)){
	   stop("Vector lengths of ErrRates and QualFreq are not equal. Exiting.")
	}
   if(round(sum(QualFreq),5) != 1) {
      warning("WARNING! Before recalibration frequencies do not sum to 1.")
   }	
	
    maxqc <- (length(QualFreq) - 1)
	temp1 <- (ErrRates - c(0:maxqc))^2
	temp2 <- temp1 * QualFreq
	flag2 <- which(is.na(temp2))
	if (length(flag2)>0){
	  FWSE <- sum(temp2[-flag2])
	} else {
	  FWSE <- sum(temp2)
	}
	flag3 <- which(QualFreq > min_freq)
	
    plot(c(0:maxqc)[flag3], ErrRates[flag3], col=col, xlab = xlabel, ylab = ylabel, 
	   xlim = lim, ylim = lim, pch = 1)
	title(main = list(main_title, cex = 1))    
	points(c(0:maxqc)[-flag3], ErrRates[-flag3], pch = 20, col = col)   
    lines(lim, lim, lty = 2)
	
	if(FWSE_out){
	   text(lim[1], lim[2], sprintf("FWSE = %.2f", FWSE), adj = c(0,1), cex = 1.5)  
       FWSE
	}
		
}