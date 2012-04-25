FWSEplot <- function(ErrRates, QualFreq, FWSE_out = TRUE, col = "blue", max_freq = 0.25,
   lim = c(0, length(QualFreq) - 1), collegend = TRUE, xlabel = "Reported Quality", 
   ylabel = "Empirical Quality", main_title = "Reported vs. Empirical Quality"){

    if(length(ErrRates) != length(QualFreq)){
	   stop("Vector lengths of ErrRates and QualFreq are not equal. Exiting.")
	}
   if(round(sum(QualFreq),5) != 1) {
      warning("WARNING! Before recalibration frequencies do not sum to 1.")
   }	
   if((col != "blue") & (col != "red")) {
      warning("WARNING! Color can only be red or blue. Resetting to blue.")
	  col = "blue"
   }	
   if(max_freq > 1) {
      warning("WARNING! max_freq > 1. Resetting to 1")
	  max_freq <- 1
   }   
   if(max_freq <= 0) {
      warning("WARNING! max_freq <= 0. Resetting to 0.01")
	  max_freq <- 0.01
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
	
	max_freq2 <- (1/max_freq)
	
    plot(c(0:maxqc), ErrRates, col=col, xlab = xlabel, ylab = ylabel, 
	   xlim = lim, ylim = lim, pch = 1)
	title(main = list(main_title, cex = 1))    
	for (i in 0:maxqc){
		if(col == "blue"){
			freqcol = rgb(max(0,1 - max_freq2*QualFreq[i+1]), max(0,1 - max_freq2*QualFreq[i+1]), 1)
		} else {
			freqcol = rgb(1, max(0,1 - max_freq2*QualFreq[i+1]), max(0,1 - max_freq2*QualFreq[i+1]))
		}
		points(i, ErrRates[i+1], pch = 20, col = freqcol)
		points(i, ErrRates[i+1], pch = 1, col = col)
	}
    lines(lim, lim, lty = 2)
	
	legcol <- c()
	if (col == "blue") {
		legcol[1] = rgb(1,1,1)
		legcol[2] = rgb(1-0.2,1-0.2,1)
		legcol[3] = rgb(1-0.4,1-0.4,1)
		legcol[4] = rgb(1-0.6,1-0.6,1)
		legcol[5] = rgb(1-0.8,1-0.8,1)
		legcol[6] = rgb(0,0,1)
	} else {
		legcol[1] = rgb(1,1,1)
		legcol[2] = rgb(1,1-0.2,1-0.2)
		legcol[3] = rgb(1,1-0.4,1-0.4)
		legcol[4] = rgb(1,1-0.6,1-0.6)
		legcol[5] = rgb(1,1-0.8,1-0.8)
		legcol[6] = rgb(1,0,0)
	}
	
	if(collegend){
    legend("bottomright", c("0",sprintf("%.2f",0.2*max_freq),sprintf("%.2f",0.4*max_freq),sprintf("%.2f",0.6*max_freq),
		sprintf("%.2f",0.8*max_freq),sprintf("> %.2f",max_freq)), title = "Frequency",cex=0.6, bty="o", fill=legcol);
	}
	if(FWSE_out){
	   text(lim[1], lim[2], sprintf("FWSE = %.2f", FWSE), adj = c(0,1), cex = 1.5)  
       FWSE
	}
		
}
