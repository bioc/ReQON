\name{FWSEplot}
\alias{FWSEplot}
\title{Plot reported vs. empirical quality.}
\description{
Plots reported vs. empirical quality scores. Also calculates and outputs Frequency-Weighted Squared Error (FWSE) and reports FWSE on the plot.
}
\usage{
FWSEplot(ErrRates, QualFreq, FWSE_out = TRUE, col = "blue", min_freq = 0.001, 
   lim = c(0, length(QualFreq) - 1), xlabel = "Reported Quality", 
   ylabel = "Empirical Quality", main_title = "Reported vs. Empirical Quality")
}
\arguments{
\item{ErrRates}{vector of empirical error rates on the Phred scale.}
\item{QualFreq}{vector of relative frequencies of quality scores.}
\item{FWSE_out}{option to output FWSE and report FWSE on the plot. Default = TRUE.}
\item{col}{color of plotted points. Default = "blue".}
\item{min_freq}{Any quality scores with relative frequency less than this vlue will be plotted with a solid point. Default = 0.001.}
\item{lim}{common axis limits for both the x-axis and y-axis. Default = c(0, length(QualFreq) - 1).}
\item{xlabel}{x-axis label. Default = "Reported Quality".}
\item{ylabel}{y-axis label. Default = "Empirical Quality".}
\item{main_title}{title.  Default = "Reported vs. Empirical Quality".}
}
\details{
FWSEplot plots the reported quality score against the empirical quality score and reports FWSE.  If the quality scores accurately reflect the probability of a sequencing error, then the points should fall close to the 45-degree line and FWSE should be close to zero.  If the input vectors are $ErrRatesBefore and $QualFreqBefore from the ReQON output, this function will create the bottom left diagnostic plot that is output from ReQON.  Similarly, if the input vectors are $ErrRatesAfter and $QualFreqAfter from ReQON output, then the bottom right diagnostic plot is created.

For more details and interpretation, see the vignette by: 
browseVignettes("ReQON").
}
\value{
In addition to the plot, FWSEplot calculates and outputs Frequency-Weighted Squared Error (FWSE), a measure of how close the points lie to the 45-degree line. 
}
\author{Christopher Cabanski \email{cabanski@email.unc.edu}}
\examples{
## Create relative frequency example data
require( stats )
after <- dnorm( c( 0:40 ), mean = 30, sd = 8 )
  after <- after / sum( after )
err_rate <- c( 0:40 ) + rnorm( 41, mean = 0, sd = 5) 
  err_rate[ which( err_rate < 0 ) ] <- (-1) * err_rate[ which( err_rate < 0 ) ]
  ## to guarantee that all values are positive

## plot and calculate FWSE
FWSEplot( err_rate, after, col = "red" )
}