\name{QualFreqPlot}
\alias{QualFreqPlot}
\title{Plot frequency distributions of quality scores.}
\description{
Plots the relative frequency distribution of quality scores before and after recalibration.
}
\usage{
QualFreqPlot(QualFreqBefore, QualFreqAfter, before_col = "blue", 
   after_col = "red", inc_legend = TRUE, xlabel = "Quality Score", 
   ylabel = "Relative Frequency", 
   main_title = "Frequency Distributions of Quality Scores")
}
\arguments{
\item{QualFreqBefore}{vector of relative frequencies of quality scores before recalibration. The first element in the vector corresponds to a quality score of zero.}
\item{QualFreqAfter}{vector of relative frequencies of quality scores after recalibration. The first element in the vector corresponds to a quality score of zero.}
\item{before_col}{color of line plotting the frequency before recalibration. Default = "blue".}
\item{after_col}{color of line plotting the frequency after recalibration. Default = "red".}
\item{inc_legend}{option for including a legend. Default = TRUE.}
\item{xlabel}{x-axis label. Default = "Quality Score".}
\item{ylabel}{y-axis label. Default = "Relative Frequency".}
\item{main_title}{title.  Default = "Frequency Distributions of Quality Scores".}
}
\details{
QualFreqPlot plots the relative frequency distribution of quality scores before and after recalibration.  If the input vectors are $QualFreqBefore and $QualFreqAfter from ReQON output, this function will create the top right diagnostic plot that is output from ReQON.

For more details and interpretation, see the vignette by: 
browseVignettes("ReQON").
}
\author{Christopher Cabanski \email{ccabansk@genome.wustl.edu}}
\examples{
## Create data of frequencies
require( stats )
before <- dpois( c( 0:40 ), 40 )
  before <- before / sum( before )
after <- dnorm( c( 0:40 ), mean = 30, sd = 8 )
  after <- after / sum( after )

## plot
QualFreqPlot( before, after )
}