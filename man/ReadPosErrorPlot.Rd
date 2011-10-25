\name{ReadPosErrorPlot}
\alias{ReadPosErrorPlot}
\title{Plot distribution of errors by read position.}
\description{
Plots the number of sequencing errors by read position.
}
\usage{
ReadPosErrorPlot(ReadPosErrors, error_col = "blue", thresh = 1.5, 
   thresh_col = "cyan", xlabel = "Read Position", ylabel = "# Errors", 
   main_title = "Distribution of Errors by Read Position")
}
\arguments{
\item{ReadPosErrors}{vector of sequencing error counts by read position.}
\item{error_col}{color of line plotting the errors counts. Default = "blue".}
\item{thresh}{Threshold for identifying read positions with large numbers of errors, plotted as a horizontal dashed line. Threshold is set as ``thresh * (average number of errors per read position)''. Default = 1.5.}
\item{thresh_col}{color of threshold line. Default = "cyan".}
\item{xlabel}{x-axis label. Default = "Read Position".}
\item{ylabel}{y-axis label. Default = "# Errors".}
\item{main_title}{title.  Default = "Distribution of Errors by Read Position".}
}
\details{
ReadPosErrorPlot plots the distribution of sequencing errors by read position.  If the input vector is $ReadPosErrors from ReQON output, this function will create the top left diagnostic plot that is output from ReQON.

For more details and interpretation, see the vignette by: 
browseVignettes("ReQON").
}
\author{Christopher Cabanski \email{cabanski@email.unc.edu}}
\examples{
## Create data of error counts 
x <- c( 1:30 )
err <- x^2 + ( 30 - x )^1.6 + rnorm(30, 0, 100)

## plot errors by read position
ReadPosErrorPlot( err )
}