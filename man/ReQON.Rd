\name{ReQON}
\alias{ReQON}
\title{Recalibrating Quality Of Nucleotides}
\description{
Recalibrate the nucleotide quality scores of either single-end or paired-end next-generation sequencing data that has been aligned.
}
\usage{
ReQON(in_bam, out_bam, region, max_train = -1, SNP = "", 
   RefSeq = "", plotname = "", temp_files = 0)
}
\arguments{
\item{in_bam}{file name of sorted BAM file of single-end or paired-end aligned sequencing data.  The corresponding index file (.bai file) must be located in 
the same directory.}
\item{out_bam}{file name for output BAM file with original quality scores replaced with recalibrated quality scores.}
\item{region}{training region for recalibration, as \dQuote{chromosome:start-end}.  Cannot span more than one chromosome. See note.
	  
  Example: "chr1:1-10000".}
\item{max_train}{maximum number of nucleotides to include in training region.  Useful if you want to train on e.g. the first 5 million bases of chromosome 
10.  Default = -1 (use all nucleotides from training region).}
\item{SNP}{file of SNP locations to remove from training set before recalibration.  Text or Rdata file (with variable name \dQuote{snp}) with no header and two columns: [1] chromosome, [2] position.  See note.
	
  Default: do not remove any nucleotides from training set.}
\item{RefSeq}{file of reference sequence for training set to identify sequencing errors (i.e, nucleotide is error if it does not match RefSeq).  Text or Rdata file (with variable name \dQuote{ref}) with no header and three columns: [1] chromosome, [2] position, [3] reference nucleotide (A,C,G,T).  See note.
	
  Default: errors are nucleotides not matching major allele(s) for coverage > 2, removing all nucleotides at positions with coverage of 2 or less.}
\item{plotname}{file name for saving recalibration plots in pdf.  If not specified, plots will not be produced.}
\item{temp_files}{option for keeping temporary files. 
	
  0: (default) remove all temporary files. 

  1: keep temporary files in working directory.}
}
\details{
ReQON uses logistic regression to recalibrate the nucleotide quality scores 
of a sorted BAM file.  The BAM file contains either single-end or paired-end 
next-generation sequencing data that has been aligned using any alignment tool.  
For help with sorting and indexing BAM files in R, see Rsamtools.
	
ReQON also has the option to output diagnostic plots which show the 
effectiveness of the recalibration on the training set.  

For a detailed description of usage, output and images, see the vignette by: 
browseVignettes("ReQON").

ReQON utilizes various java tools provided by Picard.  For more information on Picard, 
see http://picard.sourceforge.net
}
\value{
ReQON returns a BAM file, replacing the original quality scores with the 
recalibrated quality scores in the QUAL field.  

ReQON also outputs a data object of diagnostic data from the training set
that is plotted in the output diagnostic plots. The object variables are:
\item{$ReadPosErrors}{vector of error counts by read position.}
\item{$QualFreqBefore}{relative frequency of quality scores before recalibration.
  The first element in the vector corresponds to a quality score of zero.}
\item{$QualFreqAfter}{relative frequency of quality scores after recalibration.
  The first element in the vector corresponds to a quality score of zero.}
\item{$ErrRatesBefore}{vector of empirical error rates before recalibration, reported on 
the Phred scale.  The first element in the vector corresponds to a quality score of zero.}
\item{$ErrRatesAfter}{vector of empirical error rates after recalibration, reported on 
the Phred scale.  The first element in the vector corresponds to a quality score of zero.}
\item{FWSE}{vector of Frequency-Weighted Squared Error (FWSE) values.  The first element
is FWSE before recalibration and the second element is FWSE after recalibration.}
}
\note{
Be aware of how the chromosomes are referenced when specifying the training 
region.  For example, one BAM file may require specifying \dQuote{10:1-2000} 
while another may need \dQuote{chr10:1-2000}.

If providing SNP or RefSeq files, computations will speed up if your file 
only covers the positions in the training region.  For example, if you set 
region = \dQuote{chr10:1-2000}, then we recommend only having rows 
corresponding to chr10:1-2000 in the RefSeq/SNP file.
}
\author{Christopher Cabanski \email{cabanski@email.unc.edu}}
\examples{
## Read in sample data from seqbias package
library( ReQON )
library( seqbias )
library( Rsamtools )
ref_fn <- system.file( "extra/example.fa", package = "seqbias" )
ref_f <- FaFile( ref_fn )
open.FaFile( ref_f )
reads_fn <- system.file( "extra/example.bam", package = "seqbias" )

## Set up file of reference sequence
seqs <- scanFa( ref_f )
len <- length( seqs[[1]] )
ref <- matrix( nrow = len, ncol = 3 )
ref[,1] <- rep( "seq1", len )
ref[,2] <- c( 1:len )
str <- toString( subseq( seqs[[1]], 1, len ) )
s <- strsplit( str, NULL )
ref[,3] <- s[[1]]
write.table( ref, file = "ref_seq.txt", sep = "\t", quote = FALSE,
   row.names = FALSE, col.names = FALSE )

## Recalibrate File
sorted <- sortBam( reads_fn, tempfile() )
indexBam( sorted )
reg <- paste( "seq1:1-", len, sep = "" )
diagnostics <- ReQON( sorted, "Recalibrated_example.bam", reg, 
   RefSeq = "ref_seq.txt", plotname = "Recalibrated_example_plots.pdf" )

#Remove temporary file
unlink( "ref_seq.txt" )
}