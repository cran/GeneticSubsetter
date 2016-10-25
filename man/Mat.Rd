\name{Mat}
\alias{Mat}
\title{
Kinship Matrix Creator
}
\description{
This function creates a kinship matrix for a set of genotypes. This function is a simplified version of the function "A.mat" in the R package rrBLUP.
}
\usage{
Mat(genos)
}
\arguments{
  \item{genos}{
    A matrix of genotypes, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, or NA, where 0 represents a heterozygous marker, and NA represents missing data.  All data in this matrix must be numeric.
}
}
\value{
A matrix of kinship values between genotypes.
}
\references{
Endelman JB (2011) Ridge regression and other kernels for genomic selection with R package rrBLUP. The Plant Genome 4:250-255.
}
\author{
Ryan C. Graebner
}
\examples{
data(genotypes)
Mat(genotypes)
}
\keyword{ misc }
