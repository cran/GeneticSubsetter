\name{HET}
\alias{HET}
\title{
Heterozygosity Calculator
}
\description{
This function calculates the Expected Heterozygosity (HET; called PIC in earlier versions and in the paper describing this package) of a set of genotypes.
}
\usage{
HET(data)
}
\arguments{
  \item{data}{
    A matrix of genotypes, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, or NA, where 0 represents a heterozygous marker, and NA represents missing data. Note that this coding is different from the earlier PicCalc, which cannot handle heterozygous markers.  All data in this matrix must be numeric.
}
}
\value{
The mean Heterozygosity of all markers for the given set of genotypes.
}
\author{
Ryan C. Graebner
}
\note{
The ability to recogize heterozygous markers was included in HET, resulting in a slightly different genotype coding scheme than the earlier PicCalc.
}
\examples{
data("genotypes")
HET(genotypes)
}
\keyword{ misc }
