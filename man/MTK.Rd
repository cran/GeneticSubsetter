\name{MTK}
\alias{MTK}
\title{
MTK calculator
}
\description{
This function calculates the Mean of Transformed Kinships (MTK) of a set of genotypes.
}
\usage{
MTK(genos, subset, mat = NULL, power = 10)
}
\arguments{
  \item{genos}{
    A matrix of genotypes, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, or NA, where 0 represents a heterozygous marker, and NA represents missing data. Note that this coding is different from the earlier MtkCalc, which cannot handle heterozygous markers.  All data in this matrix must be numeric.
}
  \item{subset}{
    A vector of genotype names for which to calculate MTK.
}
  \item{mat}{
    A kinship matrix, if one has already been computed for the population.
}
  \item{power}{
    The transformation that should be made to the kinship matrix, if the MTK criterion is used. If power=1, the kinship matrix is not transformed, if power=2, the kinship matrix is squared, etc.  When the power is higher, this function preferentially eliminates genotypes that are closely related to other genotypes in the population.
}
}
\value{
The MTK of the specified subset.
}
\author{
Ryan C. Graebner
}
\note{
The ability to recognize heterozygous markers was included in MTK, resulting in a slightly different genotype coding scheme than the earlier MtkCalc.
}
\examples{
data(genotypes)
MTK(genotypes,subset=colnames(genotypes[1:5]))
}
\keyword{ misc }
