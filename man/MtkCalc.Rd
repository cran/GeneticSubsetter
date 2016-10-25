\name{MtkCalc}
\alias{MtkCalc}
\title{
MTK calculator (Depreciated)
}
\description{
***This function has been superseded by MTK.
}
\usage{
MtkCalc(genos, subset, power = 10)
}
\arguments{
  \item{genos}{
    A matrix of genotypes, that includes all genotypes that should be used to create the kinship matrix, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, where 0 represents missing data.  All data in this matrix must be numeric.
}
  \item{subset}{
    A vector of genotype names for which to calculate MTK.

}
  \item{power}{
    The transformation that should be made to the kinship matrix. If power=1, the kinship matrix is not transformed, if power=2, the kinship matrix is squared, etc.  When the power is higher, this function preferentially eliminates genotypes that are closely related to other specific genotypes.

}
}
\value{
The MTK of the specified subset.
}

\author{
Ryan C. Graebner
}

\examples{
data(genotypes)
MtkCalc(genotypes,subset=colnames(genotypes[1:5]))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
