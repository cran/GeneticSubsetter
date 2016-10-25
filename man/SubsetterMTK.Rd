\name{SubsetterMTK}
\alias{SubsetterMTK}
\title{Genotype Subsetting with PIC - Method One (Depreciated)}
\description{
***This function has been superseded by CoreSetter.
}
\usage{
SubsetterMTK(genos, save = NULL, power = 10, mat = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genos}{
    A matrix of genotypes, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, where 0 represents missing data.  All data in this matrix must be numeric.
}
  \item{save}{
    A list of genotype names, corresponding to the column names in the genotype matrix, that will not be eliminated.
}
  \item{power}{
    The transformation that should be made to the kinship matrix. If power=1, the kinship matrix is not transformed, if power=2, the kinship matrix is squared, etc.  When the power is higher, this function preferentially eliminates genotypes that are closely related to other specific genotypes.
}
  \item{mat}{
    A list of genotype names, corresponding to the column names in the genotype matrix, that will not be eliminated.
}
}
\value{
    Returns a matrix with four columns.  The first column is the importance of a particular genotype to the population's MTK, based on the order in which genotypes were eliminated.  The second column is the name of the genotype. The third and fourth columns are the MTK, and the mean of untransformed kinship values, respectively, of a population that includes the corresponding genotype, plus all genotypes that are more important.
}
\author{
Ryan C. Graebner
}
\examples{
data("genotypes")
SubsetterMTK(genotypes,save=colnames(genotypes)[c(1,5,9)])
}
\keyword{ misc }
