\name{SubsetterPIC}
\alias{SubsetterPIC}
\title{
Genotype Subsetting with PIC - Method One (Depreciated)
}
\description{
***This function has been superseeded by CoreSetter.
}
\usage{
SubsetterPIC(genos, save = NULL)
}
\arguments{
  \item{genos}{
    A matrix of genotypes, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, where 0 represents missing data.  All data in this matrix must be numeric.
}
  \item{save}{
    A list of genotype names, corresponding to the column names in the genotype matrix, that will not be eliminated.
}
}

\value{
Returns a matrix with three columns.  The first column is the importance of a particular genotype to the population's genetic diversity, based on the order in which genotypes were eliminated.  The second column is the name of the genotype, and the third column is the mean PIC of a population that includes the corresponding genotype, plus all genotypes that are more important.
}
\references{
Munoz-Amatrain M, Cuesta-Marcos A, Endelman JB, Comadran J, Bonman JM (2014) The USDA barley core collection: genetic diversity, population structure, and potential for genome-wide association studies. PloS One 9:e94688.
}
\author{
Ryan C. Graebner and Alfonso Cuesta-Marcos
}

\examples{
data("genotypes")
SubsetterPIC(genotypes,save=colnames(genotypes)[c(1,5,9)])
}
\keyword{ misc }
