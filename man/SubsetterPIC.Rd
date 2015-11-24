\name{SubsetterPIC}
\alias{SubsetterPIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genotype Subsetting with PIC - Method One
}
\description{
This function systematically eliminates genotypes from a large population on the basis of which genotype's removal would result in the population with the highest mean Polymorphism Information Content (PIC).  This method will typically return a subset with a lower mean PIC than the method used by the SubsetterCombinedPIC function if sufficient permutations are used for the later, but is quicker and will rank all genotypes, as opposed to returning a single, static subset. Subsetting with PIC is likely more beneficial for rare-trait discovery, and subsetting with an alternative criteria, the Mean of Transformed Kinships (MTK), is likely more beneficial for Genome-Wide Association Scanning (GWAS). This function uses the same method and criteria described in Munoz-Amatrain et al. (2014), but with a more computationally efficient approach.
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
