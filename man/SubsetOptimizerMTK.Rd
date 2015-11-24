\name{SubsetOptimizerMTK}
\alias{SubsetOptimizerMTK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subset Optimizer for MTK
}
\description{
This function works to systematically improve a subset's Mean of Transformed Kinships (MTK) via single-genotype replacements from a larger population.  This function will continue to work until no more single-genotype replacements can be made to increase the subset's mean MTK. Subsetting with MTK is likely more beneficial for Genome-Wide Association Scanning (GWAS), and subsetting with an alternative criteria, the Polymorphism Information Content (PIC), is likely more beneficial for rare-trait discovery.
}
\usage{
SubsetOptimizerMTK(genos, subset, mat = NULL, power = 10,
     save = NULL, print = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genos}{
    A matrix of genotypes, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, where 0 represents missing data.  All data in this matrix must be numeric.
}
  \item{subset}{
    The names of the genotypes in the starting subset.
}
  \item{mat}{
    A kinship matrix, if one has already been computed for the population.
}
  \item{power}{
    The transformation that should be made to the kinship matrix. If power=1, the kinship matrix is not transformed, if power=2, the kinship matrix is squared, etc.  When the power is higher, this function preferentially eliminates genotypes that are closely related to other specific genotypes.
}
  \item{save}{
    A list of genotype names, corresponding to the column names in the genotype matrix, that will not be eliminated.
}
  \item{print}{
    If TRUE, this function prints output showing the subset's PIC whenever a replacement is made.
}
}

\value{
  Returns a list of the genotype names included in the best subset found.
}

\author{
Ryan C. Graebner
}
\examples{
data("genotypes")
SubsetOptimizerMTK(genotypes,subset=colnames(genotypes)[c(1,3,5,7,8,9)],
    save=colnames(genotypes)[c(1,5,9)])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
