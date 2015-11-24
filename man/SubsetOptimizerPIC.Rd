\name{SubsetOptimizerPIC}
\alias{SubsetOptimizerPIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Subset Optimizer with PIC
}
\description{
This function works to systematically improve a subset's mean Polymorphism Information Content (PIC) via single-genotype replacements from a larger population.  This function will continue to work until no more single-genotype replacements can be made to increase the subset's mean PIC. Subsetting with PIC is likely more beneficial for rare-trait discovery, and subsetting with an alternative criteria, the Mean of Transformed Kinships (MTK), is likely more beneficial for Genome-Wide Association Scanning (GWAS).
}
\usage{
SubsetOptimizerPIC(genos, subset, save = NULL, print = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genos}{
    A matrix of genotypes, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, where 0 represents missing data.  All data in this matrix must be numeric.
}
  \item{subset}{
    The names of the genotypes in the starting subset.
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
Ryan C. Graebner and Alfonso Cuesta-Marcos
}

\examples{
data("genotypes")
SubsetOptimizerPIC(genotypes,subset=colnames(genotypes)[c(1,3,5,7,8,9)],
    save=colnames(genotypes)[c(1,5,9)])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ misc }
