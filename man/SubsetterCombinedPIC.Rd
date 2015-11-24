\name{SubsetterCombinedPIC}
\alias{SubsetterCombinedPIC}
\title{
Genotype Subsetting with PIC- Method Two
}
\description{
This function creates a series of random subsets.  Then, each of these subsets is improved using the SubsetOptimizerPIC function using a series of single-genotype replacements that result in a higher mean PIC for the subset, until no more single-genotype replacements can be made to improve the subset's mean PIC. Subsetting with PIC is likely more beneficial for rare-trait discovery, and subsetting with an alternative criteria, the Mean of Transformed Kinships (MTK), is likely more beneficial for Genome-Wide Association Scanning (GWAS).
}
\usage{
SubsetterCombinedPIC(genos, save = NULL, size = 100,
     permutations = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genos}{
    A matrix of genotypes, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, where 0 represents missing data.  All data in this matrix must be numeric.
}
  \item{save}{
    A list of genotype names, corresponding to the column names in the genotype matrix, that wil not be eliminated.
}
  \item{size}{
    The size of the subsets to be developed.
}
  \item{permutations}{
    The number of random subsets to improve.
}
}
\details{
  This method will typically return a subset with a higher mean PIC than the method used by the SubsetterPIC function if sufficient permutations are used. However, SubsetterCombinedPIC can take much more time, and will only return one set of unranked genotypes. If a sufficiently large number of permutations (i.e. 100-1000), different permutations will typically converge on one of several sets of genotypes for any given population and subset size.  However, there is no guarentee that this method will return the same result every time it is run.
}
\value{
  Returns a list of the genotype names included in the best subset found.
}

\author{
Ryan C. Graebner and Alfonso Cuesta-Marcos
}

\examples{
data("genotypes")
SubsetterCombinedPIC(genotypes,save=colnames(genotypes)[c(1,5,9)],
    size=6,permutations=20)
}
\keyword{ misc }
