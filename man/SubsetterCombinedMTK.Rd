\name{SubsetterCombinedMTK}
\alias{SubsetterCombinedMTK}
\title{Genotype Subsetting with MTK- Method Two}
\description{
This function creates a series of random subsets.  Then, each of these subsets is improved using the SubsetOptimizerMTK function using a series of single-genotype replacements that result in a higher MTK for the subset, until no more single-genotype replacements can be made to improve the subset's mean PIC. Subsetting with MTK is likely more beneficial for Genome-Wide Association Scanning (GWAS), and subsetting with an alternative criteria, the Polymorphism Information Content (PIC), is likely more beneficial for rare-trait discovery.
}
\usage{
SubsetterCombinedMTK(genos, save = NULL, size = 100,
     power = 10, permutations = 100, print = TRUE)}
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
  \item{power}{
    The transformation that should be made to the kinship matrix. If power=1, the kinship matrix is not transformed, if power=2, the kinship matrix is squared, etc.  When the power is higher, this function preferentially eliminates genotypes that are closely related to other specific genotypes.
}
  \item{permutations}{
    The number of random subsets to improve.
}
  \item{print}{
    If TRUE, this function prints output showing the subset's PIC whenever a replacement is made.
}
}
\details{
  This method will typically return a subset with a higher MTK than the method used by the SubsetterMTK function if sufficient permutations are used. However, SubsetterCombinedMTK can take much more time, and will only return one set of unranked genotypes. If a sufficiently large number of permutations (i.e. 100-1000), different permutations will typically converge on one of several sets of genotypes for any given population and subset size.  However, there is no guarentee that this method will return the same result every time it is run.
}
\value{
  Returns a list of the genotype names included in the best subset found.
}

\author{Ryan C. Graebner}

\examples{
data("genotypes")
SubsetterCombinedMTK(genotypes,save=colnames(genotypes)[c(1,5,9)],
    size=6,permutations=20)
}
\keyword{ misc }
