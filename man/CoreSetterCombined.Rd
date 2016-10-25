\name{CoreSetterCombined}
\alias{CoreSetterCombined}
\title{
Genotype Subsetting-Alternative Method
}
\description{
This function creates a series of random subsets.  Then, each of these subsets is improved using the CoreSetOptimizer function using a series of single-genotype replacements that result in a higher value for the subset, until no more single-genotype replacements can be made to improve the subset. This process is similar to a Local Search. Criteria that can be used to judge the value of subsets are Expected Heterozygosity (HET; for rare-trait discovery; called PIC in earlier versions and in the paper describing this package), and the Mean of Transformed Kinships (MTK; for GWAS). A complete comparison of these two criteria is presented in Graebner et al. (2015).
}
\usage{
CoreSetterCombined(genos=NULL, size=NULL, criterion =  "HET",
    save = NULL, power = 10, permutations = 100, print = TRUE,
    mat = NULL)
}
\arguments{
  \item{genos}{
    A matrix of genotypes, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, or NA, where 0 represents a heterozygous marker, and NA represents missing data. Note that this coding is different from the earlier SubsetterCombinedPIC and SubsetterCombinedMTK, which cannot handle heterozygous markers.  All data in this matrix must be numeric.
}
  \item{size}{
    The desired subset size.
}
  \item{criterion}{
    The criterion to be used for comparing subsets (HET or MTK).
}
  \item{save}{
    A list of genotype names, corresponding to the column names in the genotype matrix, that will not be eliminated.
}
  \item{power}{
    The transformation that should be made to the kinship matrix, if the MTK criterion is used. If power=1, the kinship matrix is not transformed, if power=2, the kinship matrix is squared, etc.  When the power is higher, this function preferentially eliminates genotypes that are closely related to other genotypes in the population.
}
  \item{permutations}{
    The number of random subsets to improve.
}
  \item{print}{
    If TRUE, this function prints the Heterozygosity or MTK of the best subset identified.
}
  \item{mat}{
    A kinship matrix, if one has already been computed for the population. If a kinship matrix is included, the "genos" argument may be left empty.
}
}
\value{
  Returns a list of the genotype names included in the best subset found.
}
\references{
Graebner RC, Hayes PM, Hagerty CH, Cuesta-Marcos A (2016) A comparison of polymorphism information content and mean of transformed kinships as criteria for selection informative subsets of barley (Hordeum vulgare L. s. l) from the USDA Barley Core Collection. Genet Resour Crop Evol 63:477-482.
}
\author{
Ryan C Graebner
}
\note{
The ability to recognize heterozygous markers was included in CoreSetterCombined, resulting in a slightly different genotype coding scheme than the depreciated functions SubsetterCombinedPIC and SubsetterCombinedMTK.
}
\examples{
data("genotypes")
CoreSetterCombined(genotypes,size=6,criterion="HET",permutations=10,
    save=colnames(genotypes)[c(1,5,9)])
CoreSetterCombined(genotypes,size=6,criterion="MTK",permutations=10,
    save=colnames(genotypes)[c(1,5,9)])
}
\keyword{ misc }
