\name{CoreSetOptimizer}
\alias{CoreSetOptimizer}
\title{
Subset Optimization
}
\description{
This function works to systematically improves a subset via single-genotype replacements from a larger population.  This function will continue to work until no more single-genotype replacements can be made to increase the subset's value. Criteria that can be used to judge the value of subsets are Expected Heterozygosity (HET; for rare-trait discovery; called PIC in earlier versions and in the paper describing this package), and the Mean of Transformed Kinships (MTK; for GWAS). A complete comparison of these two criteria is presented in Graebner et al. (2015).
}
\usage{
CoreSetOptimizer(genos, subset, criterion = c("HET","MTK"), 
    mat = NULL, save = NULL, power = 10, print = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genos}{
    A matrix of genotypes, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, or NA, where 0 represents a heterozygous marker, and NA represents missing data. Note that this coding is different from the earlier SubsetOptimizerPIC and SubsetOptimizerMTK, which cannot handle heterozygous markers.  All data in this matrix must be numeric.
}
  \item{subset}{
    The names of the genotypes in the starting subset.
}
  \item{criterion}{
    The criterion to be used for comparing subsets (HET or MTK).
}
  \item{mat}{
    A kinship matrix, if one has already been computed for the population.
}
  \item{save}{
    A list of genotype names, corresponding to the column names in the genotype matrix, that will not be eliminated.
}
  \item{power}{
    The transformation that should be made to the kinship matrix, if the MTK criterion is used. If power=1, the kinship matrix is not transformed, if power=2, the kinship matrix is squared, etc.  When the power is higher, this function preferentially eliminates genotypes that are closely related to other genotypes in the population.
}
  \item{print}{
    Whether to the value of intermediate subsets.
}
}
\value{
  Returns a list of the genotype names included in the best subset found.
}
\references{
Graebner RC, Hayes PM, Hagerty CH, Cuesta-Marcos A (2015) A comparison of polymorphism information content and mean of transformed kinships as criteria for selection informative subsets of barley (Hordeum vulgare L. s. l) from the USDA Barley Core Collection. Genet Resour Crop Evol doi:10.1007/s10722-015-0265-z.
}
\author{
Ryan C. Graebner
}
\note{
The ability to recogize heterozygous markers was included in CoreSetOptimizer, resulting in a slightly different genotype coding scheme than the depreciated functions SubsetOptimizerPIC and SubsetOptimizerMTK.
}
\examples{
data("genotypes")
CoreSetOptimizer(genotypes,subset=colnames(genotypes)[c(1,3,5,7,8,9)],
    criterion="HET",save=colnames(genotypes)[c(1,5,9)])
CoreSetOptimizer(genotypes,subset=colnames(genotypes)[c(1,3,5,7,8,9)],
    criterion="MTK",save=colnames(genotypes)[c(1,5,9)])
}
\keyword{ misc }
