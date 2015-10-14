\name{CoreSetter}
\alias{CoreSetter}
\title{
Genotype Subsetting
}
\description{
This function systematically eliminates genotypes from a large population to arrive at a favorable subset.  This method will typically return less favorable subsets than the method used by the CoreSetterCombined function if sufficient permutations are used for the later, but CoreSetter is quicker, and will rank all genotypes, as opposed to returning a single, static subset. Criteria that can be used to judge the value of subsets are Heterozygosity (HET; for rare-trait discovery; called PIC in earlier versions and in the paper describing this package), and the Mean of Transformed Kinships (MTK; for GWAS). A complete comparison of these two criteria is presented in Graebner et al. (2015).
}
\usage{
CoreSetter(genos, criterion = c("HET","MTK"), save = NULL,
    power = 10, mat = NULL)
}
\arguments{
  \item{genos}{
    A matrix of genotypes, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, or NA, where 0 represents a heterozygous marker, and NA represents missing data. Note that this coding is different from the earlier SubsetterPIC and SubsetterMTK, which cannot handle heterozygous markers.  All data in this matrix must be numeric.
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
  \item{mat}{
    A kinship matrix, if one has already been computed for the population.
}
}
\value{
    Returns a matrix with three columns.  The first column is the rank of a particular genotype to the population's MTK, based on the order in which genotypes were eliminated (genotypes with lower rank were retained longer, genotypes with rank of 1 were not eliminated).  The second column is the name of the genotype. The third column shows the value of the subset including that genotype and all genotypes with a lower rank, as judged by the specified criterion.
}
\references{
Graebner RC, Hayes PM, Hagerty CH, Cuesta-Marcos A (2015) A comparison of polymorphism information content and mean of transformed kinships as criteria for selection informative subsets of barley (Hordeum vulgare L. s. l) from the USDA Barley Core Collection. Genet Resour Crop Evol doi:10.1007/s10722-015-0265-z.
Munoz-Amatrain M, Cuesta-Marcos A, Endelman JB, Comadran J, Bonman JM (2014) The USDA barley core collection: genetic diversity, population structure, and potential for genome-wide association studies. PloS One 9:e94688.
}
\author{
Ryan C. Graebner and Alfonso Cuesta-Marcos
}
\note{
In Graebner et al. (2015), and in ther earlier funcitons SubsetterPIC and SubsetterMTK, heterozygous markers were counted as missing data, due to previous limitations if GeneticSubsetter. Please note that this has changed the required coding scheme for input genotype data. When using the HET criterion, this function uses the same method and criteria described in Munoz-Amatrain et al. (2014), but with a more computationally efficient approach.
}

\examples{
data("genotypes")
CoreSetter(genotypes,criterion="HET",save=colnames(genotypes)[c(1,5,9)])
CoreSetter(genotypes,criterion="MTK",save=colnames(genotypes)[c(1,5,9)])
}
\keyword{ misc }
