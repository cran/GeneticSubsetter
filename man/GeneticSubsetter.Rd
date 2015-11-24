\name{GeneticSubsetter}
\alias{GeneticSubsetter}
\docType{package}
\title{Genetic Subsetter}
\description{
This package contains a set of tools that can be used to select a subset from a larger population, using genetic data. Two criteria are used to identify subsets, in seperate functions: Polymorphism Information Content (PIC) and the Mean of Transformed Kinships (MTK).}

\details{
When selecting subsets of genotypes, two factors are important to consider: the criteria by which to judge subsets, and the method used to identfy the set of genotypes that best fit that criteria. Two criteria are Polymorphism Information Content (PIC) and the Mean of Transformed Kinships (MTK). Tests suggest that of these two criteria, PIC is better if the resulting subset will be used for rare-trait discovery, while MTK is better if the resulting subset will be used for genome-wide association scanning. To reach subsets with a high PIC or MTK, SubsetterPIC and SubsetterMTK systematially remove genotypes from the full set, creating a full ranking of genotype's contributions to their respective criteria. The function SubsetterPIC uses the same method and criteria described in Munoz-Amatrain et al. (2014), but with a more computationally efficient approach. Alternatively, SubsetterCombinedPIC and SubsetterCombinedMTK work to systematically improve a user-defined number of random subsets via single-genotype replacements, until no replacement can increase their respective criteria.  These later methods generally return subsets with higher PIC or MTK, but are subset-size specific, take more time to compute, and will not always return identical results.

}
\references{
Munoz-Amatrain M, Cuesta-Marcos A, Endelman JB, Comadran J, Bonman JM (2014) The USDA barley core collection: genetic diversity, population structure, and potential for genome-wide association studies. PloS One 9:e94688.
}
\author{
Ryan C. Graebner <ryan.graebner@gmail.com> and Alfonso Cuesta-Marcos
}

\keyword{ package }
\examples{
data("genotypes")
SubsetterPIC(genotypes,save=colnames(genotypes)[c(1,5,9)])
SubsetterMTK(genotypes,save=colnames(genotypes)[c(1,5,9)])
}
