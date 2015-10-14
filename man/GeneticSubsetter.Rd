\name{GeneticSubsetter}
\alias{GeneticSubsetter}
\docType{package}
\title{Genetic Subsetter}
\description{
This package contains a set of tools that can be used to select a subset from a larger population, using genetic data. Two criteria are used to identify subsets, in seperate functions: Heterozygosity (HET; called PIC in earlier versions and in the paper describing this package) and the Mean of Transformed Kinships (MTK).}

\details{
When selecting subsets of genotypes, two factors are important to consider: the criteria by which to judge subsets, and the method used to identfy the set of genotypes that best fit that criteria. Two criteria are Heterozygosity (HET) and the Mean of Transformed Kinships (MTK). Tests suggest that of these two criteria, Heterozygosity is better if the resulting subset will be used for rare-trait discovery, while MTK is better if the resulting subset will be used for genome-wide association scanning (Graebner et al. 2015). To reach subsets with a high Heterozygosity or MTK, CoreSetter systematially removes genotypes from the full set, creating a full ranking of genotype's contributions to their respective criteria. When the HET criterion is selected, CoreSetter uses the same method and criteria described in Munoz-Amatrain et al. (2014), except CoreSetter uses a more computationally efficient approach, and CoreSetter can consider heterozygous markers. Alternatively, CoreSetterCombined works to systematically improve a user-defined number of random subsets via single-genotype replacements, until no replacement can increase the selected criteria.  This later method generally returns subsets with a higher Heterozyosity or MTK, but are subset-size specific, take more time to compute, and will not always return identical results in subsequent runs.
}
\references{
Graebner RC, Hayes PM, Hagerty CH, Cuesta-Marcos A (2015) A comparison of polymorphism information content and mean of transformed kinships as criteria for selection informative subsets of barley (Hordeum vulgare L. s. l) from the USDA Barley Core Collection. Genet Resour Crop Evol doi:10.1007/s10722-015-0265-z.
Munoz-Amatrain M, Cuesta-Marcos A, Endelman JB, Comadran J, Bonman JM (2014) The USDA barley core collection: genetic diversity, population structure, and potential for genome-wide association studies. PloS One 9:e94688.
}
\author{
Ryan C. Graebner <ryan.graebner@gmail.com> and Alfonso Cuesta-Marcos
}

\keyword{ package }
\examples{
data("genotypes")
CoreSetter(genotypes,criterion="HET",save=colnames(genotypes)[c(1,5,9)])
CoreSetterCombined(genotypes,size=5,criterion="MTK")
}
