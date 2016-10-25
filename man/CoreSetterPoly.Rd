\name{CoreSetterPoly}
\alias{CoreSetterPoly}
\title{
Genotype Subsetting for Autopolyploids and Polymorphic Markers
}
\description{
This function systematically eliminates genotypes from a large population to arrive at a favorable subset, and can accommodate datasets with autopolyploids and polymorphic markers. At this time, CoreSetterPoly can only use the Expected Heterozygosity criterion to quantify the value of subsets, and Sequential Backward Selection to arrive at favorable subsets.
}
\usage{
CoreSetterPoly(genos, ploidy, save = NULL)
}
\arguments{
  \item{genos}{
    A matrix of genotypes, where each row is one individual, and each set of X columns (where X is the ploidy) is one locus. At each locus, any number of alleles can be included, where each allele is referred to by a different integer. Missing data should be represented by NA. The X cells for any genotype-locus combination are the alleles known to be present at that locus for that genotype, in the frequency that they are present. If at least one but not all of the cells for a genotype-locus combination are listed as NA, that data point is imputed based on the other alleles at that locus.
}
  \item{ploidy}{
    The ploidy of the organizm to be subsetted, in respect to the number of alleles that can be present at one locus.
}
  \item{save}{
    A list of genotype names, corresponding to the column names in the genotype matrix, that will not be eliminated.
}
}
\value{
    Returns a matrix with three columns.  The first column is the rank of a particular genotype to the population's MTK, based on the order in which genotypes were eliminated (genotypes with lower rank were retained longer, genotypes with rank of 1 were not eliminated).  The second column is the name of the genotype. The third column shows the value of the subset including that genotype and all genotypes with a lower rank, as judged by the Expected Heterozygosity criterion.
}
\author{
Ryan C. Graebner and Alfonso Cuesta-Marcos
}
\examples{
data("genotypes")
CoreSetterPoly(genotypes.poly,ploidy=2,save=rownames(genotypes.poly)[c(1,5,9)])
}
\keyword{ misc }
