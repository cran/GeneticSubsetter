\name{PicCalc}
\alias{PicCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PIC Calculator
}
\description{
This function calculates the mean Polymorphism Information Content (PIC) of a set of genotypes.
}
\usage{
PicCalc(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    A matrix of genotypes, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, where 0 represents missing data.  All data in this matrix must be numeric.
}
}

\value{
The mean PIC of all markers included in "data" for a given set of genotypes.
}

\author{
Ryan C. Graebner and Alfonso Cuesta-Marcos
}

\examples{
data("genotypes")
PicCalc(genotypes)
}
\keyword{ misc }
