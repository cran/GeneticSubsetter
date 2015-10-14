\name{PicCalc}
\alias{PicCalc}
\title{
PIC Calculator (Depreciated)
}
\description{

***This function has been superseeded by PIC.
}
\usage{
PicCalc(data)
}
\arguments{
  \item{data}{
    A matrix of genotypes, where each column is one individual, each row is one marker, and marker values are 1, 0, or -1, where 0 represents missing data.  All data in this matrix must be numeric.
}
}
\value{
The mean PIC of all markers included for the given set of genotypes.
}
\author{
Ryan C. Graebner and Alfonso Cuesta-Marcos
}
\examples{
data("genotypes")
PicCalc(genotypes)
}
\keyword{ misc }
