\name{UFA_locate_regex}
\alias{UFA_locate_regex}
\title{
UFA Locate regex
}
\description{
Locate indices of the pattern in the string
}
\usage{
UFA_locate_regex(string, pattern)
}
\arguments{
  \item{string}{
a string as character
}
  \item{pattern}{
a pattern to screen
}
}
\details{
This function returns `NA` when no matches is detected for the pattern.
}
\value{
A 2-column matrix of location indices. The first and second columns represent start postions and end positions, respectively.
}
\examples{
pattern <- "Cl"
string <- "NaCl.5HCl"
Location_Cl <- UFA_locate_regex(string, pattern)
}