\name{UFAx_workflow}
\alias{UFAx_workflow}
\title{
UFAx Workflow
}
\description{
This function runs the exhaustive chemical enumeration part of the IDSL.UFAx pipeline.
}
\usage{
UFAx_workflow(spreadsheet)
}
\arguments{
  \item{spreadsheet}{
IDSL.UFAx parameter spreadsheet
}
}
\value{
The MolecularFormulaAnnotationTable is saved in the assigned folder in the parameter spreadsheet.
}
\note{
You should load the IDSL.UFA package to run the IDSL.UFAx functions.
}
\examples{
\donttest{
library(IDSL.UFA) # You should load the IDSL.UFA package to run the IDSL.UFAx functions.
library(IDSL.UFAx)
s_path <- system.file("extdata", package = "IDSL.UFAx")
SSh1 <- paste0(s_path,"/UFAx_parameters.xlsx")
temp_wd <- tempdir() # update this address
temp_wd_zip <- paste0(temp_wd,"/003.mzML_UFAx_testfiles.zip")
spreadsheet <- readxl::read_xlsx(SSh1)
download.file(
  "https://github.com/idslme/IDSL.UFA/blob/main/003.mzML_UFAx_testfiles.zip?raw=true",
  destfile = temp_wd_zip)
unzip(temp_wd_zip, exdir = temp_wd)
spreadsheet[1, 4] <- temp_wd
spreadsheet[4, 4] <- temp_wd
spreadsheet[6, 4] <- temp_wd
spreadsheet[3, 4] <- "seq(1, 100)" # peak IDs to process
UFAx_workflow(spreadsheet)
}
}