\name{kinreport}
\Rdversion{1.1}
\alias{kinreport}
\title{
Creates a report of the kinetic fits
}
\description{
Function to create a report for a set of fitted models, passing it to the 
console as well as to a file, if specified.
}
\usage{
kinreport(kinobject, file = NA, data = TRUE, R2 = FALSE, vcov = FALSE, endpoint.digits = 1)
}
\arguments{
  \item{kinobject}{
	A list containing the following elements:
	The name of the parent compound to be output (\code{parent}),
	the type of the test system (\code{type}), 
	the name of the specific test system used for generating this dataset
        (\code{system}), 
  the input dataset (\code{data}),
	the list of fitted kinetic models (\code{fits}), as returned by
        \code{\link{kinfit}}, and the list of results (\code{results}) 
        as returned by \code{\link{kinresults}}.
	Optionally also
	the label position of the test compound (\code{label}) and
	the source of the data (\code{source}). }
  \item{file}{ The name of the file to which to write. }
  \item{data}{ Boolean specifying if the input data should be printed? }
  \item{R2}{ Boolean specifying if the coefficient of determination R2 should
    be printed. R2 is calculated from all values (not mean values) as 1 - RSS/TSS where
    RSS is the sum of squares of residuals, and TSS ist the total sum of squares, i.e. 
    the sum of squares of the difference of each observation from the overall mean.  }
  \item{vcov}{ Should the variance-covariance matrix/matrices be reported? }
  \item{endpoint.digits}{ How many digits should be reported for DT50 and DT90 values? }
}
\value{
The function is called for its side effect, namely the report being passed
to the R console as well as to a text file if a filename is specified.
}
\author{ Johannes Ranke }
\examples{
data(FOCUS_2006_A)
kinfits <- kinfit(FOCUS_2006_A)
kinobject <- list(
	parent = "Compound XY",
	type = "Degradation in the environment",
	system = "System 1",	
	data = FOCUS_2006_A,
	source = "Synthetic example data from FOCUS kinetics",
	fits = kinfits,
	results = kinresults(kinfits))	
kinreport(kinobject)
}
\keyword{ manip }
