\name{HS.solution}
\Rdversion{1.1}
\alias{HS.solution}
\title{ Hockey-Stick kinetics }
\description{
  Function describing two exponential decline functions with a break point
  between them.
}
\usage{
HS.solution(t, parent.0, k1, k2, tb)
}
\arguments{
  \item{t}{ Time. }
  \item{parent.0}{ Starting value for the response variable at time zero. }
  \item{k1}{ First kinetic constant. }
  \item{k2}{ Second kinetic constant. }
  \item{tb}{ Break point. Before this time, exponential decline according
    to \code{k1} is calculated, after this time, exponential decline proceeds
    according to \code{k2}. }
}
\value{
  The value of the response variable at time \code{t}.
}
\references{ 
  FOCUS (2006) \dQuote{Guidance Document on Estimating Persistence and
  Degradation Kinetics from Environmental Fate Studies on Pesticides in EU
  Registration} Report of the FOCUS Work Group on Degradation Kinetics,
  EC Document Reference Sanco/10058/2005 version 2.0, 434 pp,
  \url{http://focus.jrc.ec.europa.eu/dk} 
}
\examples{
  \dontrun{plot(function(x) HS.solution(x, 100, 2, 0.3, 0.5), 0, 2, ylim=c(0,100))}
}
\keyword{ manip }
