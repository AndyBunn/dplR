\encoding{UTF-8}
\name{ffcsaps}
\alias{ffcsaps}
\title{ Smoothing Spline with User-Specified Rigidity and Frequency Cutoff }
\description{
  Applies a smoothing spline to \code{\var{y}} with rigidity determined
  by two parameters: frequency response \code{\var{f}} at a wavelength
  of \code{\var{nyrs}} years.
}
\usage{
ffcsaps(y, x = seq_along(y), nyrs = length(y)/2, f = 0.5)
}
\arguments{
  \item{y}{ a \code{numeric} vector, typically a tree-ring series. }
  \item{x}{ a \code{numeric} vector with the same length as
    \code{\var{y}}, typically giving the years of \code{\var{y}}. }
  \item{nyrs}{ a number greater than 1, affecting the rigidity of the
    spline.  When \code{\var{f}} is kept constant, a larger
    \code{\var{nyrs}} produces a more rigid spline.  Defaults to
    \code{length(\var{y})/2}. }
  \item{f}{ a number between 0 and 1 giving the frequency response at a
    wavelength of \code{\var{nyrs}} years.  When \code{\var{nyrs}} is
    kept constant, a smaller \code{\var{f}} produces a more rigid
    spline: At one extreme, \code{\var{f} = 0} causes the function to
    return the least-squares straight line fit to the data.  At the
    other extreme, \code{\var{f} = 1} results in the natural spline,
    i.e. the function outputs \code{\var{y}}.  The default value is
    0.5.}
}
\details{
  This applies a smoothing spline similar to the spline applied in most
  dendrochronological software.  See references for more information.
  
  In \code{\var{dplR}} in version 1.7.3 the function \code{\link{caps}} was introduced. This essentially replaces \code{\var{ffcsaps}} for most uses. Where \code{\var{ffcsaps}} is written entirely in R, \code{\link{caps}} is a wrapper for a Fortran subroutine from ARSTAN that is thousands of times faster.

}
\value{
  A filtered vector.
}
\references{
  Cook, E. R. and Kairiukstis, L. A., editors (1990) \emph{Methods of Dendrochronology: Applications in the Environmental Sciences}. Springer.  \acronym{ISBN-13}: 978-0-7923-0586-6.
    
  Cook, E. R. and Peters, K. (1981) The Smoothing Spline: A New Approach to Standardizing Forest Interior Tree-Ring Width Series for Dendroclimatic Studies. Tree-Ring Bulletin, 41, 45-53.
  
}
\note{ DendroLab website: \url{https://dendrolaborg.wordpress.com/} }
\author{
  Code provided by DendroLab based on programming by F. Qeadan and
  F. Biondi, University of Nevada Reno, \acronym{USA} and adapted for
  dplR by Andy Bunn.  Patched and improved by Mikko Korpela.
}

\seealso{ \code{\link{caps}}, \code{\link{detrend}} }

\examples{library(graphics)
## Use first series from the Mesa Verde data set
data(co021)
series <- co021[, 1]
series <- series[!is.na(series)]
plot(series, type = "l", ylab = "Ring Width (mm)", col = "grey")
lines(ffcsaps(series, nyrs = 10), col = "red", lwd = 2)
lines(ffcsaps(series, nyrs = 100), col = "green", lwd = 2)
## nyrs defaults to 0.5*length(series) = 347
lines(ffcsaps(series), col = "blue", lwd = 2)
legend("topright",
       c("Series", "nyrs=10", "nyrs=100",
         paste("Default nyrs (", length(series) / 2, ")", sep="")),
       fill=c("grey", "red", "green", "blue"))

## Compare to caps
all.equal(caps(series),ffcsaps(series))

## Note behavior when NA is encountered and
## take appropriate measures as demonstrated above
y <- c(NA,NA,rnorm(100))
ffcsaps(y)

}

\keyword{ smooth }
