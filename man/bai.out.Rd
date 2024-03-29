\encoding{UTF-8}
\name{bai.out}
\alias{bai.out}
\title{ Basal Area Increment (Outside In) }
\description{
  Convert multiple ring-width series to basal area increment (i.e., ring
  area) going from the bark to the pith.
}
\usage{
bai.out(rwl, diam = NULL)
}
\arguments{
  \item{rwl}{ a \code{data.frame} with series as columns and years as
    rows such as that produced by \code{\link{read.rwl}} }
  \item{diam}{ an optional \code{data.frame} containing two variables.
    If present, variable one (\code{\var{series}} in the example below)
    gives the series \acronym{ID} as either \code{character}s or
    \code{factor}s.  These must exactly match
    \code{colnames(\var{rwl})}.  Variable two (\code{\var{diam}} in the
    example below) gives the diameter of the tree (in mm) up to the
    outermost ring (e.g., the diameter of the tree where the core was
    taken minus the thickness of the bark).  If \code{\var{diam}} is
    \code{NULL} then the diameter is taken as twice the sum of the
    widths for each series (column) in \code{\var{rwl}}. }
}
\details{
  This converts ring-width series (mm) to ring-area series (mm squared)
  (aka basal area increments) based on the diameter of the tree and the
  width of each ring moving towards the pith of the tree.  It is related
  to \code{\link{bai.in}}, which calculates each ring area starting from
  the inside of the tree and working outward.  Both methods assume a
  circular cross section (Biondi 1999).  See the references below for
  further details.
}
\value{
  A \code{data.frame} containing the ring areas for each series with
  column names, row names and dimensions of \code{\var{rwl}}.
}
\references{
  Biondi, F. (1999) Comparing tree-ring chronologies and repeated timber
  inventories as forest monitoring tools.  \emph{Ecological
    Applications}, \bold{9}(1), 216\enc{–}{--}227.

  Biondi, F. and Qeadan, F. (2008) A theory-driven approach to tree-ring
  standardization: Defining the biological trend from expected basal area
  increment.  \emph{Tree-Ring Research}, \bold{64}(2), 81\enc{–}{--}96.
}
\note{ DendroLab website: \url{https://dendrolaborg.wordpress.com/} }
\author{ Code by Andy Bunn based on work from DendroLab, University of
  Nevada Reno, \acronym{USA}.  Patched and improved by Mikko Korpela. }
\seealso{ \code{\link{bai.in}} }
\examples{
library(graphics)
library(utils)
\dontrun{
library(stats)
## Toy
n <- 100
## Make three fake tree-ring series to show that these funcs work on rwl objects
base.series <- 0.75 + exp(-0.2 * 1:n)
rwl <- data.frame(x1 = base.series + abs(rnorm(n, 0, 0.05)),
                  x2 = base.series + abs(rnorm(n, 0, 0.05)),
                  x3 = base.series + abs(rnorm(n, 0, 0.05)))

## The inside out method
foo <- bai.in(rwl = rwl)
## The outside in method
bar <- bai.out(rwl = rwl)

## Identical
head(bar)
head(foo)
}
## Use gp data
data(gp.rwl)
data(gp.dbh)
## dbh (minus the bark) from cm to mm 
gp.dbh2 <- gp.dbh[, 1:2]
gp.dbh2[, 2] <- (gp.dbh[, 2] - gp.dbh[, 3]) * 10
bar <- bai.out(rwl = gp.rwl, diam = gp.dbh2)
bar.crn <- chron(bar)
yrs <- time(bar.crn)
plot(yrs, bar.crn[, 1], type = "n",
     xlab = "Year", ylab = expression(mm^2))
lines(yrs, bar.crn[, 1], col = "grey", lty = "dashed")
lines(yrs, caps(bar.crn[, 1], nyrs = 32), col = "red", lwd = 2)
}
\keyword{ manip }
