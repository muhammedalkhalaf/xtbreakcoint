\name{summary.xtbreakcoint}
\alias{summary.xtbreakcoint}
\title{Summary Method for xtbreakcoint Objects}
\description{
Summary method for objects of class \code{"xtbreakcoint"}, providing detailed
output including individual unit results.
}
\usage{
\method{summary}{xtbreakcoint}(object, ...)
}
\arguments{
\item{object}{An object of class \code{"xtbreakcoint"}.}

\item{...}{Additional arguments (ignored).}
}
\value{
Invisibly returns \code{object}.
}
\seealso{
\code{\link{xtbreakcoint}}, \code{\link{print.xtbreakcoint}}
}
\examples{
# Generate example panel data
set.seed(42)
N <- 5
T <- 30

panel_data <- data.frame(
  id = rep(1:N, each = T),
  time = rep(1:T, N),
  y = rnorm(N * T),
  x = rnorm(N * T)
)

result <- xtbreakcoint(y ~ x, data = panel_data, id = "id", time = "time",
                       max_factors = 0)
summary(result)
}
