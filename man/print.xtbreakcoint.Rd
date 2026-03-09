\name{print.xtbreakcoint}
\alias{print.xtbreakcoint}
\title{Print Method for xtbreakcoint Objects}
\description{
Print method for objects of class \code{"xtbreakcoint"}.
}
\usage{
\method{print}{xtbreakcoint}(x, ...)
}
\arguments{
\item{x}{An object of class \code{"xtbreakcoint"}.}

\item{...}{Additional arguments (ignored).}
}
\value{
Invisibly returns \code{x}.
}
\seealso{
\code{\link{xtbreakcoint}}, \code{\link{summary.xtbreakcoint}}
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
print(result)
}
