% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{col2hex}
\alias{col2hex}
\title{Express colors as hexadecimal coding}
\usage{
col2hex(col)
}
\arguments{
\item{col}{Either numeric or character, specifying the color(s) to be
transformed.}
}
\value{
The color(s) used as input, but expressed as hexadecimal code(s).
}
\description{
Little function intended for internal use; will transform colors
(expressed either as numerical or characters) into hexadecimal code.
}
\examples{
plot(rnorm(n = 100), pch = 16, col = "red")
plot(rnorm(n = 100), pch = 16, col = col2hex("red"))

plot(rnorm(n = 100), pch = 16, col = 2)
plot(rnorm(n = 100), pch = 16, col = col2hex(2))
}
\keyword{internal}
