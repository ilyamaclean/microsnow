% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microclimatemodel.R
\name{snowtemphumE}
\alias{snowtemphumE}
\title{Above ground snow microclimate model}
\usage{
snowtemphumE(
  microsnow,
  reqhgt = 0.05,
  pai_a = NA,
  xyf = 1,
  zf = NA,
  slr = NA,
  apr = NA,
  hor = NA,
  wsa = NA,
  gmn = 0.3
)
}
\arguments{
\item{microsnow}{an object of class microsnowin as returned by \code{\link[=snowmodelin]{snowmodelin()}}}

\item{reqhgt}{specified height (m) above or below ground for which microclimate
outputs are required. Negative below ground}

\item{pai_a}{total plant area index above \code{reqhgt}. Estimated if not supplied (see details).}

\item{xyf}{optional input for called function \code{\link[=snowwind]{snowwind()}}}

\item{zf}{optional input for called function \code{\link[=wind]{wind()}}}

\item{slr}{an optional SpatRaster object of slope values (Radians). Calculated from
dtm if not supplied.}

\item{apr}{an optional SpatRaster object of aspect values (Radians). Calculated from
dtm if not supplied.}

\item{hor}{an optional array of the tangent of the angle to the horizon in
24 directions. Calculated from dtm if not supplied.}

\item{wsa}{an optional array of wind shelter coefficients in 8 directions.
Calculated from dtm if not supplied.}

\item{gmn}{optional minimum convective conductance value (mol/m^2/s). See \code{\link[=snowwind]{snowwind()}})}
}
\value{
see \code{\link[=runmicro_hr]{runmicro_hr()}}
}
\description{
Used by e.g. \code{\link[=runmicro_hr]{runmicro_hr()}} to model microclimate above ground
}
\details{
DETAILS
}
