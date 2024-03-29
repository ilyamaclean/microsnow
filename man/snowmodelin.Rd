% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/microclimatemodel.R
\name{snowmodelin}
\alias{snowmodelin}
\title{Create object of class microsnowin}
\usage{
snowmodelin(
  weather,
  precd,
  snow,
  STparams,
  dtm,
  pai,
  hgt,
  x,
  clump,
  windhgt = 2
)
}
\arguments{
\item{weather}{a data.frame of weather variables (see details)}

\item{precd}{a vector of daily precipitation}

\item{snow}{an object of class \code{snow} as created by \code{\link[=modelsnowdepth]{modelsnowdepth()}}}

\item{STparams}{now temperature model coefficients as derived by \code{\link[=fitsnowtemp]{fitsnowtemp()}}}

\item{dtm}{PackedSpatRaster or SpatRast object of elevations (see details)}

\item{pai}{a single numeric value, matrix, SpatRast or array of plant area idex values.}

\item{hgt}{ackedSpatRaster or SpatRast of vegetation heights}

\item{x}{optional single numeric value, SpatRast, matrix or array of leaf angle coefficients (see details)}

\item{clump}{a single numeric value or array of values between 0 and 1 indicating the fraction of radiation passing through larger gaps in the canopy,}

\item{windhgt}{height above round of wind speed measurement (m)}
}
\value{
an object of class \code{microsnowin} used by \code{\link[=runmicrosnow_hr]{runmicrosnow_hr()}}
}
\description{
The function \code{snowmodelin} creates an object of class \code{microsnowin}
which unpacks various component inputs and reformats as required
for running the model in hourly timesteps. Here it is assumed that the input
weather data are a data.frame - i.e. not spatially variable.
}
\details{
The format and and units of \code{weather} must follow that in the example
dataset \code{climdata}. The x,y and z units of \code{dtm} must be all be in
metres and the coordinate reference system must be defined. Canopy heights
are assumed time-invariant over the period the model is run, and must
be supplied as a SpatRast with coordinates and extent matching \code{dtm}.
\code{pai}, \code{x} and \code{clump} can either be assumed spatially and temporally invariant\verb{(supplied as single numeric values), temporally invariant (supplied as matrices or SpatRasts with extent and coordinate reference system matching}dtm\verb{or spatially and temporally variant - x & y dims must match}dtm\verb{z dim must of length}weather$temp`
). The coefficient x represents the ratio of of vertical to horizontal projections of
leaf foliage (see Campbell 1986 Agric For Meteorol, 36: 317-21 for details). As
not all wind measurements are at reference height, the height of the wind speed measurement
must be specified if not 2 m. To enable calculation of below-canopy wind profiles
in tall canopy, the wind speed is adjusted to give values for a height 2 m above
the maximum vegetation height, using the wind-height profile for a reference grass
surface.
}
