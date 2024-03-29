% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snowfitfunctions.R
\name{fitsnowdepth}
\alias{fitsnowdepth}
\title{Derives snow melt coefficient from empirical data}
\usage{
fitsnowdepth(
  weather,
  snowdepth,
  precd,
  snowenv = "Taiga",
  STparams,
  plotout = T,
  snowem = 0.99,
  zm = 0.002,
  umin = 0.5,
  astc = 1.5
)
}
\arguments{
\item{weather}{a data.frame of weather variables (see details).}

\item{snowdepth}{a vector fo measured snow depths (cm). Can be interpolated.}

\item{precd}{a vector of daily precipitation (mm).}

\item{snowenv}{one of \code{Alpine}, \code{Maritime}, \code{Prairie}, \code{Taiga} or \code{Tundra} (see details)}

\item{STparams}{an object of class \code{SnowTparams} as returned by \code{\link[=fitsnowtemp]{fitsnowtemp()}}}

\item{plotout}{optional logical indicating whether to plot observed and predicted data}

\item{snowem}{optionally, numeric value of snow emissivity}

\item{zm}{optionally, numeric value of roughness length for momentum transfer of snow surface (m)}

\item{umin}{optionally, numeric value indicating minimum wind speed used in conductivity calculations (m/s)}

\item{astc}{optionally, numeric value indicating the temperature at which precipitation falls as snow (deg C)}
}
\value{
a list of the following:

\code{psnowdepth} predicted snow depth (cm)

\code{RMS} Root-mean-square error of predicted snow depths (cm)

\code{meltfact} a snow melt factor for use with \code{\link[=pSnow]{pSnow()}}
}
\description{
The function \code{fitsnowdepth} estimates the snow melt coefficient from empirical data (or
the outputs of a more complex model) and returns this along with
estimates of snow depth and the RMS error
}
\details{
The format and and units of \code{weather} must follow that in the example
dataset \code{climdata}. The entries in weather need to be hourly. If hourly
snow depth data are unavailable, they should be derived by interpolation.
The paramater \code{snowenv} is used to compute snow density
following Sturm et al (2010) J Hydrometeorology 11: 1380-1393.
}
\examples{
# Fit snow depth
require(NicheMapR)
# Derive estimates of snow temperature using NicheMapR in place of observed data
nmrout<-runNMRSnow(climdata, precd, 67.367, 26.629, ALTT = 199)
# Derive Snow temperature parameters
STparams<-fitsnowtemp(climdata,precd,nmrout$SNOWDEP,nmrout$SNOWTEMP)
# Fit snow depth model and plot result (red = Observed)
meltfact<-fitsnowdepth(climdata,snowdepth$snowdepth,precd,STparams=STparams)$meltfact
}
