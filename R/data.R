#' A data frame of hourly weather
#'
#' A data frame of hourly weather (2nd Oct 2017 to 1st Apr 2019) for Sodankylä in Finland (67.367N, 26.629E)
#'
#' @format a data frame with the following elements:
#' \describe{
#'  \item{obs_time}{POSIXlt object of dates and times in UTC}
#'  \item{temp}{temperature (degrees C)}
#'  \item{relhum}{relative humidity (percentage)}
#'  \item{pres}{atmospheric press (kPa)}
#'  \item{swrad}{Total incoming shortwave radiation (W / m^2)}
#'  \item{difrad}{Diffuse radiation (W / m^2)}
#'  \item{skyem}{Sky emissivity (0-1)}
#'  \item{windspeed}{Wind speed (m/s)}
#'  \item{winddir}{Wind direction (decimal degrees)}
#' }
#' @source \url{https://en.ilmatieteenlaitos.fi/}
"climdata"
#'
#' A 25 m resolution digital terrain dataset.
#'
#' A dataset of elevations (m) for Sodankylä in Finland (xmin = 2963700, xmax = 2964950,
#' ymin = 10221300, ymax = 10222550, CRS: 3395).
#'
#' @format A RasterLayer object with 50 rows and 50 columns
#' @source \url{https://aws.amazon.com/public-datasets/terrain/}
"dtm"
#'
#' A 25 m resolution dataset of canopy height.
#'
#' A dataset of canopy heights (m) for Sodankylä in Finland (xmin = 2963700, xmax = 2964950,
#' ymin = 10221300, ymax = 10222550, CRS: 3395).
#'
#' @format A RasterLayer object with 50 rows and 50 columns
#' @source \url{https://land.copernicus.eu/}
"hgt"
#'
#' A 25 m resolution dataset of plant area index values.
#'
#' A dataset of plant area index values for Sodankylä in Finland (xmin = 2963700, xmax = 2964950,
#' ymin = 10221300, ymax = 10222550, CRS: 3395).
#'
#' @format A RasterLayer object with 50 rows and 50 columns
#' @source \url{https://land.copernicus.eu/}
"pai"
#'
#' Daily precipitation
#'
#' A vector of daily precipitation (2nd Oct 2017 to 1st Apr 2019) for Sodankylä in Finland (67.367N, 26.629E)
#'
#' @format A vector of daily precipitation (mm/day, snow water equivalent when snow)
#' @source \url{https://en.ilmatieteenlaitos.fi/}
"precd"
#'
#' A data frame of hourly snow depths
#'
#' A data frame of hourly snow depths (2nd Oct 2017 to 1st Apr 2019) for Sodankylä in Finland (67.367N, 26.629E)
#'
#' @format a data frame with the following elements:
#' \describe{
#'  \item{obs_time}{POSIXlt object of dates and times in UTC}
#'  \item{snowdepth}{snow dpeth (cm)}
#' }
#' @source \url{https://en.ilmatieteenlaitos.fi/}
"snowdepth"
