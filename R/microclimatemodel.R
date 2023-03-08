.Tbelow<-function(weather,precd,reqhgt,dtm,lat,long,STparams,gsnowdepth,snowtempG,pai,
                  hgta,x,lw,da,zma,zo,snowem,umin,zmin,merid,dst,slr=NA,apr=NA) {
  n<-length(weather$temp)
  # Compute wind fraction velocity
  hor<-.hor2(dtm,2)
  u2<-.winds(weather,dtm,hor)
  u2[u2<umin]<-umin
  uf<-(0.4*u2)/log((zo-da)/zma)
  # Compute wind speed at top of canopy
  l_m<-.mixinglength(hgta,pai,zmin)
  tp<-0.2*pai*hgta
  a<-(tp/l_m)^0.5
  uh<-suppressWarnings((uf/0.4)*log((hgta-da)/zma))
  uh[is.na(uh)]<-uf[is.na(uh)]
  uh[uh<uf]<-uf[uh<uf]
  # Compute wind speed at reqhgt
  uz<-uh*exp(a*((reqhgt/hgta)-1))
  uz[is.na(uz)]<-uf[is.na(uz)]
  uz[uz>uh]<-uh[uz>uh]
  # Compute conductivities
  # ~~ Conductivity from ref height to canopy top
  hgta2<-hgta
  hes<-da+0.2*zma
  hgta2[hgta2<hes]<-hes[hgta2<hes]
  gRH<-.gturb(uf,da,zma,zo,hgta2,0,0.1)
  # ~~ Canopy top to reqhgt
  gHz<-.gcanopy(l_m,a,hgta,uh,hgta,reqhgt,0)
  # ~~ Ref height to reqhgt
  gA<-1/(1/gRH+1/gHz)
  # ~~ Reqhgt height to ground
  g0<-.gcanopy(l_m,a,hgta,uh,reqhgt,0,0)
  # ~~ Leaf conductivity
  gha<-1.41*0.135*sqrt(uz/(lw*1.5))
  gL<-1/(1/gha+1/(uz*43.3))
  # Compute snow temperature
  prech<-rep(precd,each=24)
  snowalb<-.snowalb(weather,prech,astc)
  snowalb<-.vta(snowalb,dtm)
  hor<-.hor(dtm)
  rhv<-hgta*0+reqhgt
  suf<-pmax(gsnowdepth,rhv)
  paia<-((hgta-suf)/hgta)*pai
  paia[paia<0]<-0
  paia[is.na(paia)]<-0
  # Calculate absorbed radiation
  Rabs<-.snowrad(weather,dtm,hor,lat,long,paia,x,snowalb,snowem,merid,dst,slr,apr)$Rabsc
  # Cap conductivity
  zpd<-6.5*zmin
  zh<-0.2*zmin
  uv<-weather$windspeed
  uv[uv<0.5]<-0.5
  gx<-(0.4^2*43*uv)/(log((2-zpd)/zmin)*log((2-zpd)/zh))
  gmx<-max(gx)
  gmn<-min(gx)
  gA1<-gA
  gA1[gA1>gmx]<-gmx
  gA1[gA1<gmn]<-gmn
  tc<-.vta(weather$temp,dtm)
  snowtempC<-SnowTemp(Rabs,gA1,tc,STparams,snowem)
  # Compute temperature
  Tr<-.vta(weather$temp,dtm)
  Tz<-(g0*snowtempG+gL*snowtempC+gA*Tr)/(g0+gL+gA)
  # Relative humidity
  sel<-which(weather$temp<0)
  easR<-0.61078*exp((17.27*weather$temp)/(weather$temp+237.3))
  easR[sel]<-0.61078*exp((21.875*weather$temp[sel])/(weather$temp[sel]+265.5))
  eaR<-easR*(weather$relhum/100)
  sel<-which(snowtempG<0)
  eaSG<-0.61078*exp((17.27*snowtempG)/(snowtempG+237.3))
  eaSG[sel]<-0.61078*exp((21.875*snowtempG[sel])/(snowtempG[sel]+265.5))
  sel<-which(snowtempC<0)
  eaSC<-0.61078*exp((17.27*snowtempC)/(snowtempC+237.3))
  eaSC[sel]<-0.61078*exp((21.875*snowtempC[sel])/(snowtempC[sel]+265.5))
  eaR<-.vta(eaR,dtm)
  ea<-(g0*eaSG+gL*eaSC+gA*eaR)/(g0+gL+gA)
  return(list(Tz=Tz,ea=ea,snowtempC=snowtempC))
}
.Tabove<-function(weather,snowtempG,dtm,reqhgt,da,zma,zu) {
  dTH<-snowtempG-.vta(weather$temp,dtm)
  lnr<-suppressWarnings(log((reqhgt-da)/zma)/log((zu-da)/zma))
  lnr[is.na(lnr)]<-0
  lnr[lnr<0]<-0
  lnr[lnr>1]<-1
  dTz<-dTH*(1-lnr)
  Tz<-.vta(weather$temp,dtm)+dTz
  # Relative humidity
  sel<-which(weather$temp<0)
  easR<-0.61078*exp((17.27*weather$temp)/(weather$temp+237.3))
  easR[sel]<-0.61078*exp((21.875*weather$temp[sel])/(weather$temp[sel]+265.5))
  eaR<-easR*(weather$relhum/100)
  sel<-which(snowtempG<0)
  eaS<-0.61078*exp((17.27*snowtempG)/(snowtempG+237.3))
  eaS[sel]<-0.61078*exp((21.875*snowtempG[sel])/(snowtempG[sel]+265.5))
  deH<-eaS-.vta(eaR,dtm)
  dez<-deH*(1-lnr)
  ea<-.vta(eaR,dtm)+dez
  # Convert to relhum
  return(list(Tz=Tz,ea=ea))
}
#' Runs spatial snow microclimate model
#'
#' @description `runmicrosnow` runs the spatial snow microclimate model to derive an
#' array of temperatures and relative humidities for each time increment in `weather`
#' and over each grid cell of `dtm`. Only valid for grid cells and time periods with snow.
#'
#' @param weather a data.frame of weather variables (see details).
#' @param precd a vector of daily precipitation (mm).
#' @param reqhgt height above ground for which to return temperatures and relative humidities (m)
#' @param dtm a raster of elevations (m). the x and y dimensions of the raster must also be in metres
#' @param slr an optional raster object of slope values (Radians). Calculated from dtm if not supplied, but outer cells will be NA.
#' @param apr an optional raster object of aspect values (Radians). Calculated from dtm if not supplied, but outer cells will be NA.
#' @param pai a single numeric value, raster or array of plant area index values
#' @param plai a single numeric value, raster or array of the proprotion of plant area index values that are leaves
#' @param hgt a raster of vegetation heights
#' @param x optional single numeric value, raster, matrix or array of leaf distribution coefficients
#' @param lw single numeric value or array of leaf widths (m)
#' @param lat latitude of location (decimal degrees). Derived from `dtm` is not supplied, so coordinate reference of system of `dtm` must be defined.
#' @param long longitude of location (decimal degrees). Derived from `dtm` is not supplied, so coordinate reference of system of `dtm` must be defined.
#' @param sda list of snow depths and temperatures as returned by [modelsnowdepth()]. Snow epth model run if not supplied.
#' @param snowenv one of `Alpine`, `Maritime`, `Prairie`, `Taiga` or `Tundra` (only required if snow depth model not run)
#' @param STparams snow temperature model coefficients as derived by [fitsnowtemp()] (only required if snow depth model not run)
#' @param meltfact snow melt coefficient as returned by `fitsnowtemp` or `getmeltf`. Derived using `getmeltf` is not supplied.(only required if snow depth model not run)
#' @param tpi_radius radius for applying topographic positioning index when distributing snow as returned by [tpiradius()] (only required if snow depth model not run)
#' @param tfact single numeric value specifying sensitivity of snow distribution to topographic positioning index (only required if snow depth model not run)
#' @param snowem optionally, numeric value of snow emissivity
#' @param zmin optionally, numeric value of roughness length for momentum transfer of snow surface without vegetation protruding  (m)
#' @param umin optionally, numeric value indicating minimum wind speed used in conductivity calculations (m/s)
#' @param astc optionally, numeric value indicating the temperature at which precipitation falls as snow (deg C) (only required if snow model not run)
#' @param spatialmelt optional logical indicating whether or not snow melt varies spatially (e.g. greater on sun-facing slopes). Model takes longer to run if `TRUE`
#' @param Sh optionally, branch snow load coefficient (kg/m^2). 6.6 for pine and 5.9 kg for spruce. (only required if snow model not run)
#' @param zu height above ground of wind speed measurement in `weather` (m)
#' @param xyf number of grid cells over which to smooth vertical wind height profile
#' @param merid longitude of local timezone meridian (0 for UTC)
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for British Summer Time if `merid` = 0).
#' @param initdepth single numeric value or matrix of initial snow depths at start of model run (only required if snow model not run)
#' @return a list of the following:
#' @return `Tz` predicted air temperature at height `reqhgt` (deg C). Equivelent to snow temperature if `reqhgt` < `snowdepth`
#' @return `tleaf` predicted canopy leaf temperature (deg C) (equivalent to canpy snow temperature)
#' @return `relhum` predicted relative humidity (percentage) at height `reqhgt`
#' @return `sda` outputs from running [modelsnowdepth()]
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. The paramater `snowenv` is used to compute snow density
#' following Sturm et al (2010) J Hydrometeorology 11: 1380-1393. The leaf distribution angle
#' coefficient is the ratio of vertical to horizontal projections of leaf foliage
#' (~1 for decidious woodland). MOdel only valid for pixels and time-periods with snow. Use package
#' `microclimf` to derive microclimate without snow.

#' @export
#' @examples
#' require(NicheMapR)
#' # Derive estimates of snow melt temperature parameters using NicheMapR
#' nmrout<-runNMRSnow(climdata,precd,66.215,29.293,ALTT = 360)
#' STparams<-fitsnowtemp(climdata,precd,nmrout$SNOWDEP,nmrout$SNOWTEMP)
#' meltfact<-fitsnowdepth(climdata,nmrout$SNOWDEP,precd,STparams=STparams)$meltfact
#' microsnowout<-runmicrosnow(climdata,precd,reqhgt,dtm,pai,plai=0.3,hgt,x=1,
#'                            STparams=STparams,spatialmelt = TRUE)
#' # Plot snow temperature at mid day on 1st Jan
#' plot(raster(microsnowout$Tz[,,13]))
runmicrosnow <- function(weather, precd, reqhgt, dtm, slr = NA, apr = NA, pai, plai = 0.3, hgt, x = 1, lw = 0.05, lat = NA, long = NA,
                         sda = NA, snowenv = "Taiga", STparams, meltfact = NA, tpi_radius = 200, tfact = 10,
                         snowem = 0.99, zmin = 0.002, umin = 0.5, astc = 1.5, spatialmelt = FALSE,
                         Sh = 6.3, zu = 2, xyf = NA, merid = 0, dst = 0, initdepth = 0) {
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  dint<-24/as.numeric(tme[2]-tme[1])
  if (is.na(lat)) {
    ll<-.latlongfromraster(dtm)
    lat<-ll$lat
    long<-ll$long
  }
  if (class(sda) == "logical") {
    cat("Running snow model\n")
    sda<-modelsnowdepth(weather,precd,dtm,pai,hgt,STparams,meltfact,plai,x,lat,long,
                        snowenv,tpi_radius,tfact,snowem,zmin,umin,astc,spatialmelt,Sh,
                        zu,xyf,merid,dst,initdepth)
  }
  # Compute roughness length and zero plane dispacement heights
  cat("Computing temperature and humidity above canopy\n")
  n<-length(weather$temp)
  if (length(dim(pai)) != 3) {
    pai<-.unpackpai(pai,n,dtm)
  }
  if (dim(pai)[3] != n)  pai<-.unpackpai(pai,n,dtm)
  if (length(dim(plai)) != 3) {
    plai<-.unpackpai(plai,n,dtm)
  }
  if (dim(plai)[3] != n) plai<-.unpackpai(plai,n,dtm)
  x<-.unpackpai(x,n,dtm)
  # Correct wind
  gsnowdepth<-sda$gsnowdepth/100
  zo<-max(.is(hgt),gsnowdepth,na.rm=T)+2
  u2<-microctools::windadjust(weather$windspeed,zu,zo)
  weather$windspeed<-u2
  hgta<-.rta(hgt,n)
  pai[pai==0]<-0.001
  da<-.zeroplanedis(hgta,pai)
  zma<-.roughlength(hgta,pai)
  da<-.roughresample(da,dtm,xyf,zmin*6.5,30*dint)
  zma<-.roughresample(zma,dtm,xyf,zmin,30*dint)
  sel<-which(gsnowdepth>da)
  da[sel]<-gsnowdepth[sel]
  zma[sel]<-zmin
  snowtempG<-sda$snowtempG
  # Compute temperatures
  snp<-sda$snowparams
  Tze<-.Tabove(weather,snowtempG,dtm,reqhgt,da,zma,zu)
  Tz<-Tze$Tz
  ea<-Tze$ea
  cat("Computing temperature and humidity below canopy\n")
  Tzbe<-.Tbelow(weather,precd,reqhgt,dtm,lat,long,STparams,gsnowdepth,snowtempG,pai,hgta,
                x,lw,da,zma,zo,snowem,umin,zmin,merid,dst,slr,apr)
  Tzb<-Tzbe$Tz
  eab<-Tzbe$ea
  selB<-which(hgta>reqhgt)
  Tz[selB]<-Tzb[selB]
  ea[selB]<-eab[selB]
  # Comvert to relative humidity
  sel<-which(Tz<0)
  es<-0.61078*exp((17.27*Tz)/(Tz+237.3))
  es[sel]<-0.61078*exp((21.875*Tz[sel])/(Tz[sel]+265.5))
  rh<-(ea/es)*100
  rh[rh>100]<-100
  selS<-which(gsnowdepth>=reqhgt)
  Tz[selS]<-snowtempG[selS]
  rh[selS]<-100
  return(list(Tz=Tz,tleaf=Tzbe$snowtempC,relhum=rh,sda=sda))
}
