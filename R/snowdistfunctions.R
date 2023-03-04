#' Compute horizon angle in 24 directions
.hor<-function(dtm) {
  hor<-array(NA,dim=c(dim(dtm)[1:2],24))
  for (i in 1:24) hor[,,i]<-.horizon(dtm,(i-1)*15)
  hor
}
#' Compute horizon angle in 24 directions above a given height
.hor2<-function(dtm,hgt) {
  hor<-array(NA,dim=c(dim(dtm)[1:2],24))
  for (i in 1:24) hor[,,i]<-.horizon2(dtm,(i-1)*15,hgt)
  hor
}
#' Calculate horizon angle above a given height (for wind sheltering) for all time steps
.horizon2 <- function(dtm, azimuth, hgt) {
  reso<-res(dtm)[1]
  dtm<-.is(dtm)
  dtm[is.na(dtm)]<-0
  dtm<-dtm/reso
  hgt<-hgt/reso
  azi<-azimuth*(pi/180)
  horizon<-array(0,dim(dtm))
  dtm3<-array(0,dim(dtm)+200)
  x<-dim(dtm)[1]
  y<-dim(dtm)[2]
  dtm3[101:(x+100),101:(y+100)]<-dtm
  for (step in 1:10) {
    horizon[1:x,1:y]<-pmax(horizon[1:x,1:y],(dtm3[(101-cos(azi)*step^2):(x+100-cos(azi)*step^2),
                                                  (101+sin(azi)*step^2):(y+100+sin(azi)*step^2)]-dtm3[101:(x+100),101:(y+100)])/(step^2))
    m<-horizon[1:x,1:y]
    sel<-which(m<(hgt/step^2))
    m[sel]<-0
    horizon[1:x,1:y]<-m
  }
  horizon
}
#' Create array of wind after applying shelter coefficient
.winds<-function(weather,dtm,hor) {
  i<-round(weather$winddir/15,0)+1; i[i==25]<-1
  hora<-hor[,,i]
  wc<-1-atan(0.17*100*hora)/1.65
  u<-.vta(weather$windspeed,dtm)*wc
  u
}
# Calculate daily pai above snow surface where snow depth is a vector
.epaif <- function(pai,hgt,snowdepth) {
  n<-length(snowdepth)/24
  sde<-matrix(snowdepth,ncol=24,byrow=T)
  sde<-apply(sde,1,mean)
  epai<-((.rta(hgt,n)-.vta(sde,hgt))/(.rta(hgt,n))*pai)
  epai[epai<0]<-0
  epai[is.na(epai)]<-0
  epai<-.ehr(epai)
  epai
}
# Calculate pai above snow surface where snow depth is an array
.epaif1<-function(pai,hgt,snowd) {
  epai<-((hgt-snowd)/hgt)*pai
  epai[epai<0]<-0
  epai[is.na(epai)]<-0
  epai
}
#' Converts pai to an array if provided a single value or raster
.unpackpai <- function(pai,n, dtm) {
  if (class(pai)[1] == "RasterLayer") pai<-.is(pai)
  if (length(as.vector(pai)) == 1) pai<-array(pai,dim=dim(dtm))
  if (class(pai)[1] == "matrix") {
    pai<-array(pai,dim=c(dim(pai),1))
  }
  nx<-dim(pai)[3]
  i<-floor((nx/n)*c(0:(n-1))+1)
  pai<-pai[,,i]
  pai
}
#' resamples roughness lengths (or zero-plane displacement) by xyf to smooth
.roughresample<-function(zm,dtm,xyf,mnval,tint=NA) {
  if (is.na(tint)) tint<-30
  zm[is.na(zm)]<-mnval
  mdm<-trunc(min(dim(dtm)[1:2])/2)
  if (is.na(xyf)) {
    xyf<-trunc(pmax(10/res(dtm)[1],10))
    xyf<-ifelse(xyf>mdm,mdm,xyf)
  }
  if (xyf > mdm) {
    warning(paste0("xyf too large. Setting to ",mdm))
    xyf<-mdm
  }
  n<-dim(zm)[3]
  nx<-trunc(n/tint) # number of months
  i<-c(0:(nx-1))*tint+1
  zm2<-zm[,,i]
  zma<-zm2
  zmr<-raster(zm[,,1],template=dtm)
  # aggregate and resample each month
  for (i in 1:nx) {
    xx<-raster(zm2[,,i],template=dtm)
    xx<-aggregate(xx,xyf)
    zma[,,i]<-.is(resample(xx,zmr))
  }
  zma[zma<mnval]<-mnval
  # subtitute monthly values back in to days
  i<-floor((nx/n)*c(0:(n-1))+1)
  zma<-zma[,,i]
  zma
}
#' Calculates radiation absorbed by snow surface
.snowrad<-function(weather,dtm,hor,lat,long,epai,x,snowalb,snowem,merid=0,dst=0,slr=NA,apr=NA) {
  # Compute solar index
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  jd<-.jday(tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  alt<-.solalt(lt,lat,long,jd,merid,dst)
  k<-.cank(x,.vta(alt*(pi/180),dtm))
  salt<-.vta(alt,dtm)*(pi/180)
  azi<-.solazi(lt,lat,long,jd,merid,dst)
  si<-.solarindex(dtm,salt,azi,slr,apr)
  # Compute dni
  sif<-cos((90-alt)*(pi/180))
  dni<-(weather$swrad-weather$difrad)/sif
  dni[dni>1352]<-1352
  dni[dni<0]<-0
  dni[is.na(dni)]<-0
  dni<-.vta(dni,dtm)
  # Compute skyview and terrain shading
  i<-round(azi/15,0)+1; i[i==25]<-1
  hora<-hor[,,i]
  # === Calculate terrain shading
  shadowmask<-hora*0+1
  shadowmask[hora>tan(salt)]<-0
  si<-si*shadowmask
  # === Calculate sky view
  msl<-tan(apply(atan(hor),c(1,2),mean))
  svf<-0.5*cos(2*msl)+0.5
  sva<-.rta(raster(svf),length(jd))
  # === Calculate direct canopy transmission
  trdir<-.cantransdir(epai,k,snowalb)
  trdif<-.cantransdif(epai,snowalb)
  trlw<-.cantransdif(epai,1-snowem)
  Rswg<-(trdir*si*dni+trdif*sva*.vta(weather$difrad,dtm))*(1-snowalb)
  xx<-k*dni
  xx[xx>1352]<-1352
  Rswc<-((1-trdir)*xx*si+(1-trdif)*sva*.vta(weather$difrad,dtm))*(1-snowalb)
  Rlwg<-microctools::canlw(.vta(weather$temp,dtm),epai,(1-snowem),.vta(weather$skyem,dtm))$lwabs
  Rlwc<-.vta(weather$skyem,dtm)*sva*snowem*5.67*10^-8*(.vta(weather$temp,dtm)+273.15)^4
  Rabsg<-Rswg+Rlwg
  Rabsc<-Rswc+Rlwc
  return(list(Rabsg=Rabsg,Rabsc=Rabsc))
}
#' Calculates thermal conductivity to ground snow surface
.gHaG<-function(u2,zu,snowdepth,dtm,hgta,pai,zmin,weather) {
  snowd<-.vta(snowdepth,dtm)
  zpd<-(snowd-0.2)*zmin
  zpd[zpd<6.5*zmin]<-6.5*zmin
  # winds
  uf<-(0.4*u2)#/log((zu-zpd)/zmin)
  # Compute wind speed at top of canopy
  l_m<-.mixinglength(hgta,pai,zmin)
  tp<-0.2*pai*hgta
  tp[tp<0.0001]<-0.0001
  a<-(tp/l_m)^0.5
  uh<-suppressWarnings((uf/0.4)*log((hgta-zpd)/zmin))
  uh[is.na(uh)]<-uf[is.na(uh)]
  uh[uh<uf]<-uf[uh<uf]
  # Conductivities
  # ~~ Ref height to canopy top
  h<-hgta
  sel<-which(hgta<(zpd+zmin*0.2+0.001))
  h[sel]<-zpd[sel]+zmin*0.2+0.001
  gRH<-.gturb(uf,zpd,zmin,zu,h,0,0)
  # ~~ Canopy top to ground
  gH0<-.gcanopy(l_m,a,hgta,uh,hgta,0,0)
  gHa<-1/(1/gRH+1/gH0)
  gHa[gHa==0]<-gRH[gHa==0]
  # Check limits
  zpd<-6.5*zmin
  zh<-0.2*zmin
  u2<-weather$windspeed
  u2[u2<0.5]<-0.5
  gx<-(0.4^2*43*u2)/(log((2-zpd)/zmin)*log((2-zpd)/zh))
  gmx<-max(gx)
  gmn<-min(gx)
  gHa[gHa>gmx]<-gmx
  gHa[gHa<gmn]<-gmn
  gHa
}
#' Calculates thermal conductivity to canopy snow surface
.gHaC<-function(u2,zu,snowdepth,dtm,hgta,pai,zmin,weather,dint,xyf) {
  snowd<-.vta(snowdepth,dtm)
  zpd<-(snowd-0.2)*zmin
  zpd[zpd<6.5*zmin]<-6.5*zmin
  # winds
  uf<-(0.4*u2)#/log((zu-zpd)/zmin)
  # Roughness lengths and conductivity
  da<-.zeroplanedis(hgta,pai)
  zma<-.roughlength(hgta,pai)
  da<-.roughresample(da,dtm,xyf,zmin*6.5,30*dint)
  zma<-.roughresample(zma,dtm,xyf,zmin,30*dint)
  zma[zma<zmin]<-zmin
  da[da<zpd]<-zpd[da<zpd]
  gHa<-.gturb(uf,da,zma,zu,NA,0,0)
  # Check limits
  zpd<-6.5*zmin
  zh<-0.2*zmin
  u2<-weather$windspeed
  u2[u2<0.5]<-0.5
  gx<-(0.4^2*43*u2)/(log((2-zpd)/zmin)*log((2-zpd)/zh))
  gmx<-max(gx)
  gmn<-min(gx)
  gHa[gHa>gmx]<-gmx
  gHa[gHa<gmn]<-gmn
  gHa
}
#' Sets snow temperature to zero when it freezes to account for latent heat of freezing
.freeze<-function(snowtemp) {
  snowcur<-snowtemp*0
  snowpre<-snowtemp*0
  phazefreeze<-snowcur+1
  snowcur[,,2:dim(snowcur)[3]]<-snowtemp[,,2:dim(snowtemp)[3]]
  snowpre[,,2:dim(snowcur)[3]]<-snowtemp[,,1:(dim(snowtemp)[3]-1)]
  sel<-which(snowpre>0 & snowcur < 0)
  phazefreeze[sel]<-0
  snowtemp<-snowtemp*phazefreeze
  snowtemp
}
#' Converts array of daily values to hourly values
.htd<-function(a) {
  .htdx<-function(x) {
    d<-matrix(x,ncol=24,byrow=T)
    d<-apply(d,1,sum)
    d
  }
  ad<-apply(a,c(1,2),.htdx)
  ad<-aperm(ad,c(2,3,1))
  ad
}
#' Computes snow melt when snow is an array
.snowmelt <- function(weather,precd,dtm,lat,long,pai,x,hgta,snowdepth,snowenv="Taiga",
                      meltfact=0.115,STparams,snowem=0.99,zu=2,zmin=0.002,umin=0.5,
                      astc=1.5,xyf=10,initdepth=0,dint=24,merid=0,dst=0,slr=NA,apr=NA) {
  # Set snow depth to metres
  snowdepth<-snowdepth/100
  # Calculate snow albedo
  prech<-rep(0,length(precd)*24)
  sel<-(c(1:length(precd))-1)*24+1
  prech[sel]<-precd
  prech<-ifelse(weather$temp>astc,0,prech)
  snowalb<-.snowalb(weather,prech,astc)
  snowalb<-.vta(snowalb,dtm)
  hor<-.hor(dtm)
  # =====#
  epai<-.epaif(pai,hgt,snowdepth)
  # Calculate absorbed radiation
  Rabs<-.snowrad(weather,dtm,hor,lat,long,epai,x,snowalb,snowem,merid,dst,slr,apr)
  Rabsg<-Rabs$Rabsg
  Rabsc<-Rabs$Rabsc
  # Calculate conductivities
  u2<-.winds(weather,dtm,hor)
  u2[u2<umin]<-umin
  pai<-.ehr(pai)
  gHaG<-.gHaG(u2,zu,snowdepth,dtm,hgta,pai,zmin,weather)
  gHaC<-.gHaC(u2,zu,snowdepth,dtm,hgta,pai,zmin,weather,dint,xyf)
  # Calculate snow temperature
  tc<-.vta(weather$temp,dtm)
  snowtempG<-SnowTemp(Rabsg,gHaG,tc,STparams,snowem)
  snowtempC<-SnowTemp(Rabsc,gHaC,tc,STparams,snowem)
  # Set temp when freezing to zero
  snowtempG<-.freeze(snowtempG)
  snowtempC<-.freeze(snowtempC)
  # Calculate snow melt
  densfun<-.snowdens(snowenv)
  msnowG<-snowtempG
  msnowC<-snowtempC
  msnowG[snowtempG<0]<-0
  msnowC[snowtempC<0]<-0
  meltcmG<-msnowG*meltfact
  meltcmC<-msnowC*meltfact
  # Calculate rain melt in cm
  tcd<-matrix(weather$temp,ncol=24,byrow=T)
  tcd<-apply(tcd,1,mean)
  rmeltd<-ifelse(tcd>0,tcd*precd*0.0125,0)
  rmeltcm<-.vta(rep(rmeltd/24,each=24)/10,dtm)
  # Calculate total melt
  meltcmG<-meltcmG+rmeltcm
  meltcmC<-meltcmC+rmeltcm
  meltcmG<-.htd(meltcmG)
  meltcmC<-.htd(meltcmC)
  return(list(meltcmG=meltcmG,meltcmC=meltcmC,snowtempG=snowtempG,snowtempC=snowtempC))
}
#' Calculates average daily wind speed so that canopy snow interception can be computed
.winddaily<-function(weather) {
  u2<-weather$windspeed
  u<-u2*cos(weather$winddir*(pi/180))
  v<-u2*sin(weather$winddir*(pi/180))
  u2<-matrix(u2,ncol=24,byrow=T)
  u<-matrix(u,ncol=24,byrow=T)
  v<-matrix(v,ncol=24,byrow=T)
  u2<-apply(u2,1,mean)
  u<-apply(u,1,mean)
  v<-apply(v,1,mean)
  wd<-(atan2(v,u)*(180/pi))%%360
  weather2<-data.frame(windspeed=u2,winddir=wd)
  weather2
}
#' Calculates wind speed within canopy given vertical profile within the canopy
.meancanopywind<-function(uz,dtm,pai,hgt,d,zm,snowd=0,zu=2) {
  # Apply shelter coefficient
  ln<-log((zu-d)/zm)
  ln[ln<0.4]<-0.4
  uf<-(0.4*uz)/ln
  ehgt<-pmax(hgt,snowd)
  ln<-suppressWarnings(log((ehgt-d)/zm))
  ln[ln<0.001]<-0.001
  ln[is.na(ln)]<-0.001
  uh<-suppressWarnings((uf/0.4)*log((ehgt-d)/zm))
  sel<-which(is.na(uh))
  uh[sel]<-uf[sel]
  l_m<-.mixinglength(hgt,pai)
  a<-((0.2*pai*hgt)/l_m)^0.5
  ef<-exp(a*((snowd/hgt)-1))
  ha<-hgt/a
  int<-ha*(1-ef)
  ur<-int/(hgt-snowd)
  ur[is.na(ur)]<-1
  ur[ur>1]<-1
  u<-uh*ur
  sel<-which(is.na(u))
  u[sel]<-uf[sel]
  sel<-which(u<uf)
  u[sel]<-uf[sel]
  u
}
#' Calculates appropriate radius over which to apply topographic position index
#' NB function needs tuning
tpiradius<-function(windspeed, hgt, PAI) {
  mhgt<-mean(.is(hgt),na.rm=T)
  mpai<-mean(.is(PAI),na.rm=T)
  zm<-microctools::roughlength(mhgt, mpai)
  d<-microctools::zeroplanedis(mhgt,mpai)
  u<-mean(windspeed)
  u2<-microctools::windadjust(u,2,mhgt+2)
  zr<-mhgt+2
  uf<-(0.4*u)/log((zr-d)/zm)
  # above canopy
  if (mhgt<=0.05) {
    uz<-(uf/0.4)*log((0.05-d)/zm)
  } else {
    uh<-(uf/0.4)*log((mhgt-d)/zm)
    a<-microctools::attencoef(mhgt,mpai)
    uz<-uh*exp(a*((0.05/mhgt)-1))
  }
  rad<-round(400*uz,0)
  rad
}
#' Calculates topographic positioning index
.tpicalc<-function(af,me,dtm,tfact) {
  # Calculate topographic positioning index
  if (af < me/2) {
    dtmc<-aggregate(dtm,af)
    dtmc<-resample(dtmc,dtm)
  } else dtmc<-dtm*0+mean(.is(dtm),na.rm=T)
  tpi<-dtmc-dtm
  tpic<-exp(.is(tpi)/tfact)
  tpic[tpic<0.05]<-0.1
  tpic[tpic>10]<-10
  tpic
}
#' Calculates canopy snow interception
#'
#' @description Function `canopysnowint` calculates snow intercepted by the canopy
#' in a single time step using the model of Hedstrom and Pomeroy (1998) Hydrological
#' Processes, 12: 1611-1625.
#'
#' @param prec daily precipitation (mm)
#' @param uz wind speed at height `zu`  (m/s)
#' @param dtm a raster of elevations
#' @param gsnowd ground snow depth (cm)
#' @param Lt snow load in previous time step (mm SWE)
#' @param pai plant area index
#' @param plai proportion of `pai` that is living vegetation (i.e. leaves)
#' @param hgt canopy height (m)
#' @param d optionally, zero plane displacement height of wind profile
#' @param zm optionally, roughness length for momentum transfer
#' @param phs optionally, fresh snow density (kg / m^3)
#' @param Sh optionally, branch snow load coefficient (kg/m^2). 6.6 for pine and 5.9 kg for spruce.
#' @param zu height above ground of wind speed measurement `uz` (m)
#' @return `Int` snow water equivalent (mm) intercepted by canopy
#' @import microclimf
canopysnowint<-function(precd,uz,dtm,gsnowd,Lt,pai,plai,hgt,d=NA,zm=NA,phs=375,Sh=6.3,zu=2) {
  hgt<-.is(hgt)
  if (is.na(d)[1]) d<-0.65*hgt
  if (is.na(zm)[1]) zm<-d/6.5
  uc<-.meancanopywind(uz,dtm,pai,hgt,d,zm,gsnowd,zu)
  sde<-pai*0+gsnowd
  epai<-.epaif1(pai,hgt,sde)
  J<-1000 # mean forested canopy downwind distance. Set high
  w<-0.8 # terminal velocity of snow flake (m/s)
  LAI<-plai*epai # Leaf area index
  Cc<-1-exp(-LAI) # Snow leaf contact
  S<-Sh*(0.26+46/phs) # Branch snow load
  Lstr<-S*LAI  # Max canopy load (mm SWE)
  ehgt<-hgt-gsnowd
  ehgt[ehgt<0]<-0
  btm<-1-Cc*((uc*hgt)/(w*J))
  Cp<-Cc/btm  # maximum plan area of the snow±leaf contact per unit area of ground
  k<-Cp/Lstr
  k[is.na(k)]<-1
  Iinit<-(Lstr-Lt)*(1-exp(-k*precd))
  Int<-Iinit*0.678
  sel<-which(Int>precd)
  Int[sel]<-precd[sel]
  Int
}
#' Runs spatial snow model
#'
#' @description `modelsnowdepth` runs the spatial snow model to derive an
#' array of snow depths and temperatures for each time increment in `weather`
#' and over each grid cell of `dtm`
#'
#' @param weather a data.frame of weather variables (see details).
#' @param precd a vector of daily precipitation (mm).
#' @param dtm a raster of elevations (m). the x and y dimensions of the raster must also be in metres
#' @param slr an optional raster object of slope values (Radians). Calculated from dtm if not supplied, but outer cells will be NA.
#' @param apr an optional raster object of aspect values (Radians). Calculated from dtm if not supplied, but outer cells will be NA.
#' @param snowdepth a numeric vector of length equal to the number of timesteps (e.g. nrow(weather)) representing average snowdepth, in cm, across the scene
#' @param pai a single numeric value, raster or array of plant area index values
#' @param hgt a raster of vegetation heights
#' @param STparams snow temperature model coefficients as derived by [fitsnowtemp()]
#' @param meltfact snow melt coefficient as returned by `fitsnowtemp` or `getmeltf`. Derived using `getmeltf` is not supplied.
#' @param plai a single numeric value, raster or array of the proprotion of plant area index values that are leaves
#' @param x optional single numeric value, raster, matrix or array of leaf distribution coefficients
#' @param lat latitude of location (decimal degrees). Derived from `dtm` is not supplied, so coordinate reference of system of `dtm` must be defined.
#' @param long longitude of location (decimal degrees). Derived from `dtm` is not supplied, so coordinate reference of system of `dtm` must be defined.
#' @param snowenv one of `Alpine`, `Maritime`, `Prairie`, `Taiga` or `Tundra` (see details)
#' @param tpi_radius radius for applying topographic positioning index when distributing snow as returned by [tpiradius()]
#' @param tfact single numeric value specifying sensitivity of snow distribution to topographic positioning index
#' @param snowem optionally, numeric value of snow emissivity
#' @param zmin optionally, numeric value of roughness length for momentum transfer of snow surface without vegetation protruding  (m)
#' @param umin optionally, numeric value indicating minimum wind speed used in conductivity calculations (m/s)
#' @param astc optionally, numeric value indicating the temperature at which precipitation falls as snow (deg C)
#' @param spatialmelt optional logical indicating whether or not snow melt varies spatially (e.g. greater on sun-facing slopes). Model takes longer to run if `TRUE`
#' @param Sh optionally, branch snow load coefficient (kg/m^2). 6.6 for pine and 5.9 kg for spruce.
#' @param zu height above ground of wind speed measurement in `weather` (m)
#' @param xyf number of grid cells over which to smooth vertical wind height profile
#' @param merid longitude of local timezone meridian (0 for UTC)
#' @param dst an optional numeric value representing the time difference from the timezone meridian (hours, e.g. +1 for British Summer Time if `merid` = 0).
#' @param initdepth single numeric value or matrix of initial snow depths at start of model run
#' @param out optinal variable indicating whether to return hourly or daily snow depths (default is hourly)
#' @return a list of the following:
#' @return `gsnowdepth` array of predicted ground snow depth (cm) across the scene
#' @return `canswe` array of predicted canopy snow water equivalent (mm / m^2) across the scene
#' @return `snowtempG` array of predicted ground snow surface temperature (deg C) across the scene
#' @details #' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. The paramater `snowenv` is used to compute snow density
#' following Sturm et al (2010) J Hydrometeorology 11: 1380-1393. The leaf distribution angle
#' coefficient is the ratio of vertical to horizontal projections of leaf foliage (~1 for decidious woodland).
#' @export
#' @examples
#' require(NicheMapR)
#' # Derive estimates of snow melt temperature parameters using NicheMapR
#' nmrout<-runNMRSnow(climdata,precd,66.215,29.293,ALTT = 360)
#' STparams<-fitsnowtemp(climdata,precd,nmrout$SNOWDEP,nmrout$SNOWTEMP)
#' meltfact<-fitsnowdepth(climdata,nmrout$SNOWDEP,precd,STparams=STparams)$meltfact
#' # Run snow model with defaults (takes 20 seconds)
#' mout1<-modelsnowdepth(climdata,precd,dtm,pai,hgt,STparams,meltfact)
#' # Run snow model with defaults, but allowing spatially variable snow melt (takes ~100 seconds)
#' mout2<-modelsnowdepth(climdata,precd,dtm,pai,hgt,STparams,meltfact,spatialmelt = TRUE)
#' # Compare results
#' msnowdepth1<-apply(mout1$gsnowdepth,c(1,2),mean)
#' msnowdepth2<-apply(mout2$gsnowdepth,c(1,2),mean)
#' snowdif<-msnowdepth2-msnowdepth1
#' par(mfrow=c(2,2))
#' plot(raster(msnowdepth1))
#' plot(raster(msnowdepth2))
#' plot(raster(snowdif))
modelsnowdepth<-function(weather, precd, snowdepth, dtm, slr = NA, apr = NA,
                         pai, hgt, STparams, meltfact=NA, plai = 0.3,
                         x = 1, lat = NA, long = NA, snowenv = "Taiga", tpi_radius = 200, tfact = 10,
                         snowem=0.99, zmin=0.002, umin=0.5, astc=1.5, spatialmelt = FALSE,
                         Sh = 6.3, zu = 2, xyf = NA, merid = 0, dst = 0, initdepth = 0, out = "hourly") {
  zo<-max(.is(hgt),na.rm=T)+2
  u2<-microctools::windadjust(weather$windspeed,zu,zo)
  weather$windspeed<-u2
  # Calculate time interval
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  dint<-24/(as.numeric(tme[2]-tme[1]))
  if (is.na(lat)) {
    ll<-.latlongfromraster(dtm)
    lat<-ll$lat
    long<-ll$long
  }
  ndays<-length(weather$temp)/24
  if (ndays%%1 != 0) stop ("weather data frame must contain hourly weather for complete days")
  pai<-.unpackpai(pai,ndays,dtm)
  plai<-.unpackpai(plai,ndays,dtm)
  x<-.unpackpai(x,ndays*24,dtm)
  hgta<-.rta(hgt,length(weather$temp))
  if (is.na(meltfact)) meltfact<-getmeltf(snowenv)
  # Calculate melt
  initdepths<-mean(initdepth,na.rm=T)
  snp<-pSnow(weather,precd,snowenv,STparams,meltfact,snowem,zmin,umin,astc,initdepths)
  # Check conditions
  reso<-res(dtm)[1]
  af<-round(tpi_radius/reso,0)
  me<-min(dim(dtm)[1:2]) # extent
  if (spatialmelt) {
    melt<-.snowmelt(weather,precd,dtm,lat,long,pai,x,hgta,snowdepth,snowenv,meltfact,
                    STparams,snowem,zu,zmin,umin,astc,xyf,initdepth,dint,merid,dst,
                    slr,apr)
  } else {
    mh<-snp$snowmelt
    md<-matrix(mh,ncol=24,byrow=TRUE)
    md<-apply(md,1,sum)
    melt1<-.vta(md,dtm)
    melt<-list(meltcmG=melt1,meltcmC=melt1)
  }
  # Convert snowfall to daily and snow water equivelent (kg /m3)
  snowfall<-matrix(snp$snowfall,ncol=24,byrow=TRUE)
  snowfall<-apply(snowfall,1,sum)
  densfun<-.snowdens(snowenv)
  snowfall<-(snowfall*densfun[2]*10)
  n<-length(snowfall)
  fallswe<-.vta(snowfall,dtm)
  # Calculate topographic positioning index
  reso<-res(dtm)[1]
  af<-round(tpi_radius/reso,0)
  me<-min(dim(dtm)[1:2]) # extent
  # Calculate mean daily terrain-corrected wind speed
  hor<-.hor2(dtm,2)
  we2<-.winddaily(weather)
  u2<-.winds(we2,dtm,hor)
  # Calculate roughness lengths
  d<-.zeroplanedis(.rta(hgt,ndays),pai)
  zm<-.roughlength(.rta(hgt,ndays),pai)
  d<-.roughresample(d,dtm,xyf,zmin*6.5)
  zm<-.roughresample(zm,dtm,xyf,zmin)
  # Calculate topographic positioning index
  tpic<-.tpicalc(af,me,dtm,tfact)
  # Calculate snow water equivelent of snow melt
  meltsweG<-(melt$meltcmG*densfun[2]*10) # Ground snow melt (kg / m3)
  meltsweC<-(melt$meltcmC*densfun[2]*10) # Canopy snow melt(kg / m3)
  snowswe<-fallswe*0
  # Calculate daily snow density and snow water equivelent
  phs<-matrix(snp$snowdens,ncol=24,byrow=T)
  phs<-apply(phs,1,mean)
  dsnowdepth<-matrix(snp$snowdepth,ncol=24,byrow=T)
  dsnowdepth<-apply(dsnowdepth,1,mean)
  # Calculate initial ratio of snow depths
  rint<-suppressWarnings(canopysnowint(apply(fallswe,c(1,2),mean),
                                       mean(weather$windspeed),dtm,0,0,pai[,,1],plai[,,1],
                                       hgt,d[,,1],zm[,,1],phs[1]*1000,Sh,zu))/
    apply(fallswe,c(1,2),mean)
  rint[is.na(rint)]<-0
  twe<-.vta(dsnowdepth*(phs/10),dtm) # Total water equivelent (kg / m3)
  swe<-twe*(1-.rta(raster(rint),ndays)) # Ground water equivelent (kg / m3)
  canswe<-twe*.rta(raster(rint),ndays) # Canopy water equivelent (kg / m3)
  # Calculate steps when topographic positioning index needs recalculating
  sdepthchange<-dsnowdepth[2:length(dsnowdepth)]-dsnowdepth[1:(length(dsnowdepth)-1)]
  sdepthchange<-abs(c(0,sdepthchange))
  # Run model in daily timesteps
  ndays<-length(snowfall)
  for (i in 2:ndays) {
    # Calculate ground snow depth (needed for effective pai)
    gsnowd<-swe[,,i-1]/(phs[i-1]*10)
    # Snow intecerption by canopy (kg / m3)
    int<-suppressWarnings(canopysnowint(fallswe[,,i],u2[i],dtm,gsnowd,canswe[,,i-1],
                                        pai[,,i],plai[,,i],hgt,d[,,i],zm[,,i],phs[i]*1000,Sh))
    # Calculate spatially distributed fall less that intercepted by canopy
    gfall<-fallswe[,,i]-int
    gfall<-gfall*tpic
    ad<-mean(snowfall[i]-int)/mean(gfall,na.rm=T)
    ad[is.na(ad)]<-1
    ad[is.infinite(ad)]<-1
    gfall<-gfall*ad
    # Calculate canopy snow water equivelent
    m<-canswe[,,i-1]+int-meltsweC[,,i]
    m[m<0]<-0
    canswe[,,i]<-m
    # Calculate ground snow water equivelent
    m<-swe[,,i-1]+gfall-meltsweG[,,i]
    m[m<0]<-0
    swe[,,i]<-m
    if (sdepthchange[i]>0) {
      dsnow<-swe[,,i]/(phs[i]*10)
      dtms<-dtm+raster(dsnow,template=dtm)
      tpic<-.tpicalc(af,me,dtms,tfact)
    }
  }
  gsnowdepth<-swe/.vta(phs*10,dtm)
  if (out == "hourly") {
    gsnowdepth<-.ehr(gsnowdepth)
    canswe<-.ehr(canswe)
  }
  if (spatialmelt) {
    snowtempG<-melt$snowtempG
  } else snowtempG<-.vta(snp$snowtemp,dtm)
  return(list(gsnowdepth=gsnowdepth,canswe=canswe,snowtempG=snowtempG))
}
