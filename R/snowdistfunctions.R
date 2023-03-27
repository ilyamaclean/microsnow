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
.epaif <- function(pai,hgt,snowdepth,clump) {
  # Adjust for snow dpeth
  nd<-length(snowdepth)/24
  sde<-matrix(snowdepth,ncol=24)
  sde<-apply(sde,1,mean)
  sde<-.vta(sde,hgt)
  h<-.rta(hgt,nd)
  mu<-(h-sde)/h
  mu[mu<0.001]<-0.001
  mu[is.na(mu)]<-1
  epai<-mu*pai
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
  if (class(pai)[1] == "SpatRaster") pai<-.is(pai)
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
.roughresample<-function(zm,dtm,xyf,tint=NA) {
  if (is.na(tint)) tint<-30*24
  mdm<-trunc(min(dim(dtm)[1:2])/2)
  if (xyf > mdm | is.na(xyf)) xyf<-NA
  if (is.na(xyf) == F) {
    if (xyf > 1) {
      n<-dim(zm)[3]
      nx<-trunc(n/tint)
      i<-c(0:(nx-1))*tint+1
      # aggregate and resample each month
      zma<-array(NA,dim=c(dim(zm)[1:2],nx))
      for (ii in 1:length(i)) {
        xx<-rast(zm[,,i[ii]])
        ext(xx)<-ext(dtm)
        crs(xx)<-crs(dtm)
        xx<-aggregate(xx,10,fun=mean)
        zma[,,ii]<-.is(resample(xx,dtm))
      }
      i<-floor((nx/n)*c(0:(n-1))+1)
      zma<-zma[,,i]
    } else zma<-zm
  } else {
    zmv<-apply(zm,3,mean,na.rm=T)
    zma<-.vta(zmv,dtm)
  }
  return(zma)
}
#' Calculates radiation absorbed by snow surface
# DK: added arguments `slr` and `apr` to allow for user input
.snowrad<-function(weather,snowdepth,dtm,hor,lat,long,epai,x,snowalb,snowem,slr=NA,apr=NA,clump=0,clumpd="daily") {
  # === (1a) Calculate solar altitude
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  jd<-.jday(tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  alt<-.solalt(lt,lat,long,jd,merid=0,dst=0)
  salt<-.vta(alt,dtm)*(pi/180)
  azi<-.solazi(lt,lat,long,jd,merid=0,dst=0)
  # === (1b) Calculate solar index
  si<-.solarindex(dtm,salt,azi,slr,apr)
  # === (1c) Calculate horizon angles
  i<-round(azi/15,0)+1; i[i==25]<-1
  hora<-hor[,,i]
  # === (1d) Calculate terrain shading
  shadowmask<-hora*0+1
  shadowmask[hora>tan(salt)]<-0
  si<-si*shadowmask
  # === (1e) Get two-stream parameters
  if (clumpd == "daily")  {
    cl<-.ehr(clump)
  } else cl<-clump
  twostreamp<-.twostreamparams(epai,x,cl,snowalb,snowalb,0,salt*(180/pi),si)
  # === (1f) Calculate sky view
  msl<-tan(apply(atan(hor),c(1,2),mean))
  svf<-0.5*cos(2*msl)+0.5
  sva<-.rta(rast(svf),length(jd))
  # === (1g) Calculate albedo
  Kc<-with(twostreamp,kkd$kd/kkd$k0)
  albd<-with(twostreamp,p1+p2+cl^2*snowalb)
  albb<-with(twostreamp,p5/sig+p6+p7+cl^Kc*gref2)
  # radiation
  dv<-weather$swrad-weather$difrad
  dirr<-.vta(dv,dtm)/si
  dirr[is.na(dirr)]<-0
  difr<-.vta(weather$difrad,dtm)
  albedo<-(dirr*sin(salt)*albb+difr*albd)/(dirr*si+difr)
  albedo[is.na(albedo)]<-albd[is.na(albedo)]
  albedo[albedo>1]<-albd[albedo>1]
  albedo[albedo>0.99]<-0.99
  albedo[albedo<0.01]<-0.01
  # === (1m) Calculate ground absorbed radiation
  # Shortwave
  Rbgm<-(1-cl^Kc)*exp(-twostreamp$kkd$kd*epai)+cl^Kc
  Rdbm<-with(twostreamp,(1-cl^2)*((p8/sig)*exp(-kkd$kd*epai)+p9*exp(-h*epai)+p10*exp(h*epai)))
  Rddm<-with(twostreamp,(1-cl^2)*(p3*exp(-h*epai)+p4*exp(h*epai))+cl^2)
  radGsw<-(1-snowalb)*(dirr*si*Rbgm+dirr*sin(salt)*Rdbm+difr*sva*Rddm)
  radGsw[is.na(radGsw)]<-0
  radGsw[radGsw<0]<-0
  # Longwave
  trd<-(1-cl^2)*exp(-epai)+cl^2
  tc<-.vta(weather$temp,dtm)
  skyem<-.vta(weather$skyem,dtm)
  Rem<-0.97*5.67*10^-8*(tc+273.15)^4 # Longwave emitted
  lwsky<-skyem*Rem # Longwave radiation down from sky
  radGlw<-0.97*(trd*sva*lwsky+(1-trd)*Rem+(1-sva)*Rem)
  # === (1m) Calculate canopy and ground combined absorbed radiation
  trb<-(1-cl^Kc)*exp(-twostreamp$kkd$kd*epai)+cl^Kc
  radCsw<-(1-albedo)*(dirr*sin(salt)+sva*difr)
  radCsw[is.na(radCsw)]<-0
  radCsw[is.infinite(radCsw)]<-0
  radClw<-sva*lwsky
  Rabsg<-radGsw+radGlw
  return(list(Rabsg=Rabsg,radCsw=radCsw,radClw=radClw,twostreamp=twostreamp))
}
.calcG<-function(u2,zu,snowdepth,dtm,hgta,pai,zmin,xyf) {
  snowd<-.vta(snowdepth,dtm)
  sel<-which(snowd>=hgta)
  # zero plane displacement
  d<-.zeroplanedis(hgta,pai)
  d[sel]<-snowd[sel]
  zm<-.roughlength(hgta,pai,d)
  zm[zm<zmin]<-zmin
  zm[sel]<-zmin
  # smooth roughness lengths
  da<-.roughresample(d,dtm,xyf)
  zma<-.roughresample(zm,dtm,xyf)
  uf<-(0.4*u2)/log((zu-da)/zma)
  # to canopy heat exchange surface
  gHes<-.gturb(uf,da,zma,zu)
  # To ground
  h<-hgta
  h[sel]<-snowd[sel]
  h[h<10*zmin]<-10*zmin
  gHa<-.gturb(uf,d,zm,zu,h,0,0.05) # From ref hgt to canopy top (or snow surface)
  gHs<-.gcanopy(uf,h,d,h,snowd) # From canopy top to snow surface
  g0<-1/(1/gHa+1/gHs) # From ground to ref hgt (above canopy)
  g0[sel]<-gHa[sel] # replace with gHa if snow depth > canopy height
  return(list(gHes=gHes,g0=g0))
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
    d<-apply(d,1,sum,na.rm=T)
    d
  }
  ad<-apply(a,c(1,2),.htdx)
  ad<-aperm(ad,c(2,3,1))
  ad
}
#' Computes snow melt when snow is an array
.snowmelt <- function(weather,precd,dtm,lat,long,pai,x,hgt,hgta,snowdepth,snowenv="Taiga",
                      meltfact=0.115,STparams,snowem=0.99,zu=2,zmin=0.002,umin=0.5,
                      astc=1.5,xyf=10,initdepth=0,dint=24,slr=NA,apr=NA,clump) {
  # Set snow depth to metres
  snowdepth<-snowdepth/100
  # Calculate snow refletcance
  prech<-rep(0,length(precd)*24)
  sel<-(c(1:length(precd))-1)*24+1
  prech[sel]<-precd
  prech<-ifelse(weather$temp>astc,0,prech)
  snowalb<-.snowalb(weather,prech,astc)
  snowalb<-.vta(snowalb,dtm)
  hor<-.hor(dtm)
  # =====#
  epai<-.epaif(pai,hgt,snowdepth,clump)
  epai<-.ehr(epai)
  Rabs<-.snowrad(weather,snowdepth,dtm,hor,lat,long,epai,x,snowalb,snowem,slr,apr,clump)
  Rabsg<-Rabs$Rabsg
  Rabsc<-Rabs$radCsw+Rabs$radClw
  # Calculate conductivities
  u2<-.winds(weather,dtm,hor)
  u2[u2<umin]<-umin
  pai<-.ehr(pai)
  gHs<-.calcG(u2,zu,snowdepth,dtm,hgta,pai,zmin,xyf)
  gHaC<-gHs$gHes
  gHaG<-gHs$g0
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
  return(list(meltcmG=meltcmG,meltcmC=meltcmC,snowtempG=snowtempG,snowtempC=snowtempC,Radparams=Rabs))
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
  uf<-suppressWarnings((0.4*uz)/log((zu-d)))
  ehgt<-pmax(hgt,snowd)
  uh<-suppressWarnings((uf/0.4)*log((ehgt-d)/zm))
  s<-which(uh<uf)
  uh[s]<-uf[s]
  Be<-0.205*pai^0.445+0.1
  a<-pai/hgt
  Lc<-(0.25*a)^-1
  Lm<-2*Be^3*Lc
  a<-(Be*hgt)/Lm
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
    dtmc<-aggregate(dtm,af,na.rm=T)
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
#' @param snowdepth a numeric vector of length equal to the number of timesteps (e.g. nrow(weather)) representing average snowdepth, in cm,
#' across the scene as e.g. returned by [runNMRSnow()] or spline interpolated from measurments.
#' @param dtm a SpatRast of elevations (m). the x and y dimensions of the raster must also be in metres
#' @param pai a single numeric value, SpatRast or array of plant area index values
#' @param hgt a SpatRast of vegetation heights
#' @param STparams snow temperature model coefficients as derived by [fitsnowtemp()]
#' @param slr an optional SpatRast of slope values (Radians). Calculated from dtm if not supplied, but outer cells will be NA.
#' @param apr an optional SpatRast of aspect values (Radians). Calculated from dtm if not supplied, but outer cells will be NA.
#' @param meltfact snow melt coefficient as returned by `fitsnowtemp` or `getmeltf`. Derived using `getmeltf` is not supplied.
#' @param plai a single numeric value, SpatRast or array of the proprotion of plant area index values that are leaves
#' @param x optional single numeric value, SpatRast, matrix or array of leaf distribution coefficients
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
#' @param initdepth single numeric value or matrix of initial snow depths at start of model run
#' @param out optional variable indicating whether to return hourly or daily snow depths (default is hourly)
#' @param clump a single numeric value or array of values between 0 and 1 indicating the fraction of radiation passing through larger gaps in the canopy,
#' ( see [microclimf::clumpestimate()])
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
#' # Run snow model with defaults and inbuilt datasets and model coefficients (takes 20 seconds)
#' snd<-nowdepth$snowdepth # snow depth vector
#' mout1<-modelsnowdepth(climdata,precd,snd,dtm,pai,hgt,STparams,SDparams$meltfact)
#' # Run snow model with defaults, but allowing spatially variable snow melt (takes ~100 seconds)
#' mout2<-modelsnowdepth(climdata,precd,snd,dtm,pai,hgt,STparams,SDparams$meltfact,spatialmelt = TRUE)
#' # Compare results
#' msnowdepth1<-apply(mout1$gsnowdepth,c(1,2),mean)
#' msnowdepth2<-apply(mout2$gsnowdepth,c(1,2),mean)
#' snowdif<-msnowdepth2-msnowdepth1
#' par(mfrow=c(2,2))
#' plot(rast(msnowdepth1))
#' plot(rast(msnowdepth2))
#' plot(rast(snowdif))
modelsnowdepth<-function(weather, precd, snowdepth, dtm, pai, hgt, STparams, meltfact=NA,
                         slr = NA, apr = NA, plai = 0.3, x = 1, lat = NA, long = NA,
                         snowenv = "Taiga", tpi_radius = 200, tfact = 10, snowem=0.99,
                         zmin=0.002, umin=0.5, astc=1.5, spatialmelt = FALSE, Sh = 6.3,
                         zu = 2, xyf = NA, initdepth = 0, out = "hourly", clump = 0.2) {
  # Unpack PackedSpatRasters if necessary
  if (class(dtm)[1]=="PackedSpatRaster") dtm<-rast(dtm)
  if (class(pai)[1]=="PackedSpatRaster") pai<-rast(pai)
  if (class(hgt)[1]=="PackedSpatRaster") hgt<-rast(hgt)
  #  Adjust wind speeed to 2m above tallest vegetation
  zo<-max(.is(hgt),na.rm=T)+2
  u2<-microctools::windadjust(weather$windspeed,zu,zo)
  weather$windspeed<-u2
  # Calculate time interval
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  dint<-24/(as.numeric(tme[2]-tme[1]))
  if (is.na(lat)) {
    ll<-.latlongfromrast(dtm)
    lat<-ll$lat
    long<-ll$long
  }
  ndays<-length(weather$temp)/24
  if (ndays%%1 != 0) stop ("weather data frame must contain hourly weather for complete days")
  pai<-.unpackpai(pai,ndays,dtm)
  plai<-.unpackpai(plai,ndays,dtm)
  x<-.unpackpai(x,ndays*24,dtm)
  clump<-.unpackpai(clump,ndays,dtm)
  hgta<-.rta(hgt,ndays*24)
  if (is.na(meltfact)) meltfact<-getmeltf(snowenv)
  # Calculate melt
  initdepths<-mean(initdepth,na.rm=T)
  snp<-pSnow(weather,precd,snowenv,STparams,meltfact,snowem,zmin,umin,astc,initdepths)
  # Check conditions
  reso<-res(dtm)[1]
  af<-round(tpi_radius/reso,0)
  me<-min(dim(dtm)[1:2]) # extent
  if (spatialmelt) {
    melt<-.snowmelt(weather,precd,dtm,lat,long,pai,x,hgt,hgta,snowdepth,snowenv,meltfact,
                    STparams,snowem,zo,zmin,umin,astc,xyf,initdepth,dint,
                    slr,apr,clump)
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
  zm<-.roughlength(.rta(hgt,ndays),pai,d)
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
                                       hgt,d[,,1],zm[,,1],phs[1]*1000,Sh,zo))/apply(fallswe,c(1,2),mean)
  rint[is.na(rint)]<-0
  twe<-.vta(dsnowdepth*(phs/10),dtm) # Total water equivelent (kg / m3)
  swe<-twe*(1-.rta(raster(rint),ndays)) # Ground water equivelent (kg / m3)
  canswe<-twe*.rta(raster(rint),ndays) # Canopy water equivelent (kg / m3)
  # Calculate steps when topographic positioning index needs recalculating
  sdepthchange<-dsnowdepth[2:length(dsnowdepth)]-dsnowdepth[1:(length(dsnowdepth)-1)]
  sdepthchange<-abs(c(0,sdepthchange))
  # Run model in daily timesteps
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
      rsnow<-rast(dsnow)
      ext(rsnow)<-ext(dtm)
      crs(rsnow)<-crs(dtm)
      dtms<-dtm+rsnow
      tpic<-.tpicalc(af,me,dtms,tfact)
    }
  }
  gsnowdepth<-swe/.vta(phs*10,dtm)
  if (out == "hourly") {
    # DK: presumed typo fix, previously saving to `snowdepth` but now `gsnowdepth`
    gsnowdepth<-.ehr(gsnowdepth)
    canswe<-.ehr(canswe)
  }
  if (spatialmelt) {
    snowtempG<-melt$snowtempG
    Radparams<-melt$Radparams
  } else {
    snowtempG<-.vta(snp$snowtemp,dtm)
    Radparams<-NA
  }
  snowout<-list(gsnowdepth=gsnowdepth,canswe=canswe,snowtempG=snowtempG,Radparams=Radparams)
  class(snowout)<-"snow"
  return(snowout)
}
