#' Computes snow reflectance
.snowalb <- function(weather, prech, astc=1.5) {
  snowyn<-ifelse(prech>0,1,0)
  age<-0
  for (i in 2:length(snowyn)) {
    age[i]<-ifelse(snowyn[i]==1,0,age[i-1]+1)
  }
  aged<-matrix(age,ncol=24,byrow=T)
  aged<-apply(aged,1,mean)
  aged<-aged/24
  # Calculate snow albedo
  snowalb<-(-9.8740*log(aged)+78.3434)/100
  snowalb[snowalb>0.95]<-0.95
  snowalb<-rep(snowalb,each=24)
}
#' Computes snow energy balance
.SnowEnergyBalance<-function(tc,u2,rh,pk,Rsw,skyem,snowalb,snowem=0.99,zm=0.002,umin=0.5) {
  # Variables set in function
  zpd=6.5*zm
  zh<-0.2*zm
  sb<-5.67*10^-8
  ph<-43.3
  cp<-29.3
  lambdaS<-51135.14
  lambdaW<-44526
  lambda<-tc*0+lambdaW
  sel<-which(tc<0)
  lambda[sel]<-lambdaS
  gHa<-(0.4^2*ph*u2)/(log((2-zpd)/zm)*log((2-zpd)/zh))
  # Radiation emitted
  Rem<-snowem*sb*(tc+273.15)^4
  # Radiation absorbed
  Rlwd<-skyem*sb*(tc+273.15)^4
  Rabs<-(1-snowalb)*Rsw+snowem*Rlwd
  # Latent heat
  ests<-0.61078*exp((17.27*tc)/(tc+237.3))
  ests[sel]<-0.61078*exp((21.875*tc[sel])/(tc[sel]+265.5))
  ea<-(rh/100)*ests
  L<-((lambda*gHa)/pk)*(ests-ea)
  # Snow heat storage
  return(list(Rabs=Rabs,Rem=Rem,L=L))
}
#' Fits snow model (one step)
.FitSnowOne<-function(weather,snowdepth,precd,snowenv,STparams,meltfact=0.11,
                      snowem=0.99,zm=0.002,umin=0.5,astc=1.5) {
  predsnowd<-pSnow(weather,precd,snowenv,STparams,meltfact,snowem,zm,umin,astc)$snowdepth
  RMS<-sqrt(sum((predsnowd-snowdepth)^2)/length(snowdepth))
  RMS
}
#' Estimates snow temperature
#'
#' The function `SnowTemp` estimates snow temperatures
#'
#' @param Rabs radiation absorbed by snow surface (W/m^2)
#' @param gHa molar conductivity for heat (mol/m^2/s)
#' @param tc air temperature (deg G)
#' @param STparams an object of class `SnowTparams` as returned by [fitsnowtemp()]
#' @param snowem optionally, numeric value(s) of snow emissivity (~0.95)
#' @return snow temperature (deg C)
#' @export
SnowTemp<-function(Rabs,gHa,tc,STparams,snowem=0.99) {
  m1<-STparams$m1
  m2<-STparams$m2
  sb<-5.67*10^-8
  Rem<-sb*snowem*(tc+273.15)^4
  Rnet<-Rabs-Rem
  pdT<-m1$coef[1,1]+m1$coef[2,1]*Rnet+m1$coef[3,1]*gHa+m1$coef[4,1]*Rnet*gHa
  psnowtemp<-pdT+tc
  pdT2<-m2$coef[1,1]+m2$coef[2,1]*Rnet+m2$coef[3,1]*gHa+m2$coef[4,1]*Rnet*gHa
  psnowtemp2<-pdT2+tc
  sel<-which(psnowtemp>0)
  psnowtemp[sel]<-psnowtemp2[sel]
  psnowtemp
}
#' Estimates snow depth and other snow variables
#'
#' The function `pSnow` estimates snow depth, albedo, density and temperature
#'
#' @param weather a data.frame of weather variables (see details).
#' @param precd a vector of daily precipitation (mm).
#' @param snowenv one of `Alpine`, `Maritime`, `Prairie`, `Taiga` or `Tundra` (see details)
#' @param STparams an object of class `SnowTparams` as returned by [fitsnowtemp()]
#' @param meltfact a snow melt coefficient, as returned by [fitsnowdepth()] or [getmeltf()]
#' @param snowem optionally, snow emissivity
#' @param zm optionally, roughness length for momentum transfer of snow surface (m)
#' @param umin optionally, numeric value indicating minimum wind speed used in conductivity calculations (m/s)
#' @param astc optionally, numeric value indicating the temperature at which precipitation falls as snow (deg C)
#' @param initdepth single numeric value giving the initial snow depth (cm)
#' @return a list of the following:
#' @return `snowdepth` snow depth (cm)
#' @return `snowdens` snow density (g/cm^3)
#' @return `snowtemp` snow temperature (deg C)
#' @return `snowref` snow reflectance (range 0 to 1)
#' @return `snowfall` snow fall (cm) at each time increment
#' @return `snowmelt` snow melt (cm) at each time increment
#' @export
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. The paramater `snowenv` is used to compute snow density
#' following Sturm et al (2010) J Hydrometeorology 11: 1380-1393.
pSnow <- function(weather, precd, snowenv = "Taiga", STparams, meltfact=0.11,
                  snowem=0.99,zm=0.002,umin=0.5,astc=1.5, initdepth=0) {
  # Calculate snow albedo
  prech<-rep(0,length(precd)*24)
  sel<-(c(1:length(precd))-1)*24+1
  prech[sel]<-precd
  prech<-ifelse(weather$temp>astc,0,prech)
  snowalb<-.snowalb(weather,prech,astc)
  # Estimate snow temperature
  sb<-5.67*10^-8
  Rabs<-(1-snowalb)*weather$swrad+weather$skyem*sb*snowem*(weather$temp+273.15)^4
  u2<-weather$windspeed
  u2[u2<umin]<-umin
  zpd<-6.5*zm
  zh<-0.2*zm
  gHa<-(0.4^2*43*u2)/(log((2-zpd)/zm)*log((2-zpd)/zh))
  psnowtemp<-SnowTemp(Rabs,gHa,weather$temp,STparams,snowem)
  # Set to zero when freezing
  snowcur<-psnowtemp[2:length(psnowtemp)]
  snowprev<-psnowtemp[1:(length(psnowtemp)-1)]
  phasefreeze<-c(1,ifelse(snowprev>0 & snowcur < 0,0,1))
  psnowtemp<-psnowtemp*phasefreeze
  # Calculate snow melt
  densfun<-.snowdens(snowenv)
  msnow<-ifelse(psnowtemp<0,0,psnowtemp)
  meltcm<-msnow*meltfact
  # Calculate rain melt in cm
  tcd<-matrix(weather$temp,ncol=24,byrow=T)
  tcd<-apply(tcd,1,mean)
  rmeltd<-ifelse(tcd>0,tcd*precd*0.0125,0)
  rmeltcm<-rep(rmeltd/24,each=24)/10
  # Calculate snow fall (cm)
  fallcm<-(prech*0.997)/(densfun[2]*10)
  # Calculate depth
  predsnowd<-initdepth
  for (i in 2:(length(fallcm)+1)) {
    predsnowd[i]<-predsnowd[i-1]+fallcm[i-1]-meltcm[i-1]-rmeltcm[i-1]
    predsnowd[i]<-ifelse(predsnowd[i]<0,0,predsnowd[i])
  }
  predsnowd<-predsnowd[-1]
  # calculate snow pack age
  age<-0
  for (i in 1:length(predsnowd)) {
    age[i+1]<-ifelse(predsnowd[i]>0,age[i]+1/24,0)
  }
  age<-age[-1]
  # Calculate snow density
  pdensity<-(densfun[1]-densfun[2])*(1-exp(-densfun[3]*predsnowd-densfun[4]*age))+densfun[2]
  predsnowd<-predsnowd*(densfun[2]/pdensity)
  return(list(snowdepth=predsnowd,snowdens=pdensity,snowtemp=psnowtemp,
              snowref=snowalb,snowfall=fallcm,snowmelt=meltcm+rmeltcm))
}
#' The function `fitsnowtemp` derives model parameters for fitting the snow temperature model
#'
#' @param weather a data.frame of weather variables (see details).
#' @param precd a vector of daily precipitation (mm).
#' @param snowdepth a vector of snow depth measurements (cm) corresponding
#' to each time period for which snow temperature measurements are available (can be interpolated)
#' @param snowtemp an vector of snow temperature (deg C) measurements (see details)
#' @param snowem optionally, numeric value of snow emissivity
#' @param zm optionally, numeric value of roughness length for momentum transfer of snow surface (m)
#' @param umin optionally, numeric value indicating minimum wind speed used in conductivity calculations (m/s)
#' @param astc optionally, numeric value indicating the temperature at which precipitation falls as snow (deg C)
#' @param plotout optional logical indicating whether to plot fitted relationship
#' @return an object of class `SnowTparams` used for fitting the snow temperature model
#' @export
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. However, entries do not need to be hourly. They canshould correspond to
#' the time for which snow temperature measurements are available.
#' @examples
#' require(NicheMapR)
#' # Derive estimates of snow temperature (and depth) using NicheMapR
#' nmrout<-runNMRSnow(weather, precd, 67.367, 26.629, snowenv = "Taiga", ALTT = 199)
#' # Use NicheMapR outputs as inputs to fitsnowtemp
#' STparams<-fitsnowtemp(climdata,precd,nmrout$SNOWDEP,nmrout$SNOWTEMP,plotout=TRUE)
fitsnowtemp<-function(weather,precd,snowdepth,snowtemp,snowem=0.99,zm=0.002,umin=0.5,astc=1.5,plotout=FALSE) {
  # set constants
  sb<-5.67*10^-8
  ph<-43
  zpd=6.5*zm
  zh<-0.2*zm
  # Compute net radiation
  prech<-rep(0,length(precd)*24)
  sel<-(c(1:length(precd))-1)*24+1
  prech[sel]<-precd
  prech<-ifelse(weather$temp>astc,0,prech)
  snowalb<-.snowalb(weather,prech,astc)
  Rabs<-(1-snowalb)*weather$swrad+weather$skyem*sb*snowem*(weather$temp+273.15)^4
  Rem<-sb*snowem*(weather$temp+273.15)^4
  Rnet<-Rabs-Rem
  # Compute conductivity
  u2<-weather$windspeed
  u2[u2<umin]<-umin
  gHa<-(0.4^2*ph*u2)/(log((2-zpd)/zm)*log((2-zpd)/zh))
  # Model
  dT<-snowtemp-weather$temp
  sel<-which(snowdepth>0 & snowtemp <= 0)
  m1<-summary(lm(dT[sel]~Rnet[sel]*gHa[sel]))
  sel<-which(snowdepth>0 & snowtemp > 0)
  m2<-summary(lm(dT[sel]~Rnet[sel]*gHa[sel]))
  STparams<-list(m1=m1,m2=m2)
  class(STparams)<-"SnowTparams"
  if (plotout) {
    pT<-m1$coef[1,1]+m1$coef[2,1]*Rnet+m1$coef[3,1]*gHa+m1$coef[4,1]*Rnet*gHa+weather$temp
    pT2<-m2$coef[1,1]+m2$coef[2,1]*Rnet+m2$coef[3,1]*gHa+m2$coef[4,1]*Rnet*gHa+weather$temp
    sel<-which(snowdepth>0 & snowtemp >0.5 | snowtemp < -0.5)
    plot(snowtemp[sel]~pT[sel],pch=15,cex=0.5,ylab="Predicted",xlab="Observed")
  }
  return(STparams)
}
#' Derives snow melt coefficient from empirical data
#'
#' The function `fitsnowdepth` estimates the snow melt coefficient from empirical data (or
#' the outputs of a more complex model) and returns this along with
#' estimates of snow depth and the RMS error
#'
#' @param weather a data.frame of weather variables (see details).
#' @param snowdepth a vector fo measured snow depths (cm). Can be interpolated.
#' @param precd a vector of daily precipitation (mm).
#' @param snowenv one of `Alpine`, `Maritime`, `Prairie`, `Taiga` or `Tundra` (see details)
#' @param STparams an object of class `SnowTparams` as returned by [fitsnowtemp()]
#' @param plotout optional logical indicating whether to plot observed and predicted data
#' @param snowem optionally, numeric value of snow emissivity
#' @param zm optionally, numeric value of roughness length for momentum transfer of snow surface (m)
#' @param umin optionally, numeric value indicating minimum wind speed used in conductivity calculations (m/s)
#' @param astc optionally, numeric value indicating the temperature at which precipitation falls as snow (deg C)
#' @return a list of the following:
#' @return `psnowdepth` predicted snow depth (cm)
#' @return `RMS` Root-mean-square error of predicted snow depths (cm)
#' @return `meltfact` a snow melt factor for use with [pSnow()]
#' @export
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. The entries in weather need to be hourly. If hourly
#' snow depth data are unavailable, they should be derived by interpolation.
#' The paramater `snowenv` is used to compute snow density
#' following Sturm et al (2010) J Hydrometeorology 11: 1380-1393.
#' @examples
#' # Fit snow depth
#' require(NicheMapR)
#' # Derive estimates of snow temperature using NicheMapR in place of observed data
#' nmrout<-runNMRSnow(climdata, precd, 67.367, 26.629, ALTT = 199)
#' # Derive Snow temperature parameters
#' STparams<-fitsnowtemp(climdata,precd,nmrout$SNOWDEP,nmrout$SNOWTEMP)
#' # Fit snow depth model and plot result (red = Observed)
#' meltfact<-fitsnowdepth(climdata,snowdepth$snowdepth,precd,STparams=STparams)$meltfact
fitsnowdepth<-function(weather,snowdepth,precd,snowenv="Taiga",STparams,plotout=T,
                       snowem=0.99,zm=0.002,umin=0.5,astc=1.5) {
  rms<-0
  for (i in 1:20) {
    mf<-i*10
    rms[i]<-.FitSnowOne(weather,snowdepth,precd,snowenv,STparams,mf,snowem,zm,umin,astc)
  }
  sel0<-which(rms==min(rms))[1]
  for (i in 1:20) {
    mf<-sel0*10-10+i
    rms[i]<-.FitSnowOne(weather,snowdepth,precd,snowenv,STparams,mf,snowem,zm,umin,astc)
  }
  sel1<-which(rms==min(rms))[1]
  for (i in 1:20) {
    mf<-sel0*10-10+sel1-1+i/10
    rms[i]<-.FitSnowOne(weather,snowdepth,precd,snowenv,STparams,mf,snowem,zm,umin,astc)
  }
  sel2<-which(rms==min(rms))[1]
  for (i in 1:20) {
    mf<-sel0*10-10+sel1-1+sel2/10-1/10+i/100
    rms[i]<-.FitSnowOne(weather,snowdepth,precd,snowenv,STparams,mf,snowem,zm,umin,astc)
  }
  sel3<-which(rms==min(rms))[1]
  for (i in 1:20) {
    mf<-sel0*10-10+sel1-1+sel2/10-1/10+sel3/100-1/100+i/1000
    rms[i]<-.FitSnowOne(weather,snowdepth,precd,snowenv,STparams,mf,snowem,zm,umin,astc)
  }
  sel4<-which(rms==min(rms))[1]
  mf<-sel0*10-10+sel1-1+sel2/10-1/10+sel3/100-1/100+sel4/1000
  # Generate snow predictions
  psnowd<-pSnow(weather,precd,snowenv,STparams,mf,snowem,zm,umin,astc)$snowdepth
  if (plotout) {
    mx<-max(psnowd,snowdepth)
    plot(psnowd,type="l",col=rgb(0,0,0,0.5),ylim=c(0,mx),xlab="Hour of Year",ylab="Snow depth",main="")
    par(new=T)
    plot(snowdepth,col=rgb(1,0,0,0.5),type="l",ylim=c(0,mx),xlab="",ylab="",main=paste("RMS:",round(rms[sel4],3)))
  }
  out<-list(psnowdepth=psnowd,RMS=rms[sel4],meltfact=mf)
  class(out)<-"SnowDparams"
  return(out)
}
.snowdens<-function(snowenv="Tundra") {
  densfun<-c(0.5975,0.2237,0.0012,0.0038)
  if (snowenv == "Maritime") densfun<-c(0.5979,0.2578,0.001,0.0038)
  if (snowenv == "Prairie") densfun<-c(0.594,0.2332,0.0016,0.0031)
  if (snowenv == "Tundra") densfun<-c(0.363,0.2425,0.0029,0.0049)
  if (snowenv == "Taiga") densfun<-c(0.217,0.217,0,0)
  densfun
}
#' Wrapper function for running NicheMapR snow model
#'
#' The function `runNMRSnow` runs the multi-layer snow model in `NicheMapR`
#' (Kearney and POrter 2017, Ecography 40: 664-674) and can be used
#' to derive hourly snow depths and temperatures for fitting the simpler model emulator
#' in `microsnow`. Requires `NicheMapR` to be installed (https://github.com/mrke/NicheMapR)
#'
#' @param weather a data.frame of weather variables (see details).
#' @param precd a vector of daily precipitation (mm).
#' @param lat latitude of location for which snow depth and temperature are required (decimal degrees)
#' @param long longitude of location for which snow depth and temperature are required (decimal degrees)
#' @param snowenv one of `Alpine`, `Maritime`, `Prairie`, `Taiga` or `Tundra` (see details)
#' @param zm optionally, numeric value of roughness length for momentum transfer of snow surface (m)
#' @param snowem optionally, snow emissivity (~0.99) speed used in conductivity calculations (m/s)
#' @param ALTT elevation of location (m)
#' @return a list of the following:
#' @return `SNOWDEP` vector of predicted snow depth (cm)
#' @return `SNOWTEMP` vector of predicted snow surface temperature (deg C)
#' @import NicheMapR
#' @export
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. The paramater `snowenv` is used to compute snow density
#' following Sturm et al (2010) J Hydrometeorology 11: 1380-1393.
runNMRSnow <- function(weather, precd, lat, long, snowenv = "Taiga", zm = 0.002, snowem = 0.99, ALTT = 0) {
  loc<-c(long,lat)
  tmehr<-as.POSIXlt(weather$obs_time,tz="UTC")
  nyears<-length(unique(tmehr$year))
  fail<-nyears*24*365
  ystart<-tmehr$year[1]+1900
  yfinish<-tmehr$year[length(tmehr)]+1900
  yearlist<-seq(ystart,(ystart+(nyears-1)),1)
  doy <- unique(as.numeric(strftime(tmehr, format = "%j")))
  ndays <- unique(paste(as.numeric(strftime(tmehr, format = "%j")),
                        as.numeric(strftime(tmehr, format = "%y"))))
  ndays <- length(ndays)
  doynum<-ndays
  ida<-ndays
  microdaily<-1
  daystart<-1
  # Set root properties
  L<-c(0,0,8.2,8,7.8,7.4,7.1,6.4,5.8,4.8,4,1.8,0.9,0.6,0.8,0.4,0.4,0,0)*10000
  R1<-0.001
  RW<-2.5e+10
  RL<-2e+06
  PC<- -1500
  SP<-10
  IM<-1e-06
  # Set snow properties
  snowtemp<-1.5
  densfun<-.snowdens(snowenv)
  snowdens<-densfun[1]
  snowmelt<-1
  undercatch<-1
  rainmelt<-0.0125
  grasshade<-0
  MAXSHADES<-rep(100,ndays)
  MINSHADES<-rep(0,ndays)
  intercept<-0 # Snow interception by canopy. Could allow for this in NMR
  x<-t(as.matrix(as.numeric(c(loc[1],loc[2]))))
  ALREF<-abs(trunc(x[1]))
  HEMIS<-ifelse(x[2]<0,2,1)
  ALAT<-abs(trunc(x[2]))
  AMINUT<-(abs(x[2])-ALAT)*60
  ALONG<- abs(trunc(x[1]))
  ALMINT<-(abs(x[1])-ALONG)*60
  azmuth<-0
  lat<-as.numeric(loc[2])
  long<-as.numeric(loc[1])
  VIEWF<-1
  # Set soil properties
  Density<-2.56
  Thcond<-2.5
  SpecHeat<-870
  PE<-rep(2.6,19)
  BB<-rep(5.2,19)
  BD<-rep(1.529643,19)
  KS<-rep(6.4e-05,19)
  BulkDensity<-BD[seq(1,19,2)]
  #########################
  lt<-tmehr$hour+tmehr$min/60+tmehr$sec/3600
  jd<-microctools::jday(tme=tmehr)
  sa<-microctools::solalt(lt,lat,long,jd,0)
  ZENhr<-90-sa
  ZENhr[ZENhr>90]<-90
  TAIRhr<-weather$temp
  SOLRhr<-weather$swrad*VIEWF
  sb<-5.67*10^-8
  IRDhr<-weather$skyem*sb*(weather$temp+273.15)^4
  RHhr<-weather$relhum
  RHhr[RHhr>100]<-100
  RHhr[RHhr<0]<-0
  e0<-microctools::satvap(TAIRhr)
  ea<-e0*(RHhr/100)
  eo<-1.24*(10*ea/(TAIRhr+273.15))^(1/7)
  CLDhr<-((weather$skyem-eo)/(1-eo))*100
  CLDhr[CLDhr<0]<-0
  CLDhr[CLDhr>100]<-100
  WNhr<-weather$windspeed
  WNhr[is.na(WNhr)]<-0.1
  PRESShr<-weather$pres*1000
  RAINFALL<-precd
  #RAINFALL[RAINFALL<0.1]<-0
  ZENhr2<-ZENhr
  ZENhr2[ZENhr2!=90]<-0
  dmaxmin<-function(x,fun) {
    dx <- t(matrix(x, nrow = 24))
    apply(dx, 1, fun)
  }
  TMAXX<-dmaxmin(TAIRhr,max)
  TMINN<-dmaxmin(TAIRhr,min)
  CCMAXX<-dmaxmin(CLDhr,max)
  CCMINN<-dmaxmin(CLDhr,min)
  RHMAXX<-dmaxmin(RHhr,max)
  RHMINN<-dmaxmin(RHhr,min)
  WNMAXX<-dmaxmin(WNhr,max)
  WNMINN<-dmaxmin(WNhr,min)
  PRESS<-dmaxmin(PRESShr,min)
  ###
  slope <- 0
  azmuth <- 0
  relhum <- 1
  optdep.summer<-as.data.frame(rungads(loc[2],loc[1],relhum, 0))
  optdep.winter<-as.data.frame(rungads(loc[2],loc[1],relhum, 1))
  optdep<-cbind(optdep.winter[,1],rowMeans(cbind(optdep.summer[,2],optdep.winter[,2])))
  optdep<-as.data.frame(optdep)
  colnames(optdep)<-c("LAMBDA","OPTDEPTH")
  a<-lm(OPTDEPTH~poly(LAMBDA,6,raw=TRUE),data=optdep)
  LAMBDA<-c(290,295,300,305,310,315,320,330,340,350,360,370,380,390,400,420,440,460,480,
            500,520,540,560,580,600,620,640,660,680,700,720,740,760,780,800,820,840,860,
            880,900,920,940,960,980,1000,1020,1080,1100,1120,1140,1160,1180,1200,1220,
            1240,1260,1280,1300,1320,1380,1400,1420,1440,1460,1480,1500,1540,1580,1600,
            1620,1640,1660,1700,1720,1780,1800,1860,1900,1950,2000,2020,2050,2100,2120,
            2150,2200,2260,2300,2320,2350,2380,2400,2420,2450,2490,2500,2600,2700,2800,
            2900,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000)
  TAI<-predict(a,data.frame(LAMBDA))
  RAINFALL<-RAINFALL
  ALLMINTEMPS<-TMINN
  ALLMAXTEMPS<-TMAXX
  ALLTEMPS<-cbind(ALLMAXTEMPS,ALLMINTEMPS)
  WNMAXX<-WNMAXX
  WNMINN<-WNMINN
  WNhr<-WNhr
  REFLS<-rep(0.15,ndays)
  PCTWET <-rep(0,ndays)
  soilwet<-RAINFALL
  soilwet[soilwet<=1.5]<-0
  soilwet[soilwet>0]<-90
  if (ndays < 1) PCTWET<-pmax(soilwet,PCTWET)
  Intrvls<-rep(0,ndays)
  Intrvls[1]<-1
  Numtyps<-10
  Nodes<-matrix(data=0,nrow=10,ncol=ndays)
  Nodes[1:10,]<-c(1:10)
  ALREF<-abs(trunc(x[1]))
  HEMIS<-ifelse(x[2]<0,2,1)
  ALAT<-abs(trunc(x[2]))
  AMINUT<-(abs(x[2])-ALAT)*60
  ALONG<-abs(trunc(x[1]))
  ALMINT<-(abs(x[1])-ALONG)*60
  avetemp<-(sum(TMAXX)+sum(TMINN))/(length(TMAXX)*2)
  soilinit<-rep(avetemp,20)
  tannul<-mean(unlist(ALLTEMPS))
  deepsoil<-rep(mean(TAIRhr),ndays)
  SLES<-matrix(nrow=ndays,data=0)
  SLES<-SLES+0.99
  SoilMoist_Init <- c(0.1,0.12,0.15,0.2,0.25,0.3,0.3,0.3,0.3,0.3)
  moists2<-matrix(nrow=10,ncol=ndays,data=0)
  moists2[1:10,]<-SoilMoist_Init
  moists<-moists2
  soilprops<-matrix(data=0,nrow=10,ncol=5)
  soilprops[,1]<-BulkDensity
  soilprops[,2]<-min(0.26,1-BulkDensity/Density)
  soilprops[,3]<-Thcond
  soilprops[,4]<-SpecHeat
  soilprops[,5]<-Density
  soilprops[1:2,3] <- 0.2
  soilprops[1:2,4] <- 1920
  hourly<-1
  if (length(precd) == length(TAIRhr)) {
    rainhourly<-1
    RAINhr<-precd
    RAINhr[RAINhr<0.1]<-0
    raintest<-RAINhr
  } else if (length(precd) == length(TAIRhr)/24) {
    rainhourly<-0
    RAINhr<-rep(precd/24,each=24)
  } else stop("Rainfall must be daily or hourly")
  snowmodel<-1
  D0<-0.0
  microinput<-c(ndays,zm,1.5,0.05,2,Numtyps,0,0,0,0,1,ida,
                HEMIS,ALAT,AMINUT,ALONG,ALMINT,ALREF,slope,azmuth,ALTT,1,
                microdaily,tannul,0.0167238,VIEWF,snowtemp,snowdens,snowmelt,undercatch,
                1,runshade=0,1,1000,0,snowmodel,rainmelt,
                0,densfun,hourly,rainhourly,0,0,RW,PC,RL,SP,R1,IM,
                500,0,0,fail,0,intercept,grasshade,0,0,D0)
  doy1<-matrix(data=0,nrow=ndays,ncol=1)
  SLES1<-matrix(data=0,nrow=ndays,ncol=1)
  MAXSHADES1<-matrix(data=0,nrow=ndays,ncol=1)
  MINSHADES1<-matrix(data=0,nrow=ndays,ncol=1)
  TMAXX1<-matrix(data=0,nrow=ndays,ncol=1)
  TMINN1<-matrix(data=0,nrow=ndays,ncol=1)
  CCMAXX1<-matrix(data=0,nrow=ndays,ncol=1)
  CCMINN1<-matrix(data=0,nrow=ndays,ncol=1)
  RHMAXX1<-matrix(data=0,nrow=ndays,ncol=1)
  RHMINN1<-matrix(data=0,nrow=ndays,ncol=1)
  WNMAXX1<-matrix(data=0,nrow=ndays,ncol=1)
  WNMINN1<-matrix(data=0,nrow=ndays,ncol=1)
  REFLS1<-matrix(data= 0,nrow=ndays,ncol=1)
  PCTWET1<-matrix(data=0,nrow=ndays, ncol=1)
  RAINFALL1<-matrix(data=0,nrow=ndays,ncol=1)
  tannul1<-matrix(data=0,nrow=ndays,ncol=1)
  moists1<-matrix(data=0,nrow=10,ncol=ndays)
  SLES1[1:ndays]<-SLES
  MAXSHADES1[1:ndays]<-MAXSHADES
  MINSHADES1[1:ndays]<-MINSHADES
  TMAXX1[1:ndays]<-TMAXX
  TMINN1[1:ndays]<-TMINN
  CCMAXX1[1:ndays]<-CCMAXX
  CCMINN1[1:ndays]<-CCMINN
  RHMAXX1[1:ndays]<-RHMAXX
  RHMINN1[1:ndays]<-RHMINN
  WNMAXX1[1:ndays]<-WNMAXX
  WNMINN1[1:ndays]<-WNMINN
  REFLS1[1:ndays]<-REFLS
  PCTWET1[1:ndays]<-PCTWET
  raind<-matrix(RAINhr,ncol=24,byrow=TRUE)
  raind<-apply(raind,1,sum)
  RAINFALL1[1:ndays]<-raind
  tannul1[1:ndays]<-tannul
  moists1[1:10, 1:ndays] <- moists
  tides<-matrix(data=0,nrow=24*ndays,ncol=3)
  TIMAXS<-c(1,1,0,0)
  TIMINS<-c(0,0,1,1)
  LAI<-rep(0,ndays)
  DEP <- c(0,2.5,5,10,15,20,30,50,100,200)
  DD <- rep(2.65,19)
  hori <-rep(0, 36)
  micro<-list(tides=tides,microinput=microinput,doy=doy,SLES=SLES1,DEP=DEP,Nodes=Nodes,
              MAXSHADES=MAXSHADES,MINSHADES=MINSHADES,TIMAXS=TIMAXS,TIMINS=TIMINS,TMAXX=TMAXX1,
              TMINN=TMINN1,RHMAXX=RHMAXX1,RHMINN=RHMINN1,CCMAXX=CCMAXX1,CCMINN=CCMINN1,
              WNMAXX=WNMAXX1,WNMINN=WNMINN1,TAIRhr=TAIRhr,RHhr=RHhr,WNhr=WNhr,CLDhr=CLDhr,
              SOLRhr=SOLRhr,RAINhr=RAINhr,ZENhr=ZENhr,IRDhr=IRDhr,REFLS=REFLS1,PCTWET=PCTWET1,
              soilinit=soilinit,hori=hori,TAI=TAI,soilprops=soilprops,moists=moists1,
              RAINFALL=RAINFALL1,tannulrun=deepsoil,PE=PE,KS=KS,BB=BB,BD=BD,DD=DD,L=L,LAI=LAI)
  microut<-microclimate(micro)
  metout<-as.data.frame(microut$metout)
  SNOWDEP<-metout$SNOWDEP
  snow<-microut$sunsnow
  snowtemp<-snow[,3]
  return(list(SNOWDEP=SNOWDEP,SNOWTEMP=snowtemp))
}
#' Returns inbuilt snow melt factors for different snow environments
#'
#' The function `getmeltf` returns pre-calculated snow melt parameters,
#' derived using [fitsnowdepth()] for different snow environments.
#'
#' @param snowenv one of `Alpine`, `Maritime`, `Prairie`, `Taiga` or `Tundra`
#' @return a single numeric value indicating the melt factor
#' @export
getmeltf <- function(snowenv="Taiga") {
  meltf<-0.2
  if (snowenv=="Maritime") meltf<-0.133
  if (snowenv=="Prairie") meltf<-0.091
  if (snowenv=="Taiga") meltf<-0.051
  if (snowenv=="Tundra") meltf<-0.042
  meltf
}
