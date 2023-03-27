#' @title Create object of class microsnowin
#' @description The function `snowmodelin` creates an object of class `microsnowin`
#' which unpacks various component inputs and reformats as required
#' for running the model in hourly timesteps. Here it is assumed that the input
#' weather data are a data.frame - i.e. not spatially variable.
#' @param weather a data.frame of weather variables (see details)
#' @param precd a vector of daily precipitation
#' @param snow an object of class `snow` as created by [modelsnowdepth()]
#' @param STparams now temperature model coefficients as derived by [fitsnowtemp()]
#' @param dtm  PackedSpatRaster or SpatRast object of elevations (see details)
#' @param pai a single numeric value, matrix, SpatRast or array of plant area idex values.
#' @param hgt ackedSpatRaster or SpatRast of vegetation heights
#' @param x optional single numeric value, SpatRast, matrix or array of leaf angle coefficients (see details)
#' @param clump a single numeric value or array of values between 0 and 1 indicating the fraction of radiation passing through larger gaps in the canopy,
#' @param windhgt height above round of wind speed measurement (m)
#' @return an object of class `microsnowin` used by [runmicrosnow_hr()]
#' @details The format and and units of `weather` must follow that in the example
#' dataset `climdata`. The x,y and z units of `dtm` must be all be in
#' metres and the coordinate reference system must be defined. Canopy heights
#' are assumed time-invariant over the period the model is run, and must
#' be supplied as a SpatRast with coordinates and extent matching `dtm`.
#' `pai`, `x` and `clump` can either be assumed spatially and temporally invariant`
#' (supplied as single numeric values), temporally invariant (supplied as matrices or
#' SpatRasts with extent and coordinate reference system matching `dtm` or spatially
#' and temporally variant - x & y dims must match `dtm` z dim must of length `weather$temp`
#' ). The coefficient x represents the ratio of of vertical to horizontal projections of
#' leaf foliage (see Campbell 1986 Agric For Meteorol, 36: 317-21 for details). As
#' not all wind measurements are at reference height, the height of the wind speed measurement
#' must be specified if not 2 m. To enable calculation of below-canopy wind profiles
#' in tall canopy, the wind speed is adjusted to give values for a height 2 m above
#' the maximum vegetation height, using the wind-height profile for a reference grass
#' surface.
#' @rdname snowmodelin
#' @export
snowmodelin <- function(weather,precd,snow,STparams,dtm,pai,hgt,x,clump,windhgt=2) {
  # Unpack PackedSpatRasters if necessary
  if (class(dtm)[1]=="PackedSpatRaster") dtm<-rast(dtm)
  if (class(pai)[1]=="PackedSpatRaster") pai<-rast(pai)
  if (class(hgt)[1]=="PackedSpatRaster") hgt<-rast(hgt)
  soiltype<-7
  n<-length(weather$temp)
  pai<-.unpackpai(pai,n,dtm)
  x<-.unpackpai(x,n,dtm)
  clump<-.unpackpai(clump,n,dtm)
  hgt<-.rta(hgt,n)
  # correct wind height
  mxhgt<-max(.is(hgt),na.rm=T)+2
  weather$windspeed<-.windcorrect(weather$windspeed,windhgt,mxhgt)
  r<-dtm
  ll<-.latlongfromrast(r)
  tme<-as.POSIXlt(weather$obs_time,tz="UTC")
  tc<-.vta(weather$temp,r)
  difr<-.vta(weather$difrad,r)
  di<-weather$swrad-weather$difrad
  di[di<0]<-0
  jd<-.jday(tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  sa<-.solalt(lt,ll$lat,ll$long,jd)
  ze<-90-sa
  si<-cos(ze*(pi/180))
  si[si<0]<-0
  dni<-di/si
  dni[is.na(dni)]<-0
  dni[dni>1352]<-1352
  dirr<-.vta(dni,r)
  dp<-weather$difrad/weather$swrad
  dp[is.na(dp)]<-0.5; dp[dp<0]<-0; dp[dp>1]<-1
  dp<-.vta(dp,dtm)
  skyem<-.vta(weather$skyem,r)
  # == (0b) Derived variables
  estl<-.satvap(weather$temp)
  ea<-(weather$relhum/100)*estl
  tdew<-.dewpoint(ea,weather$temp)
  estl<-.vta(estl,r)
  ea<-.vta(ea,r)
  tdew<-.vta(tdew,r)
  # Calculate snow refletcance
  astc=1.5
  prech<-rep(0,n)
  sel<-(c(1:length(precd))-1)*24+1
  prech[sel]<-precd
  prech<-ifelse(weather$temp>astc,0,prech)
  snowalb<-.snowalb(weather,prech,astc)
  lref<-.vta(snowalb,dtm)
  ltra<-lref*0
  pk<-.vta(weather$pres,r)
  gref<-lref
  soiltype<-.unpackpai(soiltype,n,dtm)[,,1]
  soilp<-.soilinit(soiltype)
  out<-list(tme=tme,tc=tc,difr=difr,dirr=dirr,dp=dp,skyem=skyem,
            estl=estl,ea=ea,tdew=tdew,pk=pk,pai=pai,vegx=x,lref=lref,ltra=ltra,veghgt=hgt,
            gsmax=1000,clump=clump,gref=gref,rho=soilp$rho,Vm=soilp$Vm,leafd=0.2,
            Vq=soilp$Vq,Mc=soilp$Mc,soilb=soilp$soilb,psi_e=soilp$psi_e,Smax=soilp$Smax,
            dtm=dtm,lat=ll$lat,long=ll$long,climdata=weather,
            prec=precd,snow=snow,maxhgt=mxhgt,progress=0)
  class(out) <-"microsnowin"
  out
}
#' @title Create object of class microsnowin with weather data as an array
#' @description The function `snowmodelina` creates an object of class `microsnowin`
#' which unpacks various component inputs and reformats as required
#' for running the model in hourly timesteps. Here it is assumed that the input
#' weather data are as arrays - i.e. variable in space
#' @param climarray a list of arrays of weather variables (see details). See also [nctoarray()]
#' @param precarray an array of daily precipitation (see details)
#' @param tme an object of class POSIXlt giving the dates and times for each weather variable stroed in the array
#' @param r a SpatRaster object giving the resolution, spatial extent, and projection of the weather data (see details)
#' @param altcorrect a single numeric value indicating whether to apply an elevational lapse rate correction to temperatures (0 = no correction, 1 = fixed lapse rate correction, 2 = humidity-dependent variable lapse rate correction, see details)
#' @param snow an object of class `snow` as created by [modelsnowdepth()]
#' @param STparams now temperature model coefficients as derived by [fitsnowtemp()]
#' @param dtm a SpatRaster object of elevations (see details)
#' @param pai a single numeric value, matrix, SpatRast or array of plant area idex values.
#' @param hgt ackedSpatRaster or SpatRast of vegetation heights
#' @param x optional single numeric value, SpatRast, matrix or array of leaf angle coefficients (see details)
#' @param clump a single numeric value or array of values between 0 and 1 indicating the fraction of radiation passing through larger gaps in the canopy,
#' @param windhgt height above round of wind speed measurement (m)
#' @return an object of class `microsnowin` used by [runmicrosnow_hr()]
#' @details The units of `climarray` must follow those in the dataset `climdata`.
#' It must be a list with each component of the list an array, named using the same
#' names as the column headers in weather (e.g. temp for temperature), excluding `obs_time`.
#' Dimensions 1 and 2 of the array must be the same as `r` and dimension 3 must have
#' the same length as `tme`. If `r` has a different resolution to `dtm` the climate
#' data are resampled to match the resolution of `dtm`. The array of Plant Area index values in `vegp` must
#' of the same x and y dims as `dtm` but can contain any number of repeated
#' measures up to the number of entries in `tme`. Data are interpolated to the
#' time increment of `tme`. Other vegetation paramaters, including vegetation
#' height are assumed time-invarient. The SpatRaster datasets in `soilc` must have
#' the same x and y dims as `dtm`. The x,y and z units of `dtm` must be all be in
#' meters and the coordinate reference system must be defined. If `altcorrect`>0,
#' and the dimensions of `r` are not identical to those of `dtm`, the elevation
#' difference between each pixel of the dtm and the dtm coarsened to the resolution of
#' `r` is calculated and an elevational lapse rate correction is applied to the
#' temperature data to accoutn for these elevation differences. If `altcorrect`=1,
#' a fixed lapse rate of 5 degrees per 100m is applied. If `altcorrect`=2, humidity-dependent
#' lapse rates are calculate and applied.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  require(terra)
#'  require(NicheMapR)
#'  tme<-as.POSIXlt(climdata$obs_time,tz="UTC")
#' sel<-which(tme$mon+1 == 4 & tme$year+1900==2018)
#' weather<-climdata[sel,]
#' tmed<-as.POSIXlt(c(0:546)*24*3600,origin="2017-10-02",tz="UTC")
#' sel<-which(tmed$mon+1 == 4 & tmed$year+1900==2018)
#' precapr<-microsnow::precd[sel]
#' # Derive estimates of snow melt temperature parameters using NicheMapR
#' nmrout<-runNMRSnow(weather,precapr,66.215,29.293,ALTT = 360)
#' STparams<-fitsnowtemp(weather,precapr,nmrout$SNOWDEP,nmrout$SNOWTEMP)
#' meltfact<-fitsnowdepth(weather,nmrout$SNOWDEP,precapr,STparams=STparams)$meltfact
#' # Run snow model with defaults (takes 20 seconds)
#' snowdepth<-nmrout$SNOWDEP
#' snow<-modelsnowdepth(weather,precapr,snowdepth,dtm,pai,hgt,STparams,meltfact,spatial=T)
#' # ======== Create dummy array datasets ========= #
#' # ~~ Function to turn vector into 5 x 5 array ~~ #
#' toarray<-function(x,xdim=5,ydim=5) {
#'   xl<-rep(x,each=ydim*xdim)
#'   xl<-array(xl,dim=c(ydim,xdim,length(x)))
#'   xl
#' }
#' # ~~ Create list of climate variable arrays ~~ #
#' climarray<-list(temp=toarray(weather$temp),
#'                 relhum=toarray(weather$relhum),
#'                 pres=toarray(weather$pres),
#'                 swrad=toarray(weather$swrad),
#'                 difrad=toarray(weather$difrad),
#'                 skyem=toarray(weather$skyem),
#'                 windspeed=toarray(weather$windspeed),
#'                 winddir=toarray(weather$winddir))
#' # ~~ Create precipitation array ~~ #
#' precarray<-toarray(precapr)
#' # ~~ Create other variables ~~ #
#' tme<-as.POSIXlt(weather$obs_time)
#' microsnow<-snowmodelina(climarray,precarray,tme,r,altcorrect=0,snow,STparams,dtm,pai,hgt,x,clump,windhgt=2)
#'  }
#' }
#' @rdname snowmodelina
#' @export
snowmodelina<-function(climarray,precarray,tme,r,altcorrect=0,snow,STparams,dtm,pai,hgt,x,clump,windhgt=2) {
  # Unpack PackedSpatRasters if necessary
  if (class(dtm)[1]=="PackedSpatRaster") dtm<-rast(dtm)
  if (class(pai)[1]=="PackedSpatRaster") pai<-rast(pai)
  if (class(hgt)[1]=="PackedSpatRaster") hgt<-rast(hgt)
  soiltype<-7
  n<-length(weather$temp)
  pai<-.unpackpai(pai,n,dtm)
  x<-.unpackpai(x,n,dtm)
  clump<-.unpackpai(clump,n,dtm)
  hgt<-.rta(hgt,n)
  # correct wind height
  # Correct wind speed
  mxhgt<-max(.is(hgt),na.rm=T)+2
  climarray$windspeed<-.windcorrect(climarray$windspeed,windhgt,mxhgt)
  # Create weather and precipitation dataset
  weather<-.catoweather(climarray,tme)
  precip<-apply(precarray,3,mean,na.rm=T)
  # Resample climdata
  tc<-suppressWarnings(.resa(climarray$temp,r,dtm))
  difr<-suppressWarnings(.resa(climarray$difrad,r,dtm))
  # ~~ Calculate dni
  di<-climarray$swrad-climarray$difrad
  jd<-.jday(tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  ll<-.latslonsfromr(r)
  n<-length(tme)
  sa<-.solalt(.vta(lt,r),.rta(rast(ll$lats),n),.rta(rast(ll$lons),n),.vta(jd,r))
  ze<-90-sa
  si<-cos(ze*(pi/180))
  si[si<0]<-0
  dirr<-suppressWarnings(.resa(di/si,r,dtm))
  dirr[is.na(dirr)]<-0
  dirr[dirr>1352]<-1352; dirr[dirr<0]<-0
  dp<-climarray$difrad/climarray$swrad
  dp[is.na(dp)]<-0.5; dp[dp<0]<-0; dp[dp>1]<-1
  dp<-suppressWarnings(.resa(dp,r,dtm))
  skyem<-suppressWarnings(.resa(climarray$skyem,r,dtm))
  pk<-suppressWarnings(.resa(climarray$pres,r,dtm))
  # ~~ Derived variables
  estl<-.satvap(climarray$temp)
  ea<-(climarray$relhum/100)*estl
  tdew<-.dewpoint(ea,climarray$temp)
  estl<-suppressWarnings(.resa(estl,r,dtm))
  ea<-suppressWarnings(.resa(ea,r,dtm))
  tdew<-suppressWarnings(.resa(tdew,r,dtm))
  # Calculate snow refletcance
  astc=1.5
  prech<-rep(0,n)
  sel<-(c(1:length(precd))-1)*24+1
  prech[sel]<-precd
  prech<-ifelse(weather$temp>astc,0,prech)
  snowalb<-.snowalb(weather,prech,astc)
  lref<-.vta(snowalb,dtm)
  ltra<-lref*0
  gref<-lref
  soiltype<-.unpackpai(soiltype,n,dtm)[,,1]
  soilp<-.soilinit(soiltype)
  # Elevation correction
  # ~~ Fix lapse rate
  if (altcorrect>0) {
    dtmc<-resample(dtm,r)
    dtmc[is.na(dtmc)]<-0
    dtmc<-resample(dtmc,dtm)
    elevd<-dtmc-dtm
  }
  if (altcorrect==1) {
    tcdif<-elevd*(5/1000)
    tc<-.rta(tcdif,n)+tc
  }
  if (altcorrect==2) {
    lr<-.lapserate(climarray$temp,climarray$relhum,climarray$pres)
    lr<-suppressWarnings(.resa(lr,r,dtm))
    tcdif<-.rta(elevd,n)*lr
    tc<-tcdif+tc
  }
  ll<-.latlongfromrast(dtm)
  out<-list(tme=tme,tc=tc,difr=difr,dirr=dirr,dp=dp,skyem=skyem,
            estl=estl,ea=ea,tdew=tdew,pk=pk,pai=pai,vegx=x,lref=lref,ltra=ltra,veghgt=hgt,
            gsmax=1000,clump=clump,gref=gref,rho=soilp$rho,Vm=soilp$Vm,leafd=0.2,
            Vq=soilp$Vq,Mc=soilp$Mc,soilb=soilp$soilb,psi_e=soilp$psi_e,Smax=soilp$Smax,
            dtm=dtm,lat=ll$lat,long=ll$long,climdata=weather,prec=precip,snow=snow,
            maxhgt=mxhgt,progress=0)
  class(out) <-"microsnowin"
  out
}
#' Snow two-stream radiation function
.snowtwostream<-function(microsnow, reqhgt = 0.05, pai_a = NA, slr = NA, apr = NA, hor = NA) {
  # Calculate horizon angle of NA
  dtm<-microsnow$dtm
  dtm[is.na(dtm)]<-0
  if (class(hor)[1]=="logical") {
    hor<-array(NA,dim=c(dim(dtm)[1:2],24))
    for (i in 1:24) hor[,,i]<-.horizon(dtm,(i-1)*15)
  }
  # Run snow rad if not run
  if (class(microsnow$snow$Radparams) == "logical") {
    snowdepth<-apply(microsnow$snow$gsnowdepth,3,mean,na.rm=T)
    epai<-.epaif1(microsnow$pai,microsnow$veghgt,microsnow$snow$gsnowdepth)
    prech<-rep(0,length(microsnow$prec)*24)
    sel<-(c(1:length(precd))-1)*24+1
    prech[sel]<-precd
    prech<-ifelse(microsnow$climdata$temp>1.5,0,prech)
    snowalb<-.snowalb(microsnow$climdata,prech,1.5)
    snowalb<-.vta(snowalb,dtm)
    Radparams<-.snowrad(microsnow$climdata,snowdepth,microsnow$dtm,hor,
                        microsnow$lat,microsnow$long,epai,microsnow$vegx,snowalb,
                        0.99,slr,apr,microsnow$clump,clumpd="hourly")
  } else Radparams<-microsnow$snow$Radparams
  # Assign canopy radiation
  microsnow$radCsw<-Radparams$radCsw
  microsnow$radClw<-Radparams$radClw
  # Calculate additional canopy variables needed
  if (class(pai_a)[1]=="logical") {
    fd<-foliageden(reqhgt,microsnow$veghgt,microsnow$pai)
    microsnow$pai_a<-fd$pai_a
    microsnow$leafden<-fd$leafden
  }
  twostreamp<-Radparams$twostreamp
  n<-(microsnow$veghgt-reqhgt)/microsnow$veghgt
  n[n<0]<-0
  pai_a<-microsnow$pai_a/(1-microsnow$clump*n)
  pai_t<-microsnow$pai/(1-microsnow$clump)
  Kc<-with(twostreamp,kkd$kd/kkd$k0)
  # Calculate additional terrain variables needed
  tme<-microsnow$tme
  jd<-.jday(tme)
  lt<-tme$hour+tme$min/60+tme$sec/3600
  alt<-.solalt(lt,microsnow$lat,microsnow$long,jd)
  salt<-.vta(alt,dtm)*(pi/180)
  azi<-.solazi(lt,microsnow$lat,microsnow$long,jd)
  si<-.solarindex(dtm,salt,azi,slr,apr)
  i<-round(azi/15,0)+1; i[i==25]<-1
  hora<-hor[,,i]
  shadowmask<-hora*0+1
  shadowmask[hora>tan(salt)]<-0
  si<-si*shadowmask
  msl<-tan(apply(atan(hor),c(1,2),mean))
  svf<-0.5*cos(2 *msl)+0.5
  svfa<-.rta(rast(svf),length(jd))
  # === Calculate up and downstream at reqhgt
  # Downward
  Rbgm<-with(twostreamp,(1-(microsnow$clump^(Kc*n)))*exp(-kkd$kd*pai_a)+(microsnow$clump^(Kc*n)))
  Rdbm<-with(twostreamp,(1-(microsnow$clump^(2*n)))*((p8/sig)*exp(-kkd$kd*pai_a)+p9*exp(-h*pai_a)+p10*exp(h*pai_a)))
  Rdbm[salt<=0]<-0
  Rddm<-with(twostreamp,(1-microsnow$clump^(2*n))*(p3*exp(-h*pai_a)+p4*exp(h*pai_t))+microsnow$clump^(2*n))
  microsnow$Rbdown<-microsnow$dirr*si*Rbgm
  microsnow$Rddown<-microsnow$difr*svfa*Rddm+microsnow$dirr*sin(salt)*Rdbm
  # Upward
  Rdbm<-with(twostreamp,(1-(microsnow$clump^(2*n)))*((p5/sig)*exp(-kkd$kd*pai_a)+p6*exp(-h*pai_a)+p7*exp(h*pai_a)))
  Rdbm[salt<=0]<-0
  Rddm<-with(twostreamp,(1-(microsnow$clump^(2*n)))*(p1*exp(-h*pai_a)+p2*exp(h*pai_t))+microsnow$clump^(2*n))
  microsnow$Rdup<-microsnow$difr*svfa*Rddm+microsnow$dirr*si*Rdbm
  # === (1n) Calculate leaf swabs
  microsnow$radLsw<-with(twostreamp,0.5*(1-microsnow$lref)*(microsnow$Rddown+microsnow$Rdup+kkd$k*si*microsnow$Rbdown))
  microsnow$k<-twostreamp$kkd$k
  microsnow$kc<-Kc
  microsnow$progress<-1
  # Clean micro
  microsnow$lat<-NULL
  microsnow$long<-NULL
  return(microsnow)
}
#' Snow hourly surface tmep function
.snowsoiltemp_hr  <- function(microsnow, reqhgt = 0.05, pai_a = NA, slr = NA, apr = NA,
                             hor = NA, wsa = NA) {
  # run two-stream and wind functions if not run
  if (microsnow$progress<1) {
    microsnow<-.snowtwostream(microsnow,reqhgt,pai_a,slr,apr,hor)
  }
  # Calculate net radiation at ground
  T0<-microsnow$snow$snowtempG
  # Calculate soil conductivity and specific heat capacity
  cs<-(2400*microsnow$rho/2.64+4180*microsnow$Smax)
  phs<-(microsnow$rho*(1-microsnow$Smax)+microsnow$Smax)*1000
  pcs<-cs*phs # Total specific heat capacity
  frs<-microsnow$Vm+microsnow$Vq
  c1<-(0.57+1.73*microsnow$Vq+0.93*microsnow$Vm)/(1-0.74*microsnow$Vq-0.49*microsnow$Vm)-2.8*frs*(1-frs)
  c2<-1.06*microsnow$rho*microsnow$Smax
  c3<-1+2.6*microsnow$Mc^-0.5
  c4<-0.03+0.7*frs^2
  k<-c1+c2*microsnow$Smax-(c1-c4)*exp(-(c3*microsnow$Smax)^4)
  # Calculate damping depth etc
  om<-(2*pi)/(24*3600)
  ka<-k/pcs
  DD<-sqrt((2*ka)/om)
  # Calculate G
  A0<-aperm(apply(T0,c(1,2),.A0f),c(2,3,1))
  r<-rast(A0[,,1])
  cda<-microsnow$climdata
  t0<-.vta(.t0f(cda$temp),r)
  tt<-(c(1:length(cda$temp))-1)%%24
  tt<-.vta(tt*3600,r)
  hiy<-length(microsnow$tme)
  G<-(sqrt(2)*A0*.rta(rast(k),hiy)*sin(om*(tt-t0)+(pi/4)))/.rta(rast(DD),hiy)
  # Save elements to micro
  microsnow$T0<-T0
  microsnow$ka<-ka
  microsnow$DD<-DD
  microsnow$G<-G
  microsnow$progress<-2
  # Clean micro
  microsnow$prec<-NULL
  microsnow$rho<-NULL
  microsnow$Vm<-NULL
  microsnow$Vq<-NULL
  microsnow$Mc<-NULL
  #micro$radGsw<-NULL
  #micro$radGlw<-NULL
  return(microsnow)
}
#' Snow wind function
.snowwind <- function(microsnow, reqhgt = 0.05, pai_a = NA, xyf = 1, zf = NA, slr = NA, apr = NA,
                     hor = NA, wsa = NA, gmn = 0.3) {
  # Calculate zfactor
  tme<-microsnow$tme
  ed<-as.numeric(tme[length(tme)])
  st<-as.numeric(tme[1])
  ti<-round(((ed-st)/3600)/(length(tme)-1))
  if (microsnow$progress<2) {
    if (ti==1) {
      microsnow<-.snowsoiltemp_hr(microsnow,reqhgt,pai_a,soilinit,tfact,slr,apr,hor,wsa,twi,soilmcoefs,soiltcoefs)
    } else stop("Need to run soiltemp_dy")
  }
  # Calculate wind shelter coefficient
  vhgt<-.rast(microsnow$veghgt[,,1],microsnow$dtm)
  dsm<-microsnow$dtm+vhgt
  s<-1
  if (res(microsnow$dtm)[1]<=100) s<-10
  ws<-.windshelter(microsnow$climdata$winddir,dsm,2,s,wsa)
  if (is.na(zf)) zf<-ifelse(ti==1,30*24,30)
  # Calculate roughness lengths  etc
  snowd<-microsnow$snow$gsnowdepth/100
  hgta<-microsnow$veghgt
  sel<-which(snowd>=hgta)
  # zero plane displacement
  d<-.zeroplanedis(hgta,microsnow$pai)
  d[sel]<-snowd[sel]
  zm<-.roughlength(hgta,microsnow$pai,d)
  zm[zm<0.003]<-0.003
  zm[sel]<-0.003
  # smooth roughness lengths
  dm<-.roughresample(d,microsnow$dtm,xyf)
  zma<-.roughresample(zm,microsnow$dtm,xyf)
  uref<-.vta(microsnow$climdata$windspeed,microsnow$dtm)
  uf<-suppressWarnings((0.4*uref)/log((microsnow$maxhgt-dm)/zma))*ws
  # Calculate wind speed
  s1<-which(reqhgt>=hgta & reqhgt>snowd) # reqhgt above canopy and snow
  s2<-which(reqhgt<hgta & reqhgt>snowd) # reqhgt below canopy but above snow
  uz<-uf*0 # # reqhgt at or below snowdepth
  uz[s1]<-(uf[s1]/0.4)*suppressWarnings(log((reqhgt-dm[s1])/zma[s1])) # above canopy
  sel<-which(is.na(uz[s1]))
  uz[s1[sel]]<-uf[s1[sel]]
  uh<-(uf[s2]/0.4)*log((hgta[s2]-d[s2])/zm[s2]) # at canopy top
  sel<-which(is.na(uh))
  uh[sel]<-uf[s2[sel]]
  Be<-0.205*microsnow$pai[s2]^0.445+0.1
  a<-microsnow$pai[s2]/hgta[s2]
  Lc<-(0.25*a)^-1
  Lm<-2*Be^3*Lc
  uz[s2]<-uh*exp(Be/Lm*(reqhgt-hgta[s2]))
  wsp<-.vta(microsnow$climdata$windspeed,microsnow$dtm)
  sel<-which(uz>wsp)
  uz[sel]<-wsp[sel]
  sel<-which(uz<uf)
  uz[sel]<-uf[sel]
  # Calculate convective conductance to heat exchange surface
  gHa<-.gturb(uf,d,zm,microsnow$maxhgt)
  Test<-microsnow$snow$snowtempG
  sel<-which(Test<microsnow$tdew)
  Test[sel]<-microsnow$tdew[sel]
  gHa[gHa>4]<-4
  # Save outputs
  microsnow$uf<-uf
  microsnow$uz<-uz
  microsnow$gHa<-gHa
  microsnow$Test<-Test
  microsnow$d<-d
  microsnow$zm<-zm
  microsnow$progress<-3
  return(microsnow)
}
#' Snow Penman-Monteith
.snowPenMont <- function(tc,pk,ea,radabs,gHa,G,es=NA,tdew=NA,T_est=NA,allout=TRUE) {
  if (is.na(es[1])) es<-.satvap(tc)
  if (is.na(tdew[1])) tdew<-.dewpoint(ea,tc)
  if (class(T_est)=="logical") T_est<-tc
  Tav<-(tc+T_est)/2
  delta<-.delta(Tav)
  sb<-5.67*10^-8
  gHr<-gHa+(4*0.97*sb*(Tav+273.15)^3)/29.3
  Rem<-0.97*sb*(tc+273.15)^4
  m<-44526*(gHa/pk)
  L<-m*(es-ea)
  dTmx<- -0.6273*max(tc,na.rm=T)+49.79
  dT<-(radabs-Rem-L-G)/(29.3*gHr+m*delta)
  # Limit
  dT[dT>dTmx]<-dTmx
  tcan<-tc+dT
  tcan<-.lim(tcan,tdew)
  tcan[tcan>72]<-72
  H<-29.3*gHa*(tcan-tc)
  HR<-H/(radabs-Rem)
  HR[HR<0]<-0
  HR[HR>1]<-1
  if (allout) {
    estl<-.satvap(tcan)
    L<-m*(estl-ea)
    L<-.lim(L,0)
    out<-list(tcan=tcan,H=H,L=L,G=G,HR=HR)
  } else out<-list(tcan=tcan,H=H,HR=HR)
  return(out)
}
#' @title Above ground snow microclimate model
#' @description Used by e.g. [runmicro_hr()] to model microclimate above ground
#' @param microsnow an object of class microsnowin as returned by [snowmodelin()]
#' @param reqhgt specified height (m) above or below ground for which microclimate
#' outputs are required. Negative below ground
#' @param pai_a total plant area index above `reqhgt`. Estimated if not supplied (see details).
#' @param xyf optional input for called function [snowwind()]
#' @param zf optional input for called function [wind()]
#' @param slr an optional SpatRaster object of slope values (Radians). Calculated from
#' dtm if not supplied.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' dtm if not supplied.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from dtm if not supplied.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from dtm if not supplied.
#' @param gmn optional minimum convective conductance value (mol/m^2/s). See [snowwind()])
#' @return see [runmicro_hr()]
#' @details DETAILS
#' @rdname snowtemphumE
#' @export
snowtemphumE<-function(microsnow, reqhgt = 0.05, pai_a = NA, xyf = 1, zf = NA, slr = NA,
                       apr = NA, hor = NA, wsa = NA, gmn = 0.3) {
  # run soiltemp function if not run
  if (microsnow$progress<3) {
    microsnow<-.snowwind(microsnow,reqhgt,pai_a,xyf,zf,slr,apr,hor,wsa,gmn)
  }
  # Calculate vapour conductivity
  #climdata<-micro$climdata
  # Calculate T0
  canabs<-microsnow$radCsw+microsnow$radClw
  # Calculate Latent heat flux
  TH<-.snowPenMont(microsnow$tc,microsnow$pk,microsnow$ea,canabs,microsnow$gHa,microsnow$G,
                  microsnow$estl,microsnow$tdew,microsnow$Test)
  # Calculate temperature and vapour pressure above canopy, setting height to canopy top of reghgt < hgt
  z<-microsnow$veghgt
  sel<-which(reqhgt>=z) # above canopy
  z[sel]<-reqhgt
  Tzv<-.TVabove(TH,microsnow,z)
  microsnow$Tz<-Tzv$Tz
  ez<-Tzv$ez
  # Computes below canopy temperatures
  # NB Automatically replaces with above canopy if needed
  Tzb<-.LangrangianSimT(reqhgt,microsnow,TH)
  microsnow$Tz<-Tzb$To
  # Compute leaf temperature
  microsnow<-.leaftemp(microsnow,reqhgt,TH$tcan,Tzb$leafabs)
  # Compute relative humidity
  rh<-.LangrangianSimV(reqhgt,microsnow,ez)
  # Return values needed
  soilm<-.rta(rast(microsnow$Smax),dim(rh)[3])
  microsnow<-list(Tz=microsnow$Tz,tleaf=microsnow$tleaf,T0=microsnow$T0,soilm=soilm,
                  relhum=rh,windspeed=microsnow$uz,Rdirdown=microsnow$Rbdown,
                  Rdifdown=microsnow$Rddown,Rlwdown=Tzb$Rldown,
                  Rswup=microsnow$Rdup,Rlwup=Tzb$Rlwup)
  class(microsnow)<-"microout"
  return(microsnow)
}
#' @title Below ground microclimate model
#' @description Used by e.g. [runmicro_hr()] to model microclimate below ground
#' @param microsnow an object of class microsnowin as returned by [snowmodelin()]
#' @param reqhgt specified height (m) below ground for which microclimate
#' outputs are required. Must be negative.
#' @param pai_a total plant area index above `reqhgt`. Estimated if not supplied (see details).
#' @param slr an optional SpatRaster object of slope values (Radians). Calculated from
#' dtm if not supplied.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' dtm if not supplied.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from dtm if not supplied.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from dtm if not supplied.
#' @return an array of soil temperatures (deg C)
#' @rdname snowbelow_hr
#' @export
#' @importFrom stats filter
snowbelow_hr<-function(microsnow, reqhgt, pai_a = NA, slr = NA, apr = NA, hor = NA, wsa = NA) {
  xx <- function(x,y) as.numeric(stats::filter(c(x,x),rep(1/y,y), sides = 2))
  if (is.null(microsnow$T0[1])) {
    microsnow<-.snowsoiltemp_hr(microsnow,reqhgt,pai_a,slr,apr,hor,wsa)
  }
  hiy<-length(microsnow$tme)
  Tav<-apply(microsnow$T0,c(1,2),mean)
  Tan<-microsnow$T0-.rta(rast(Tav),hiy)
  hr_sm<-(2*pi*abs(reqhgt)^2)/(7200*mean(microsnow$ka,na.rm=TRUE)) # calcaulate effective periodicity of soil fluctuations
  hr_ps<-round((24*abs(reqhgt))/(2*pi*mean(microsnow$DD,na.rm=TRUE)),0)
  if (hr_sm < hiy) {
    Tz<-aperm(apply(Tan,c(1,2),xx,hr_sm),c(2,3,1))
    Tz<-pmax(Tz[,,1:hiy],Tz[,,(hiy+1):dim(Tz)[3]],na.rm=T)
    if (hr_ps > 0) Tz<-abind::abind(Tz[,,(dim(Tz)[3]-(hr_ps-1)):dim(Tz)[3]],Tz)
    Tz<-Tz[,,1:hiy]+.rta(rast(Tav),hiy)
  } else Tz<-.rta(rast(Tav),hiy)
  # Return values needed
  return(Tz)
}
#' @title Run snow microclimate model (hourly)
#' @description  Run snow microclimate model in hourly time increments
#' @param microsnow an object of class microsnowin as returned by [snowmodelin()]
#' @param reqhgt specified height (m) above or below ground for which microclimate
#' outputs are required. Negative below ground
#' @param pai_a total plant area index above `reqhgt`. Estimated if not supplied (see details).
#' @param xyf optional input for called function [snowwind()]
#' @param zf optional input for called function [wind()]
#' @param slr an optional SpatRaster object of slope values (Radians). Calculated from
#' dtm if not supplied.
#' @param apr an optional SpatRaster object of aspect values (Radians). Calculated from
#' dtm if not supplied.
#' @param hor an optional array of the tangent of the angle to the horizon in
#' 24 directions. Calculated from dtm if not supplied.
#' @param wsa an optional array of wind shelter coefficients in 8 directions.
#' Calculated from dtm if not supplied.
#' @param gmn optional minimum convective conductance value (mol/m^2/s). See [snowwind()])
#' @return an object of class microsnowout with the following components:
#' @return `Tz` Array of air temperatures at height `reqhgt` (deg C). Identical to `T0`
#' if `reqhgt` = 0.
#' @return `tleaf` Array of leaf temperatures at height `reqhgt` (deg C).
#' NA if `reqhgt` greater than canopy height or `reqhgt` <= 0.
#' @return `T0` Array of ground surface temperatures (deg C)
#' @return `relhum` Array of relative humidities at height `reqhgt` (percentage).
#' NA if `reqhgt` <= 0.
#' @return `windspeed` Array of wind speeds at height `reqhgt` (m/s).
#' NA if `reqhgt` <= 0.
#' @return `Rdirdown` Array of downward direct shortwave radiation incident on
#' horizontal surface (W/m^2)
#' @return `Rdifdown` Array of downward diffuse shortwave radiation incident on
#' horizontal surface (W/m^2)
#' @return `Rlwdown` Array of downward longwave radiation incident on horizontal
#' surface (W/m^2)
#' @return `Rswup` Array of upward shortwave radiation (assumed diffuse) incident
#' on underside of horizontal surface (W/m^2)
#' @return `Rlwup` Array of upward longwave radiation incident on underside of
#' horizontal surface (W/m^2)
#' @details `pai_a` is used to calculate the radiation intercepted by leaves at `reqhgt` if
#' below canopy. If not supplied it is calculated from total plant area index by
#' assuming a realistic shape to the vertical profile foliage within the canopy. Wind speed and
#' radiation values are only returned when `reqhgt > 0`. To derive radiation values when
#' `reqhgt = 0`, set `pai_a` to `microsnow$pai`. If supplied, `pai_a` must have the
#' same dimensions as `microsnow$pai`. I.e. with the same x and y dims as the the
#' supplied dtm and values for each hour as the z dimension.
#' @examples
#' \dontrun{
#' if(interactive()){
#'  require(NicheMapR)
#'  require(terra)
#'  require(abind)
#'  # Derive estimates of snow melt temperature parameters using NicheMapR
#'  nmrout<-runNMRSnow(climdata,precd,66.215,29.293,ALTT = 360)
#'  STparams<-fitsnowtemp(climdata,precd,nmrout$SNOWDEP,nmrout$SNOWTEMP)
#'  meltfact<-fitsnowdepth(climdata,nmrout$SNOWDEP,precd,STparams=STparams)$meltfact
#'  # Run snow model in spatial mode (takes ~3 mins)
#'  snowdepth<-nmrout$SNOWDEP
#'  snow<-modelsnowdepth(climdata,precd,snowdepth,dtm,pai,hgt,STparams,meltfact,spatial=T)
#'  # Prepare microclimate model inputs
#'  microsnow<-snowmodelin(climdata,precd,snow,STparams,dtm,pai,hgt,x=1,clump=0.2,soiltype=7,windhgt=2)
#'  mout<-runmicrosnow_hr(microsnow, reqhgt = 0.05)
#'  plot(rast(mout$Tz[,,500])) # Air temperature, 500th hour
#'  }
#' }
#' @rdname runmicrosnow_hr
#' @export
runmicrosnow_hr <- function(microsnow, reqhgt, pai_a = NA, xyf = 1, zf = NA, slr = NA,
                            apr = NA, hor = NA, wsa = NA, gmn = 0.3) {
  # Calculate soil surface temperature and soil moisture
  microsnow<-snowsoiltemp_hr(microsnow,reqhgt,pai_a,slr,apr,hor,wsa)
  # Run above ground
  if (reqhgt > 0) {
    mout<-snowtemphumE(microsnow,reqhgt,pai_a,xyf,zf,slr,apr,hor,wsa,gmn)
  }
  # Run at ground level
  if (reqhgt == 0) {
    n<-length(microsnow$climdata$temp)
    soilm<-.rta(rast(microsnow$Smax),n)
    mout<-list(Tz=microsnow$T0,tleaf=NA,T0=microsnow$T0,soilm=soilm,
               relhum=NA,windspeed=NA,Rdirdown=NA,Rdifdown=NA,Rlwdown=NA,
               Rswup=NA,Rlwup=NA)
    class(mout)<-"microout"
  }
  # Run below ground
  if (reqhgt < 0) {
    Tz<-below_hr(microsnow,reqhgt,pai_a,slr,apr,hor,wsa)
    n<-length(microsnow$climdata$temp)
    soilm<-.rta(rast(microsnow$Smax),n)
    mout<-list(Tz=Tz,tleaf=NA,T0=microsnow$T0,soilm=soilm,
               relhum=NA,windspeed=NA,Rdirdown=NA,Rdifdown=NA,Rlwdown=NA,
               Rswup=NA,Rlwup=NA)
    class(mout)<-"microsnowout"
  }
  return(mout)
}
