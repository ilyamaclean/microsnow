#' Create SpatRaster object using a template
#' @import terra
.rast <- function(m,tem) {
  r<-rast(m)
  ext(r)<-ext(tem)
  crs(r)<-crs(tem)
  r
}
#' Corrects wind profile
.windcorrect <- function(uz, windhgt, maxhgt) {
  d<-0.08
  zm<-0.01476
  uf<-(0.4*uz)/log((windhgt-d)/zm)
  uo<-(uf/0.4)*log((maxhgt+2-d)/zm)
  uo
}
#' Calculates wind shelter coefficient in specified direction
.windcoef <- function(dsm, direction, hgt = 1, res = 1) {
  reso<-res(dsm)[1]
  dsm<-.is(dsm)
  dsm[is.na(dsm)]<-0
  dtm<-dsm/reso
  hgt<-hgt/reso
  azi<-direction*(pi/180)
  horizon <- array(0, dim(dtm))
  dtm3 <- array(0, dim(dtm) + 200)
  x <- dim(dtm)[1]
  y <- dim(dtm)[2]
  dtm3[101:(x + 100), 101:(y + 100)] <- dtm
  for (step in 1:10) {
    horizon[1:x,1:y]<-pmax(horizon[1:x,1:y],(dtm3[(101-cos(azi)*step^2):(x+100-cos(azi)*step^2),
                                                  (101+sin(azi)*step^2):(y+100+sin(azi)*step^2)]-dtm3[101:(x+100),101:(y+100)])/(step^2))
    horizon[1:x,1:y]<-ifelse(horizon[1:x,1:y]<(hgt/step^2),0,horizon[1:x,1:y])
  }
  index<-1-atan(0.17*100*horizon)/1.65
  index
}
#' Calculates array of wind shelter coefficients in wdct directions
.windsheltera<-function(dtm, whgt, s) {
  if (is.na(s)) {
    s<-min(dim(dtm)[1:2])
    s<-ifelse(s>10,10,s)
  }
  dct2<-16
  drs<-c(0:(dct2-1))*360/dct2
  a<-array(NA,dim=c(dim(dtm)[1:2],dct2))
  for (i in 1:dct2) {
    r<-.rast(.windcoef(dtm,drs[i],whgt,res(dtm)[1]),dtm)
    r<-resample(aggregate(r,fact=s,fun="mean"),dtm)
    a[,,i]<-as.matrix(r,wide=T)
  }
  a2<-array(NA,dim=c(dim(dtm)[1:2],8))
  for (i in 1:8) {
    if (i == 1) {
      m<-0.5*a[,,1]+0.25*a[,,2]+0.25*a[,,dct2]
    } else m<-0.5*a[,,(i*2-1)]+0.25*a[,,i*2]+0.25*a[,,(i*2)-2]
    a2[,,i]<-m
  }
  a2
}
#' Calculate shelter coefficient for all wind directions
.windshelter <- function(wdir,dtm,whgt,s,wsa = NA) {
  # Create array of wind shelter coefficients
  if (class(wsa)[1] == "logical") wsa<-.windsheltera(dtm,whgt,s)
  # Apply array of wind shelter coefficients
  i<-360/dim(wsa)[3]
  sel<-(round(wdir/i,0))+1
  sel[sel==(dim(wsa)[3])+1]<-1
  sel[sel==0]<-dim(wsa)[3]
  wa<-wsa[,,sel]
  wa
}
#' Calculate saturated vapour pressure
.satvap <- function(tc) {
  e0<-(tc<0)*610.78/1000+(tc>=0)*611.2/1000
  L <- (tc<0)*2.834*10^6+(tc>=0)*((2.501*10^6)-(2340*tc))
  T0<-(tc<0)*273.15+(tc>=0)*273.15
  estl<-e0*exp((L/461.5)*(1/T0-1/(tc+273.15)))
  estl
}
#' Calculate dewpoint
.dewpoint <- function(ea, tc) {
  e0<-611.2/1000
  L<-(2.501*10^6)-(2340*tc)
  it<-1/273.15-(461.5/L)*log(ea/e0)
  Tdew<-1/it-273.15
  e0<-610.78/1000
  L<-2.834*10^6
  it<-1/273.15-(461.5/L)*log(ea/e0)
  Tfrost<-1/it-273.15
  sel<-which(Tdew<0)
  Tdew[sel]<-Tfrost[sel]
  Tdew
}
#' Derive soil parameters form soil type
.soilinit <- function(soiltype) {
  u<-unique(as.vector(.is(soiltype)))
  u<-u[is.na(u)==F]
  rho<-array(NA,dim=dim(soiltype)[1:2])
  Vm<-rho; Vq<-rho; Mc<-rho;
  Smax<-rho; Smax<-rho
  psi_e<-rho; soilb<-rho
  for (i in 1:length(u)) {
    sel<-which(microclimf::soilparameters$Number==u[i])
    sop<-microclimf::soilparameters[sel,]
    sel<-which(.is(soiltype)==u[i])
    rho[sel]<-sop$rho
    Vm[sel]<-sop$Vm
    Vq[sel]<-sop$Vq
    Mc[sel]<-sop$Mc
    psi_e[sel]<-sop$psi_e
    soilb[sel]<-sop$b
    Smax[sel]<-sop$Smax
  }
  return(list(rho=rho,Vm=Vm,Vq=Vq,Mc=Mc,psi_e=psi_e,soilb=soilb,Smax=Smax))
}
#' Calculates diurnal temperature fluctuation
.A0f<-function(tc) {
  tc<-matrix(tc,ncol=24,byrow=T)
  mn<-apply(tc,1,min)
  mx<-apply(tc,1,max)
  A0<-(mx-mn)/2
  rep(A0,each=24)
}
#' Calculates reference time for phase of diurnal temperature fluctuation
.t0f<-function(tc) {
  tc<-matrix(tc,ncol=24,byrow=T)
  tmx<-apply(tc,1,which.max)-1
  t0<-(tmx-6)%%24
  rep(t0*3600,each=24)
}
#' Calculates reference time for phase of diurnal temperature fluctuation (daily)
.t0fd<-function(tc) {
  tc<-matrix(tc,ncol=24,byrow=T)
  tmx<-apply(tc,1,which.max)-1
  t0<-(tmx-6)%%24
  t0<-t0*3600
  t0
}
#' Calculates slope of the saturated vapour pressure curve
.delta <- function(tc) {
  es1<-.satvap(tc-0.5)
  es2<-.satvap(tc+0.5)
  delta<-es2-es1
  delta
}
#' Cap values
.lim<-function(x,l,up=FALSE) {
  if (length(l)==1) l<-x*0+l
  if (up) {
    s<-which(x>l)
    x[s]<-l[s]
  } else {
    s<-which(x<l)
    x[s]<-l[s]
  }
  x
}
#' Calculates temperature and vapour above canopy
.TVabove<-function(TH,microsnow,z) {
  # Temperature
  hgt<-microsnow$veghgt
  hgt[hgt<0.001]<-0.001
  d<-.zeroplanedis(hgt,microsnow$pai)
  zm<-.roughlength(hgt,microsnow$pai,d)
  lnr<-log((z-d)/zm)/log((microsnow$maxhgt-d)/zm)
  dH<-TH$tcan-microsnow$tc
  dZ<-dH*(1-lnr)
  Tz<-dZ+microsnow$tc
  # Vapour pressure
  dH<-.satvap(TH$tcan)-microsnow$ea
  dZ<-dH*(1-lnr)
  ez<-microsnow$ea+dZ
  return(list(Tz=Tz,ez=ez))
}
#' Calculates source concentration
.SourceD<-function(reqhgt,swabs,tcan,tground,tref,skyem,pai,pai_a,clump,Hratio,folden) {
  paie<-pai/(1-clump)
  pai_ae<-pai_a/(1-clump)
  # Lw down
  sb<-5.67*10^-8
  trd<-exp(-pai_ae)+clump
  lwdown<-trd*skyem*sb*(tref+273.15)^4+
    (1-trd)*0.97*sb*(tcan+273.15)^4
  # Lw up
  tru<-exp(-(paie-pai_ae))+clump
  lwup<-tru*sb*0.97*(tground+273.15)^4+
    (1-tru)*0.97*sb*(tcan+273.15)^4
  lwabs<-0.97*0.5*(lwdown+lwup)
  # H
  Rem<-0.97*sb*(tcan+273.15)^4
  H<-Hratio*(swabs+lwabs-Rem)
  S<-folden*H
  return(list(S=S,radabs=swabs+lwabs,Rldown=lwdown,Rlwup=lwup))
}
#' Computes near field concentration
.Cnear<-function(S,uf,hgt) {
  p1<-sqrt(uf)*(hgt^(2/3))
  slp<-0.42458*p1^0.85443
  Tn<-S/(slp*1000)
  Tn
}
#' Performs Langrangian simulation (temperature)
.LangrangianSimT<-function(reqhgt,microsnow,TH) {
  # Calculate d and zh of ground-layer
  ghgt<-microsnow$snow$gsnowdepth/100
  #ghgt[is.na(ghgt)]<-0
  d<-0.075*microsnow$veghgt+ghgt
  zh<-0.001506*microsnow$veghgt
  lnr<-suppressWarnings(log((reqhgt-d)/zh)/log((microsnow$veghgt-d)/zh))
  sel<-which(reqhgt<(d+zh))
  # Calculate log-linear gradient
  Tzfo<-microsnow$T0-lnr*(microsnow$T0-microsnow$Tz)
  # Calculate linear gradient
  Tzfl<-microsnow$T0-(reqhgt/microsnow$veghgt)*
        (microsnow$T0-microsnow$Tz)
  Tzfo[sel]<-Tzfl[sel]
  # Partition according to canopy cover
  tr<-exp(-microsnow$pai)
  Tf<-tr*Tzfo+(1-tr)*Tzfl
  Tmx<-pmax(microsnow$T0,microsnow$Tz)
  Tmn<-pmin(microsnow$T0,microsnow$Tz)
  s1<-which(Tf>Tmx)
  s2<-which(Tf<Tmn)
  Tf[s1]<-Tmx[s1]
  Tf[s2]<-Tmn[s2]
  # Compute Source Density
  skyem<-.vta(microsnow$climdata$skyem,microsnow$dtm)
  SR<-.SourceD(reqhgt,microsnow$radLs,TH$tcan,microsnow$T0,microsnow$tc,skyem,
               microsnow$pai,microsnow$pai_a,microsnow$clump,TH$HR,microsnow$leafden)
  S<-SR$S
  leafabs<-SR$radabs
  # Compute near field temperature
  Tn<-.Cnear(S,microsnow$uf,microsnow$veghgt)
  To<-Tf+Tn
  # Replace with above canopy temperature if reqhgt > canopy height
  sel<-which(microsnow$veghgt<=reqhgt)
  To[sel]<-microsnow$Tz[sel]
  # Replace with snow temperature if reqhgt <= snowdepth
  sel<-which(reqhgt<=ghgt)
  To[sel]<-microsnow$snow$snowtempG[sel]
  # Set limits
  dT<-To-microsnow$tc
  dTmx<- -0.6273*max(microsnow$tc,na.rm=T)+49.79
  dT[dT>dTmx]<-dTmx
  To<-microsnow$tc+dT
  To[To>72]<-72
  To<-.lim(To,microsnow$tdew)
  return(list(To=To,leafabs=leafabs,Rldown=SR$Rldown,Rlwup=SR$Rlwup))
}
# Calculate free convection (used for calculating minimum conductivity)
.gfree<-function(leafd,H) {
  d<-0.71*leafd
  dT<-0.7045388*(d*H^4)^(1/5)
  gha<-0.0375*(dT/d)^0.25
  gha[gha<0.1]<-0.1
  gha
}
# Calculates leaf temperature
.leaftemp<-function(microsnow,reqhgt,tcan,leafabs) {
  Rem<-0.97*5.67*10^-8*(microsnow$tc+273.15^4)
  # Calculate conductivity
  n<-dim(Rem)[3]
  ld<-0.71*microsnow$leafd
  gh<-0.023256*sqrt(microsnow$uz/ld)
  He<-0.7*(leafabs-Rem)
  gmn<-.gfree(ld/0.71,abs(He))
  sel<-which(gh<gmn)
  gh[sel]<-gmn[sel]
  TH<-snowPenMont(microsnow$Tz,microsnow$pk,microsnow$ea,leafabs,gh,
              0,NA,microsnow$tdew,1,T_est=tcan)
  microsnow$tleaf<-TH$tcan
  sel<-which(reqhgt>microsnow$veghgt)
  microsnow$tleaf[sel]<-tcan[sel]
  microsnow$gh<-gh
  return(microsnow)
}
.LangrangianSimV<-function(reqhgt,microsnow,ez) {
  # Calculate soil surface effective vapour pressure
  n<-dim(ez)[3]
  e0<-.satvap(microsnow$T0)
  # Calculate d and zh of ground-layer
  ghgt<-microsnow$snow$gsnowdepth/100
  d<-0.075*microsnow$veghgt+ghgt
  zh<-0.001506*microsnow$veghgt
  # Calculate far field vapour pressure
  lnr<-suppressWarnings(log((reqhgt-d)/zh)/log((microsnow$veghgt-d)/zh))
  sel<-which(reqhgt<(d+zh))
  # Calculate log-linear gradient
  efo<-e0-lnr*(e0-ez)
  # Calculate linear gradient
  efl<-e0-(reqhgt/microsnow$veghgt)*(e0-ez)
  efo[sel]<-efl[sel]
  tr<-exp(-microsnow$pai)
  ef<-tr*efo+(1-tr)*efl
  # Compute near field vapour pressure
  es<-.satvap(microsnow$tleaf)
  L<-44526*microsnow$gh/microsnow$pk*(es-microsnow$ea)
  S<-L*microsnow$leafden
  Cn<-.Cnear(S,microsnow$uf,microsnow$veghgt)*29.3*43
  Cn[Cn>2500]<-2500
  # Vapour pressure
  eo<-ef+microsnow$pk/(44526*43)*Cn
  sel<-which(microsnow$veghgt<=reqhgt)
  eo[sel]<-ez[sel]
  eas<-.satvap(microsnow$Tz)
  # Relative humidity
  rh<-(eo/eas)*100
  # Set to 100% if below snow
  sel<-which(reqhgt<=ghgt)
  rh[sel]<-100
  rh[rh>100]<-100
  rh[rh<30]<-30
  return(rh)
}
#' convert climate array to data.frame of weather
.catoweather<-function(climarray,tme) {
  # ~~~~ Wind
  uwind<-climarray$windspeed*cos(climarray$winddir*(pi/180))
  vwind<-climarray$windspeed*sin(climarray$winddir*(pi/180))
  uw<-apply(uwind,3,mean,na.rm=T)
  vw<-apply(vwind,3,mean,na.rm=T)
  ws=sqrt(uw^2+vw^2)
  wd=(atan2(vw,uw)*(180/pi))%%360
  weather<-data.frame(obs_time=tme,
                      temp=apply(climarray$temp,3,mean,na.rm=T),
                      relhum=apply(climarray$relhum,3,mean,na.rm=T),
                      pres=apply(climarray$pres,3,mean,na.rm=T),
                      swrad=apply(climarray$swrad,3,mean,na.rm=T),
                      difrad=apply(climarray$difrad,3,mean,na.rm=T),
                      skyem=apply(climarray$skyem,3,mean,na.rm=T),
                      windspeed=ws,
                      winddir=wd)
  return(weather)
}
#' Resample array
#' @import terra
.resa<-function(varin,r,dtm) {
  if (dim(varin)[1]!=dim(dtm)[1] | dim(varin)[2]!=dim(dtm)[2]) {
    b<-rast(varin)
    ext(b)<-ext(r)
    crs(b) <- crs(r)
    if (as.character(crs(b)) != as.character(crs(dtm))) {
      b <- project(b,dtm, method = 'near')
    }
    b<-resample(b,dtm)
    a<-as.array(b)
    a<-round(a,3)
  } else a<-varin
  a
}
#' Latitudes from SpatRaster object
.latsfromr <- function(r) {
  e <- ext(r)
  lts <- rep(seq(e$ymax - res(r)[2] / 2, e$ymin + res(r)[2] / 2, length.out = dim(r)[1]), dim(r)[2])
  lts <- array(lts, dim = dim(r)[1:2])
  lts
}
#' Longitudes from SpatRaster object
.lonsfromr <- function(r) {
  e <- ext(r)
  lns <- rep(seq(e$xmin + res(r)[1] / 2, e$xmax - res(r)[1] / 2, length.out = dim(r)[2]), dim(r)[1])
  lns <- lns[order(lns)]
  lns <- array(lns, dim = dim(r)[1:2])
  lns
}
#' Lats and longs from SpatRaster object, including reprojection
.latslonsfromr <- function(r) {
  lats<-.latsfromr(r)
  lons<-.lonsfromr(r)
  xy<-data.frame(x=as.vector(lons),y=as.vector(lats))
  xy <- sf::st_as_sf(xy, coords = c('x', 'y'), crs = crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(lat = sf::st_coordinates(ll)[,2],
                   long = sf::st_coordinates(ll)[,1])
  lons<-array(ll$long,dim=dim(lons))
  lats<-array(ll$lat,dim=dim(lats))
  return(list(lats=lats,lons=lons))
}
#' Lape rates
.lapserate <- function(tc, rh, pk) {
  e0 <- 0.6108 * exp(17.27 * tc / (tc + 237.3))
  ws <- 0.622 * e0 / pk
  ea <- e0 * (rh / 100)
  rv <- 0.622 * ea / (pk - ea)
  lr <- 9.8076 * (1 + (2501000 * rv) / (287 * (tc + 273.15))) / (1003.5 +
                                                                   (0.622 * 2501000 ^ 2 * rv) / (287 * (tc + 273.15) ^ 2))
  lr
}
