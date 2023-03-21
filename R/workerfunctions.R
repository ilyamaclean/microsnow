#' Check if input is a raster and convert to matrix if it is
#' @import raster
.is <- function(r) {
  getValues(r,format="matrix")
}
#' Convert matrix to array
.rta <- function(r,n) {
  m<-.is(r)
  a<-array(rep(m,n),dim=c(dim(r)[1:2],n))
  a
}
#' Convert vector to array
.vta <- function(v,r) {
  m<-.is(r)
  va<-rep(v,each=dim(m)[1]*dim(m)[2])
  a<-array(va,dim=c(dim(m),length(v)))
  a
}
#' Calculates latitude and longitude from raster
#' @import raster
.latlongfromraster<-function (r) {
  e <- extent(r)
  xy <- data.frame(x = (e@xmin + e@xmax)/2, y = (e@ymin + e@ymax)/2)
  coordinates(xy) = ~x + y
  proj4string(xy) = crs(r)
  ll <- as.data.frame(spTransform(xy, CRS("+init=epsg:4326")))
  ll <- data.frame(lat = ll$y, long = ll$x)
  ll
}
#' Calculate zero plane displacement
.zeroplanedis <- function(hgt, pai) {
  d<-(1-(1-exp(-sqrt(7.5*pai)))/sqrt(7.5*pai))*hgt
  d
}
#' Calculate roughness length governing momentum transfer
.roughlength <- function(hgt, pai, zm0 = 0.003) {
  d<-.zeroplanedis(hgt,pai)
  ur<-sqrt(zm0+(0.3*pai)/2)
  ur[ur>0.3]<-0.3
  zm<-(hgt-d)*exp(-0.4*ur-0.193)
  zm[zm<zm0]<-zm0
  zm
}
#' Calculates below canopy mixing length
.mixinglength <- function(hgt, pai, zm0 = 0.003, mnlm = 0.04) {
  d<-.zeroplanedis(hgt,pai)
  zm<-.roughlength(hgt,pai,zm0)
  l_m<-(0.32*(hgt-d))/log((hgt-d)/zm)
  l_m[l_m<mnlm]<-mnlm
  l_m
}
#' expand daily array to hourly array
.ehr<-function(a) {
  n<-dim(a)[1]*dim(a)[2]
  o1<-rep(c(1:n),24*dim(a)[3])
  o2<-rep(c(1:dim(a)[3]),each=24*n)-1
  o2<-o2*max(o1)
  o<-o1+o2
  ah<-rep(a,24)
  ah<-ah[o]
  ah<-array(ah,dim=c(dim(a)[1:2],dim(a)[3]*24))
  ah
}
#' Calculates tangent of horizon angle
.horizon <- function(dtm, azimuth) {
  reso<-res(dtm)[1]
  dtm<-.is(dtm)
  dtm[is.na(dtm)]<-0
  dtm<-dtm/reso
  azi<-azimuth*(pi/180)
  horizon<-array(0,dim(dtm))
  dtm3<-array(0,dim(dtm)+200)
  x<-dim(dtm)[1]
  y<-dim(dtm)[2]
  dtm3[101:(x+100),101:(y+100)]<-dtm
  for (step in 1:10) {
    horizon[1:x,1:y]<-pmax(horizon[1:x,1:y],(dtm3[(101-cos(azi)*step^2):(x+100-cos(azi)*step^2),
                                                  (101+sin(azi)*step^2):(y+100+sin(azi)*step^2)]-dtm3[101:(x+100),101:(y+100)])/(step^2))
  }
  horizon
}
#' Calculates the astronomical Julian day
.jday <- function(tme) {
  yr<-tme$year+1900
  mth<-tme$mon+1
  dd<-tme$mday+(tme$hour+(tme$min+tme$sec/60)/60)/24
  madj<-mth+(mth<3)*12
  yadj<-yr+(mth<3)*-1
  jd<-trunc(365.25*(yadj+4716))+trunc(30.6001*(madj+1))+dd-1524.5
  b<-(2-trunc(yadj/100)+trunc(trunc(yadj/100)/4))
  jd<-jd+(jd>2299160)*b
  jd
}
#' Calculates solar time
.soltime <- function(localtime, long, jd, merid = 0, dst = 0) {
  m<-6.24004077+0.01720197*(jd-2451545)
  eot<- -7.659*sin(m)+9.863*sin(2*m+3.5932)
  st<-localtime+(4*(long-merid)+eot)/60-dst
  st
}
#' Calculates the solar altitude
.solalt <- function(localtime, lat, long, jd, merid = 0, dst = 0) {
  st<-.soltime(localtime,long,jd,merid,dst)
  tt<-0.261799*(st-12)
  d<-(pi*23.5/180)*cos(2*pi*((jd-171)/365.25))
  sh<-sin(d)*sin(lat*pi/180)+cos(d)*cos(lat*pi/180)*cos(tt)
  sa<-(180*atan(sh/sqrt(1-sh^2)))/pi
  sa
}
#' Calculates the solar azimuth
.solazi <- function(localtime, lat, long, jd, merid = 0, dst = 0) {
  st<-.soltime(localtime,long,jd,merid,dst)
  tt<-0.261799*(st-12)
  d<-(pi*23.5/180)*cos(2*pi*((jd-171)/365.25))
  sh<-sin(d)*sin(lat*pi/180)+cos(d)*cos(lat*pi/180)*cos(tt)
  hh<-(atan(sh/sqrt(1-sh^2)))
  sazi<-cos(d)*sin(tt)/cos(hh)
  cazi<-(sin(lat*pi/180)*cos(d)*cos(tt)-cos(pi*lat/180)*sin(d))/
    sqrt((cos(d)*sin(tt))^2+(sin(pi*lat/180)*cos(d)*cos(tt)-cos(pi*lat/180)*sin(d))^2)
  sqt <- 1-sazi^2
  sqt[sqt<0]<-0
  solz<-180+(180*atan(sazi/sqrt(sqt)))/pi
  solz[cazi<0 & sazi<0]<-180-solz[cazi<0 & sazi<0]
  solz[cazi<0 & sazi>=0]<-540-solz[cazi<0 & sazi>=0]
  solz
}
#' Calculates the solar coefficient
.solarindex<- function(dtm,alt,azi,slr=NA,apr=NA) {
  if (class(slr)[1]=="logical") slr<-terrain(dtm,opt="slope")
  if (class(apr)[1]=="logical") apr<-terrain(dtm,opt="aspect")
  sl<-.rta(slr,length(azi))
  ap<-.rta(apr,length(azi))
  zen<-pi/2-alt
  sazi<-.vta(azi,dtm)*(pi/180)
  i<-cos(zen)*cos(sl)+sin(zen)*sin(sl)*cos(sazi-ap)
  i[i<0]<-0
  i
}
#' Calculates radiation extinction coefficient for canopy
.cank <- function(x,sa) {
  zen<-pi/2-sa
  k<-sqrt((x^2+(tan(zen)^2)))/(x+1.774*(x+1.182)^(-0.733))
  k
}
#' Calculates direct radiation transmission through canopy
.cantransdir <- function(l, k, ref = 0.23, clump = 0) {
  f<-1/(1-clump)
  s<-sqrt(1-ref)
  ks<-k*s
  tr<-exp(-ks*l*f)
  tr<-(1-clump)*tr+clump
  tr
}
#' Calculates diffuse radiation transmission through canopy
.cantransdif <- function(l, ref = 0.23, clump = 0) {
  f<-1/(1-clump)
  s<-sqrt(1-ref)
  tr<-exp(-s*l*f)
  tr<-(1-clump)*tr+clump
  tr
}
#' Calculates turbulent molar conductivity above canopy
.gturb<-function(uf,d,zm,z1,z0=NA,psi_h=0,gmin) {
  zh<-0.2*zm
  if (is.na(z0)[1]) {
    z0<-d+zh
  }
  xx<-(z1-d)/(z0-d)
  xx[xx<0.001]<-0.001
  ln<-log(xx)
  lnr<-ln*log((z1-d)/zh)
  psx<-lnr*psi_h
  g<-(0.4*43*uf)/(ln+psx)
  sel<-which(g<gmin)
  # DK: Originally was g[sel]<-gmin[sel]. But gmin is just a single constant
  # (here, 0(. So I don't think we're supposed to select inside gmin
  g[sel]<-gmin
  sel<-which(g>1e10)
  g[sel]<-1e10
  g
}
#' Calculates turbulent molar conductivity within canopy
.gcanopy <- function(l_m,a,hgt,uh,z1,z0,gmin) {
  e0<-exp(-a*(z0/hgt-1))
  e1<-exp(-a*(z1/hgt-1))
  g<-(l_m*21.5*uh*a)/(e0-e1)
  # Set minimum
  sel<-which(g<gmin)
  g[sel]<-gmin[sel]
  g[is.na(g)]<-1/mean(1/gmin)
  sel<-which(g>1e10)
  g[sel]<-1e10
  g
}
