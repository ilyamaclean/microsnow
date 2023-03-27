#' Check if input is a raster and convert to matrix if it is
#' @import raster
.is <- function(r) {
  as.matrix(r,wide=T)
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
#' @import terra
.latlongfromrast<-function (r) {
  e <- ext(r)
  xy <- data.frame(x=(e$xmin+e$xmax)/2,y=(e$ymin+e$ymax)/2)
  xy <- sf::st_as_sf(xy, coords = c("x", "y"),
                     crs = crs(r))
  ll <- sf::st_transform(xy, 4326)
  ll <- data.frame(lat = sf::st_coordinates(ll)[2], long = sf::st_coordinates(ll)[1])
  return(ll)
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
.solarindex<- function(dtm,salt,azi,slr=NA,apr=NA) {
  if (class(slr)[1]=="logical") slr<-terrain(dtm,v="slope")
  if (class(apr)[1]=="logical") apr<-terrain(dtm,v="aspect")
  sl<-.rta(slr*pi/180,length(azi))
  ap<-.rta(apr*pi/180,length(azi))
  sl[is.na(sl)]<-0
  ap[is.na(ap)]<-0
  sazi<-.vta(azi,dtm)*(pi/180)
  i<-sin(salt)*cos(sl)+cos(salt)*sin(sl)*cos(sazi-ap)
  i[i<0]<-0
  i
}
#' Calculates radiation extinction coefficient for canopy
.cank <- function(x,sa,si) {
  # Raw k
  sa[sa<0]<-0
  zen<-pi/2-sa
  k<-sqrt((x^2+(tan(zen)^2)))/(x+1.774*(x+1.182)^(-0.733))
  k0<-sqrt(x^2)/(x+1.774*(x+1.182)^(-0.733))
  # k dash
  kd<-k*cos(zen)/si
  sel<-which(si==0)
  kd[sel]<-1
  return(list(k=k,kd=kd,k0=k0))
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
.gturb<-function(uf,d,zm,z1,z0=NA,psi_h=0,gmin=0.05) {
  zh<-0.2*zm
  if (class(z0)=="logical") z0<-d+zh
  g<-(0.4*43*uf)/log((z1-d)/(z0-d))
  sel<-which(g<gmin)
  g[g<gmin]<-gmin
  g
}
#' Calculates turbulent molar conductivity within canopy
.gcanopy <- function(uf,h,d,z1,z0=0.0122,phi_h=1,gmin=0.05) {
  sel<-which(z0<0.0122)
  z0[sel]<-0.0122
  g<-(0.4*43*uf*(h-d)*phi_h)/log(z1/z0)
  g[g<gmin]<-gmin
  g
}
#' Calculate zero plane displacement
.zeroplanedis<-function(h,pai) {
  pai[pai<0.001]<-0.001
  d<-(1-(1-exp(-sqrt(7.5*pai)))/sqrt(7.5*pai))*h
  d
}
#' Calculate roughness length
.roughlength<-function(h,pai,d=NA,psi_h=0) {
  if (class(d)=="logical") d<-.zeroplanedis(h,pai)
  Be<-sqrt(0.003+(0.2*pai)/2)
  zm<-(h-d)*exp(-0.4/Be)*exp(psi_h)
  zm
}
#' Calculate parameters of two-stream radiation model
.twostreamparams<-function(pait,x,clump,lref,gref,ltra,alt,si) {
  si[si<0]<-0
  alt[alt<0]<-0
  si[alt<=0]<-0
  # === (1a) Calculate canopy k
  kkd<-.cank(x,alt*pi/180,si)
  # === (1c) Adjust paramaters for gap fraction and inclined surface
  pai_t<-pait/(1-clump)
  sk<-(1-(1-gref)*si)/(1-(1-gref)*sin(alt*pi/180))
  gref2<-gref*sk
  # === (1d) Calculate two stream base parameters
  om<-lref+ltra
  a<-1-om
  del<-lref-ltra
  mla<-(9.65*(3+x)^(-1.65))
  mla[mla>pi/2]<-pi/2
  J<-cos(mla)^2
  gma<-0.5*(om+J*del)
  s<-0.5*(om+J*del/kkd$k)*kkd$k
  sstr<-om*kkd$k-s
  # === (1e) Calculate two stream base parameters
  h<-sqrt(a^2+2*a*gma)
  sig<-kkd$kd^2+gma^2-(a+gma)^2
  S1<-exp(-h*pai_t)
  S2<-exp(-kkd$kd*pai_t)
  u1<-a+gma*(1-1/gref)
  u2<-a+gma*(1-gref)
  D1<-(a+gma+h)*(u1-h)*1/S1-(a+gma-h)*(u1+h)*S1
  D2<-(u2+h)*1/S1-(u2-h)*S1
  # === (1f) Calculate Diffuse radiation parameters
  p1<-(gma/(D1*S1))*(u1-h)
  p2<-(-gma*S1/D1)*(u1+h)
  p3<-(1/(D2*S1))*(u2+h)
  p4<-(-S1/D2)*(u2-h)
  # === (1f) Calculate Direct radiation parameters
  u1<-a+gma*(1-1/gref2)
  u2<-a+gma*(1-gref2)
  D1<-(a+gma+h)*(u1-h)*1/S1-(a+gma-h)*(u1+h)*S1
  D2<-(u2+h)*1/S1-(u2-h)*S1
  p5<- -s*(a+gma-kkd$kd)-gma*sstr
  v1<-s-(p5*(a+gma+kkd$kd))/sig
  v2<-s-gma-(p5/sig)*(u1+kkd$kd)
  p6<-(1/D1)*((v1/S1)*(u1-h)-(a+gma-h)*S2*v2)
  p7<-(-1/D1)*((v1*S1)*(u1+h)-(a+gma+h)*S2*v2)
  p8<-sstr*(a+gma+kkd$kd)-gma*s
  v3<-(sstr+gma*gref2-(p8/sig)*(u2-kkd$kd))*S2
  p9<-(-1/D2)*((p8/(sig*S1))*(u2+h)+v3)
  p10<-(1/D2)*(((p8*S1)/sig)*(u2-h)+v3)
  # Tests
  return(list(p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6,p7=p7,p8=p8,p9=p9,p10=p10,
              h=h,sig=sig,gref2=gref2,kkd=kkd))
}
