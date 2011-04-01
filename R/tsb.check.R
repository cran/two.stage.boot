#####################
# check inputs of tsb
#####################

tsb.check <- function(y,N,M.0,psu,psu.size,g,B,conf.level){
  if(!class(y)=='data.frame') stop('y must be of class data.frame',call.=FALSE)
  if(!(class(N)=='numeric'|class(N)=='integer')) stop('N must be a number')
  if(!(round(N)==N&length(N)==1)) stop('N must be a single integer')
  if(!(class(M.0)=='numeric'|class(M.0)=='integer')) stop('M.0 must be a number')
  if(!(round(M.0)==M.0&length(M.0)==1)) stop('M.0 must be a single integer')
  if(!length(psu)==nrow(y)) stop('length(psu) must equal nrow(y)')
  if(!length(psu.size)==nrow(y)) stop('length(psu.size) must equal nrow(y)')  
  if(!(class(psu.size)=='numeric'|class(psu.size)=='integer')) stop('psu.size must be numbers')
  if(sum(!round(psu.size)==psu.size)>0) stop('psu.size must be integers')
  if(!class(g)=='function') stop('g must be of class function')
  if(!(class(B)=='numeric'|class(B)=='integer')) stop('B must be a number')
  if(!(round(B)==B&length(B)==1)) stop('B must be a single integer')
  if(!(class(conf.level)=='numeric')) stop('conf.level must be a number')
  if(!length(conf.level)==1) stop('conf.level must be a single integer')
  if(!(0<conf.level&conf.level<1)) stop('conf.level must be within (0,1)')
  if(!sum(tapply(psu.size,psu,var)>0)==0) stop('psu.size must be constant within a PSU')
}
