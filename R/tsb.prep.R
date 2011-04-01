##########################
# summarize sample for tsb
##########################

tsb.prep <- function(y,N,M.0,psu,psu.size,g,B){
  # design stuff
  P <- ncol(y) # number of observed variables
  n <- length(unique(psu)) # number of sampled PSUs
  f.1 <- n/N # primary sampling fraction
  m <- tapply(psu,psu,length) # secondary sample sizes
  m.partial.sum <- 0; i <-2 # how many observations before this PSU's observations
                            # used to tell C how far down sample vector to go
  while(i <= n) {m.partial.sum <- c(m.partial.sum,sum(m[1:(i-1)])); i <- i+1}
  m.sum <- sum(m) # total sample size
  M <- tapply(psu.size,psu,mean) # PSU sizes
  M.bar <- sum(M)/N # average PSU size
  lambda.1 <- sqrt(n/(n-1) * (1 - f.1)) # weird function from Rao & Wu

  # variable stuff
  y.totals <- apply(y,2,function(x)tapply(x,psu,sum)*M/m) # matrix nxP
  y.bar <- N/n*apply(y.totals,2,sum)/M.0 # vector of length P with HT estimated pop means
  y.in <- y[,1]; i <- 2; while(i <= P) {y.in <- c(y.in, y[,i]);i <- i+1} # concatenated columns of y

  out <- list(P,n,f.1,m,m.partial.sum,m.sum,M,M.bar,lambda.1,y.totals,y.bar,y.in)
  names(out) <- c('P','n','f.1','m','m.partial.sum','m.sum','M','M.bar','lambda.1',
                  'y.totals','y.bar','y.in')
  return(out)
}
