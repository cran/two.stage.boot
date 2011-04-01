###############################################################
# compare two.stage.boot to a version implemented entirely in R
###############################################################

# write checking function
tsb.check <- function(y,       # dataframe whose columns are arguments to gp
                     N,       # number of PSUs in population
                     M.0,     # number of SSUs in entire population
                     psu,     # PSU index
                     psu.size,# the number of SSUs in a given PSU
                     g,       # non-linear function of interest
                     B,       # bootstrap size
                     conf.level=0.95){      
  
  # check if data has proper nested structure
  # check if psu.size is constant within a PSU

  ##############
  # sample stats
  ##############
  # design stuff
  P <- ncol(y) # number of observed variables
  n <- length(unique(psu)) # number of sampled PSUs
  f.1 <- n/N # primary sampling fraction
  m <- tapply(psu,psu,length) # secondary sample sizes
  m.partial.sum <- 0; i <-2 # how many observations before this PSU's observations
                            # used to tell C how far down sample vector to go
  while(i <= n) {m.partial.sum <- c(m.partial.sum,sum(m[1:(i-1)])); i <- i+1}
  m.sum <- sum(m) # total sample size
  star.y.space <- max(m) # size of the largest secondary sample
                         # used to allocate memory for resampling in C
  M <- tapply(psu.size,psu,mean) # PSU sizes
  M.bar <- sum(M)/N # average PSU size
  lam.1 <- sqrt(n/(n-1) * (1 - f.1)) # weird function from Rao & Wu

  # variable stuff
  y.totals <- apply(y,2,function(x)tapply(x,psu,sum)*M/m) # matrix nxP
  y.bar <- N/n*apply(y.totals,2,sum)/M.0 # vector of length P with HT estimated pop means
  y.in <- y[,1]; i <- 2; while(i <= P) {y.in <- c(y.in, y[,i]);i <- i+1} # concatenated columns of y

  ###################
  # boot calculations
  ###################
  # compute resampled pop mean estimates
  tilde.bar <- rep(0,P*B)
  check <- rep(0,n)
  for(b in 1:B){
    for(i in 1:n){
      #pick PSU and evaluate primary sampling stuff
      u <- runif(1)
      star_psu <- ceiling(u*n)
      check[i] <- star_psu
      star.m <- m[star_psu]
      star.m.partial.sum <- m.partial.sum[star_psu]
      star.M <- M[star_psu]
      lam.2 <- sqrt(f.1*(1-star.m/star.M)*star.m/(star.m-1))
      
      #pick SSUs and evaluate secondary sampling stuff
      v <- rep(0,star.m)
      j <- 1
      while(j <= star.m){
        u <- runif(1)
        v[j] <- ceiling(u*star.m)
        j <- j+1
      }

      for(p in 1:P){
        star.y <- y.in[(p-1)*m.sum + star.m.partial.sum + v]
        syt <- sum(star.y)*star.M/star.m
        tilde.y <- y.bar[p]+lam.1*(syt/M.bar-y.bar[p])+lam.2*(star.M*star.y/M.bar-syt/M.bar)
        tilde.bar[(b-1)*P+p] <- tilde.bar[(b-1)*P+p] + sum(tilde.y)/n/star.m
      }
    }
  }
    
  tilde.mat <- matrix(tilde.bar,,ncol=P,byrow=T) # put estimates in a nice matrix
  g.boot <- apply(tilde.mat,1,g) # compute bootstrap estimate of g

  ################
  # return results
  ################
  out <- list(g.boot,
              mean(g.boot),
              sd(g.boot),
              c(quantile(g.boot,(1-conf.level)/2),quantile(g.boot,1-(1-conf.level)/2)),
              check)
  names(out) <- c('boot sample','boot mean','boot SE','boot CI','check')
  return(out)
}

library(two.stage.boot)
data(tiny.eg)

set.seed(1)
tsb.out <- tsb(y=tiny.eg$y,N=tiny.eg$N,M.0=tiny.eg$M.0,psu=tiny.eg$psu,
               psu.size=tiny.eg$psu.size,g=tiny.eg$g,B=tiny.eg$B)

set.seed(1)
check.out <- tsb.check(y=tiny.eg$y,N=tiny.eg$N,M.0=tiny.eg$M.0,psu=tiny.eg$psu,
                     psu.size=tiny.eg$psu.size,g=tiny.eg$g,B=tiny.eg$B)

all.equal(tsb.out[[1]],check.out[[1]])

