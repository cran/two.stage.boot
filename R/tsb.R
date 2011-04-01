##############
# tsb function
##############

tsb <- function(y,       # dataframe whose columns are arguments to gp
                N,       # number of PSUs in population
                M.0,     # number of SSUs in entire population
                psu,     # PSU index
                psu.size,# the number of SSUs in a given PSU
                g,       # non-linear function of interest
                B,       # bootstrap size
                conf.level=0.95){
  # check to make sure that arguments are the right thing
  tsb.check(y,N,M.0,psu,psu.size,g,B,conf.level)
  
  # sample stats
  prep <- tsb.prep(y,N,M.0,psu,psu.size,g,B)


  # boot calculations
  out.C <- .C('work',
              Bin = as.integer(B),
              Pin = as.integer(prep$P),
              nin = as.integer(prep$n),
              f_1in = as.double(prep$f.1),
              m = as.integer(prep$m),
              m_partial_sum = as.integer(prep$m.partial.sum),
              m_sumin = as.integer(prep$m.sum),
              m_max = as.integer(max(prep$m)),
              M = as.integer(prep$M),
              M_bar_in = as.double(prep$M.bar),
              lam1_in = as.double(prep$lambda.1),
              y_bar = as.double(prep$y.bar),
              y = as.double(prep$y.in),
              tilde_bar = as.double(rep(0,prep$P*B)))
  
  tilde.mat <- matrix(out.C$tilde_bar,,ncol=prep$P,byrow=T) # put estimates in a nice matrix
  g.boot <- apply(tilde.mat,1,g) # compute bootstrap estimate of g

  # return results
  out <- list(g.boot,
              mean(g.boot),
              sd(g.boot),
              c(quantile(g.boot,(1-conf.level)/2),quantile(g.boot,1-(1-conf.level)/2)))
  names(out) <- c('values','mean','SE','CI')
  return(out)
}

  
