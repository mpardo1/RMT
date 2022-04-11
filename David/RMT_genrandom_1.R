####### RANDOM MATRICES FOR METAPOPULATION MODELS ##

library(tidyverse)
library(ggplot2)
library(ggforce)
library(copula)
library(tmvtnorm)

####### GENERATE MATRICES ##########################

### Index
#rand_mat: simple random matrix, with given mu and sigma and choice of distribution
#rand_mat_ell: as above, allowing also for elliptic correlated entries
#rand_mat_comm: matrix for the commuting model (ie multiplied already by the betas
#and subtracted the gammas, which are allowed to be random). no correlation
#rand_mat_cor_norm: random matrix with general correlations as in Baron et al; Gaussian
#rand_mat_cor_norm_N: as above with different rescaling in N (experiment)
#rand_mat_cor_beta: random matrix with general correlations as in Baron et al; beta dist

#compute parameters for a beta distribution
#so that mean and sigma are as given
alphabeta <- function(mu,sig) {
  alphabeta <- ((1-mu)/sig^2 -1/mu)*mu^2
}
betabeta <- function(mu,sig) {
  betabeta <- alphabeta(mu,sig)*(1/mu-1)
}

#compute parameters for a gamma distribution
#so that mean and sigma are as given
alphagamma <- function(mu,sig) {
  alphagamma <- (mu/sig)^2
}
betagamma <- function(mu,sig) {
  betagamma <- mu/(sig^2)
}

# mu1 <- .6
# mu2 <- .7
# sig1 <- sig2 <- .2
# rho <- .25

bivbeta <- function(N,mu1,sig1,mu2,sig2,rho) {
  
  cors <- vector()
  maxalpha <- max(c((1-mu1)*(mu1*(1-mu1)/(sig1^2)-1),(1-mu2)*(mu2*(1-mu2)/(sig2^2)-1)))
  for (alpha5 in seq(0.01,maxalpha,by = 0.01)) {
  
  alpha1 <- mu1*(mu1*(1-mu1)/(sig1^2)-1) - (1-mu2)*(mu2*(1-mu2)/(sig2^2)-1) + alpha5
  alpha2 <- mu2*(mu2*(1-mu2)/(sig2^2)-1) - (1-mu1)*(mu1*(1-mu1)/(sig1^2)-1) + alpha5
  alpha3 <- (1-mu2)*(mu2*(1-mu2)/(sig2^2)-1) - alpha5
  alpha4 <- (1-mu1)*(mu1*(1-mu1)/(sig1^2)-1) - alpha5
  
  U1 <- rgamma(N, shape = alpha1, scale = 1)
  U2 <- rgamma(N, shape = alpha2, scale = 1)
  U3 <- rgamma(N, shape = alpha3, scale = 1)
  U4 <- rgamma(N, shape = alpha4, scale = 1)
  U5 <- rgamma(N, shape = alpha5, scale = 1)
  
  X <- (U1+U3)/(U1+U3+U4+U5)
  Y <- (U2+U4)/(U2+U3+U4+U5)
  
  cors[as.character(alpha5)] <- cor(X,Y)
  
  }

  alpha5 <- 0.01 + 0.01*which.min(abs(cors - rho))
  alpha1 <- mu1*(mu1*(1-mu1)/(sig1^2)-1) - (1-mu2)*(mu2*(1-mu2)/(sig2^2)-1) + alpha5
  alpha2 <- mu2*(mu2*(1-mu2)/(sig2^2)-1) - (1-mu1)*(mu1*(1-mu1)/(sig1^2)-1) + alpha5
  alpha3 <- (1-mu2)*(mu2*(1-mu2)/(sig2^2)-1) - alpha5
  alpha4 <- (1-mu1)*(mu1*(1-mu1)/(sig1^2)-1) - alpha5

  U1 <- rgamma(N, shape = alpha1, scale = 1)
  U2 <- rgamma(N, shape = alpha2, scale = 1)
  U3 <- rgamma(N, shape = alpha3, scale = 1)
  U4 <- rgamma(N, shape = alpha4, scale = 1)
  U5 <- rgamma(N, shape = alpha5, scale = 1)
  
  X <- (U1+U3)/(U1+U3+U4+U5)
  Y <- (U2+U4)/(U2+U3+U4+U5)
  
  return(cbind(X,Y))
  
}
  
# generate a random matrix with entries coming from the specified distribution
# distrib can be set to gamma or lognormal; more distributions can be included
# as long as they depend on two parameters

#simple random matrix
rand_mat <- function(N,mu,sig,distrib = "beta"){
  
  muln <- log(mu^2/sqrt(mu^2 + sig^2))
  sdln <- sqrt(log(1+sig^2/mu^2))
  
  rmatrix <- dplyr::case_when(
    distrib == "gamma" ~ matrix(rgamma(N^2,shape = (mu/sig)^2,rate = mu/(sig^2)), nrow = N),
    distrib == "lnorm" ~ matrix(rlnorm(N^2,meanlog = muln,sdlog = sdln), nrow = N),
    distrib == "beta" ~  matrix(rbeta(N^2,alphabeta(mu,sig),betabeta(mu,sig)), nrow = N),
    TRUE ~ matrix(rep(0, N^2), nrow = N)
  )
  
  rmatrix <- matrix(rmatrix, nrow = N)
  diag(rmatrix) <- rep(0,N)
  
  return(rmatrix)
}

#random matrix with elliptic correlated entries
rand_mat_ell <- function(N,mu,sig,rho,distrib = "beta"){
  
  param1 <- dplyr::case_when(
    distrib == "gamma" ~ (mu/sig)^2,
    distrib == "lnorm" ~ log(mu^2/sqrt(mu^2 + sig^2)),
    distrib == "beta" ~ alphabeta(mu,sig),
    TRUE ~ 0
  )
  
  param2 <- dplyr::case_when(
    distrib == "gamma" ~ mu/(sig^2),
    distrib == "lnorm" ~ sqrt(log(1+sig^2/mu^2)),
    distrib == "beta" ~ betabeta(mu,sig),
    TRUE ~ 0
  )
  
  cop <- normalCopula(param=rho, dim = 2, dispstr = "un")
  bivdist <- mvdc(copula = cop, margins=c(distrib,distrib),
                  paramMargins=list(list(param1, param2),
                                    list(param1, param2)))
  bivdistvalues <- rMvdc(N^2, bivdist)
  
  rmatrix <- matrix(0,N,N)
  for(i in c(1:N)){
    for(j in c(i:N)){
      ind = j + (i-1)*N
      rmatrix[i,j] = bivdistvalues[ind,1]
      rmatrix[j,i] = bivdistvalues[ind,2]
    }
  }
  
  diag(rmatrix) <- rep(0,N)
  return(rmatrix)
}

#simple random matrix for the commuting model
#random betas and gammas are allowed
#output is the whole jacobian! already
#coupled with the betas
rand_mat_comm <- function(N,mu,sig,mb,sb,mg,sg,distrib){
  
  muln <- log(mu^2/sqrt(mu^2 + sig^2))
  sdln <- sqrt(log(1+sig^2/mu^2))
  
  rmatrix <- dplyr::case_when(
    distrib == "gamma" ~ matrix(rgamma(N^2,shape = (mu/sig)^2,rate = mu/(sig^2)), nrow = N),
    distrib == "lnorm" ~ matrix(rlnorm(N^2,meanlog = muln,sdlog = sdln), nrow = N),
    distrib == "beta" ~  matrix(rbeta(N^2,alphabeta(mu,sig),betabeta(mu,sig)), nrow = N),
    TRUE ~ matrix(rep(0, N^2), nrow = N)
  )
  rmatrix <- matrix(rmatrix, nrow = N)
  
  diag(rmatrix) <- rep(1,N)
  
  betas <- matrix(rep(0,N^2), N)
  diag(betas) <- rnorm(N, mb, sb)
  
  rmatrix <- rmatrix %*% betas
  diag(rmatrix) <- diag(rmatrix) - rnorm(N,mg,sg)
  
  return(rmatrix)
}

#random matrix with general correlated entries
#let us see what happens with gaussian noise
#so that the mean is isolated from the noise
rand_mat_cor_norm <- function(N,mu,sig,rho,G,r,c){

  copfc <- normalCopula(param=G/sqrt(r*c), dim = 2, dispstr = "un")
  bivfc <- mvdc(copula = copfc, margins=c("norm","norm"),
                paramMargins=list(list(0,sqrt((r/N))*sig/N),
                                  list(0,sqrt((c/N))*sig/N)))
  
  valfc <- rMvdc(N, bivfc)
  
  fils <- valfc[,1] #betas of baron et al
  cols <- valfc[,2] #kappas of baron et al
  
  copmob <- normalCopula(param= (N*rho-2*G)/(N-(r+c)), dim = 2, dispstr = "un")
  bivmob <- mvdc(copula = copmob, margins=c("norm","norm"),
                 paramMargins=list(list(0, sqrt(sig^2*(N-r-c)/(N^3))),
                                   list(0, sqrt(sig^2*(N-r-c)/(N^3)))))

  valmob <- rMvdc(N*(N-1)/2, bivmob)
  
  rmatrix <- matrix(0,N,N)
  ind <- 1
  for(i in c(1:(N-1))){
    for(j in c((i+1):N)){
      rmatrix[i,j] = mu/N + valmob[ind,1] + fils[i] + cols[j]
      rmatrix[j,i] = mu/N + valmob[ind,2] + fils[j] + cols[i]
      ind <- ind+1
    }
  }
  diag(rmatrix) <- mu/N - 1
  return(rmatrix)
}

#### cor check
# prueba <- vector()
# for (j in c(1:100)){
#   valmob <- rMvdc(N*(N-1)/2, bivmob)
#   prueba[j]<- cor(valmob[,1],valmob[,2])
# }
# mean(prueba)
####

#random matrix with general correlated entries
#and re-scaled in N so that all the cors
#are in the same scale (N^0)
rand_mat_cor_norm_N <- function(N,mu,sig,rho,G,r,c){
  
  copfc <- normalCopula(param=G/sqrt(r*c), dim = 2, dispstr = "un")
  bivfc <- mvdc(copula = copfc, margins=c("norm","norm"),
                paramMargins=list(list(0,sig*r),
                                  list(0,sig*c)))
  valfc <- rMvdc(N, bivfc)
  
  fils <- valfc[,1] #betas of baron et al
  cols <- valfc[,2] #kappas of baron et al
  
  copmob <- normalCopula(param= (rho-2*G)/(1-(r+c)), dim = 2, dispstr = "un")
  bivmob <- mvdc(copula = copmob, margins=c("norm","norm"),
                 paramMargins=list(list(0, sig*(1-r-c)),
                                   list(0, sig*(1-r-c))))
  valmob <- rMvdc(N*(N-1)/2, bivmob)
  
  rmatrix <- matrix(0,N,N)
  ind <- 1
  for(i in c(1:(N-1))){
    for(j in c((i+1):N)){
      rmatrix[i,j] = mu + valmob[ind,1] + fils[i] + cols[j]
      rmatrix[j,i] = mu + valmob[ind,2] + fils[j] + cols[i]
      ind <- ind+1
    }
  }
  diag(rmatrix) <- mu - 1
  return(rmatrix)
}

#random matrix with general correlated entries
#now with noise coming from beta distribution
#recall mu\in(0,1), sig^2\in(0,mu*(1-mu))
rand_mat_cor_beta <- function(N,mu,sig,rho,G,r,c){
  
  copfc <- normalCopula(param=G/sqrt(r*c), dim = 2, dispstr = "un")
  bivfc <- mvdc(copula = copfc, margins=c("beta","beta"),
                paramMargins=list(list(alphabeta(mu/(3*N),sqrt(r)*sig/N),betabeta(mu/(3*N),sqrt(r)*sig/N)),
                                  list(alphabeta(mu/(3*N),sqrt(c)*sig/N),betabeta(mu/(3*N),sqrt(c)*sig/N))))
  valfc <- rMvdc(N, bivfc)
  
  fils <- valfc[,1] #betas of baron et al
  cols <- valfc[,2] #kappas of baron et al
  
  copmob <- normalCopula(param= (rho-2*G/N)/(1-(r+c)/N), dim = 2, dispstr = "un")
  bivmob <- mvdc(copula = copmob, margins=c("beta","beta"),
                 paramMargins=list(list(alphabeta(mu/(3*N), sig*sqrt(N-r-c)/N),betabeta(mu/(3*N), sig*sqrt(N-r-c)/N)),
                                   list(alphabeta(mu/(3*N), sig*sqrt(N-r-c)/N),betabeta(mu/(3*N), sig*sqrt(N-r-c)/N))))
  valmob <- rMvdc(N*(N-1)/2, bivmob)
  
  # prueba <- vector()
  # for (j in c(1:100)){
  #   valmob <- rMvdc(N*(N-1)/2, bivmob)
  #   prueba[j]<- cor(valmob[,1],valmob[,2])
  # }
  # mean(prueba)

  rmatrix <- matrix(0,N,N)
  ind <- 1
  for(i in c(1:(N-1))){
    for(j in c((i+1):N)){
      rmatrix[i,j] = valmob[ind,1] + fils[i] + cols[j]
      rmatrix[j,i] = valmob[ind,2] + fils[j] + cols[i]
      ind <- ind+1
    }
  }
  diag(rmatrix) <- rep(0,N)
  
  paramdf <- data.frame(wantedN = c(mu/N,sig/N,rho,G/N,r/N,c/N),
                        simulated = c(mean(fils)+mean(cols)+mean(valmob),
                                      sqrt(var(fils)+var(cols)+var(valmob[,1])),
                                      (2*cov(fils,cols)+cov(valmob[,1],valmob[,2]))/(var(fils)+var(cols)+var(valmob[,1])),
                                      cov(fils,cols)/(var(fils)+var(cols)+var(valmob[,1])),
                                      var(fils)/(var(fils)+var(cols)+var(valmob[,1])),
                                      var(cols)/(var(fils)+var(cols)+var(valmob[,1])))) %>%
    mutate(simulated = round(simulated,6))
  rownames(paramdf) <- c("mu","sigma","rho","Gamma","r","c")
  
  return(list(rmatrix,paramdf))
}

#### SAVE ANTES DE ENREDAR CON BARON GALLA

rand_mat_cor_beta <- function(N,mu,sig,rho,G,r,c){

  copfc <- normalCopula(param=G/sqrt(r*c), dim = 2, dispstr = "un")
  bivfc <- mvdc(copula = copfc, margins=c("beta","beta"),
                paramMargins=list(list(alphabeta(mu/(3*N),sqrt((r/N))*sig/N),betabeta(mu/(3*N),sqrt((r/N))*sig/N)),
                                  list(alphabeta(mu/(3*N),sqrt((c/N))*sig/N),betabeta(mu/(3*N),sqrt((c/N))*sig/N))))
  valfc <- rMvdc(N, bivfc)

  fils <- valfc[,1] #betas of baron et al
  cols <- valfc[,2] #kappas of baron et al

  copmob <- normalCopula(param= (N*rho-2*G)/(N-(r+c)), dim = 2, dispstr = "un")
  bivmob <- mvdc(copula = copmob, margins=c("beta","beta"),
                 paramMargins=list(list(alphabeta(mu/(3*N), sqrt(sig^2*(N-r-c)/(N^3))),betabeta(mu/(3*N), sqrt(sig^2*(N-r-c)/(N^3)))),
                                   list(alphabeta(mu/(3*N), sqrt(sig^2*(N-r-c)/(N^3))),betabeta(mu/(3*N), sqrt(sig^2*(N-r-c)/(N^3))))))
  valmob <- rMvdc(N*(N-1)/2, bivmob)

  # prueba <- vector()
  # for (j in c(1:100)){
  #   valmob <- rMvdc(N*(N-1)/2, bivmob)
  #   prueba[j]<- cor(valmob[,1],valmob[,2])
  # }
  # mean(prueba)

  rmatrix <- matrix(0,N,N)
  ind <- 1
  for(i in c(1:(N-1))){
    for(j in c((i+1):N)){
      rmatrix[i,j] = valmob[ind,1] + fils[i] + cols[j]
      rmatrix[j,i] = valmob[ind,2] + fils[j] + cols[i]
      ind <- ind+1
    }
  }
  diag(rmatrix) <- rep(0,N)

  paramdf <- data.frame(wantedN = c(mu/N,sig/N,rho,G/N,r/N,c/N),
                        simulated = c(mean(fils)+mean(cols)+mean(valmob),
                                      sqrt(var(fils)+var(cols)+var(valmob[,1])),
                                      (2*cov(fils,cols)+cov(valmob[,1],valmob[,2]))/(var(fils)+var(cols)+var(valmob[,1])),
                                      cov(fils,cols)/(var(fils)+var(cols)+var(valmob[,1])),
                                      var(fils)/(var(fils)+var(cols)+var(valmob[,1])),
                                      var(cols)/(var(fils)+var(cols)+var(valmob[,1])))) %>%
    mutate(simulated = round(simulated,6))
  rownames(paramdf) <- c("mu","sigma","rho","Gamma","r","c")

  return(list(rmatrix,paramdf))
}
####

rand_mat_cor_beta2 <- function(N,mu,sig,rho,G,r,c){
  
  copfc <- normalCopula(param=G/sqrt(r*c), dim = 2, dispstr = "un")
  bivfc <- mvdc(copula = copfc, margins=c("beta","beta"),
                paramMargins=list(list(alphabeta(mu/(3*N),sqrt((r/N))*sig/N),betabeta(mu/(3*N),sqrt((r/N))*sig/N)),
                                  list(alphabeta(mu/(3*N),sqrt((c/N))*sig/N),betabeta(mu/(3*N),sqrt((c/N))*sig/N))))
  valfc <- rMvdc(N, bivfc)
  
  fils <- valfc[,1] #betas of baron et al
  cols <- valfc[,2] #kappas of baron et al
  
  copmob <- normalCopula(param= (N*rho-2*G)/(N-(r+c)), dim = 2, dispstr = "un")
  bivmob <- mvdc(copula = copmob, margins=c("beta","beta"),
                 paramMargins=list(list(alphabeta(mu/(3*N), sqrt(sig^2*(N-r-c)/(N^3))),betabeta(mu/(3*N), sqrt(sig^2*(N-r-c)/(N^3)))),
                                   list(alphabeta(mu/(3*N), sqrt(sig^2*(N-r-c)/(N^3))),betabeta(mu/(3*N), sqrt(sig^2*(N-r-c)/(N^3))))))
  valmob <- rMvdc(N*(N-1)/2, bivmob)
  
  # prueba <- vector()
  # for (j in c(1:100)){
  #   valmob <- rMvdc(N*(N-1)/2, bivmob)
  #   prueba[j]<- cor(valmob[,1],valmob[,2])
  # }
  # mean(prueba)
  
  rmatrix <- matrix(0,N,N)
  ind <- 1
  for(i in c(1:(N-1))){
    for(j in c((i+1):N)){
      rmatrix[i,j] = valmob[ind,1] + fils[i] + cols[j]
      rmatrix[j,i] = valmob[ind,2] + fils[j] + cols[i]
      ind <- ind+1
    }
  }
  diag(rmatrix) <- mu/N - 1
  
  return(rmatrix)
}

#random matrix with general correlated entries
#and without copulas; using instead multivariate
#truncated normal distributions
rand_mat_cor_trun <- function(N,mu,sig,rho,G,r,c){
  
  covfc <- (sqrt((r/N))*sig/N)*(sqrt((c/N))*sig/N)*(G/sqrt(r*c))
  valfc <- rtmvnorm(N, mean = rep(mu/(3*N),2),
                    sigma = matrix(c(sqrt((r/N))*sig/N, covfc, covfc, sqrt((c/N))*sig/N), nrow = 2),
                    lower = c(0,0), upper = c(1,1))
  
  fils <- valfc[,1] #betas of baron et al
  cols <- valfc[,2] #kappas of baron et al
  
  covmob <- (sqrt(sig^2*(N-r-c)/(N^3)))*(sqrt(sig^2*(N-r-c)/(N^3)))/(N-(r+c))
  valmob <- rtmvnorm(N*(N-1)/2, mean = rep(mu/(3*N),2),
                     sigma = matrix(c(sqrt(sig^2*(N-r-c)/(N^3)), covmob, covmob, sqrt(sig^2*(N-r-c)/(N^3))), nrow = 2),
                     lower = c(0,0), upper = c(1,1))
  
  rmatrix <- matrix(0,N,N)
  ind <- 1
  for(i in c(1:(N-1))){
    for(j in c((i+1):N)){
      rmatrix[i,j] = valmob[ind,1] + fils[i] + cols[j]
      rmatrix[j,i] = valmob[ind,2] + fils[j] + cols[i]
      ind <- ind+1
    }
  }
  diag(rmatrix) <- rep(0,N)
  return(rmatrix)
}


####### PLOT EIGENVALUES ###########################

# save eigenvalues of matrix
eigen_mat <- function(mat){
  eigen_m <- as.complex(eigen(mat, only.values = TRUE)$values)
  df <- data.frame(re = Re(eigen_m), im = Im(eigen_m))
  return(df)
}

# plot eigenvalues
plot_eigen <- function(mymat){
  eigmat <- eigen_mat(mymat)
  plot_eigen <- ggplot(eigmat) + geom_point(aes(re,im), size = 0.05) 
}

# plot eigenvalues and predictions
plot_eigen_rmt <- function(jacobian,
                           N,mub,mug,
                           muw,sw,rhow,Gammaw,
                           muc,sc,rhoc,Gammac,
                           tau) {
  
  sigma <- sqrt(mub^2*sw^2 + 2*mub*tau + sc^2)
  rho <- (mub^2*rhow*sw^2 + rhoc*sc^2)/(sigma^2)
  Gamma <- (mub^2*Gammaw*(sw^2)/N + Gammac*(sc^2)/N)/(sigma^2)
  
  #center <- mub*(1-muw) - mug - (N-1)*(muc*(N^2-2*N+1)/(N^2-N)+mucstar*(N-1)/(N^2-N))-muc
  
  center <- mub*(1-muw) - mug - N*muc
  radius <- sigma*sqrt(N)
  outlier <- ifelse(((Gammaw == 0) & (Gammac == 0)),
                    mub - mug + mub*muw*(N-1),
                    mub - mug + mub*muw*(N-1) +
                      (1/2)*(N*mub*muw+N*muc)*(1+rho/Gamma)*(sqrt(1+(4*Gamma*sigma^2)/((mub*muw + muc)^2))-1))
  
  #outlier <- -1 + muw + (muw/2)*(1+rho/G)*(sqrt(1+(4*G*(sw/sqrt(N))^2)/((muw)^2))-1) # Baron et al
  #outlier <- -1 + N*muw + (N*muw/2)*(1+rho/G)*(sqrt(1+(4*G*sw^2)/(N*muw^2))-1)# direct rescaling
  #outlier <- -1 + N*muw + (N*muw/2)*(1+rho/(G*N))*(sqrt(1+(4*N*G*sw^2)/(N*muw^2))-1)# retuned rescaling
  #outlier <- mu*N*(3/2+rho/(2*G*N))*(sqrt(1+4*G*N*sw^2/muw^2)-1)
  
  eigmat <- eigen_mat(jacobian)
  
  plot_eigen_rmt <- ggplot(eigmat) + 
    coord_fixed() +
    geom_vline(xintercept = 0, color = "blue") +
    geom_point(aes(x = outlier, y = 0), color = "red") +
    geom_ellipse(aes(x0 = center, y0 = 0, a = (1+rho)*radius, b = (1-rho)*radius, angle = 0), color = "red") +
    geom_point(aes(re,im), size = 0.05) +
    labs(x="re", y="im", subtitle = paste0("Real outlier = ",max(eigmat$re),". RMT prediction = ",mub - mug + mub*muw*(N-1),
                                           ". Baron-Galla prediction = ",mub - mug + mub*muw*(N-1) +
                                             (1/2)*(N*mub*muw+N*muc)*(1+rho/Gamma)*(sqrt(1+(4*Gamma*sigma^2)/((mub*muw + muc)^2))-1))) +
    theme_bw()
  
  print(plot_eigen_rmt)
}


