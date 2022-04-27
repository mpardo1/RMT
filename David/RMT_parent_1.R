####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### parent script
###
### generate, plot and integrate metapopulation
### epidemiological models
### 
#########################################################

setwd("~/RMT/David/")

source("~/RMT/David/RMT_genrandom_1.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")

####### GENERATE JACOBIAN ###############################

# number of patches
N <- 300

# epidemiological
#all rates must lie in (0,1) except for betas

Deltas <- rep(0.1, N) # birth rate
mub <- 0.1
sb <- 0.001
betas <- rep(mub, N) # transmission rates
#betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.1, N) # loss of immunity rates
mud <- 0.1
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.6
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
#gammas = deaths + alphas + deltas

# mobility
#commuting and migration networks

muw <- 0.2
sw <- 0.05
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- .15 #gamma of baron et al
rw <- .1
cw <- .3

(Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)
com_mat <- rand_mat_cor_norm_MPA(N,muw,sw,rhow,Gammaw,rw,cw)
com_mat_d <- rand_mat_cor_norm_N(N,muw,sw,rhow,Gammaw,rw,cw)

muc <- 0.001
sc <- 0.0005
rhoc <- .001
Gammac <- .004
rc <- .03
cc <- .06

(Gammac/sqrt(rc*cc) < 1) & ((N*rhoc-2*Gammac)/(N-(rc+cc)) < 1)

# COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
# COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
# COMMUTING <- rand_mat_cor_norm(N, muw, sw, rhow, Gammaw, rw, cw)
# COMMUTING <- rand_mat_cor_trun(N, muw, sw, rhow, Gammaw, rw, cw)
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)[[1]]
# COMMUTING[sample.int(N^2, round(p*N^2))] <- 0

# MIGRATION <- matrix(rep(0, N^2), nrow = N)
# MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
# MIGRATION <- rand_mat_ell(N, muc, sc, rhoc, distrib = "beta")
# MIGRATION <- rand_mat_cor_norm(N, muc, sc, rhoc, Gammac, rc, cc)
# MIGRATION <- rand_mat_cor_trun(N, muc, sc, rhoc, Gammac, rc, cc)
MIGRATION <- rand_mat_cor_beta(N, muc*N, sc*N, rhoc, Gammac, rc, cc)[[1]]
sum(colSums(MIGRATION) > 1)

# jacobian

diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)

####### INTEGRATE SYSTEM ################################

# initial populations
# for constant populations, set deltas = 0, Deltas = deaths

sus_init <- rep(100000, N) # initial susceptibles
inf_init <- rep(100, N)    # initial infecteds

end_time <- 25

# integro el sistema con condiciones iniciales 
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

####### PLOTS ###########################################

# plot the mobility network
#legend for plotmobility2 can be found in RMT_plotmobility
plotmobility(COMMUTING)
plotmobility(MIGRATION)
plotmobility2(MIGRATION, COMMUTING)

# plot the eigenvalues of the system

plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0) #+ xlim(c(-215,-210))

# plot SUS, INF, REC or TOT population
plot_int(N, sol, state = "INF")

sol_df <-  as.data.frame(sol)
for(i in c(1:N)){
  colnames(sol_df)[i+1] <-  paste0("S",i)
  colnames(sol_df)[N+i+1] <-  paste0("I",i)
  colnames(sol_df)[2*N+i+1] <-  paste0("R",i)
}

####### TEST 1: ROLE OF THE 2 SCALES OF MOBILITY ########

# MU

mutilde <- 0.181
mobplots <- infplots <- eigplots <- list()

for (lambda in seq(0.1,0.9,by = 0.1)) {
  
  muw <- (1-lambda)*mutilde/mub
  muc <- lambda*mutilde
  
  COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
  MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
  
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
             COMMUTING,MIGRATION,
             sus_init,inf_init,end_time)
  
  mobplots[[as.character(lambda)]] <- plotmobility2(MIGRATION, COMMUTING)
  infplots[[as.character(lambda)]] <- plot_int(N, sol, state = "INF")
  eigplots[[as.character(lambda)]] <- plot_eigen_rmt(jacobian,
                                                     N,mub,mug = mud + mua + mudel,
                                                     muw,sw,rhow,Gammaw,
                                                     muc,sc,rhoc,Gammac,
                                                     tau = 0)
}

stilde <- 0.045
mobplots <- infplots <- eigplots <- list()

for (lambda in seq(0,1,by = 0.1)) {
  
  sw <- sqrt(1-lambda)*stilde/mub
  sc <- sqrt(lambda)*stilde
  
  COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
  MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
  
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
             COMMUTING,MIGRATION,
             sus_init,inf_init,end_time)
  
  mobplots[[as.character(lambda)]] <- plotmobility2(MIGRATION, COMMUTING)
  infplots[[as.character(lambda)]] <- plot_int(N, sol, state = "INF")
  eigplots[[as.character(lambda)]] <- plot_eigen_rmt(jacobian,
                                                     N,mub,mug = mud + mua + mudel,
                                                     muw,sw,rhow,Gammaw,
                                                     muc,sc,rhoc,Gammac,
                                                     tau = 0)
}

# SIGMA

stilde <- 0.045
mobplots <- infplots <- eigplots <- list()
  
for (lambda in seq(0,1,by = 0.1)) {
  
  sw <- sqrt(1-lambda)*stilde/mub
  sc <- sqrt(lambda)*stilde
  
  COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
  MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
  
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
             COMMUTING,MIGRATION,
             sus_init,inf_init,end_time)
  
  mobplots[[as.character(lambda)]] <- plotmobility2(MIGRATION, COMMUTING)
  infplots[[as.character(lambda)]] <- plot_int(N, sol, state = "INF")
  eigplots[[as.character(lambda)]] <- plot_eigen_rmt(jacobian,
                                                     N,mub,mug = mud + mua + mudel,
                                                     muw,sw,rhow,Gammaw,
                                                     muc,sc,rhoc,Gammac,
                                                     tau = 0)
}

####### TEST 2: RANDOM VS DIRECTED CONTROL ##############

MIGRATION <- matrix(rep(0,N^2),nrow = N)
MIGRATION <- rand_mat_ell(N, muc, sc, rhoc, distrib = "beta")
COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)

COMMUTINGD <- COMMUTINGR <- COMMUTING

muwstar <- 0.8

COMMUTINGD[1,2:N] <- muwstar
COMMUTINGR[sample(c(1:N^2),N-1)] <- muwstar

MIGRATIOND <- MIGRATIONR <- MIGRATION

mucstar <- 0.2

MIGRATIOND[1,2:N] <- mucstar
MIGRATIONR[sample(c(1:N^2),N-1)] <- mucstar

jacobiand <- (COMMUTINGD + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
jacobianr <- (COMMUTINGR + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)

jacobiand <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATIOND - diag(colSums(MIGRATIOND) + deaths + alphas + deltas)
jacobianr <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATIONR - diag(colSums(MIGRATIONR) + deaths + alphas + deltas)

plotmobility(COMMUTINGD)
plotmobility(COMMUTINGR)
plotmobility(MIGRATIOND)
plotmobility(MIGRATIONR)
plotmobility2(MIGRATIOND, COMMUTING)
plotmobility2(MIGRATIONR, COMMUTING)

sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
plot_int(N, sol, state = "INF")

sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATIOND,
           sus_init,inf_init,end_time)
plot_int(N, sol, state = "INF")

plot_eigen_rmt(jacobiand,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0)
eigen_dfd <- eigen_mat(jacobiand)

plot_eigen_rmt(jacobianr,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0)
eigen_dfr <- eigen_mat(jacobianr)

#[1,2:N] and N-1 seem to give the same result for muw* = 0
#looking at the expression this seems to make sense
#for muw* = .8 this seems to be the case as well, although
#the dynamics of the epidemics are quite different: the
#directed perturbation causes a single patch to have a
#very high number of infecteds (mirar cuando se une, si se une),
#while the random one causes
#a more homogeneous and milder rise
#random case seems to always give a slightly larger outlier
#the same seems to happen with [2:N,1], for muw* = 0. for
#muw* = .8 we have a single patch with a now much lower
#number of infecteds in the directed case

#for migration:
#[1,2:N] and mucstar = 0 means no difference, as expected. mucstar = .2
#turns a stable into an asymptotically system into an asymptotically
#unstable one. No particular patches show a qualitatively different
#behaviour, and the directed modification results in a steeper growth
#[2:N,1], mucstar = .2 again turns a stable system into an unstable one,
#and the directed perturbation does result in a single patch with a high
#number of infecteds. esto no tiene sentido porque deber?a ser al rev?s.
#revisar por qu? es

#three nodes - qualitatively same results with 3 nodes instead of one
COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
MIGRATION <- rand_mat_ell(N, muc, sc, rhoc, distrib = "beta")
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)

COMMUTINGD <- COMMUTINGR <- COMMUTING

muwstar <- 0.8

COMMUTINGD[1,2:N] <- COMMUTINGD[2,3:N] <- COMMUTINGD[2,1] <- COMMUTINGD[3,1:2] <- COMMUTINGD[3,4:N] <- muwstar
COMMUTINGD[2:N,1] <- COMMUTINGD[3:N,2] <- COMMUTINGD[1,2] <- COMMUTINGD[1:2,3] <- COMMUTINGD[4:N,3] <- muwstar
COMMUTINGR[sample(c(1:N^2),3*(N-1))] <- muwstar

MIGRATIOND <- MIGRATIONR <- MIGRATION

mucstar <- 0.2

MIGRATIOND[1,2:N] <- MIGRATIOND[2,3:N] <- MIGRATIOND[2,1] <- MIGRATIOND[3,1:2] <- MIGRATIOND[3,4:N] <- muwstar
MIGRATIOND[2:N,1] <- MIGRATIOND[3:N,2] <- MIGRATIOND[1,2] <- MIGRATIOND[1:2,3] <- MIGRATIOND[4:N,3] <- muwstar
MIGRATIONR[sample(c(1:N^2),3*(N-1))] <- mucstar

jacobiand <- (COMMUTINGD + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
jacobianr <- (COMMUTINGR + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)

jacobiand <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATIOND - diag(colSums(MIGRATIOND) + deaths + alphas + deltas)
jacobianr <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATIONR - diag(colSums(MIGRATIONR) + deaths + alphas + deltas)

plotmobility(COMMUTINGD)
plotmobility(COMMUTINGR)

plotmobility(MIGRATIOND)
plotmobility(MIGRATIONR)

sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATIOND,
           sus_init,inf_init,end_time)
plot_int(N, sol, state = "INF")

sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATIONR,
           sus_init,inf_init,end_time)
plot_int(N, sol, state = "INF")

plot_eigen_rmt(jacobiand,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0)
eigen_dfd <- eigen_mat(jacobiand)

plot_eigen_rmt(jacobianr,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0)
eigen_dfr <- eigen_mat(jacobianr)

####### TEST 2.5: PERTURBATIONS ###################################

# improve

errores <- vector()
for (j in c(1:1000)) {
  
  N <- 300
  beta <- 2
  alpha <- 0
  gamma <- 4
  muw <- 3
  muws <- 20
  sw <- 3
  
  betas <- matrix(rep(0,N^2), nrow = N)
  diag(betas) <- rep(beta,N)
  betas[1,1] <- beta+alpha
  
  BIGT <- rand_mat(N, muw, sw, distrib = "gamma")
  diag(BIGT) <- rep(1,N)
  #BIGT[1,2] <- muws # single edge, single direction
  #BIGT[1,2] <- BIGT[2,1] <- muws # single edge, both directions
  BIGT[1,c(2:N)] <- BIGT[2,1] <- BIGT[2,c(3:N)] <- BIGT[2,1]+muws-muw # two rows
  #BIGT[1,c(2:N)] <- BIGT[2,1] <- BIGT[2,c(3:N)] <- BIGT[3,c(1:2)] <- BIGT[3,c(4:N)] <- BIGT[2,1]+muws-muw # three edges, all directions
  BIGT <- BIGT%*%betas
  
  BIGS <- matrix(rep(0,N^2), nrow = N)
  diag(BIGS) <- rep(-gamma, N)
  
  jacobian <- BIGT+BIGS
  ngm <- -BIGT%*%solve(BIGS)
  
  outlier <- beta*(1-muw +(muw*(N-1))/2 +muws/2 + sqrt(muw^2*((N-3)^2)/4 - (muw*muws*(3*N-5))/2 + (muws^2)/4) ) - gamma
  #outlier <- beta*(1-muw +muw*N/2 + (1/2)*(alpha/beta) + sqrt((muw^2)*(N^2)/4+(N-1)*muw^2*alpha/beta-(N/2-1)*muw*alpha/beta+((1/2)*(alpha/beta))^2))-gamma # single edge, single direction
  
  plot_eigen(jacobian) + coord_fixed() + theme_bw() +
    geom_vline(xintercept = 0, color = "blue") +
    #geom_circle(aes(x0 = beta-gamma, y0 = 0, r = beta*(N-1)*muw), color = "green") +
    #geom_circle(aes(x0 = beta-gamma, y0 = 0, r = beta*((N-2)*muw+muws)), color = "green") +
    geom_point(aes(x = trueout, y = 0), color = "red") +
    geom_ellipse(aes(x0 = beta*(1-muw)-gamma, y0 = 0, a = sqrt(N)*beta*sw, b = sqrt(N)*beta*sw, angle = 0), color = "red")
  
  #outlier
  #eigen_mat(jacobian)$re[1]
  
  errores[j] <- outlier-eigen_mat(jacobian)$re[1]
  
}
###### TEST 3: CLT ON EIGENVALUE DISTRIBUTION ###########

MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")

mub*(1-muw)-(mua+mud+mudel)-N*muc 
mub*(1-muw)-(mua+mud+mudel)-(colSums(MIGRATION) + muc)
mean(mub*(1-muw)-(mua+mud+mudel)-(colSums(MIGRATION) + muc))
sd(mub*(1-muw)-(mua+mud+mudel)-(colSums(MIGRATION) + muc))

sc*sqrt(N)
colSums(MIGRATION) - (N-1)*muc
mean(colSums(MIGRATION) - (N-1)*muc)
sd(colSums(MIGRATION) - (N-1)*muc)

errores <- vector()
for (j in c(1:100)) {
  COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
  MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
  
  diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  errores[j] <- eigen_mat(jacobian)[100,1] - (mub - mua-mud-mudel + mub*muw*(N-1))
}

mean(errores)
sd(errores)

###### TEST 4: ACCURACY OF OUR BARON-GALLA MATRICES ##################

# let's do first several iterations with the same parameters and compare the
# estimates to discard statistical variations

outliers <- data.frame(real = vector(), rmt = vector(), bg = vector())
reports <- list()
for (j in c(1:1000)) {
  
  COM <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)
  MIG <- rand_mat_cor_beta(N, muc*N, sc*N, rhoc, Gammac, rc, cc)
  
  COMMUTING <- COM[[1]]
  MIGRATION <- MIG[[1]]
  
  tau <- 0
  mug <- mud + mua + mudel
  sigma <- sqrt(mub^2*sw^2 + 2*mub*tau + sc^2)
  rho <- (mub^2*rhow*sw^2 + rhoc*sc^2)/(sigma^2)
  Gamma <- (mub^2*Gammaw*(sw^2)/N + Gammac*(sc^2)/N)/(sigma^2)
  
  diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  outliers[j,] <- c(eigen_mat(jacobian)[1,1],mub - mug + mub*muw*(N-1),
                    mub - mug + mub*muw*(N-1) +
                    (1/2)*(N*mub*muw+N*muc)*(1+rho/Gamma)*(sqrt(1+(4*Gamma*sigma^2)/((mub*muw + muc)^2))-1))
  
  reports[[j]] <- list(COM[[2]],MIG[[2]])
  
  if(j%%20 == 0){print(j)}
}
print(c(mean(outliers[,1]),outliers[1,2:3]))

# different parameters now

mubs <- c(.1,1,1.5)
muas <- c(.1,.4,.9)
muws <- c(0.05, 0.2, 0.5)
sws <- c(0.03, 0.1, 0.3)
mucs <- c(0.001, 0.1, 0.3)
scs <- c(0.0005, 0.05, 0.1)
rhos <- Gammas <- rs <- cs <- c(.1,.3,.5,.8)

rhoc <- 0
Gammac <- 0
rc <- 0
cc <- 0

parameters <- expand.grid(mub = mubs, mua = muas,
                          muw = muws, sw = sws,
                          muc = mucs, sc = scs,
                          rho = rhos, Gamma = Gammas,
                          r = rs, c = cs) %>%
  mutate(real = rep(0,nrow(.)), rmtpred = rep(0,nrow(.)), bgpred = rep(0,nrow(.)))

reports <- list()

for (i in c(1:nrow(parameters))) {
  
  mub <- parameters$mub[i]
  mua <- parameters$mua[i]
  
  muw <- parameters$muw[i]
  sw <- parameters$sw[i]
  rhow <- parameters$rho[i]
  Gammaw <- parameters$Gamma[i]
  rw <- parameters$r[i]
  cw <- parameters$c[i]
  muc <- parameters$muc[i]
  sc <- parameters$sc[i]
  
  if ((Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)) {
  
  genbeta <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)
  
  COMMUTING <- genbeta[[1]]
  MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
  diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  sigma <- sqrt(mub^2*sw^2 + 2*mub*tau + sc^2)
  rho <- (mub^2*rhow*sw^2 + rhoc*sc^2)/(sigma^2)
  Gamma <- (mub^2*Gammaw*(sw^2)/N + Gammac*(sc^2)/N)/(sigma^2)
  
  parameters$real[i] <- eigen_mat(jacobian)[1,1]
  parameters$rmtpred[i] <- mub - mua - mud - mudel + mub*muw*(N-1)
  parameters$bgpred[i] <- mub - mua - mud - mudel + mub*muw*(N-1) +
    (1/2)*(N*mub*muw+N*muc)*(1+rho/Gamma)*(sqrt(1+(4*Gamma*sigma^2)/((mub*muw + muc)^2))-1)
  
  reports[[i]] <- genbeta[[2]]
  
  }
}


####### TEST 5: DIFFERENTLY CORRELATED MOBILITY NETWORKS ##############

# "base" network
muw <- 0.3
sw <- 0.15
rhow <- 0.1 #original rho (Gamma of baron et al)
Gammaw <- 0 #gamma of baron et al
rw <- 0.1
cw <- 0.1
MIGRATION <- matrix(rep(0, N^2), nrow = N)
muc <- sc <- rhoc <- Gammac <- rc <- cc <- 0

(abs(Gammaw/sqrt(rw*cw)) < 1) & (abs((N*rhow-2*Gammaw)/(N-(rw+cw))) < 1)

rmat <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)
rmat[[2]]
COMMUTING <- rmat[[1]] 
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
plotmobility(COMMUTING)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
plot_int(N, sol, state = "INF")
plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0) #+ xlim(c(-215,-210))

# move sw
#lower/higher values account for an almost constant/very high localized flows (caution)
#and lower/higher quantitative variability across patches, as expected. same outlier
sw <- 0.01
sw <- 0.3
(abs(Gammaw/sqrt(rw*cw)) < 1) & (abs((N*rhow-2*Gammaw)/(N-(rw+cw))) < 1)
rmat <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)
rmat[[2]]
COMMUTING <- rmat[[1]]
sum(COMMUTING > 1)
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
plotmobility(COMMUTING)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
plot_int(N, sol, state = "INF")
plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0) #+ xlim(c(-215,-210))

# move rhow
#the simulated rho misses (specially for negative values) but a wide range is available
#positive/negative values yield a more symmetric/antisymmetric network
#and lower/higher quantitative variability across patches, as expected. same outlier
sw <- .15
rhow <- .6
rhow <- -.6
(abs(Gammaw/sqrt(rw*cw)) < 1) & (abs((N*rhow-2*Gammaw)/(N-(rw+cw))) < 1)
rmat <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)
rmat[[2]]
COMMUTING <- rmat[[1]]
sum(COMMUTING > 1)
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
plotmobility(COMMUTING)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
plot_int(N, sol, state = "INF")
plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc = 0,sc = 0,rhoc = 0,Gammac = 0,
               tau = 0) #+ xlim(c(-215,-210))

# move Gammaw
#we need high r and c for high Gamma. These miss a little bit (more than Gamma)
#positive/negative values yield a less/more structured/? network? the effect is not
#very noticeable (recall that Gamma is in a 1/N scale). positive/negative values
#yield lower/higher quantitative variability across patches.
#Baron-Galla's prediction for the outlier is noticeably worse than RMT. it improves
#after using the real values but still worse
rhow <- .1
Gammaw <- .5
Gammaw <- -.4
rw <- cw <- .6
(abs(Gammaw/sqrt(rw*cw)) < 1) & (abs((N*rhow-2*Gammaw)/(N-(rw+cw))) < 1)
rmat <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)
rmat[[2]]
COMMUTING <- rmat[[1]]
sum(COMMUTING > 1)
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
plotmobility(COMMUTING)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
plot_int(N, sol, state = "INF")
plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc = 0,sc = 0,rhoc = 0,Gammac = 0,
               tau = 0) #+ xlim(c(-215,-210))

# move rw and cw
#we need high r and c for high Gamma. These miss a little bit (more than Gamma)
#positive/negative values yield a less/more structured/? network? the effect is not
#very noticeable (recall that Gamma is in a 1/N scale). positive/negative values
#yield lower/higher quantitative variability across patches.
#Baron-Galla's prediction for the outlier is noticeably worse than RMT. it improves
#after using the real values but still worse
Gammaw <- 0
rw <- .7
cw <- .1
(abs(Gammaw/sqrt(rw*cw)) < 1) & (abs((N*rhow-2*Gammaw)/(N-(rw+cw))) < 1)
rmat <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)
rmat[[2]]
COMMUTING <- rmat[[1]]
sum(COMMUTING > 1)
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
plotmobility(COMMUTING)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
plot_int(N, sol, state = "INF")
plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc = 0,sc = 0,rhoc = 0,Gammac = 0,
               tau = 0) #+ xlim(c(-215,-210))

####### TEST 6: STABILITY REGIONS (MARTA) #############################

N <- 100

step <- 0.005
beta_vec <- seq(0.01,0.9,step)
muw_vec <- seq(0.01,0.9,step)
mug <- mud + mua + mudel

df_sol <- data.frame(beta = 0, gamma = 0, N = 0, muw = 0, state = FALSE)
for(i in c(1:length(beta_vec))){
  print(paste0(" i : ", i,"/",length(beta_vec)))
  for(j in c(1:length(muw_vec))){
    if(is.na(muw_vec[j]) | muw_vec[j] > 1 | muw_vec[j] <0  ){
      print(paste0("Problem with muw: ", muw_vec[j]))
    }else{
      # Computed by the numerically computed eigevalues;
      
      # COMMUTING <- rand_mat(N, muw_vec[j], sw, distrib = "beta")
      # diag(COMMUTING) <- 0
      # MIGRATION <- rand_mat(N, muc, sc, distrib = "beta")
      # diag(MIGRATION) <- 0      
      
      # EPI param:
      # betas <- rep(beta_vec[i], N)
      # deltas <- rep(mudel, N)
      # deaths <- rep(mud, N)
      # alphas <- rep(mua, N)
      
      # jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
      #   diag(deaths + alphas + deltas + colSums(MIGRATION))
      # eigen <- eigen_mat(jacobian)
      # max_eig <- max(eigen$re)
      
      # Computed by the prediction from RMT:
      state <- ifelse(((N-1)*muw_vec[j]) < ((mug/beta_vec[i]) - 1), TRUE, FALSE)
      df_sol <- rbind(df_sol, c(beta_vec[i], mug, N, muw_vec[j], state))
    }
  }
}

df_sol <- df_sol[-1,] %>% mutate(Stability = ifelse(.$state == TRUE, "Stable", "Unstable"))

library(latex2exp)

# Values for the points in the area graph:
stab_par <- 0.07
unstab_par <- 0.2

# Create annotate for labels at each point in the area graph:
annotation <- data.frame(
  x = c(stab_par + 0.5,stab_par + 0.5, unstab_par + 0.04),
  y = c(stab_par + 0.3,unstab_par + 0.02 , stab_par + 0.3),
  label = c("c", "d", "e"))

library(ggstar)
plot_area <- ggplot(df_sol) +
  geom_point(aes(beta,muw, colour = Stability)) + theme_bw()  +
  scale_color_manual(values=c("#3066BE", "#A63446")) +
  ylab(TeX("$\\mu_w$")) +
  xlab(TeX("$\\beta$")) +
  # ggtitle(""*gamma/beta~": 4")
  ggtitle(paste0("N: ",N)) +
  coord_fixed() +
  geom_point(aes(stab_par,stab_par), colour= "#ADA544", size = 3) +
  geom_point(aes(stab_par,unstab_par), colour= "#ADA544", size = 3) +
  geom_point(aes(unstab_par,stab_par), colour= "#ADA544", size = 3) +
  geom_text(data=annotation, aes(x=x, y=y, label=label),
            color="#ADA544", size=9, angle=0, fontface="bold") +
  theme(text = element_text(size = 30), legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=5)))

plot_area <- plot_area + labs(title = "b")
plot_area

#######################################################################

