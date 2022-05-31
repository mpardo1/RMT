####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### parent script
###
### generate, plot and integrate metapopulation
### epidemiological models
### 
#########################################################
rm(list=ls()) 

source("~/RMT/David/RMT_genrandom_1.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")

####### GENERATE JACOBIAN ###############################

# number of patches
N <- 400

muw <- 1
sw <- 0.7
rhow <- -0.4 #original rho (Gamma of baron et al)
Gammaw <- -0.5 #gamma of baron et al
rw <- 1.4
cw <- 1.6

(Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)
# Outliers:
outl_BG <- -1 + muw + (muw/2)*(1+rhow/Gammaw)*(sqrt(1+(4*Gammaw*sw^2/muw^2))-1)
outl <- -1 + muw + rhow*sw^2/muw

outl_BG_abs <- -1 + muw + (muw/2)*(1+rhow/Gammaw)*(sqrt(1+(4*abs(Gammaw)*sw^2/muw^2))-1)

# Matrix MPA form:
com_mat <- rand_mat_cor_norm_MPA(N,muw,sw,rhow,Gammaw,rw,cw)
plot_eig <- plot_eigen(com_mat[[1]])
com_mat[[2]]

plot_eig

plot_eig + 
  geom_point(aes(-outl_BG_abs,0), colour = "green") 

eig_MPA <- eigen_mat(com_mat[[1]])
max_eig <- max(eig_MPA$re)
eig_MPA_circ <- eig_MPA[-which(eig_MPA$re == max(eig_MPA$re)),]

library(latex2exp)
plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green") +
  geom_point(aes(outl,0), colour = "blue") +
  theme_bw()

eig_MPA <- eigen_mat(com_mat[[1]])
eig_MPA_circ <- eig_MPA[-which(eig_MPA$re == max(eig_MPA$re)),]

center <-  -1
radius <- sw
ggplot(eig_MPA) + 
  geom_point(aes(re, im), size = 0.1) + 
  geom_point(aes(outl_BG,0), colour = "red") +
  geom_point(aes(outl,0), colour = "green") +
  geom_ellipse(aes(x0 = center, y0 = 0, a = (1+rhow)*radius,
                   b = (1-rhow)*radius, angle = 0), color = "red") 

# Matrix DGG form:
com_mat_D <- rand_mat_cor_norm_DGG(N,muw,sw,rhow,Gammaw,rw,cw)
plot_eig <- plot_eigen(com_mat_D[[1]])
com_mat_D[[2]]

plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green")

plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green") +
  geom_point(aes(outl,0), colour = "blue", size= 0.2) +
  theme_bw() 

eig_DD <- eigen_mat(com_mat_D[[1]])
eig_DD_circ <- eig_DD[-which(eig_DD$re == max(eig_DD$re)),]

center <-  -1
radius <- sw/sqrt(N)
ggplot(eig_DD_circ) + 
  geom_point(aes(re, im), size = 0.1) + 
  geom_ellipse(aes(x0 = center, y0 = 0, a = (1+rhow)*radius,
                   b = (1-rhow)*radius, angle = 0), color = "red") 
  
# Reescale variables:
N = 300
muw <- 0.05
sw <- 0.3
rhow <- 0.45
Gammaw <- -0.5
rw <- 0.75
cw <- 0.6


df_err <- data.frame(it = 0, out_re = 0)
count=0
max_it <- 200
while(count < max_it){
  count = count +1
  com_mat_resc <- rand_mat_cor_norm_MPA_resc1(N,muw,sw,rhow,Gammaw,rw,cw)
  eig_resc <- eigen_mat(com_mat_resc[[1]])
  max_eig <- max(eig_resc$re)
  df_err[nrow(df_err) + 1, ] <-  c(count, max_eig)
}

outl_BG <- -1 + N*muw + ((N*muw/2)*(1+(rhow/(Gammaw)))*(sqrt(1+(4*Gammaw*sw^2/(N*muw^2)))-1))
outl <- -1 + N*muw 

df_err$out_RMT <- outl
df_err$out_BG <- outl_BG

mean(df_err$out_BG - df_err$out_re)
mean(df_err$out_RMT - df_err$out_re)

N = 300
muw <- 0.05
sw <- 0.3
rhow <- 0.4
Gammaw <- -0.5
rw <- 0.75
cw <- 0.6

# outl_BG <- -1 + N*muw + ((N*muw/2)*(1+(rhow/(N*Gammaw)))*(sqrt(1+(4*Gammaw*sw^2/muw^2))-1))
outl_BG <- -1 + N*muw + ((N*muw/2)*(1+(rhow/Gammaw))*(sqrt(1+(4*Gammaw*sw^2/(N*muw^2)))-1))
outl <- -1 + N*muw 

ifelse(((4*Gammaw*sw^2)< (-muw^2)),FALSE,TRUE)

com_mat_resc <- rand_mat_cor_norm_MPA_resc1(N,muw,sw,rhow,Gammaw,rw,cw)
plot_eig <- plot_eigen(com_mat_resc[[1]])
com_mat_resc[[2]]

plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green") +
  geom_point(aes(outl,0), colour = "blue") +
  theme_bw()

plot_eig + 
  geom_point(aes(outl,0), colour = "blue") +
  theme_bw()

com_mat_resc[[2]]
eig_MPA_resc <- eigen_mat(com_mat_resc[[1]])
eig_MPA_circ <- eig_MPA[-which(eig_MPA_resc$re == max(eig_MPA_resc$re)),]


center <-  -1
radius <- sqrt(N)*sw
plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green") +
  geom_point(aes(outl,0), colour = "blue") +
  geom_ellipse(aes(x0 = center, y0 = 0, a = (1+rhow)*radius,
                   b = (1-rhow)*radius, angle = 0), color = "red")

# Area of stability for rho and Gamma:
N = 100
muw <- 0.009
sw <- 0.1
rhow <- 0.4
Gammaw <- -0.5
rw <- 0.75
cw <- 0.6


outl_BG <- -1 + N*muw + 
  ((N*muw/2)*(1+(rhow/Gammaw))*(sqrt(1+(4*Gammaw*sw^2/(N*muw^2)))-1))

rho_vec <- seq(-1,1,0.01) 
Gamma_vec <- seq(-1,1,0.01) 
df_stab <- data.frame(rho = 0 , gamma = 0, stab = 0, pob = 0)
for(i in c(1:length(rho_vec))){
  for(j in c(1:length(Gamma_vec))){
    cond1 <- ifelse(abs(Gammaw) > sqrt(rw*cw), FALSE, TRUE)
    cond2 <- ifelse(abs(N*rhow - (2*Gammaw)) > (N-((rw+cw))), FALSE, TRUE)
    
   if(cond1 == FALSE){
     pob_s <- 1
   }else if( cond2 == FALSE ){
     pob_s <- 2
   }else{
     pob_s <- 0
    }
    
    rhow <- rho_vec[i]
    Gammaw <- Gamma_vec[j]
    if(Gammaw == 0 | Gammaw < (-N*muw^2)/(4*sw^2)){
      stab = -1
    }else{
      outl_BG <- -1 + N*muw + 
        ((N*muw/2)*(1+(rhow/Gammaw))*(sqrt(1+(4*Gammaw*sw^2/(N*muw^2)))-1))
      if(outl_BG>0){
        stab = 0
      }else if(outl_BG<0){
        stab = 1
      }
    }
   
   
    df_stab[nrow(df_stab)+1,] <- c(rhow, Gammaw, stab, pob_s)
  }
}

df_stab$stab[which(is.na(df_stab$stab) == TRUE)] <- -1

library("latex2exp")
ggplot(df_stab) + 
  geom_point(aes(rho, gamma, colour = factor(stab)))  + 
  scale_colour_manual(values = c("#464D77", "#36827F", "#F9DB6D")) + 
  scale_fill_discrete(name = "Stability", labels = c("NA", "FALSE", "TRUE")) +
  ylab(TeX("$\\Gamma$")) +
  xlab(TeX("$\\rho$")) +
  theme_bw()

ggplot(df_stab) + 
 geom_point(aes(rho, gamma, colour = factor(pob)), size = 0.1)

Gammaw <- 0.5
rhow <- -0.55
# outl_BG <- -1 + N*muw + ((N*muw/2)*(1+(rhow/(N*Gammaw)))*(sqrt(1+(4*Gammaw*sw^2/muw^2))-1))
outl_BG <- -1 + N*muw + ((N*muw/2)*(1+(rhow/Gammaw))*(sqrt(1+(4*Gammaw*sw^2/(N*muw^2)))-1))
outl <- -1 + N*muw 

ifelse(((4*Gammaw*sw^2)< (-muw^2)),FALSE,TRUE)

com_mat_resc <- rand_mat_cor_norm_MPA_resc1(N,muw,sw,rhow,Gammaw,rw,cw)
plot_eig <- plot_eigen(com_mat_resc[[1]])
com_mat_resc[[2]]

max_im <- max(eigen_mat(com_mat_resc[[1]])$im)

plot_eig + 
  geom_segment(aes(x = 0, y = -max_im, xend = 0, yend = max_im), color = "red") +
  theme_bw()

plot_eig + 
  geom_point(aes(outl_BG,0), colour = "green") +
  geom_point(aes(outl,0), colour = "blue") +
  geom_segment(aes(x = 0, y = -max_im, xend = 0, yend = max_im), color = "red") +
  theme_bw()

##### JACOBIAN ###
# Area of stability for rho and Gamma:
N = 50
Deltas <- rep(0.3, N) # birth rate
mub <- 1
sb <- 0.001
betas <- rep(mub, N) # transmission rates
# betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.1, N) # loss of immunity rates
mud <- 0.3
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.65
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas
mug <- gammas[1]
tau <- 0
# Mobility variables:
muw <- 0.01
sw <- 0.09
rhow <- - 0.4
Gammaw <- - 0.8
rw <- 0.9
cw <- 0.9

muc <- 0.01
sc <- sw
rhoc <- rhow
Gammac <- Gammaw
rc <- rw
cc <- cw

COMMUTING <- rand_mat_cor_norm_MPA_resc1(N,muw,sw,rhow,Gammaw,rw,cw)[[1]]
diag(COMMUTING) <- 1
MIGRATION <- rand_mat_cor_norm_MPA_resc1(N,muc,sc,rhoc,Gammac,rc,cc)[[1]]
diag(MIGRATION) <- 1

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

plot_eig <- plot_eigen(jacobian)
plot_eig

sh <- sqrt(mub^2*sw^2 + sc^2)
rhoh <- (rhow*mub^2*sw^2 + rhoc*sc^2)/sh^2
Gammah <- (Gammaw*mub^2*sw^2 + Gammac*sc^2)/(sh^2)
ch <- (cw*mub^2*sw^2 + cc*sc^2)/(sh^2)
rh <- (rw*mub^2*sw^2 + rc*sc^2)/(sh^2)

sq <- sqrt(1+((4*Gammah*sh^2)/N*(mub*muw+muc)^2))-1
outl_BG_h <- (N-1)*mub*muw + mub - mug + (N*(mub*muw + muc)/2)*(1+(Gammah/rhoh))*sq

out_RMT_h <- mub - mug + mub*muw*(N-1)

plot_eig +
  geom_point(aes(outl_BG_h,0), colour = "green") +
  geom_point(aes(out_RMT_h,0), colour = "blue", size = 0.1) +
  theme_bw()


rho_vec <- seq(-1,1,0.1)
Gamma_vec <- seq(-1,1,0.1)
df_out <- data.frame(rho = 0, gamma = 0, inc = 0)
for(i in c(1:length(rho_vec))){
  for(j in c(1:length(Gamma_vec))){
    rhoh <- rho_vec[i]
    Gammah <- Gamma_vec[j]
    
    inc <- (mub/2)*(mub*muw + muc)*(1+(rhoh/Gammah))*(sqrt(1+(4*Gammah*sh^2)/(mub*muw + muc)^2)-1)
    df_out[nrow(df_out)+1,] <- c(rhoh,Gammah, inc) 
  }
}

sw <- 0.05
sc <- 0.05
sh <- sqrt(mub^2*sw^2 + sc^2)
mub <- 0.01
muc <-  0.1
muw <-  0.1
sh <- 0.1
Gammah <- 0.8
rhoh <- 0
sq <- sqrt(1+((4*Gammah*sh^2)/((mub*muw+muc)^2)))-1
outl_BG_h <- (N-1)*mub*muw + mub - mug + (N*(mub*muw + muc)/2)*(1+(rhoh/(N*Gammah)))*sq
outl_BG_h - ((N-1)*mub*muw + mub - mug)

Gamma_vec <- seq(-1,1,0.01)
df_out <- data.frame(gamma = 0, outl_BG = 0, outl_RMT = 0)
outl_RMT <- (N-1)*mub*muw + mub - mug
for(j in c(1:length(Gamma_vec))){
  Gammah <- Gamma_vec[j]
  
  if((1+((4*Gammah*sh^2)/((mub*muw+muc)^2)) > 0) | Gammah == 0){
    sq <- sqrt(1+((4*Gammah*sh^2)/((mub*muw+muc)^2)))-1
    outl_BG_h <- (N-1)*mub*muw + mub - mug + (N*(mub*muw + muc)/2)*(1+(rhoh/(N*Gammah)))*sq
    outl_RMT <- (N-1)*mub*muw + mub - mug
  }else{
    outl_BG_h = 0
  }
 
  df_out[nrow(df_out)+1,] <- c(Gammah, outl_BG_h, outl_RMT) 
}

df_out <- df_out[-which(df_out$gamma == 0),]
ggplot(df_out) + 
  geom_line(aes(gamma, outl_BG), size = 0.4, colour = "green") + 
  geom_point(aes(gamma, outl_BG), size = 0.8, colour = "blue") + 
  geom_line(aes(gamma, outl_RMT), size = 0.4, colour = "red") + 
  theme_bw()
