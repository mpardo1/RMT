####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### article plots
###
### plots for the article
### works over "RMT_parent.R"
###
rm(list = ls())
source("~/RMT/David/RMT_genrandom.R")
source("~/RMT/David/RMT_plotmobility.R")
source("~/RMT/David/d_functions_eigen_int.R")
library("viridis")

####### COMMUTING: PANEL ################################
library("copula")
library("copula")
library("viridis")
library("ggsci")
library("grid")
library("reshape")
library("ggpubr")
library("cowplot")
library("ggforce")
library("latex2exp")
#----------------------------------------------------------------------------#
####### Parameter values #####

mub*muw*(N-1)+mub-mua-mud-mudel

col_unstab <- "#0D0C18"
colA <- "#F69A79"
colB <- "#F26430"
colC <- "#D2430F"
colD <- "#85DCFF"
colE <- "#0ABAFF"
colF <- "#0084B8"

# N <- 100
N <- 50
Deltas <- rep(0.4, N) # birth rate
mub <- 0.15
# sb <- 0.001
betas <- rep(mub, N) # transmission rates
#betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.1, N) # loss of immunity rates
mud <- 0.4
deaths <- rep(mud, N) # not disease-related death rates

mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates

muw <- 0.08
sw <- 0.095          # original rho (Gamma of baron et al)

muc <- 0.001
sc <- 0.0005

#### UNSTABLE ####
mua <- 0.2 
alphas <- rep(mua, N) # recovery rates
gammas = deaths + alphas + deltas
mug <- gammas[1]

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <-  rep(0,N)

jacobian <- (COMMUTING + diag(N)) %*% diag(betas)  -
  diag(deaths + alphas + deltas )
Kmat <- (COMMUTING + diag(N)) %*% diag(betas)%*%diag(1/(deaths + alphas + deltas))


eig_J <- eigen_mat(jacobian)
eig_J$type <- "Jacobian"
eig_K <- eigen_mat(Kmat)
eig_K$type <- "NGM"
eig_t <- rbind(eig_J, eig_K)
maxJ1 <- max(eig_J$re)
maxK1 <- max(abs(eig_K$re))

colA <- "#2E294E"
colB <- "#16BAC5"

ggunst <- ggplot() + 
  geom_circle(aes(x0 = 0, y0 = 0, r = 1), 
              size =  0.3, linetype = "dashed", color = "red",
              alpha = 0.5) + 
  geom_vline(xintercept = 0, 
             size =  0.4, linetype = "dashed", color = "red") +
  geom_point(data = eig_t, aes(re, im, color = type), size = 0.3) + 
  coord_fixed() + 
  theme_bw() +
  scale_color_manual(values = c(colA, colB),
                     name = NULL ) +
  geom_point(aes(maxJ1,0), size = 1, color = colA) +
  geom_point(aes(maxK1,0), size = 1, color = colB) +
  theme_bw() + 
  theme(text = element_text(size = 15),legend.position = "left",
        legend.text.align = 0) +
  guides(color = guide_legend(override.aes = list(size = 3))) 
ggunst
### Stable
mua <- 0.5 # stable
alphas <- rep(mua, N) # recovery rates
gammas = deaths + alphas + deltas
mug <- gammas[1]

muw <- 0.1
sw <- 0.095          # original rho (Gamma of baron et al)

muc <- 0.001
sc <- 0.0005

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
diag(COMMUTING) <-  rep(0,N)

jacobian <- (COMMUTING + diag(N)) %*% diag(betas)  -
  diag(deaths + alphas + deltas )
Kmat <- (COMMUTING + diag(N)) %*% diag(betas)%*%diag(1/(deaths + alphas + deltas))


eig_J <- eigen_mat(jacobian)
eig_J$type <- "Jacobian"
eig_K <- eigen_mat(Kmat)
eig_K$type <- "NGM"
eig_t <- rbind(eig_J, eig_K)
maxJ <- max(eig_J$re)
maxK <- max(abs(eig_K$re))

colA <- "#2E294E"
colB <- "#16BAC5"

ggst <- ggplot() + 
  geom_circle(aes(x0 = 0, y0 = 0, r = 1), 
              size =  0.3, linetype = "dashed", color = "red",
              alpha = 0.5) + 
  geom_vline(xintercept = 0, 
             size =  0.4, linetype = "dashed", color = "red") +
  geom_point(data = eig_t, aes(re, im, color = type), size = 0.3) + 
  coord_fixed() + 
  theme_bw() +
  scale_color_manual(values = c(colA, colB),
                     name = NULL ) +
  geom_point(aes(maxJ,0), size = 1, color = colA) +
  geom_point(aes(maxK,0), size = 1, color = colB) +
  theme_bw() + 
  theme(text = element_text(size = 15),legend.position = "left",
        legend.text.align = 0) +
  guides(color = guide_legend(override.aes = list(size = 3))) 
ggst

ggarrange( ggst + xlab("") + ylab("") ,ggunst + xlab("") + ylab("") , common.legend = TRUE, widths  = c(0.93,1))

