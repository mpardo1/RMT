####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### article plots
###
### plots for the article
### works over "RMT_parent.R"
###
# rm(list = ls())
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
### FUNCTION FOR PLOTS

if( getwd() == "/home/marta"){
  Path <- "~/Documentos/PHD/2022/RMT_SIR/Plots/panel4/"
}else{
  Path <- "~/Documents/PHD/2022/RMT_SIR/Plots/panel4/"
}

plot_panelc <- function(N, sol, COMMUTING, colormean, color_low = "#FFFFFF", color_high, IMIN = "N", IMAX = "N",CMMIN = "N", CMMAX = "N",
                        size_text = 15) {
  
  comm_flows <- data.frame(variable = paste0(rep("S",N),as.character(c(1:N))), 
                           incoming = rowSums(COMMUTING)/(N-1),
                           outgoing = colSums(COMMUTING)/(N-1)) %>%
    mutate(both = (incoming+outgoing)/2)
  
  sol_df <- as.data.frame(sol) %>% select(c(1,(N+2):(2*N+1)))
  names(sol_df) <- c("time", paste0(rep("S",N),as.character(c(1:N))))
  sol_df <- mutate(sol_df, Smean = rowSums(sol_df[,-1])/N) %>%
    reshape2::melt(id.vars = c("time")) %>%
    left_join(comm_flows, by = "variable")
  
  IMIN <- ifelse(IMIN == "N", min(sol_df$value), IMIN)
  IMAX <- ifelse(IMAX == "N", max(sol_df$value), IMAX)
  CMMIN <- ifelse(CMMIN == "N", min(sol_df$both), CMMIN)
  CMMAX <- ifelse(CMMAX == "N", max(sol_df$both), CMMAX)
  
  ggplot(filter(sol_df, variable != "Smean"), aes(time, value, group = variable)) + 
    geom_line(aes(colour = both), size = .5) +
    scale_color_gradient(low = color_low, high = color_high,
                         na.value = "white", limits = c(CMMIN,CMMAX)) +
    # scale_fill_gradient(low="grey90", high="purple") +
    geom_line(data = filter(sol_df, variable == "Smean"), aes(y = value), 
              color = colormean, size = 0.8, linetype = "dashed") +
    ylab("Infected individuals") +
    ylim(c(IMIN,IMAX)) +
    theme(legend.position = "none") +
    guides(color = "none") +
    theme_bw() + 
    theme(text = element_text(size = 15), legend.position = "bottom") 
  
}

### NETWORK
plotmobilityMPA <- function(mat, high_col, low_col = "#FFFFFF"){
  diag(mat) <- NA
  longData <- melt(mat)
  colnames(longData) <- c("X1", "X2", "value")
  ggplot(longData, aes(x = X2, y = X1)) + 
    geom_raster(aes(fill=value)) + 
    # scale_color_brewer(palette = "Dark2")+
    scale_fill_gradient(low=low_col, high=high_col, na.value = "white") +
    theme_bw() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position = "none") + 
    theme_void()  +
    coord_fixed() +
    ylim(c(N,0)) 
}

# N <- 100
N <- 40
# para los plots
sus_init <- rep(100000, N) # initial susceptibles
inf_init <- runif(N, min = 50,100)  # initial infecteds
end_time <- 200
end_time_rand <- 20

Deltas <- rep(0.1, N) # birth rate
mub <- 0.1
sb <- 0.001
betas <- rep(mub, N) # transmission rates
#betas <- rgamma(N, shape = (mub/sb)^2, rate = mub/(sb^2))
thetas <- rep(0.1, N) # loss of immunity rates
mud <- 0.1
deaths <- rep(mud, N) # not disease-related death rates
mua <- 0.45
alphas <- rep(mua, N) # recovery rates
mudel <- 0
deltas <- rep(mudel, N) # disease-related death rates
gammas = deaths + alphas + deltas

muw <- 0.12
sw <- 0.07/3
rhow <- 0 #original rho (Gamma of baron et al)

muc <- 0.0001
sc <- 0.00001
rhoc <- .001

mub*muw*(N-1)+mub-mua-mud-mudel

col_unstab = "#0D0C18"   # negro
col_stabA = "#06D622"  # verde
col_stabB = "#E96F1D"  # naranja
col_stabC = "#F8DC22"  # amarillo
col_stabD = "#B62F2F"  # rojo
col_stabE = "#2C5CE1"  # azul
col_stabF = "#AA2FB5"  # violeta

COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
MIGRATION <- rand_mat_ell(N, muc, sc, rhoc, distrib = "beta")
# MIGRATION <- matrix(0,N,N)
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

col_low <- "#90F0CE"
col_high <- "#3F83E9"
col_mean <- "#000205"
plot_stab_mob <- plotmobility(COMMUTING, col_low, col_high) + 
  theme(legend.position = "left") 
ggsave(file=paste0(Path,"unstnet.svg"))
plot_stab_int <- plot_panelc(N, sol, COMMUTING,col_mean, col_low, col_high )
ggsave(file=paste0(Path,"unstsol.svg"))

gg_stable <- ggarrange(plot_stab_mob,
                       plot_stab_int, 
                       nrow  =2, labels = c("Unstable scenario",""),
                       vjust = 1.1)

plot_grid(plot_stab_mob + ggtitle("Unstable Scenario"),
          plot_stab_int, 
          nrow  =2 ,
          vjust = 1.1)
# gg_stable<- annotate_figure(gg_stable, 
#                               top = textGrob("Stable scenario:",
#                                              x = 0.13,
#                                              gp = gpar(cex = 1.3)))
gg_stable
ggsave(file=paste0(Path,"unstable.svg"))


#### Eigenvalues
mug <- gammas[1]
tau <- 0
c_unstab <- mub*(1-muw) - mug - N*muc
r_unstab <- sqrt(mub^2*sw^2 + 2*mub*tau + sc^2)*sqrt(N)
o_unstab <- mub - mug + muw*mub*(N-1)

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))
plot_eigen(jacobian) + 
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = r_unstab,
                   b = r_unstab, 
                   angle = 0), color = col_circ) + geom_point(aes(o_unstab,0), 
                                          color = col_unstab, shape = 21) +
  geom_vline(xintercept = 0) +
  theme_bw()

stragAB <- function(k, nu){
  (N/2)*(mub*muw+muc) + ((k-1)/2)*(mub*nu) + mub*(1-muw) - mug - N*muc + 
    (1/2)*sqrt(N^2*(mub*muw + muc)^2 + ((k-1)*mub*nu)^2 + 2*(N+N*k-4*k)*(mub*muw + muc)*mub*nu)
  # N*(mub*muw+muc)/2 + (k-1)*mub*nu/2 +
  #   sqrt((N*(mub*muw+muc))^2+((k-1)*mub*nu)^2+2*(mub*nu)*(mub*muw+muc)*(N+N*k-4*k))/2  +
  #   mub*(1-muw) - mug - N*muc
} 

stragC <- function(k, nu){
  N*(mub*muw+muc)/2 + (k/2-1)*mub*nu/2 +
    sqrt((N*(mub*muw+muc))^2+mub*nu*(2*mub*muw+muc)*(N+3*N*k/2-2*(k/2+1))+
           ((mub*nu)^2)*(2*k*(N-k/2)+(k/2-1)^2))/2 +
    mub*(1-muw) - mug - N*muc
} 

stragDEF <- function(k, muwstar){
  mub*(muw*(N-k)+muwstar*k/N)*(N-1)/N + mub - mug
} 

 ### CONTROL
nodes <- 4
muwstar <- 0
COMMUTINGA <- COMMUTINGB <- COMMUTINGC <- COMMUTINGD <- COMMUTINGE <- COMMUTINGF <- COMMUTING

# STRATEGY A
fs <- sample(c(1:N), nodes)
COMMUTINGA[fs,] <- muwstar
solA <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTINGA,MIGRATION,
           sus_init,inf_init,end_time)

jacobianA <- (COMMUTINGA + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

nu <- muwstar - muw
outlA <- stragAB(nodes,nu)

plot_eigen(jacobianA) + 
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = r_unstab,
                   b = r_unstab, 
                   angle = 0), color = col_circ) +
  geom_point(aes(outlA,0), 
             color = col_stabA, shape = 21) +
  theme_bw()

  # STRATEGY B
cs <- sample(c(1:N), nodes)
COMMUTINGB[,fs] <- muwstar
solB <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGB,MIGRATION,
            sus_init,inf_init,end_time)

jacobianB <- (COMMUTINGB + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

nu <- muwstar - muw
outlB <- stragAB(nodes,nu)

plot_eigen(jacobianB) + 
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = r_unstab,
                   b = r_unstab, 
                   angle = 0), color = col_circ) +
  geom_point(aes(outlB,0), 
             color = col_stabB, shape = 21) +
  theme_bw()


# STRATEGY C
fcs <- sample(c(1:N), nodes/2)
fcs <- fs[c(1:(nodes/2))]
COMMUTINGC[fcs,] <- COMMUTINGC[,fcs] <- muwstar
solC <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGC,MIGRATION,
            sus_init,inf_init,end_time)

jacobianC <- (COMMUTINGC + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

nu <- muwstar - muw
outlC <- stragC(nodes,nu)

plot_eigen(jacobianC) + 
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = r_unstab,
                   b = r_unstab, 
                   angle = 0), color = col_circ) +
  geom_point(aes(outlC,0), 
             color = col_stabC, shape = 21) +
  theme_bw()

# RAND STRATEGIES #
end_time <- 200

# STRATEGY D
inds <- sample(c(1:N^2), nodes*N)
COMMUTINGD[inds] <- muwstar

solD <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGD,MIGRATION,
            sus_init,inf_init,end_time)

jacobianD <- (COMMUTINGD + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

nu <- muwstar - muw
outlD <- stragDEF(nodes,nu)

sig <- sqrt(mub^2*sw^2  + sc^2)
p <- (nodes*N)/(N^2)
sigma_spa <- sqrt((p*(1-p)*sig^2)+(p*(1-p)*(mub*muw + muc)^2) + sig^2*p^2)
radius_spa <- sqrt(N)*sigma_spa

plot_eigen(jacobianD) + 
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = radius_spa,
                   b = radius_spa, 
                   angle = 0), color = col_circ) +
  geom_point(aes(outlD,0), 
             color = col_stabD, shape = 21) +
  theme_bw() + coord_fixed()

# STRATEGY E
rinds <- inds[c(1:(nodes*N/2))]
COMMUTINGE[rinds] <- muwstar
COMMUTINGE <- t(COMMUTINGE)
COMMUTINGE[rinds] <- muwstar
COMMUTINGE <- t(COMMUTINGE)

solE <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGE,MIGRATION,
            sus_init,inf_init,end_time)

jacobianE <- (COMMUTINGE + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

nu <- muwstar - muw
outlE <- stragDEF(nodes,nu)

# correlation:
MAT_CORR <- jacobianE
vec1 <- MAT_CORR[upper.tri(MAT_CORR)]
vec2 <- MAT_CORR[lower.tri(MAT_CORR)]
rhoa <- cor(vec1,vec2)

MAT_CORR <-  COMMUTINGE
vec1 <- MAT_CORR[upper.tri(MAT_CORR)]
vec2 <- MAT_CORR[lower.tri(MAT_CORR)]
rhoa <- cor(vec1,vec2)
# rhoa <- length(rinds)/(N^2) - N
rhoa <- (mub^2*rhoa*sw^2)/(mub^2*sw+sc)^2

ar <- 0.05
br <-  0.008
plot_eigen(jacobianE) + 
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = ar,#sig*sqrt(N)*(1-rhoa),
                   b = br,#sig*(1+rhoa), 
                   angle = 0), color = col_circ) +
  geom_point(aes(outlE,0), 
             color = col_stabE, shape = 21) +
  theme_bw()
# STRATEGY F
COMMUTINGF <- COMMUTINGF*(1-nodes/N)
solF <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGF,MIGRATION,
            sus_init,inf_init,end_time)

jacobianF <- (COMMUTINGF + diag(N)) %*% diag(betas) + MIGRATION -
  diag(deaths + alphas + deltas + colSums(MIGRATION))

nu <- muwstar - muw
outlF <- stragDEF(nodes,nu)

plot_eigen(jacobianF) + 
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = r_unstab,
                   b = r_unstab, 
                   angle = 0), color = col_circ) +
  geom_point(aes(outlF,0), 
             color = col_stabF, shape = 21) +
  theme_bw()

####EIGENFULL###
eigstab <- eigen_mat(jacobian)
eigstab$str <- ""
eigstab$type <- "" 
eigenA <- eigen_mat(jacobianA)
eigenA$str <- "A"
eigenA$type <- "not-rand" 
eigenB <- eigen_mat(jacobianB)
eigenB$str <- "B"
eigenB$type <- "not-rand" 
eigenC <- eigen_mat(jacobianC)
eigenC$str <- "C"
eigenC$type <- "not-rand" 
eigenD <- eigen_mat(jacobianD)
eigenD$str <- "D"
eigenD$type <- "rand" 
eigenE <- eigen_mat(jacobianE)
eigenE$str <- "E"
eigenE$type <- "rand" 
eigenF <- eigen_mat(jacobianF)
eigenF$str <- "F"
 eigenF$type <- "rand" 

eigenT <- rbind(eigstab,eigenA,eigenB,eigenC,eigenD,eigenE,eigenF)

col_unstab = "#0D0C18"   # negro
col_stabA = "#06D622"  # verde
col_stabB = "#E96F1D"  # naranja
col_stabC = "#F8DC22"  # amarillo
col_stabD = "#B62F2F"  # rojo
col_stabE = "#2C5CE1"  # azul
col_stabF = "#AA2FB5"  # violeta

supsmall = 0.2
small = 0.5
mid = 1
big = 1.5

gg_eigen_full <- ggplot(eigenT) + 
  geom_point(aes(re, im), size = 0.2) +
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = r_unstab,
                   b = r_unstab, 
                   angle = 0), color = col_unstab, size = 2.5) + 
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = r_unstab,
                   b = r_unstab, 
                   angle = 0), color = col_stabA, size = 1.8) + 
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = r_unstab,
                   b = r_unstab, 
                   angle = 0), 
               color = col_stabB, size = 1.3) +
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = r_unstab,
                   b = r_unstab, 
                   angle = 0), color = col_stabC, size = 0.8) +
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = radius_spa,
                   b = radius_spa, 
                   angle = 0), color = col_stabD, size = 0.3)  +
  geom_ellipse(aes(x0 = c_unstab, y0 = 0, 
                   a = ar,#sig*sqrt(N)*(1-rhoa),
                   b = br,#sig*(1+rhoa), 
                   angle = 0), color = col_stabE, size = 0.3) +
  geom_ellipse(aes(x0 = c_unstab, y0 = 0 , 
                   a = r_unstab,
                   b = r_unstab, 
                   angle = 0), color = col_stabF, size = 0.3) +
  geom_point(aes(o_unstab,0), 
             color = col_unstab, shape = 21, size = 1.5) +
  geom_point(aes(outlA,0), 
             color = col_stabA, shape = 21, size = 1) +
  geom_point(aes(outlB,0), 
             color = col_stabB, shape = 21, size = 2.5) +
  geom_point(aes(outlC,0), 
             color = col_stabC, shape = 21, size = 1.5) +
  geom_point(aes(outlD,0), 
             color = col_stabD, shape = 23, size = 2.5) +
  geom_point(aes(outlE,0), 
             color = col_stabE, shape = 23, size = 2) +
  geom_point(aes(outlF,0), 
             color = col_stabF, shape = 23, size = 1.5) +
  geom_vline(xintercept = 0, color = "blue",  linetype = "dashed") +
  theme_bw()

gg_eigen_full
 #######

#### MOB HEAT MAP########

###########################################################

plotmobility <- function(mob, color1 = "#0000FF", cmin = "N", cmax = "N"){
  
  # cmin <- ifelse(cmin == "N", min(mob), cmin)
  # cmax <- ifelse(cmax == "N", max(mob), cmax)
  cmin = min(mob)
  cmax = max(mob)
  collow <- "white"
  
  N <- nrow(mob)
  if (is.matrix(mob) & (nrow(mob) == ncol(mob))) {
    
    mob <- t(mob)
    
    diag(mob) <- rep(0,N)
    
    #mob <- mob[seq(N,1,-1),]
    
    mob <- as.data.frame(mob)
    names(mob) <- c(1:N)
    mob$x <- c(N:1)
    mob <- mob %>% pivot_longer(as.character(c(1:N)), names_to = "y", values_to = "mob")
    
    mob$x <- factor(mob$x, c(N:1))
    mob$y <- factor(mob$y, c(N:1))
    
    ggplot(mob, aes(x,y)) +
      geom_tile(aes(fill = mob)) +
      scale_fill_gradient(low = collow, high = color1, na.value = "yellow", limits = c(cmin,cmax)) +
      theme_void() +
      theme(legend.position="none", panel.background=element_blank()) +
      coord_fixed()
    
  } else {
    print("mob needs to be a square matrix")
  }
}

CMAXDET <- max(c(COMMUTINGA,COMMUTINGB,COMMUTINGC))
CMAXRAND <- max(c(COMMUTINGD,COMMUTINGE, COMMUTINGF))

IMAX <- max(unlist(c(select(as.data.frame(solA), c((N+2):(2*N+1))),
            select(as.data.frame(solB), c((N+2):(2*N+1))),
            select(as.data.frame(solC), c((N+2):(2*N+1))),
            select(as.data.frame(solD), c((N+2):(2*N+1))),
            select(as.data.frame(solE), c((N+2):(2*N+1))),
            select(as.data.frame(solF), c((N+2):(2*N+1)))))) + 50
IMIN <- min(unlist(c(select(as.data.frame(solA), c((N+2):(2*N+1))),
              select(as.data.frame(solB), c((N+2):(2*N+1))),
              select(as.data.frame(solC), c((N+2):(2*N+1))),
              select(as.data.frame(solD), c((N+2):(2*N+1))),
              select(as.data.frame(solE), c((N+2):(2*N+1))),
              select(as.data.frame(solF), c((N+2):(2*N+1)))))) - 50
CMMAX <- max(c(rowSums(COMMUTINGA)/(N-1),colSums(COMMUTINGA)/(N-1),
               rowSums(COMMUTINGB)/(N-1),colSums(COMMUTINGB)/(N-1),
               rowSums(COMMUTINGC)/(N-1),colSums(COMMUTINGC)/(N-1),
               rowSums(COMMUTINGD)/(N-1),colSums(COMMUTINGD)/(N-1),
               rowSums(COMMUTINGE)/(N-1),colSums(COMMUTINGE)/(N-1),
               rowSums(COMMUTINGF)/(N-1),colSums(COMMUTINGF)/(N-1)))
CMMIN <- min(c(rowSums(COMMUTINGA)/(N-1),colSums(COMMUTINGA)/(N-1),
               rowSums(COMMUTINGB)/(N-1),colSums(COMMUTINGB)/(N-1),
               rowSums(COMMUTINGC)/(N-1),colSums(COMMUTINGC)/(N-1),
               rowSums(COMMUTINGD)/(N-1),colSums(COMMUTINGD)/(N-1),
               rowSums(COMMUTINGE)/(N-1),colSums(COMMUTINGE)/(N-1),
               rowSums(COMMUTINGF)/(N-1),colSums(COMMUTINGF)/(N-1)))              
# 
# commpanel <-ggarrange(plotmobility(COMMUTINGA, cmax = CMAX), plot_panelc(N, solA, COMMUTINGA, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
#                      plotmobility(COMMUTINGB, cmax = CMAX), plot_panelc(N, solB, COMMUTINGB, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
#                      plotmobility(COMMUTINGC, cmax = CMAX), plot_panelc(N, solC, COMMUTINGC, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
#                      plotmobility(COMMUTINGD, cmax = CMAX), plot_panelc(N, solD, COMMUTINGD, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
#                      plotmobility(COMMUTINGE, cmax = CMAX), plot_panelc(N, solE, COMMUTINGE, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
#                      plotmobility(COMMUTINGF, cmax = CMAX), plot_panelc(N, solF, COMMUTINGF, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
#                      ncol = 2, nrow = 6)
col_lstab <- "#ACA8C4"
col_stab <- "#615E71"
col_ldet <- "#B1C2B6"
col_det <- "#40B565"
col_lrand <- "#7F9BD5"
col_rand <- "#567DD1"
reescale = 1
mobstab <- plotmobility(reescale*COMMUTING, cmax = CMAXDET, color1 =col_stab )
mobstab
mobA <- plotmobility(reescale*COMMUTINGA , cmax = CMAXDET, color1 =col_det )
# mobA <- plotmobilityMPA(reescale*COMMUTINGA - 1, col_det )
mobA

mobB <- plotmobility(COMMUTINGB, cmax = CMAXDET, color1 =col_det )
mobB
mobC <- plotmobility(COMMUTINGC, cmax = CMAXDET, color1 =col_det )
mobC
mobD <- plotmobility(COMMUTINGD, cmax = CMAXRAND, color1 =col_rand )
mobD
mobE <- plotmobility(COMMUTINGE, cmax = CMAXRAND, color1 =col_rand )
mobE
mobF <- plotmobility(COMMUTINGF, cmax = CMAXRAND, color1 =col_rand )
mobF

####INTEGRATION####
int <- plot_panelc(N, sol, COMMUTING,col_mean, col_lstab, col_stab )
int
intA <- plot_panelc(N, solA, COMMUTINGA,col_mean, col_ldet, col_det, IMIN = IMIN, 
            IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + ylim(c(0,100))
intA
intB <- plot_panelc(N, solB, COMMUTINGB,col_mean, col_ldet, col_det, IMIN = IMIN, 
            IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + ylim(c(0,100))
intB
intC <- plot_panelc(N, solC, COMMUTINGC,col_mean, col_ldet, col_det, IMIN = IMIN, 
            IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX)  + ylim(c(0,100))
intC
intD <- plot_panelc(N, solD, COMMUTINGD,col_mean, col_lrand, col_rand, IMIN = IMIN, 
                    IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + ylim(c(0,100))
intD
intE <- plot_panelc(N, solE, COMMUTINGE,col_mean, col_lrand, col_rand, IMIN = IMIN, 
                    IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + ylim(c(0,100))
intE
intF <- plot_panelc(N, solF, COMMUTINGF,col_mean, col_lrand, col_rand, IMIN = IMIN, 
                    IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX)  + ylim(c(0,100))
intF
text_tit = 15
gstab <- plot_grid(mobstab + ggtitle("Unstable Scenario"),
          int ,
          plot_eig_lines + ggtitle("a") + theme(legend.position = "bottom"),
          nrow = 3, rel_heights = c(0.8,1,1.5),
          align = "v")

                
gdet <- plot_grid(mobA + ylab("") + ggtitle("Scenario A \n"),
                  mobB + ggtitle("Scenario B \n"),
                  mobC + ggtitle("Scenario C \n"),
                  intA + rremove("ylab") + rremove("xlab"),
                  intB + rremove("ylab") + rremove("xlab"),
                  intC + rremove("ylab") + rremove("xlab"), 
                  ncol = 3, nrow = 2,rel_heights = c(1,1),
                  align = "v")

gdet <- annotate_figure(gdet,
                left = textGrob("Infected individuals",
                                rot = 90, vjust = 1, y = 0.35, gp = gpar(cex = 1.3)),
                bottom = textGrob("Time", gp = gpar(cex = 1.3)))

grand <- plot_grid(mobD + ylab("") + ggtitle("Scenario D \n"),
                  mobE + ggtitle("Scenario E \n"),
                  mobF + ggtitle("Scenario F \n"),
                  intD + rremove("ylab") + rremove("xlab"),
                  intE + rremove("ylab") + rremove("xlab"),
                  intF + rremove("ylab") + rremove("xlab"), 
                  ncol = 3, nrow = 2,rel_heights = c(1,1),
                  align = "v")

grand <- annotate_figure(grand,
                        left = textGrob("Infected individuals",
                                        rot = 90, vjust = 1, y = 0.35, gp = gpar(cex = 1.3)),
                        bottom = textGrob("Time", gp = gpar(cex = 1.3)))

plotr <- plot_grid(gdet,
          grand,
          gg_eigen_full + ggtitle("b"),
          ncol = 1, nrow = 3,
          rel_heights = c(2,2,1))

plotl <- plot_grid(mobstab + ggtitle("Unstable Scenario\n"),
                   int ,
                   plot_eig_lines + ggtitle("a") + theme(legend.position = "none"),
                   nrow = 3, rel_heights = c(0.8,1,2),
                   align = "v")

plot_grid(plotl, plotr, ncol = 2, rel_widths = c(1,3))

legend <- get_legend(
  # create some space to the left of the legend
  plot_eig_lines + theme(legend.box.margin = margin(0, 0, 0, 12))
)
# commpanel11 <-ggarrange(plotmobility(COMMUTINGA, col_low, col_high),
#                         plot_panelc(N, solA, COMMUTINGA,col_mean, col_low, col_high, IMIN = IMIN, 
#                                     IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
#                           rremove("xlab") +  rremove("ylab"),
#                         ncol = 1, nrow = 2, legend = "none")
# 
# commpanel11<- annotate_figure(commpanel11, 
#                               top = textGrob("Scenario A:", x = 0.1, gp = gpar(cex = 1.3)),
#                               left = textGrob("Infected individuals",
#                                               rot = 90, vjust = 1, y = 0.3, gp = gpar(cex = 1.3)),
#                               bottom = textGrob("Time", gp = gpar(cex = 1.3)))
# commpanel11
# ggsave(file=paste0(Path,"scenarioA.svg"))
# 
# #--------#
# commpanel12 <-ggarrange(plotmobility(COMMUTINGB, col_low, col_high),
#                         plot_panelc(N, solB, COMMUTINGB,col_mean, col_low, col_high, IMIN = IMIN, 
#                                     IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
#                           rremove("xlab") +  rremove("ylab"),
#                         ncol = 1, nrow = 2, legend = "none")
# 
# commpanel12<- annotate_figure(commpanel12, 
#                               top = textGrob("Scenario B:", x = 0.1, gp = gpar(cex = 1.3)),
#                               left = textGrob("",
#                                               rot = 90, vjust = 1, y = 0.3, gp = gpar(cex = 1.3)),
#                               bottom = textGrob("Time", gp = gpar(cex = 1.3)))
# commpanel12
# ggsave(file=paste0(Path,"scenarioB.svg"))
# 
# #--------#
# commpanel13 <-ggarrange(plotmobility(COMMUTINGC, col_low, col_high),
#                         plot_panelc(N, solC, COMMUTINGC,col_mean, col_low, col_high, IMIN = IMIN, 
#                                     IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
#                           rremove("xlab") +  rremove("ylab") ,
#                         ncol = 1, nrow = 2, legend = "none")
# 
# commpanel13<- annotate_figure(commpanel13, 
#                               top = textGrob("Scenario C:", x = 0.1, gp = gpar(cex = 1.3)),
#                               left = textGrob("Infected individuals",
#                                               rot = 90, vjust = 1, y = 0.3, gp = gpar(cex = 1.3)),
#                               bottom = textGrob("Time", gp = gpar(cex = 1.3)))
# commpanel13
# ggsave(file=paste0(Path,"scenarioC.svg"))
# 
# #--------#
# commpanel21 <-ggarrange(plotmobility(COMMUTINGD, col_low, col_high),
#                         plot_panelc(N, solD, COMMUTINGD,col_mean, col_low, col_high, IMIN = IMIN, 
#                                     IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
#                           rremove("xlab") +  rremove("ylab"),
#                         ncol = 1, nrow = 2, legend = "none")
# 
# commpanel21<- annotate_figure(commpanel21, 
#                               top = textGrob("Scenario D:", x = 0.1, gp = gpar(cex = 1.3)),
#                               left = textGrob("",
#                                               rot = 90, vjust = 1, y = 0.3, gp = gpar(cex = 1.3)),
#                               bottom = textGrob("Time", gp = gpar(cex = 1.3)))
# commpanel21
# ggsave(file=paste0(Path,"scenarioD.svg"))
# 
# #--------#
# commpanel22 <-ggarrange(plotmobility(COMMUTINGE, col_low, col_high),
#                         plot_panelc(N, solE, COMMUTINGE,col_mean, col_low, col_high, IMIN = IMIN, 
#                                     IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
#                           rremove("xlab") +  rremove("ylab"),
#                         ncol = 1, nrow = 2, legend = "none")
# 
# commpanel22<- annotate_figure(commpanel22, 
#                               top = textGrob("Scenario E:", x = 0.1, gp = gpar(cex = 1.3)),
#                               left = textGrob("Infected individuals",
#                                               rot = 90, vjust = 1, y = 0.3, gp = gpar(cex = 1.3)),
#                               bottom = textGrob("Time", gp = gpar(cex = 1.3)))
# commpanel22
# ggsave(file=paste0(Path,"scenarioE.svg"))
# 
# #--------#
# commpanel23 <-ggarrange(plotmobility(COMMUTINGF, col_low, col_high),
#                         plot_panelc(N, solF, COMMUTINGF,col_mean, col_low, col_high, IMIN = IMIN, 
#                                     IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
#                           rremove("xlab") +  rremove("ylab"),
#                         ncol = 1, nrow = 2, legend = "none")
# 
# commpanel23<- annotate_figure(commpanel23, 
#                               top = textGrob("Scenario F:", x = 0.1, gp = gpar(cex = 1.3)),
#                               left = textGrob("",
#                                               rot = 90, vjust = 1, y = 0.3, gp = gpar(cex = 1.3)),
#                               bottom = textGrob("Time", gp = gpar(cex = 1.3)))
# commpanel23
# ggsave(file=paste0(Path,"scenarioF.svg"))


####
# plotmob <- ggarrange(plotmobility(COMMUTINGA, col_low, col_high),
#                   plotmobility(COMMUTINGB, col_low, col_high),
#                   plotmobility(COMMUTINGC, col_low, col_high),
#                   plotmobility(COMMUTINGD, col_low, col_high),
#                   plotmobility(COMMUTINGE, col_low, col_high),
#                   plotmobility(COMMUTINGF, col_low, col_high),
#                   nrow = 2, ncol = 3, labels = c("Strategy A","Strategy B",
#                                                  "Strategy C","Strategy D",
#                                                  "Strategy E","Strategy F"),
#                   legend = "none")
# 
# plotmob
# ggsave(file=paste0(Path,"mob_plots.svg"))
# 
# plotint <- ggarrange(plot_panelc(N, solA, COMMUTINGA,col_mean, col_low, col_high, IMIN = IMIN, 
#                                  IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
#                        rremove("xlab") +  rremove("ylab") ,
#                      plot_panelc(N, solB, COMMUTINGB,col_mean, col_low, col_high, IMIN = IMIN, 
#                                  IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
#                        rremove("xlab") +  rremove("ylab") ,
#                      plot_panelc(N, solC, COMMUTINGC,col_mean, col_low, col_high, IMIN = IMIN, 
#                                  IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
#                        rremove("xlab") +  rremove("ylab") ,
#                      plot_panelc(N, solD, COMMUTINGD,col_mean, col_low, col_high, IMIN = IMIN, 
#                                  IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
#                        rremove("xlab") +  rremove("ylab") ,
#                      plot_panelc(N, solE, COMMUTINGE,col_mean, col_low, col_high, IMIN = IMIN, 
#                                  IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
#                        rremove("xlab") +  rremove("ylab") ,
#                      plot_panelc(N, solF, COMMUTINGF,col_mean, col_low, col_high, IMIN = IMIN, 
#                                  IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
#                        rremove("xlab") +  rremove("ylab") ,
#                   nrow = 2, ncol = 3, labels = c("Strategy A","Strategy B",
#                                                  "Strategy C","Strategy D",
#                                                  "Strategy E","Strategy F"),
#                   legend = "none", vjust = 1.3,
#                   heights = c(0.8,0.8,0.8,0.8,0.8,0.8) ,
#                   align = "h")
# 
# plotint<- annotate_figure(plotint, 
#                               left = textGrob("Infected individuals",
#                                               rot = 90, vjust = 1, y = 0.5, gp = gpar(cex = 1.3)),
#                               bottom = textGrob("Time", gp = gpar(cex = 1.3)))
# 
# plotint
# ggsave(file=paste0(Path,"int_plots.svg"))

### COMPARISON

eig_comp <- data.frame(k = integer(0),
                       RMT = numeric(0),
                       LRP_in = numeric(0),
                       LRP_out = numeric(0),
                       LRP_inout = numeric(0),
                       RMT_allr = numeric(0),
                       Real_in = numeric(0),
                       Real_out = numeric(0),
                       Real_inout = numeric(0),
                       Real_sparse = numeric(0),
                       Real_corr = numeric(0),
                       Real_dec = numeric(0)
                       )

nu <- muwstar - muw
mug <- mud + mua + mudel
max_nodes <- floor(N/4)
max_nodes <- 20
for (nodes in seq(0,max_nodes,by = 2)) {
  
  COMMUTINGA <- COMMUTINGB <- COMMUTINGC <- COMMUTINGD <- COMMUTINGE <- COMMUTINGF <- COMMUTING
  fs <- sample(c(1:N), nodes)
  COMMUTINGA[fs,] <- muwstar
  COMMUTINGB[,fs] <- muwstar
  fcs <- sample(c(1:N), nodes/2)
  fcs <- fs[c(1:(nodes/2))]
  COMMUTINGC[fcs,] <- COMMUTINGC[,fcs] <- muwstar
  inds <- sample(c(1:N^2), nodes*N)
  COMMUTINGD[inds] <- muwstar
  rinds <- inds[c(1:(nodes*N/2))]
  COMMUTINGE[rinds] <- muwstar
  COMMUTINGE <- t(COMMUTINGE)
  COMMUTINGE[rinds] <- muwstar
  COMMUTINGE <- t(COMMUTINGE)
  COMMUTINGF <- COMMUTINGF*(1-nodes/N)
  
  jacA <- (COMMUTINGA + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacB <- (COMMUTINGB + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacC <- (COMMUTINGC + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacD <- (COMMUTINGD + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacE <- (COMMUTINGE + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacF <- (COMMUTINGF + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  k <- nodes
  eig_comp[k+1,] <- c(k,
                      mub - mug + muw*mub*(N-1),
                    N*(mub*muw+muc)/2 + (k-1)*mub*nu/2 +
                      sqrt((N*(mub*muw+muc))^2+((k-1)*mub*nu)^2+2*(mub*nu)*(mub*muw+muc)*(N+N*k-4*k))/2  +
                      mub*(1-muw) - mug - N*muc,
                    N*(mub*muw+muc)/2 + (k-1)*mub*nu/2 +
                      sqrt((N*(mub*muw+muc))^2+((k-1)*mub*nu)^2+2*(mub*nu)*(mub*muw+muc)*(N+N*k-4*k))/2  +
                      mub*(1-muw) - mug - N*muc,
                    N*(mub*muw+muc)/2 + (k/2-1)*mub*nu/2 +
                      sqrt((N*(mub*muw+muc))^2+mub*nu*(2*mub*muw+muc)*(N+3*N*k/2-2*(k/2+1))+
                             ((mub*nu)^2)*(2*k*(N-k/2)+(k/2-1)^2))/2 +
                      mub*(1-muw) - mug - N*muc,
                    mub*(muw*(N-k)+muwstar*k/N)*(N-1)/N + mub - mug,
                    max(eigen_mat(jacA)$re),
                    max(eigen_mat(jacB)$re),
                    max(eigen_mat(jacC)$re),
                    max(eigen_mat(jacD)$re),
                    max(eigen_mat(jacE)$re),
                    max(eigen_mat(jacF)$re)
                    )
}

eig_comp <- filter(eig_comp, !is.na(k))

eig_comp_real <- select(eig_comp, k, Real_in, Real_out, Real_inout, Real_sparse, Real_corr, Real_dec) %>%
  pivot_longer(cols = (2:ncol(.)), names_to = "prediction", values_to = "eigenvalue")
eig_comp_pred <- select(eig_comp, k, LRP_in, LRP_inout, RMT_allr) %>%
  pivot_longer(cols = (2:ncol(.)), names_to = "prediction", values_to = "eigenvalue")

# plot_eig_lines <- 
  eig_comp_pred %>% ggplot(aes(x = k, y = eigenvalue, color = prediction)) +
  geom_line(size = 0.8) + 
  geom_point(data = eig_comp_real, 
             aes(x = k, y = eigenvalue, color = prediction)) +
  scale_color_manual(values=c("#0000FF","#FF0000","#63B159",
                              "#0000FF","#0000FF","#FF0000",
                              "#63B159","#63B159","#63B159")) + 
  theme_bw() + theme(text = element_text(size = 15))
ggsave(file=paste0(Path,"commparison.svg"))

####### MIGRATION: DIRECTED VS RANDOM EXODUS ############

# base parameters + reduced muw (to induce stability with control)
muc <- 0.001
sc <- 0.0005

mucstar <- 0.005
MIGRATION <- rand_mat_cor_beta(N, muc*N, sc*N, rhoc, Gammac, rc, cc)[[1]]

# control: two rows
f1 <- 30
f2 <- 55
MIGRATION[f1,c(1:(f1-1),(f1+1):N)] <- MIGRATION[f2,c(1:(f2-1),(f2+1):N)] <- mucstar

# control: two columns
c1 <- 45
c2 <- 80
MIGRATION[c(1:(c1-1),(c1+1):N),c1] <- MIGRATION[c(1:(c2-1),(c2+1):N),c2] <- mucstar

# control: one row and one column
fc <- 65
MIGRATION[fc,c(1:(fc-1),(fc+1):N)] <- MIGRATION[c(1:(fc-1),(fc+1):N),fc] <- mucstar

# control: random
MIGRATION[sample(c(1:N^2),2*(N-1))] <- mucstar


diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
sum(colSums(MIGRATION) > 1)
mobplot <- plotmobility(MIGRATION)
mobplot
#print(mobplot, vp = viewport(angle = -90))
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)

plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0) #+ xlim(c(-215,-210))

sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
plot_int(N, sol, state = "INF")

####### CORRELATIONS: MONOTONICITY IN GAMMA #############

N <- 500

muw <- .05
sw <- .01
rhow <- .1
rw <- cw <- .7
COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")

muc <- 0.001
sc <- 0.0005
rhoc <- 0
Gammac <- 0
MIGRATION <- rand_mat_ell(N, muc, sc, rhoc, distrib = "beta")

outlier_df <- data.frame("Gammaw" = vector(), "realGammaw" = vector(), "bgpred" = vector(), "rmtpred" = vector(), "real" = vector())
for (Gammaw in seq(0.01,0.99,by = .005)) {
  
  if ((Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)) {
  
  gCOMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)  
  
  COMMUTING <- gCOMMUTING[[1]]
  diag(COMMUTING) <- rep(0,N)
  
  diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
  sum(colSums(MIGRATION) > 1)
  
  jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  sigma <- sqrt(mub^2*sw^2 + sc^2)
  rho <- (mub^2*rhow*sw^2 + rhoc*sc^2)/(sigma^2)
  Gamma <- (mub^2*Gammaw*(sw^2)/N + Gammac*(sc^2)/N)/(sigma^2)
  
  outlier_df <- rbind(outlier_df, c(Gammaw,
                                    gCOMMUTING[[2]][4,2],
                                    mub - mud - mua - mudel + mub*muw*(N-1) +
                                      (1/2)*(N*mub*muw+N*muc)*(1+rho/Gamma)*(sqrt(1+(4*Gamma*sigma^2)/((mub*muw + muc)^2))-1),
                                    mub - mud - mua - mudel + mub*muw*(N-1),
                                    max(eigen_mat(jacobian)$re)))
  print(Gammaw)
  }
}

names(outlier_df) <- c("Gammaw", "realGammaw", "bgpred", "rmtpred", "real")

outlier_df <- pivot_longer(outlier_df, cols = (2:ncol(outlier_df)) , names_to = "prediction", values_to = "outlier")

filter(outlier_df, prediction %in% c("realGammaw")) %>%
  ggplot(aes(x = Gammaw, y = outlier, color = prediction)) +
  geom_line(size = 1)

filter(outlier_df, prediction %in% c("bgpred", "rmtpred", "real")) %>%
  ggplot(aes(x = Gammaw, y = outlier, color = prediction)) +
  geom_line(size = 1)

outlier_df <- pivot_wider(outlier_df, names_from = prediction, values_from = outlier) %>%
  mutate(bgpredreal = mub - mud - mua - mudel + mub*muw*(N-1) +
           (1/2)*(N*mub*muw+N*muc)*(1+rho/realGammaw)*(sqrt(1+(4*realGammaw*sigma^2)/((mub*muw + muc)^2))-1))

outlier_df %>%
  ggplot(aes(x = Gammaw, y = N*realGammaw)) +
  geom_point() +
  geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) +
  geom_line(aes(x = Gammaw, y = rmtpred)) +
  geom_line(aes(x = Gammaw, y = bgpred))

filter(outlier_df, prediction %in% c("bgpred", "rmtpred", "real", "bgpredreal")) %>%
  ggplot(aes(x = Gammaw, y = outlier, color = prediction)) +
  geom_line(size = 1)

####### CORRELATIONS: INSTABILITY #######################

N <- 200

mub <- 0.02
betas <- rep(mub, N) # transmission rates
mud <- 0.1282
deaths <- rep(mud, N)
mua <- 0.25
alphas <- rep(mua, N)
mudel <- 0
deltas <- rep(mudel, N)

# re-scaling only mu and sigma: RMT and BG are indistinguishable
# muw <- 0.09
# sw <- 0.06
# rhow <- 0.45
# Gammaw <- 0.3
# rw <- 0.4
# cw <- 0.4

# re-scaling all parameters: ellipse is lost, outlier is still accurate
muw <- 0.09
sw <- 0.06
rhow <- 0.45
Gammaw <- 0.3*N
rw <- 0.4*N
cw <- 0.4*N

(Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

muc <- muw/100
sc <- sw/120
rhoc <- rhow/100
Gammac <- Gammaw/100
rc <- rw/100
cc <- cw/100

(Gammac/sqrt(rc*cc) < 1) & ((N*rhoc-2*Gammac)/(N-(rc+cc)) < 1)

COMMUTING <- rand_mat_cor_beta(N, N*muw, sqrt(N)*sqrt(N)*sw, rhow, Gammaw, rw, cw)[[1]]
MIGRATION <- rand_mat_cor_norm(N, N*muc, sqrt(N)*sqrt(N)*sc, rhoc, Gammac, rc, cc)[[1]]
# COM[[2]]
# MIG[[2]]
# COMMUTING <- COM[[1]]
# MIGRATION <- MIG[[1]]
#MIGRATION <- matrix(rep(0,N^2), nrow = N)
#muc <- sc <- rhoc <- Gammac <- 0
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigmayor <- eigen_mat(jacobian)
redmayor <- plotmobility(COMMUTING[1:80,1:80])

COMMUTING <- rand_mat_cor_norm(N, N*muw, sqrt(N)*sqrt(N)*sw, -rhow, -Gammaw, rw, cw)[[1]]
MIGRATION <- rand_mat_cor_norm(N, N*muc, sqrt(N)*sqrt(N)*sc, -rhoc, -Gammac, rc, cc)[[1]]
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigmenor <- eigen_mat(jacobian)
redmenor <- plotmobility(COMMUTING[1:80,1:80])

COMMUTING <- rand_mat_cor_norm(N, N*muw, sqrt(N)*sqrt(N)*sw, rhow, 0.001, 0.001, 0.001)[[1]]
MIGRATION <- rand_mat_ell(N, muc, sc, rhoc, distrib = "beta")
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigRMT <- eigen_mat(jacobian)
reduncor <- plotmobility(COMMUTING[1:80,1:80])

sig <- sqrt((mub*sw)^2 + sc^2)
rho <- (rhow*((mub*sw)^2) + rhoc*sc^2)/sig^2
Gamma <- (Gammaw*((mub*sw)^2) + Gammac*sc^2)/sig^2

outlierREmayor <- max(eigmayor$re)
outlierREmenor <- max(eigmenor$re)
outlierRMT <- (N-1)*mub*muw + mub - mud - mua - mudel
outlierBGmayor <- (N-1)*mub*muw + mub - mud - mua - mudel +
  (N/2)*(mub*muw+muc)*(1+rho/Gamma)*(sqrt(1+(4*Gamma*sig^2/(N*(mub*muw+muc)^2)))-1)
outlierBGmenor <- (N-1)*mub*muw + mub - mud - mua - mudel +
  (N/2)*(mub*muw+muc)*(1+rho/Gamma)*(sqrt(1-(4*Gamma*sig^2/(N*(mub*muw+muc)^2)))-1)

center <- mub*(1-muw)-mua-mud-mudel-N*muc
radius <- sqrt(N)*sig

color0 <- "black"
color1 <- "red"
color2 <- "blue"
color3 <- "green"
color4 <- "purple"

BGvsRMT <- ggplot(data = eigmenor) + 
  geom_ellipse(aes(x0 = center, y0 = 0, a = (1+rho)*radius,
                   b = (1-rho)*radius, angle = 0), color = color1) +
  geom_vline(xintercept = 0, color = color2, linetype = "dashed") +
  geom_point(aes(x = outlierRMT, y = 0), color = color1) +
  geom_point(aes(x = outlierBGmenor, y = 0), color = color3) +
  geom_point(aes(x = outlierBGmayor, y = 0), color = color4) +
  geom_point(aes(re,im), size = 0.05, color = color3, alpha = .5) +
  geom_point(data = eigmayor, aes(re,im), size = 0.05, color = color4, alpha = .5) +
  geom_point(data = eigRMT, aes(re,im), size = 0.05, color = color0, alpha = .5) +#+ xlim(c(-0.7,-0.6))
  labs(x="Real part", y="Imaginary part") +
  theme_bw()

ggsave(file="BGvsRMT.svg", plot=BGvsRMT, width=10, height=8)

# prueba <- data.frame(param = numeric(0), eigen_dist = numeric(0))
# ind <- 0
# for (mub in seq(0,1,by = .05)) {
#   if ((Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)) {
#   
#   ind <- ind+1
#   sig <- sqrt((mub*sw)^2 + sc^2)
#   rho <- (rhow*((mub*sw)^2) + rhoc*sc^2)/sig^2
#   Gamma <- (Gammaw*((mub*sw)^2) + Gammac*sc^2)/sig^2
#   
#   outlierRMT <- (N-1)*mub*muw + mub - mud - mua - mudel
#   outlierBG <- (N-1)*mub*muw + mub - mud - mua - mudel +
#     (N/2)*(mub*muw+muc)*(1+rho/Gamma)*(sqrt(1+(4*Gamma*sig^2/(N*(mub*muw+muc)^2)))-1)
#   
#   prueba[ind,] <- c(mub,abs(outlierRMT-outlierBG))
#   }
# }
# ggplot(data = prueba, aes(x = param, y = eigen_dist)) + geom_point() #+ coord_fixed()

####### CORRELATIONS: MOBILITY NETWORKS #################

N <- 30

# positive rhow vs negative rhow
rhow <- .85
rhow <- -.85
COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING)

ggsave(file="mobnet_negrho.svg")

# high rw/cw (exchange values)
muw <- 0.2
sw <- 0.06
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- .15 #gamma of baron et al
rw <- .7*N
cw <- .1*N
(Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)[[1]]
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING)

ggsave(file="mobnet_posr.svg")

# positive vs negative Gammaw
muw <- 0.25
sw <- 0.1
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- .22*N #gamma of baron et al
Gammaw <- -.22*N #gamma of baron et al
rw <- .25*N
cw <- .25*N
(Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)[[1]]
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING)

ggsave(file="mobnet_negG.svg")

# "high" Gammaw and rhow

muw <- 0.25
sw <- 0.1
rhow <- 0.85 #original rho (Gamma of baron et al)
Gammaw <- .4*N #gamma of baron et al
rw <- .45*N #positive
cw <- .45*N #positive
(abs(Gammaw)/sqrt(rw*cw) < 1) & (rw+cw<N) & (abs(N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

# negative values generate a network with the same properties

muw <- 0.25
sw <- 0.1
rhow <- -0.85 #original rho (Gamma of baron et al)
Gammaw <- -.4*N #gamma of baron et al
rw <- .45*N #positive
cw <- .45*N #positive
(abs(Gammaw)/sqrt(rw*cw) < 1) & (rw+cw<N) & (abs(N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)[[1]]
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING)

# preventing disease

muw <- 0.25
sw <- 0.1
rhow <- 0.35 #original rho (Gamma of baron et al)
Gammaw <- -.1*N #gamma of baron et al
rw <- .15*N #positive
cw <- .15*N #positive
(abs(Gammaw)/sqrt(rw*cw) < 1) & (rw+cw<N) & (abs(N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

COMMUTING <- rand_mat(N, muw, sw, distrib = "beta")
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)[[1]]
diag(COMMUTING) <- rep(0,N)
plotmobility(COMMUTING)

####### CORRELATIONS: PANEL #################

N <- 200

mub <- 0.02
betas <- rep(mub, N) # transmission rates
mud <- 0.1282
deaths <- rep(mud, N)
mua <- 0.25
alphas <- rep(mua, N)
mudel <- 0
deltas <- rep(mudel, N)

# re-scaling all parameters: ellipse is lost, outlier is still accurate
muw <- 0.09
sw <- 0.06
rhow <- 0.45
Gammaw <- 0.3*N
rw <- 0.4*N
cw <- 0.4*N

muc <- muw/100
sc <- sw/120
rhoc <- rhow/100
Gammac <- Gammaw/100
rc <- rw/100
cc <- cw/100

COMMUTING <- rand_mat_cor_beta(N, N*muw, sqrt(N)*sqrt(N)*sw, rhow, Gammaw, rw, cw)[[1]]
MIGRATION <- rand_mat_cor_beta(N, N*muc, sqrt(N)*sqrt(N)*sc, rhoc, Gammac, rc, cc)[[1]]
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigmayor <- eigen_mat(jacobian)
redmayor <- plotmobility(COMMUTING[1:60,1:60])

COMMUTING <- rand_mat_cor_beta(N, N*muw, sqrt(N)*sqrt(N)*sw, -rhow, -Gammaw, rw, cw)[[1]]
MIGRATION <- rand_mat_cor_beta(N, N*muc, sqrt(N)*sqrt(N)*sc, -rhoc, -Gammac, rc, cc)[[1]]
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigmenor <- eigen_mat(jacobian)
redmenor <- plotmobility(COMMUTING[1:60,1:60])

COMMUTING <- rand_mat_cor_beta(N, N*muw, sqrt(N)*sqrt(N)*sw, rhow, 0.001, 0.001, 0.001)[[1]]
MIGRATION <- rand_mat_cor_beta(N, N*muc, sqrt(N)*sqrt(N)*sc, rhoc, 0.00001, 0.00001, 0.00001)[[1]]
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
eigRMT <- eigen_mat(jacobian)
reduncor <- plotmobility(COMMUTING[1:60,1:60])

sig <- sqrt((mub*sw)^2 + sc^2)
rho <- (rhow*((mub*sw)^2) + rhoc*sc^2)/sig^2
Gamma <- (Gammaw*((mub*sw)^2) + Gammac*sc^2)/sig^2

outlierREmayor <- max(eigmayor$re)
outlierREmenor <- max(eigmenor$re)
outlierRMT <- (N-1)*mub*muw + mub - mud - mua - mudel
outlierBGmayor <- (N-1)*mub*muw + mub - mud - mua - mudel +
  (N/2)*(mub*muw+muc)*(1+rho/Gamma)*(sqrt(1+(4*Gamma*sig^2/(N*(mub*muw+muc)^2)))-1)
outlierBGmenor <- (N-1)*mub*muw + mub - mud - mua - mudel +
  (N/2)*(mub*muw+muc)*(1+rho/Gamma)*(sqrt(1-(4*Gamma*sig^2/(N*(mub*muw+muc)^2)))-1)

center <- mub*(1-muw)-mua-mud-mudel-N*muc
radius <- sqrt(N)*sig

color0 <- "black"
color1 <- "red"
color2 <- "blue"
color3 <- "green"
color4 <- "purple"

BGvsRMT <- ggplot(data = eigmenor) + 
  geom_ellipse(aes(x0 = center, y0 = 0, a = (1+rho)*radius,
                   b = (1-rho)*radius, angle = 0), color = color2) +
  geom_vline(xintercept = 0, color = color2, linetype = "dashed") +
  geom_point(aes(x = outlierRMT, y = 0), size = 3, color = color1, shape = 15) +
  geom_point(aes(x = outlierBGmenor, y = 0), size = 3, color = color3, shape = 16) +
  geom_point(aes(x = outlierBGmayor, y = 0), size = 3, color = color4, shape = 17) +
  geom_point(aes(re,im), size = 1.7, color = color3, shape = 16) +  #, alpha = .4
  geom_point(data = eigmayor, aes(re,im), size = 1.7, color = color4, shape = 17) +
  geom_point(data = eigRMT, aes(re,im), size = 1.7, color = color1, shape = 15) +#+ xlim(c(-0.7,-0.6))
  labs(x="Real part", y="Imaginary part") +
  theme_bw()
BGvsRMT
ggsave(file="BGvsRMT.svg", plot=BGvsRMT)

ggsave(file="redmenor.svg", plot=redmenor)
ggsave(file="redmayor.svg", plot=redmayor)
ggsave(file="reduncor.svg", plot=reduncor)

corpanel <-ggarrange(ggarrange(redmenor, reduncor, redmayor, labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1),
          BGvsRMT,
          ncol = 1, nrow = 2)
corpanel

ggsave(file="corpanel.svg", plot=corpanel)
ggsave(corpanel, file="corpanel.eps", device = cairo_ps, fallback_resolution = 600)

