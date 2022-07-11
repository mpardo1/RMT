####### RANDOM MATRICES FOR METAPOPULATION MODELS #######
### 
### article plots
###
### plots for the article
### works over "RMT_parent.R"
###

####### BASE PARAMETERS #################################

# number of patches
N <- 100

sus_init <- rep(100000, N) # initial susceptibles
inf_init <- rep(100, N)    # initial infecteds
end_time <- 100

# epidemiological

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

muw <- 0.2
sw <- 0.05
rhow <- 0 #original rho (Gamma of baron et al)
Gammaw <- .15 #gamma of baron et al
rw <- .1
cw <- .3

(Gammaw/sqrt(rw*cw) < 1) & ((N*rhow-2*Gammaw)/(N-(rw+cw)) < 1)

muc <- 0.001
sc <- 0.0005
rhoc <- .001
Gammac <- .004
rc <- .03
cc <- .06

(Gammac/sqrt(rc*cc) < 1) & ((N*rhoc-2*Gammac)/(N-(rc+cc)) < 1)

COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)[[1]]
MIGRATION <- rand_mat_cor_beta(N, muc*N, sc*N, rhoc, Gammac, rc, cc)[[1]]
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
sum(colSums(MIGRATION) > 1)

# jacobian

jacobian <- (COMMUTING + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

####### PERTURBED TRANSMISSION AT SEVERAL PATCHES #######

mub <- 0.1
betas <- rep(mub, N)
alpha <- 3*mub
mug <- mud + mua + mudel

eig_comp <- data.frame(k = integer(0),
                       RMT = numeric(0),
                       LRPk = numeric(0),
                       LRP1 = numeric(0),
                       Real = numeric(0))
                  
for (k in seq(1,26,by = 1)) {
  
  betask <- betask1 <- betas
  betask[1:k] <- mub + alpha
  betask1[1] <- mub + k*alpha
  
  jacobian <- (COMMUTING + diag(N)) %*% diag(betask) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacobian1 <- (COMMUTING + diag(N)) %*% diag(betask1) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  eig_comp[k,] <- c(k,
                    (mub+alpha*k/N)*(muw*(N-1)+1)-mug,
                    N*(mub*muw+muc)/2 + alpha*((k-1)*muw+1)/2 +
                      sqrt((N*(mub*muw+muc))^2+(alpha^2)*((k-1)*muw+1)^2+2*alpha*(mub*muw+muc)*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1)))/2 +
                      mub*(1-muw) - mua - mud - mudel - N*muc,
                    max(eigen_mat(jacobian)$re),
                    max(eigen_mat(jacobian1)$re))
}

eig_comp %>% pivot_longer(cols = (2:ncol(eig_comp)), names_to = "prediction", values_to = "eigenvalue") %>%
  ggplot(aes(x = k, y = eigenvalue, color = prediction)) +
  geom_line(size = 1)

####### PERTURBED TRANSMISSION SCALING WITH N ###########

mub
eig_comp <- data.frame(k = numeric(0),
                       RMT = numeric(0),
                       LRP = numeric(0),
                       real = numeric(0))
alpha <- 0.05

#replace beta at a node by alpha = kbeta, for k = -1,0,1,2,...N

for (k in seq(1,2000)) {
  
  betasalpha <- betas
  betasalpha[1] <- mub + k*alpha
  
  jacobian <- (COMMUTING + diag(N)) %*% diag(betasalpha) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  #jacobiank <- (COMMUTING + diag(N)) %*% diag(betasalphak) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  eig_comp[k,] <- c(k*alpha,
                    mub*(1+k*alpha/N)*(muw*(N-1)+1)-mud-mua-mudel,
                    N*(mub*muw+muc)/2 + k*alpha/2 +
                      sqrt((N*(mub*muw+muc))^2 + (k^2)*alpha^2 + 2*k*alpha*(mub*muw+muc)*(N*1 + 2*(N-1)*(muw-1)))/2 +
                      mub*(1-muw) - mua - mud - mudel - N*muc,
                    #(N/2)*(mub*muw + muc) + (alpha/2)*(1 + (k-1)*muw) + mub*(1-muw) -
                    #  mua - mua - mudel - N*muc + (1/2)*sqrt(N^2*(mub*muw + muc)^2 + alpha^2*(1 + (k-1)*muw)^2 +
                    #                             2*alpha*(mub*muw + muc)*(N*(1+(k-1)*muw) + 2*(N-k)*(muw-1))),
                    max(eigen_mat(jacobian)$re))
}

eig_comp %>% select(-RMT) %>% pivot_longer(cols = (2:ncol(.)), names_to = "prediction", values_to = "eigenvalue") %>%
  ggplot(aes(x = k, y = eigenvalue, color = prediction)) +
  geom_line(size = 1) +
  labs(subtitle = paste0("beta = ",mub))

#con mub = 1 se nota mas, aun asi poco
ggplot(data = eig_comp, aes(x = k, y = LRP-real)) + geom_point()

####### MOBILITY: GENERALITIES ##########################

# pruebas

plotmobility(COMMUTING)
plotmobility(MIGRATION)
plotmobility2(MIGRATION, COMMUTING)

plot_eigen_rmt(jacobian,
               N,mub,mug = mud + mua + mudel,
               muw,sw,rhow,Gammaw,
               muc,sc,rhoc,Gammac,
               tau = 0) #+ xlim(c(-215,-210))

plot_int(N, sol, state = "INF")

# base
muw <- .2
sw <- .05
muc <- .001
sc <- .0005

# muc
muc <- .01
sc <- .005

# sigma
sw <- .1

COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)[[1]]
MIGRATION <- rand_mat_cor_beta(N, muc*N, sc*N, rhoc, Gammac, rc, cc)[[1]]
sum(colSums(MIGRATION) > 1)

sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)
plot_int(N, sol, state = "INF")

sol_df <-  as.data.frame(sol)
for(i in c(1:N)){
  colnames(sol_df)[i+1] <-  paste0("S",i)
  colnames(sol_df)[N+i+1] <-  paste0("I",i)
  colnames(sol_df)[2*N+i+1] <-  paste0("R",i)}

sol_df_I <- select(sol_df,starts_with("I"))

plot_int(N, sol, state = "INF") +labs(subtitle = paste0("total = ",round(max(rowSums(sol_df_I))),
                                                        " max = ",round(max(sol_df_I))))

# grafica areas estabilidad N vs muw

library(latex2exp)
library(ggstar)

step <- 0.005
N_vec <- seq(1,200,by = 1)
muw_vec <- seq(0.01,0.9,step)
mug <- mud + mua + mudel

mub <- .1
mub <- .5

df_sol <- data.frame(beta = 0, gamma = 0, N = 0, muw = 0, state = FALSE)
for(i in c(1:length(N_vec))){
  print(paste0("i : ", i,"/",length(N_vec)))
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
      state <- ifelse(((N_vec[i]-1)*muw_vec[j]) < ((mug/mub) - 1), TRUE, FALSE)
      df_sol <- rbind(df_sol, c(mub, mug, N_vec[i], muw_vec[j], state))
    }
  }
}

df_sol <- df_sol[-1,] %>% mutate(Stability = ifelse(.$state == TRUE, "Stable", "Unstable"))

plot_area <- ggplot(df_sol) +
  geom_point(aes(N,muw,colour = Stability)) + theme_bw() +
  scale_color_manual(values=c("#3066BE", "#A63446")) +
  ylab(TeX("$\\mu_w$")) +
  xlab(TeX("$N$")) +
  labs(title = paste0("beta = ",mub)) +
  #coord_fixed() +
  theme(text = element_text(size = 30), legend.position = "bottom") +
  theme(
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=5)))

plot_area

####### COMMUTING: DIRECTED VS RANDOM CONTROL ###########

# base parameters + reduced muw (to induce stability with control)
muw <- .061
sw <- .02

muwstar <- 0
COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)[[1]]

# control: two rows
f1 <- 30
f2 <- 55
COMMUTING[f1,c(1:(f1-1),(f1+1):N)] <- COMMUTING[f2,c(1:(f2-1),(f2+1):N)] <- muwstar

# control: two columns
c1 <- 45
c2 <- 80
COMMUTING[c(1:(c1-1),(c1+1):N),c1] <- COMMUTING[c(1:(c2-1),(c2+1):N),c2] <- muwstar

# control: one row and one column
fc <- 65
COMMUTING[fc,c(1:(fc-1),(fc+1):N)] <- COMMUTING[c(1:(fc-1),(fc+1):N),fc] <- muwstar

# control: random
COMMUTING[sample(c(1:N^2),2*(N-1))] <- muwstar


diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
plotmobility(COMMUTING)
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

# now in a loop for the graph

r0df <- data.frame("muwstar" = vector(),
                   "r0_ff" = vector, "r0_cc" = vector(), "r0_fc" = vector(), "r0_r" = vector(),
                   "r0_ffn" = vector, "r0_ccn" = vector(), "r0_fcn" = vector(), "r0_rn" = vector())

for (muwstar in seq(0.01,0.99,by = .005)) {
  
  COMMUTING <- rand_mat_cor_beta(N, muw*N, sw*N, rhow, Gammaw, rw, cw)[[1]]
  diag(COMMUTING) <- rep(0,N)
  COMFF <- COMCC <- COMFC <- COMR <- COMMUTING
  
  # control: two rows
  f1 <- sample(c(2:(N-1)),1)
  f2 <- sample(c(2:(N-1)),1)
  COMFF[f1,c(1:(f1-1),(f1+1):N)] <- COMFF[f2,c(1:(f2-1),(f2+1):N)] <- muwstar
  
  # control: two columns
  c1 <- sample(c(2:(N-1)),1)
  c2 <- sample(c(2:(N-1)),1)
  COMCC[c(1:(c1-1),(c1+1):N),c1] <- COMCC[c(1:(c2-1),(c2+1):N),c2] <- muwstar
  
  # control: one row and one column
  fc <- sample(c(2:(N-1)),1)
  COMFC[fc,c(1:(fc-1),(fc+1):N)] <- COMFC[c(1:(fc-1),(fc+1):N),fc] <- muwstar
  
  # control: random
  COMR[sample(c(1:N^2),2*(N-1))] <- muwstar
  
  jacff <- (COMFF + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jaccc <- (COMCC + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacfc <- (COMFC + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacr <- (COMR + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  muwr <- ((N^2-N-2*(N-1))*muw + 2*(N-1)*muwstar)/(N^2-N) 
  
  r0df <- rbind(r0df,c(muwstar,
                       mub*(muw*(N-1)/2 + muwstar/2 +sqrt((muw^2*((N-3)^2))/4 + muw*muwstar*(3*N-5)/2 + muwstar^2/4) + 1-muw) - mud - mua - mudel,
                       mub*(muw*(N-1)/2 + muwstar/2 +sqrt((muw^2*((N-3)^2))/4 + muw*muwstar*(3*N-5)/2 + muwstar^2/4) + 1-muw) - mud - mua - mudel,
                       max(eigen_mat(jacfc)$re),
                       mub*(1+muwr*(N-1)) - mud - mua - mudel,
                       max(eigen_mat(jacff)$re), max(eigen_mat(jaccc)$re),
                       max(eigen_mat(jacfc)$re), max(eigen_mat(jacr)$re)))
  
  print(muwstar)
}

names(r0df) <- c("muwstar", "r0_ff", "r0_cc", "r0_fc", "r0_r", "r0_ffn", "r0_ccn", "r0_fcn", "r0_rn")
  
r0df <- pivot_longer(r0df, cols = (2:ncol(r0df)) , names_to = "prediction", values_to = "r0")

r0df %>%
  ggplot(aes(x = muwstar, y = r0, color = prediction)) +
    geom_line(size = 1) +
    coord_fixed()

filter(r0df, prediction %in% c("r0_ff","r0_cc","r0_ffn","r0_ccn")) %>%
  ggplot(aes(x = muwstar, y = r0, color = prediction)) +
  geom_line()

filter(r0df, prediction %in% c("r0_ff","r0_fc","r0_ffn","r0_fcn")) %>%
  ggplot(aes(x = muwstar, y = r0, color = prediction)) +
  geom_line()

filter(r0df, prediction %in% c("r0_ff","r0_ffn","r0_r","r0_rn")) %>%
  ggplot(aes(x = muwstar, y = r0, color = prediction)) +
  geom_line()

####### COMMUTING: DIRECTED PERTURBATIONS ###############

nu <- -muw
mug <- mud + mua + mudel

eig_comp <- data.frame(k = integer(0),
                       RMT = numeric(0),
                       LRP_in = numeric(0),
                       LRP_inout = numeric(0),
                       Real_r = numeric(0),
                       Real_in = numeric(0),
                       Real_inout = numeric(0))

for (k in seq(2,26,by = 2)) {
  
  COMM_IN <- COMM_INOUT <- COMM_R <- COMMUTING
  COMM_IN[,1:k] <- muw + nu
  COMM_INOUT[,1:(k/2)] <- COMM_INOUT[1:(k/2),] <- muw + nu
  COMM_R[sample(c(1:N^2),k*(N-1))] <- muw + nu
  diag(COMM_IN) <- diag(COMM_INOUT) <- diag(COMM_R) <- rep(0,N)
  
  jacobian_in <- (COMM_IN + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacobian_inout <- (COMM_INOUT + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacobian_r <- (COMM_R + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  eig_comp[k,] <- c(k,
                    mub*(N-1)*(muw+nu*k/N)+mub-mug,
                    N*(mub*muw+muc)/2 + (k-1)*mub*nu/2 +
                      sqrt((N*(mub*muw+muc))^2+((k-1)*mub*nu)^2+2*(mub*nu)*(mub*muw+muc)*(N+N*k-4*k))/2  +
                      mub*(1-muw) - mug - N*muc,
                    N*(mub*muw+muc)/2 + (k/2-1)*mub*nu/2 +
                      sqrt((N*(mub*muw+muc))^2+mub*nu*(2*mub*muw+muc)*(N+3*N*k/2-2*(k/2+1))+((mub*nu)^2)*(2*k*(N-k/2)+(k/2-1)^2))/2 +
                      mub*(1-muw) - mug - N*muc,
                    max(eigen_mat(jacobian_r)$re),
                    max(eigen_mat(jacobian_in)$re),
                    max(eigen_mat(jacobian_inout)$re))
}

eig_comp <- filter(eig_comp, !is.na(k))
eig_comp_real <- select(eig_comp, k, Real_r, Real_in, Real_inout) %>%
  pivot_longer(cols = (2:ncol(.)), names_to = "prediction", values_to = "eigenvalue")
eig_comp_pred <- select(eig_comp, k, RMT, LRP_in, LRP_inout) %>%
  pivot_longer(cols = (2:ncol(.)), names_to = "prediction", values_to = "eigenvalue")

eig_comp_pred %>% ggplot(aes(x = k, y = eigenvalue, color = prediction)) +
  geom_line(size = 1) + 
  geom_point(data = eig_comp_real, aes(x = k, y = eigenvalue, color = prediction)) +
  scale_color_manual(values=c("#0000FF","#FF0000","#0000FF","#FF0000","#63B159","#63B159"))

####### COMMUTING: DIRECTED CORRELATED PERTURBATIONS ####

nu <- -muw
mug <- mud + mua + mudel

eig_comp <- data.frame(k = integer(0),
                       random = numeric(0),
                       correlated = numeric(0),
                       real_r = numeric(0),
                       real_c = numeric(0))


for (k in seq(2,(N^2)/10,by = 2)) {
  
  COMM_R <- COMM_C <- COMMUTING
  COMM_R[sample(c(1:N^2),k)] <- muw + nu
  
  ind <- vector()
  for (K in c(2:N)) {
    ind <- c(ind, seq(((K-2)*N+K),(K-1)*N))
  }
  cor_ind <- sample(ind, k/2)
  COMM_C[cor_ind] <- muw + nu
  COMM_C <- t(COMM_C)
  COMM_C[cor_ind] <- muw + nu
  COMM_C <- t(COMM_C)
  
  jacobian_r <- (COMM_R + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  jacobian_c <- (COMM_C + diag(N)) %*% diag(betas) + MIGRATION - diag(colSums(MIGRATION) + deaths + alphas + deltas)
  
  eig_comp[k,] <- c(k,
                    mub*(N-1)*(muw+nu*k/(N^2-N))+mub-mug + 1/(N*mug),
                    mub*(N-1)*(muw+nu*k/(N^2-N))+mub-mug,
                    max(eigen_mat(jacobian_r)$re),
                    max(eigen_mat(jacobian_c)$re))
}

eig_comp <- filter(eig_comp, !is.na(k))
eig_comp %>%
  pivot_longer(cols = (2:ncol(.)), names_to = "prediction", values_to = "eigenvalue") %>% 
  ggplot(aes(x = k, y = eigenvalue, color = prediction)) +
  geom_line(size = 1)
  #scale_color_manual(values=c("#0000FF","#FF0000","#0000FF","#FF0000","#63B159","#63B159"))

####### COMMUTING: PANEL ################################
library("copula")
N <- 60

# para los plots
sus_init <- rep(100000, N) # initial susceptibles
inf_init <- runif(N, min = 1400,1600)  # initial infecteds
end_time <- 25

Deltas <- rep(0.1, N) # birth rate
mub <- 0.15
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
#gammas = deaths + alphas + deltas

muw <- 0.05
sw <- 0.05/3
rhow <- 0 #original rho (Gamma of baron et al)

muc <- 0.05
sc <- 0.05/4
rhoc <- .001

mub*muw*(N-1)+mub-mua-mud-mudel

COMMUTING <- rand_mat_ell(N, muw, sw, rhow, distrib = "beta")
MIGRATION <- rand_mat_ell(N, muc, sc, rhoc, distrib = "beta")
diag(COMMUTING) <- diag(MIGRATION) <- rep(0,N)
sol <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
           COMMUTING,MIGRATION,
           sus_init,inf_init,end_time)

### FUNCTION FOR PLOTS

plot_panelc <- function(N, sol, COMMUTING, colormean, color_low, color_high, IMIN = "N", IMAX = "N",CMMIN = "N", CMMAX = "N",
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
              color = colormean, size = 1.2, linetype = "dashed") +
    ylab("No infected individuals") +
    ylim(c(IMIN,IMAX)) +
    theme(legend.position = "none") +
    guides(color = "none") +
    theme_bw() + 
    theme(text = element_text(size = 15), legend.position = "bottom") 

}

### UNSTABLE NETWORK
library("copula")
library("viridis")
library("ggsci")
library("grid")
library("reshape")

plotmobility <- function(mat, low_col, high_col){
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
    ylim(c(60,0)) 
}

col_low <- "#90F0CE"
col_high <- "#063F95"
col_mean <- "#3C1642"
plot_stab_mob <- plotmobility(COMMUTING, col_low, col_high) + 
  theme(legend.position = "left") 
ggsave(file="unstnet.svg")
plot_stab_int <- plot_panelc(N, sol, COMMUTING,col_mean, col_low, col_high )
ggsave(file="unstsol.svg")

gg_stable <- ggarrange(plot_stab_mob,
          plot_stab_int)

gg_stable<- annotate_figure(gg_stable, 
                              top = textGrob("Estable scenario:",
                                             x = 0.13,
                                             gp = gpar(cex = 1.3)))
gg_stable
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

# STRATEGY B
cs <- sample(c(1:N), nodes)
COMMUTINGB[,fs] <- muwstar
solB <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGB,MIGRATION,
            sus_init,inf_init,end_time)

# STRATEGY C
fcs <- sample(c(1:N), nodes/2)
fcs <- fs[c(1:nodes/2)]
COMMUTINGC[fcs,] <- COMMUTINGC[,fcs] <- muwstar
solC <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGC,MIGRATION,
            sus_init,inf_init,end_time)

# STRATEGY D
inds <- sample(c(1:N^2), nodes*(N-1))
COMMUTINGD[inds] <- muwstar
solD <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGD,MIGRATION,
            sus_init,inf_init,end_time)

# STRATEGY E
rinds <- inds[c(1:(N^2)*nodes/(N-1)/2)]
COMMUTINGE[rinds] <- muwstar
COMMUTINGE <- t(COMMUTINGE)
COMMUTINGE[rinds] <- muwstar
COMMUTINGE <- t(COMMUTINGE)
solE <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGE,MIGRATION,
            sus_init,inf_init,end_time)

# STRATEGY F
COMMUTINGF <- COMMUTINGF*(1-nodes/N)
solF <- int(N, Deltas,betas,deaths,thetas,alphas,deltas,
            COMMUTINGF,MIGRATION,
            sus_init,inf_init,end_time)

CMAX <- max(c(COMMUTINGA,COMMUTINGB,COMMUTINGC,COMMUTINGD,COMMUTINGE, COMMUTINGF))
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

commpanel11 <-ggarrange(plotmobility(COMMUTINGA, col_low, col_high),
                        plot_panelc(N, solA, COMMUTINGA,col_mean, col_low, col_high, IMIN = IMIN, 
                                    IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
                          rremove("xlab") +  rremove("ylab"),
                        ncol = 1, nrow = 2, legend = "none")

commpanel11<- annotate_figure(commpanel11, 
                              top = textGrob("Scenario A:", x = 0.1, gp = gpar(cex = 1.3)),
                              left = textGrob("No Infected individuals",
                                              rot = 90, vjust = 1, y = 0.3, gp = gpar(cex = 1.3)),
                              bottom = textGrob("Time", gp = gpar(cex = 1.3)))
commpanel11

#--------#
commpanel12 <-ggarrange(plotmobility(COMMUTINGB, col_low, col_high),
                        plot_panelc(N, solB, COMMUTINGB,col_mean, col_low, col_high, IMIN = IMIN, 
                                    IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
                          rremove("xlab") +  rremove("ylab"),
                        ncol = 1, nrow = 2, legend = "none")

commpanel12<- annotate_figure(commpanel12, 
                              top = textGrob("Scenario B:", x = 0.1, gp = gpar(cex = 1.3)),
                              left = textGrob("No Infected individuals",
                                              rot = 90, vjust = 1, y = 0.3, gp = gpar(cex = 1.3)),
                              bottom = textGrob("Time", gp = gpar(cex = 1.3)))
commpanel12

#--------#
commpanel13 <-ggarrange(plotmobility(COMMUTINGC, col_low, col_high),
                        plot_panelc(N, solC, COMMUTINGC,col_mean, col_low, col_high, IMIN = IMIN, 
                                    IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
                          rremove("xlab") +  rremove("ylab"),
                        ncol = 1, nrow = 2, legend = "none")

commpanel13<- annotate_figure(commpanel13, 
                              top = textGrob("Scenario C:", x = 0.1, gp = gpar(cex = 1.3)),
                              left = textGrob("No Infected individuals",
                                              rot = 90, vjust = 1, y = 0.3, gp = gpar(cex = 1.3)),
                              bottom = textGrob("Time", gp = gpar(cex = 1.3)))
commpanel13

#--------#
commpanel21 <-ggarrange(plotmobility(COMMUTINGD, col_low, col_high),
                        plot_panelc(N, solD, COMMUTINGD,col_mean, col_low, col_high, IMIN = IMIN, 
                                    IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
                          rremove("xlab") +  rremove("ylab"),
                        ncol = 1, nrow = 2, legend = "none")

commpanel21<- annotate_figure(commpanel21, 
                              top = textGrob("Scenario D:", x = 0.1, gp = gpar(cex = 1.3)),
                              left = textGrob("No Infected individuals",
                                              rot = 90, vjust = 1, y = 0.3, gp = gpar(cex = 1.3)),
                              bottom = textGrob("Time", gp = gpar(cex = 1.3)))
commpanel21

#--------#
commpanel22 <-ggarrange(plotmobility(COMMUTINGE, col_low, col_high),
                        plot_panelc(N, solE, COMMUTINGE,col_mean, col_low, col_high, IMIN = IMIN, 
                                    IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
                          rremove("xlab") +  rremove("ylab"),
                        ncol = 1, nrow = 2, legend = "none")

commpanel22<- annotate_figure(commpanel22, 
                              top = textGrob("Scenario E:", x = 0.1, gp = gpar(cex = 1.3)),
                              left = textGrob("No Infected individuals",
                                              rot = 90, vjust = 1, y = 0.3, gp = gpar(cex = 1.3)),
                              bottom = textGrob("Time", gp = gpar(cex = 1.3)))
commpanel22

#--------#
commpanel23 <-ggarrange(plotmobility(COMMUTINGF, col_low, col_high),
                        plot_panelc(N, solF, COMMUTINGF,col_mean, col_low, col_high, IMIN = IMIN, 
                                    IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
                          rremove("xlab") +  rremove("ylab"),
                        ncol = 1, nrow = 2, legend = "none")

commpanel23<- annotate_figure(commpanel23, 
                              top = textGrob("Scenario F:", x = 0.1, gp = gpar(cex = 1.3)),
                              left = textGrob("No Infected individuals",
                                              rot = 90, vjust = 1, y = 0.3, gp = gpar(cex = 1.3)),
                              bottom = textGrob("Time", gp = gpar(cex = 1.3)))
commpanel23

#--------#
commpanel11 <-ggarrange(plotmobility(COMMUTINGA, col_low, col_high) + ggtitle("Scenario A"),                       plotmobility(COMMUTINGB, col_low, col_high) + ggtitle("Scenario B"), 
                      plotmobility(COMMUTINGC, col_low, col_high) + ggtitle("Scenario C"),
                      plot_panelc(N, solA, COMMUTINGA,col_mean, col_low, col_high, IMIN = IMIN, 
                                                           IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
                        rremove("xlab") +  rremove("ylab"),
                      plot_panelc(N, solB, COMMUTINGB,col_mean, col_low, col_high, IMIN = IMIN,
                                  IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX)+ 
                        rremove("xlab")+  rremove("ylab"),
                      plot_panelc(N, solC, COMMUTINGC,col_mean, col_low, col_high, IMIN = IMIN, 
                                  IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX)+ 
                        rremove("xlab")+  rremove("ylab"),
                      ncol = 3, nrow = 2, legend = "none")


commpanel11<- annotate_figure(commpanel11, 
                              left = textGrob("No Infected individuals",
                                              rot = 90, vjust = 1, y = 0.3,gp = gpar(cex = 1.3)),
                              bottom = textGrob("Time", gp = gpar(cex = 1.3)))

commpanel12 <-ggarrange(plotmobility(COMMUTINGD, col_low, col_high), 
                       plotmobility(COMMUTINGE, col_low, col_high), 
                       plotmobility(COMMUTINGF, col_low, col_high),
                       plot_panelc(N, solA, COMMUTINGA,col_mean, col_low, col_high, IMIN = IMIN, 
                                   IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX) + 
                         rremove("xlab") +  rremove("ylab"),
                       plot_panelc(N, solB, COMMUTINGB,col_mean, col_low, col_high, IMIN = IMIN,
                                   IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX)+ 
                         rremove("xlab")+  rremove("ylab"),
                       plot_panelc(N, solC, COMMUTINGC,col_mean, col_low, col_high, IMIN = IMIN, 
                                   IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX)+ 
                         rremove("xlab")+  rremove("ylab"),
                       ncol = 3, nrow = 2, legend = "none")

commpanel12<- annotate_figure(commpanel12, 
                              left = textGrob("No Infected individuals",
                                              rot = 90, vjust = 1, y = 0.3,gp = gpar(cex = 1.3)),
                              bottom = textGrob("Time", gp = gpar(cex = 1.3)))


commpanel <- ggarrange(commpanel11, commpanel12,
                       nrow = 2, ncol = 1)


commpanel <-ggarrange(plotmobility(COMMUTINGA, cmax = CMAX),
                      plotmobility(COMMUTINGB, cmax = CMAX),
                      plotmobility(COMMUTINGC, cmax = CMAX),
                      plotmobility(COMMUTINGD, cmax = CMAX),
                      plotmobility(COMMUTINGE, cmax = CMAX),
                      plot_panelc(N, solA, COMMUTINGA, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
                      plot_panelc(N, solB, COMMUTINGB, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
                      plot_panelc(N, solC, COMMUTINGC, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
                      plot_panelc(N, solD, COMMUTINGD, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
                      plot_panelc(N, solE, COMMUTINGE, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
                      ncol = 6, nrow = 2)

commpanel <-ggarrange(plotmobility(COMMUTINGA, cmax = CMAX),
                      plotmobility(COMMUTINGB, cmax = CMAX),
                      plotmobility(COMMUTINGC, cmax = CMAX),
                      plot_panelc(N, solA, COMMUTINGA, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
                      plot_panelc(N, solB, COMMUTINGB, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
                      plot_panelc(N, solC, COMMUTINGC, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
                      plotmobility(COMMUTINGD, cmax = CMAX),
                      plotmobility(COMMUTINGE, cmax = CMAX),
                      plotmobility(COMMUTINGF, cmax = CMAX),
                      plot_panelc(N, solD, COMMUTINGD, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
                      plot_panelc(N, solE, COMMUTINGE, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
                      plot_panelc(N, solF, COMMUTINGF, IMIN = IMIN, IMAX = IMAX, CMMIN = CMMIN, CMMAX = CMMAX),
                      ncol = 3, nrow = 4)

commpanel
ggsave(file="commpanel.svg", plot=commpanel)

### COMPARISON

eig_comp <- data.frame(k = integer(0),
                       LRP_in = numeric(0),
                       LRP_inout = numeric(0),
                       RMT_allr = numeric(0),
                       Real_in = numeric(0),
                       Real_out = numeric(0),
                       Real_inout = numeric(0),
                       Real_sparse = numeric(0),
                       Real_corr = numeric(0),
                       Real_dec = numeric(0))

nu <- muwstar - muw
mug <- mud + mua + mudel

for (nodes in seq(0,26,by = 2)) {
  
  COMMUTINGA <- COMMUTINGB <- COMMUTINGC <- COMMUTINGD <- COMMUTINGE <- COMMUTINGF <- COMMUTING
  fs <- sample(c(1:N), nodes)
  COMMUTINGA[fs,] <- muwstar
  COMMUTINGB[,fs] <- muwstar
  fcs <- sample(c(1:N), nodes/2)
  fcs <- fs[c(1:nodes/2)]
  COMMUTINGC[fcs,] <- COMMUTINGC[,fcs] <- muwstar
  inds <- sample(c(1:N^2), nodes*(N-1))
  COMMUTINGD[inds] <- muwstar
  rinds <- inds[c(1:(N^2)*nodes/(N-1)/2)]
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
                    N*(mub*muw+muc)/2 + (k-1)*mub*nu/2 +
                      sqrt((N*(mub*muw+muc))^2+((k-1)*mub*nu)^2+2*(mub*nu)*(mub*muw+muc)*(N+N*k-4*k))/2  +
                      mub*(1-muw) - mug - N*muc,
                    N*(mub*muw+muc)/2 + (k/2-1)*mub*nu/2 +
                      sqrt((N*(mub*muw+muc))^2+mub*nu*(2*mub*muw+muc)*(N+3*N*k/2-2*(k/2+1))+((mub*nu)^2)*(2*k*(N-k/2)+(k/2-1)^2))/2 +
                      mub*(1-muw) - mug - N*muc,
                    mub*(muw*(N-k)+muwstar*k/N)*(N-1)/N + mub - mug,
                    max(eigen_mat(jacA)$re),
                    max(eigen_mat(jacB)$re),
                    max(eigen_mat(jacC)$re),
                    max(eigen_mat(jacD)$re),
                    max(eigen_mat(jacE)$re),
                    max(eigen_mat(jacF)$re))
}

eig_comp <- filter(eig_comp, !is.na(k))

eig_comp_real <- select(eig_comp, k, Real_in, Real_out, Real_inout, Real_sparse, Real_corr, Real_dec) %>%
  pivot_longer(cols = (2:ncol(.)), names_to = "prediction", values_to = "eigenvalue")
eig_comp_pred <- select(eig_comp, k, LRP_in, LRP_inout, RMT_allr) %>%
  pivot_longer(cols = (2:ncol(.)), names_to = "prediction", values_to = "eigenvalue")

eig_comp_pred %>% ggplot(aes(x = k, y = eigenvalue, color = prediction)) +
  geom_line(size = 1) + 
  geom_point(data = eig_comp_real, aes(x = k, y = eigenvalue, color = prediction)) +
  scale_color_manual(values=c("#0000FF","#FF0000","#63B159","#0000FF","#0000FF","#FF0000","#63B159","#63B159","#63B159"))
ggsave(file="commparison.svg")

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

