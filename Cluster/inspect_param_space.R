rm(list = ls())
library("tidyverse")
library("parallel")
#--------------------FUNCTIONS------------------------------------#
cond1 <- function(N, mu_w, gam, bet){
  cond <- (N-1)*mu_w - (gam/bet) + 1
  log <- TRUE
  if(cond > 0 ){
    log <- FALSE
  }
  return(log)
}

cond2 <- function(N, mu_w, mu_c, sig_w, sig_c, tau, gam, bet){
  sig_hat <- sqrt(sig_w^2 + 2*(tau/bet) + (sig_c/bet)^2 )
  cond <- mu_w + sqrt(N)*sig_hat - (gam/bet) + 1 - N*(mu_c/bet)
  log <- TRUE
  if(cond > 0 ){
    log <- FALSE
  }
  return(log)
}


#-------------------------STAB COND -------------------------------------
# Randomly generate parameters and check the stability conditions.
# If cond_1 or cond_2 = 0 then Unstable, stable otherwise
df.param <- data.frame(N = 0, mu_w= 0, mu_c= 0, sig_w= 0,
                      sig_c= 0, tau= 0, gam= 0, bet= 0,
                      cond_1 = 0, cond_2 = 0)

it <- 0
num_seed <- 10
mat.param <- matrix(0, ncol = 3, nrow = 1)
count = 1
N <- 50
mu_w <- 0.3
while( nrow(mat.param) < num_seed){
  rand_mu_c <- runif(1,0,0.3)
  rand_sig_c <- runif(1,0,0.25)
  rand_sig_w <- runif(1,0,0.25)
  if(rand_mu_c > rand_sig_c & mu_w > rand_sig_w ){
    vec_par <- c(rand_mu_c, rand_sig_c, rand_sig_w)
    mat.param <- rbind(mat.param, vec_par)
  }
}

Cores <- 1 #parallel::detectCores()

  rand_param <- function(k){
    N <- 50
    mu_w <- 0.3
    mu_c_rand <- mat.param[k,1]
    sig_w_rand <- mat.param[k,2]
    sig_c_rand <- mat.param[k,3]
    gam <- 1
    bet <- 0.02
    cond.stab.1 <- ifelse(cond1(N, mu_w, gam, bet) == TRUE, 1, 0)
    cond.stab.2 <- ifelse(cond2(N, mu_w, mu_c_rand, sig_w_rand,
                                sig_c_rand, 0, gam, bet) == TRUE, 1, 0)
    if(cond.stab.1 == 1 & cond.stab.2 == 0){
      vec <- c(N, mu_w, mu_c_rand,sig_w_rand,sig_c_rand, 0,  gam,
               bet, cond.stab.1, cond.stab.2)
      it <- it + 1
      print("it":it)
    }else{
      vec = rep(0, 10)
    }
    vec
  }

long_df <- 10
while(nrow(df.param) < long_df){
  parall <- do.call("rbind", mclapply(seq_len(num_seed), rand_param,
                             mc.cores = Cores, mc.preschedule = F))
  vec_sum <-rowSums(parall) 
  if(sum(vec_sum) > 0){
    print("Matrix with good param:")
    print(mat)
    mat <- parall[which(rowSums(parall) > 0), ]
    df.param <- rbind(df.param, mat)
  }
}

path <- paste0("~/RMT/Integration/param_space_",Sys.Date(),".csv")
write.csv(df.param,path, row.names = TRUE)