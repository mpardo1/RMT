##### PLOT MOBILITY NETWORKS WITH TWO MOVEMENT SCALES #####

library(tidyverse)

###########################################################

plotmobility <- function(mob, color1 = "#0000FF"){
  
  N <- nrow(mob)
  if (is.matrix(mob) & (nrow(mob) == ncol(mob))) {
    
    diag(mob) <- rep(0,N)
    
    mob <- mob[seq(N,1,-1),]
    
    mob <- as.data.frame(mob)
    names(mob) <- c(1:N)
    mob$x <- c(1:N)
    mob <- mob %>% pivot_longer(as.character(c(1:N)), names_to = "y", values_to = "mob")
    
    ggplot(mob, aes(x,y, fill=mob)) +
      geom_tile() +
      scale_alpha(range = c(0.1, 1)) +
      scale_fill_gradient(low = "white", high = color1, na.value = "gray90") +
      theme_void() +
      theme(legend.position="none", panel.background=element_blank()) +
      coord_fixed()
    
  } else {
    print("mob needs to be a square matrix")
  }
}

plotmobility2 <- function(mig, comm, color1 = "#0000FF", color2 = "#FF0000"){
  
  if (is.matrix(mig) & is.matrix(comm) & (nrow(mig) == ncol(mig)) &
      (ncol(mig) == nrow(comm)) & (nrow(comm) == ncol(comm))) {
    
    mig <- mig[seq(N,1,-1),]
    comm <- comm[seq(N,1,-1),]
    
    diag(mig) <- diag(comm) <- rep(0,N)
    meanmig <- sum(mig)/(N^2-N)
    meancomm <- sum(comm)/(N^2-N)
    diag(mig) <- rep(meanmig,N)
    diag(comm) <- rep(meancomm,N)
    
    mig <- as.data.frame(mig)
    comm <- as.data.frame(comm)
    names(mig) <- names(comm) <- c(1:N)
    mig$x <- comm$x <- c(1:N)
    mig <- pivot_longer(mig, c(1:N), names_to = "y", values_to = "mig")
    comm <- pivot_longer(comm, c(1:N), names_to = "y", values_to = "comm")
    mob <- inner_join(mig, comm, by = c("x" = "x", "y" = "y")) %>%
      mutate(x = as.numeric(x), y = as.numeric(y))
    
    ggplot(mob, aes(x,y, fill=mig, alpha=comm)) +
      geom_tile() +
      scale_alpha(range = c(0.1, 1)) +
      scale_fill_gradientn(colours = c(color1, color2), na.value = "gray90") +
      theme_void() +
      theme(legend.position="none", panel.background=element_blank()) +
      coord_fixed()
    
  } else {
    print("mig and comm need to be square matrices of the same size")
  }
}

##### LEGEND ##############################################

color1 <- "#0000FF"
color2 <- "#FF0000"

#continuous scale
legend <-expand.grid(x=1:100,y=1:100)
ggplot(legend, aes(x,y, fill = atan(y/x), alpha = (x+y-2)/4)) +
  geom_tile() + 
  scale_fill_gradientn(colours = c(color1, color2), na.value = "gray90") +
  scale_alpha(range = c(0.1,1)) +
  theme_void() +
  theme(legend.position="none", panel.background=element_blank()) +
  coord_fixed()

#discrete scale
legend <- data.frame(x = c(rep(1,3),rep(2,3),rep(3,3)), y = rep(c(1,2,3),3),
                     mig = c(rep(0.25,3),rep(0.5,3),rep(0.75,3)),
                     comm = rep(c(0.25,0.5,0.75),3))
ggplot(legend, aes(x,y, fill = atan(y/x), alpha = (x+y-2)/4)) +
  geom_tile() + 
  scale_fill_gradientn(colours = c(color1, color2), na.value = "gray90") +
  scale_alpha(range = c(0.1,1)) +
  theme_void() +
  theme(legend.position="none", panel.background=element_blank()) +
  coord_fixed()
