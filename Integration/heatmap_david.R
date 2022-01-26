
library(tidyverse)
N <- 60
# create migration matrix (to be replaced by our random matrices)
mig <- matrix(runif(N^2, min = 0, max = 1), nrow = N)
diag(mig) <- rep(0,N)
mig <- as.data.frame(mig)
names(mig) <- c(1:N)
mig$x <- c(1:N)
mig <- pivot_longer(mig, c(1:N), names_to = "y", values_to = "mig")
# create commuting matrix (to be replaced by our random matrices)
comm <- matrix(runif(N^2, min = 0, max = 1), nrow = N)
diag(comm) <- rep(0,N) 
comm <- as.data.frame(comm)
names(comm) <- c(1:N)
comm$x <- c(1:N)
comm <- pivot_longer(comm, c(1:N), names_to = "y", values_to = "comm")
# store both matrices in a single dataframe for plotting
mob <- inner_join(mig, comm, by = c("x" = "x", "y" = "y")) %>%
  mutate(x = as.numeric(x), y = as.numeric(y))
# plot grid
ggplot(mob, aes(x,y, fill=mig, alpha=comm)) +
  geom_tile() +
  scale_alpha(range = c(0.1, 1)) +
  scale_fill_gradientn(colours = c("#0000FF", "#FF0000"), na.value = "gray90") +
  theme_void() +
  theme(legend.position="none", panel.background=element_blank())
# plot legend
legend <-expand.grid(x=1:100,y=1:100)
ggplot(legend, aes(x,y, fill = atan(y/x), alpha = (x+y-2)/4)) +
  geom_tile() + 
  scale_fill_gradientn(colours = c("#0000FF", "#FF0000"), na.value = "gray90") +
  scale_alpha(range = c(0.1,1)) +
  theme_void() +
  theme(legend.position="none", panel.background=element_blank())
# plot legend #2
legend <- data.frame(x = c(rep(1,3),rep(2,3),rep(3,3)), y = rep(c(1,2,3),3),
                     mig = c(rep(0.25,3),rep(0.5,3),rep(0.75,3)),
                     comm = rep(c(0.25,0.5,0.75),3))
ggplot(legend, aes(x,y, fill = atan(y/x), alpha = (x+y-2)/4)) +
  geom_tile() + 
  scale_fill_gradientn(colours = c("#0000FF", "#FF0000"), na.value = "gray90") +
  scale_alpha(range = c(0.1,1)) +
  theme_void() +
  theme(legend.position="none", panel.background=element_blank())
Contraer

