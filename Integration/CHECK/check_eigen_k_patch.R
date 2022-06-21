rm(list=ls())

# Check outlier k

N <- 300
alp <- 1
muw <- 0.3
muc <- 0.01
mub <- 0.1
a <- alp
b <- alp*muw
c <- mub*muw + muc
k <-  2
mat <- matrix(c,N,N)
mat[,1:k] <- rep(b + c,N)
diag(mat) <- a + c

eig <- data.frame(re = Re(eigen(mat)$values), im=Im(eigen(mat)$values))

outl <- (1/2)*(2*a + (k-1)*b + N*c  + 
                 sqrt((k-1)^2*b^2 + (2*(k+1)*N - 4*k)*b*c + N^2*c^2))


ggplot(eig) + 
  geom_point(aes(outl,0), color = "red") +
  geom_point(aes(re,im), size = 0.1) 
   
