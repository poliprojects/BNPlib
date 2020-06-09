library("HistogramTools")
require("VGAM")
require("metRology")
require("mvtnorm")
require("rgl")
library("plotly")
library("HistogramTools")
library("plot3D")

dens1 <- read.csv(file = 'src/python/test_res/dens1.csv', header = FALSE)
dens2 <- read.csv(file = 'src/python/test_res/dens2.csv', header = FALSE)
dens3 <- read.csv(file = 'src/python/test_res/dens3.csv', header = FALSE)
dens4 <- read.csv(file = 'src/python/test_res/dens4.csv', header = FALSE)
dens5 <- read.csv(file = 'src/python/test_res/dens5.csv', header = FALSE)


true_dens1 <- 0.5*dnorm(dens1[,1], mean=-3, sd=1)+0.5*dnorm(dens1[,1], mean=3, sd=1)
true_dens2<-0.9*dnorm(dens2[,1], mean=-5, sd=1)+0.1*dnorm(dens2[,1], mean=5, sd=1)
true_dens3<-0.3*dnorm(dens3[,1], mean=-2, sd=0.8)+0.3*dnorm(dens3[,1], mean=0, sd=0.8)+0.4*dnorm(dens3[,1], mean=2, sd=1)
true_dens4<-0.5*dt.scaled(dens4[,1],5,-5,1)+0.5*dskewnorm(dens4[,1],5,1,2)
true_dens5 = 0.5*dmvnorm(dens5[,c(1,2)], rep(-3,2),  diag(2))+0.5*dmvnorm(dens5[,c(1,2)], rep(3,2),  diag(2))


data1 <- as.matrix(read.csv(file = 'csv/test/data1.csv', header = FALSE))
data2 <- as.matrix(read.csv(file = 'csv/test/data2.csv', header = FALSE))
data3 <- as.matrix(read.csv(file = 'csv/test/data3.csv', header = FALSE))
data4 <- as.matrix(read.csv(file = 'csv/test/data4.csv', header = FALSE))
data5 <- as.matrix(read.csv(file = 'csv/test/data5.csv', header = FALSE, sep=" "))


# sd=beta/(alfa-1)=2
prior1 <- dnorm(dens1[,1], mean=mean(data1), sd=sqrt(2))
prior2<- dnorm(dens2[,1], mean=mean(data2), sd=sqrt(2))
prior3<- dnorm(dens3[,1], mean=mean(data3), sd=sqrt(2))
prior4<- dnorm(dens4[,1], mean=mean(data4), sd=sqrt(2))

############################
#TEST1
h=hist(data1)
pdf(file = "R/density_tests/test1.pdf", width = 6, height = 6)
PlotRelativeFrequency(h,main="Posterior Estimate",
                      xlab="data",ylim=c(0,0.4))
lines(dens1[,1],dens1[,2],  lwd=2)
lines(dens1[,1],true_dens1,  lwd=2,col="red")
lines(dens1[,1],prior1,  lwd=1,col="gray")
legend("topright", 
       c("Prior Density","Posterior Density", "True Density"), 
       lty=c(1, 1, 1), 
       col=c("gray","black","red"),
       bty = "n")
dev.off()



############################
#TEST2

h=hist(data2, breaks=20)
pdf(file = "R/density_tests/test2.pdf", width = 6, height = 6)
PlotRelativeFrequency(h,main="Posterior Estimate",
                      xlab="data", ylim=c(0,0.4),xlim=c(-10,10))
lines(dens2[,1],dens2[,2],  lwd=2)
lines(dens2[,1],true_dens2,  lwd=2,col="red")
lines(dens2[,1],prior2,  lwd=1,col="gray")
legend("topright", 
       c("Prior Density","Posterior Density", "True Density"), 
       lty=c(1, 1, 1), 
       col=c("gray","black","red"),
       bty = "n")
dev.off()


############################
#TEST3

h=hist(data3)
pdf(file = "R/density_tests/test3.pdf", width = 6, height = 6)
PlotRelativeFrequency(h,main="Posterior Estimate",
                      xlab="data",ylim=c(0,0.4), xlim=c(-6,6))
lines(dens3[,1],dens3[,2],  lwd=2)
lines(dens3[,1],true_dens3,  lwd=2,col="red")
lines(dens3[,1],prior3,  lwd=1,col="gray")
legend("topright", 
       c("Prior Density","Posterior Density", "True Density"), 
       lty=c(1, 1, 1), 
       col=c("gray","black","red"),
       bty = "n")
dev.off()


############################
#TEST4

h=hist(data4, breaks=20)
pdf(file = "R/density_tests/test4.pdf", width = 6, height = 6)
PlotRelativeFrequency(h,main="Posterior Estimate",
                      xlab="data", ylim=c(0,0.4), xlim=c(-10,10))
lines(dens4[,1],dens4[,2],  lwd=2)
lines(dens4[,1],true_dens4,  lwd=2,col="red")
lines(dens4[,1],prior4,  lwd=1,col="gray")
legend("topright", 
       c("Prior Density","Posterior Density", "True Density"), 
       lty=c(1, 1, 1), 
       col=c("gray","black","red"),
       bty = "n")
dev.off()



############################
#TEST5


# Histogram
x_c <- cut(data5[,1],50)
y_c <- cut(data5[,2],50)
z <- table(x_c, y_c)

x11()
hist3D(z=z, border="black")


# Posterior
grid=seq(-7,7.1, length.out=50)

open3d()
rgl::plot3d(x = dens5[,1], y = dens5[, 2], z =dens5[, 3], col="blue", xlab = "X1", ylab = "X2", zlab = "Y")
rgl::surface3d(x=grid, y =grid,
               z = matrix(dens5[,3], nrow = 50, ncol = 50),
               col = "green", alpha = 0.25, lit = FALSE)


axx <- list(
        title = "X-data"
)

axy <- list(
        title = "Y-data"
)

axz <- list(
        title = "Relative Frequency"
)

p<-plot_ly() %>% add_surface(x = grid, y = grid, z =  matrix(dens5[,3], nrow = 50, ncol = 50)) %>% layout(title = "Posterior Estimate",scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
        