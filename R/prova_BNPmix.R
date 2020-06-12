library(devtools)
install_github("rcorradin/BNPmix")

data1 <- as.matrix(read.csv(file = 'csv/test/data1.csv', header = FALSE))
dens1 <- read.csv(file = 'src/python/test_res/dens1.csv', header = FALSE)
data5 <- as.matrix(read.csv(file = 'csv/test/data5.csv', header = FALSE, sep=" "))
dens5 <- read.csv(file = 'src/python/test_res/dens5.csv', header = FALSE)


est_model <- BNPmix::MAR(data = data1, niter = 10, nburn = 1,grid = as.matrix(dens1[,1]), m0=0.0, k0=0.1,mass=1,a0=2.0,b0=2.0, hyper=FALSE,m1=0.0, s21=1.0, tau1=1.0, tau2=1.0, a1=1.0,b1=1.0)