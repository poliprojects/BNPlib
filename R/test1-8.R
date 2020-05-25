set.seed(20200229)
test1 = 0.5*rnorm(200, -3, 1) + 0.5*rnorm(200, +3, 1)
write.table(test1, "csv/test/data1.csv", row.names = F, col.names = F)

test2 = 0.9*rnorm(1000, -5, 1) + 0.1*rnorm(1000, +5, 1)
write.table(test2, "csv/test/data2.csv", row.names = F, col.names = F)

test3 = 0.3*rnorm(200, -2, 0.8^2) + 0.3*rnorm(200, 0, 0.8^2) +
	0.4*rnorm(200, +2, 1)
write.table(test3, "csv/test/data3.csv", row.names = F, col.names = F)


if (!require("VGAM")) install.packages("VGAM")
if (!require("metRology")) install.packages("metRology")
if (!require("mvtnorm")) install.packages("mvtnorm")

require("VGAM")
require("metRology")
require("mvtnorm")

test4 = 0.5*rt.scaled(400, 5,-5,1) + 0.5*rskewnorm(400, 5, 1,2)
write.table(test4, "csv/test/data4.csv", row.names = F, col.names = F)


d=c(2, 5, 10, 20)
for(i in seq(length(d))){
  test=0.5*rmvnorm(400, rep(-3,d[i]),  diag(d[i]))+ 0.5*rmvnorm(400, rep(3,d[i]),  diag(d[i]))
  write.table(test, paste("csv/test/data",i+4,".csv", sep=""), row.names = F, col.names = F)
}

