set.seed(20200229)

# Test 1
test1 = c(rnorm(100, -3, 1), rnorm(100, +3, 1))
true_clust1 = c(rep(0,100), rep(1,100))
write.table(test1, "csv/test/data1.csv", row.names = F, col.names = F)
write.table(true_clust1, "csv/test/true_clust1.csv", row.names = F,
    col.names = F)

# Test 2
test2 = c(rnorm(900, -5, 1), rnorm(100, +5, 1))
true_clust2 = c(rep(0,900), rep(1,100))
write.table(test2, "csv/test/data2.csv", row.names = F, col.names = F)
write.table(true_clust2, "csv/test/true_clust2.csv", row.names = F,
    col.names = F)

# Test 3
test3 = c(rnorm(60, -2, 0.8^2), rnorm(60, 0, 0.8^2), rnorm(80, +2, 1))
true_clust3 = c(rep(0,60), rep(1,60), rep(2,80))
write.table(test3, "csv/test/data3.csv", row.names = F, col.names = F)
write.table(true_clust3, "csv/test/true_clust3.csv", row.names = F,
    col.names = F)

# Test 4
## Include required libraries
if(!require("VGAM"))      install.packages("VGAM")
if(!require("metRology")) install.packages("metRology")
if(!require("mvtnorm"))   install.packages("mvtnorm")
require("VGAM")
require("metRology")
require("mvtnorm")
## Generate data
test4 = c(rt.scaled(200,5,-5,1), rskewnorm(200,5,1,2))
true_clust4 = c(rep(0,200), rep(1,200))
write.table(test4, "csv/test/data4.csv", row.names = F, col.names = F)
write.table(true_clust4, "csv/test/true_clust4.csv", row.names = F,
    col.names = F)

# Test 5-6
#d = c(2, 5, 10, 20)
d = c(2, 5)
for(i in seq(length(d))){
    test = rbind( rmvnorm(200, rep(-3,d[i]),  diag(d[i])),
        rmvnorm(200, rep(3,d[i]),  diag(d[i])) )
    write.table( test, paste("csv/test/data", i+4, ".csv", sep=""),
        row.names = F, col.names = F )
    true_clust = c(rep(0,200), rep(1,200))
    write.table(true_clust, paste("csv/test/true_clust", i+4, ".csv", sep=""),
        row.names = F, col.names = F)
}
