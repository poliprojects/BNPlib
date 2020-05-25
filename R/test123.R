set.seed(20200229)
test1 = 0.5*rnorm(200, -3, 1) + 0.5*rnorm(200, +3, 1)
write.table(test1, "csv/test/data1.csv", row.names = F, col.names = F)

test2 = 0.9*rnorm(1000, -5, 1) + 0.1*rnorm(1000, +5, 1)
write.table(test2, "csv/test/data2.csv", row.names = F, col.names = F)

test3 = 0.3*rnorm(200, -2, 0.8^2) + 0.3*rnorm(200, 0, 0.8^2) +
	0.4*rnorm(200, +2, 1)
write.table(test3, "csv/test/data3.csv", row.names = F, col.names = F)
