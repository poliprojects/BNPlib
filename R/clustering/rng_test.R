half = 50
data = read.csv("data.csv", header=F)
data = as.numeric(data)
hist(data)
hist(data[1:half])
hist(data[(half+1):(2*half)])
