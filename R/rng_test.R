half = 50
data = read.csv("data.csv", header=F)
data = as.numeric(data)
hist(data)
hist(data[0:(half-1)])
hist(data[half:(2*half-1)])
