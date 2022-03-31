data <- matrix(scan("actg181.txt", quiet=TRUE), nc=5, byr=TRUE)

x1 <- rep(data[,1], data[,5])
x2 <- rep(data[,2], data[,5])
y1 <- rep(data[,3], data[,5])
y2 <- rep(data[,4], data[,5])

actg181 <- cbind(x1,x2,y1,y2)
dimnames(actg181) <- list(NULL, c("x1","x2","y1","y2"))

rm(data,x1,x2,y1,y2)

