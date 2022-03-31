data <- matrix(scan("actg181.txt", quiet=TRUE), nc=5, byr=TRUE)

x1 <- rep(data[,1], data[,5])
x2 <- rep(data[,2], data[,5])
y1 <- rep(data[,3], data[,5])
y2 <- rep(data[,4], data[,5])

actg181Mod <- cbind(x1-0.5,x2+0.5,y1-0.5,y2+0.5)
dimnames(actg181Mod) <- list(NULL, c("x1","x2","y1","y2"))

rm(data,x1,x2,y1,y2)

