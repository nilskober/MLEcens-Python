data <- matrix(scan("menopause.txt", quiet=TRUE), nc=5, byr=TRUE)
# first column:  age at which women were observed
# second column: total number of women observed at that age
# third column:  number of women observed to have had operative menopause
# fourth column: number of women observed to have had natural menopause
# fifth column:  number of women observed to have had no menopause yet

# I order rows according to operative menopause, natural menopause, no menopause
x1 <- c(rep(rep(0,nrow(data)), data[,3]),
        rep(rep(0,nrow(data)), data[,4]),
        rep(data[,1], data[,5]))
x2 <- c(rep(data[,1], data[,3]), 
        rep(data[,1], data[,4]), 
        rep(rep(100,nrow(data)), data[,5]))
y1 <- c(rep(rep(.75,nrow(data)), data[,3]),
        rep(rep(1.75,nrow(data)), data[,4]),
        rep(rep(.75,nrow(data)), data[,5]))
y2 <- c(rep(rep(1.25,nrow(data)), data[,3]),
        rep(rep(2.25,nrow(data)), data[,4]),
        rep(rep(2.25,nrow(data)), data[,5]))
menopauseMod <- cbind(x1,x2,y1,y2)
dimnames(menopauseMod) <- list(NULL, c("x1","x2","y1","y2"))

rm(data,x1,x2,y1,y2)

