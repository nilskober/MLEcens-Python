cosmesis <- matrix(scan("cosmesis.txt", quiet=TRUE), nc=3, byr=TRUE)
dimnames(cosmesis) <- list(NULL, c("x1","x2","tr"))

