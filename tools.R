acc <- function(x) {
    a <- 0
    res <- rep(0, length(x))

    for(i in 1:length(x)) {
       a <- a + x[i]
       res[i] <- a
    }

    return(res)
}

