# author: jakob fischer (jakob@automorph.info)
# date: 9. July 2015
# description: 
# little helper functions that implement general functions that are not sufficiently
# available through the language


# function "acc" calculates the cummulative sum of numeric vectors

acc <- function(x) {
    a <- 0
    res <- rep(0, length(x))

    for(i in 1:length(x)) {
       a <- a + x[i]
       res[i] <- a
    }

    return(res)
}

