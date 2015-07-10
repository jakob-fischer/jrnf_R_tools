# author: jakob fischer (jakob@automorph.info)
# date: 9. July 2015
# description: 
# little helper functions that implement general functions that are not sufficiently
# available through the language


# 'acc' was a self implemented version with the functionality of cumsum before
# TODO If all code using "acc" has been changed to "cumsum", remove

acc <- cumsum


# This function removes all the duplicates / permutations if one checks for 
# subisomorphisms with a directed ring with the igraph method 
# "graph.get.subisomorphisms.vf2(g, graph.ring(n, directed=T))".
# This is done by simply requiring that the first edge in the sumisomorphism
# is the one with the lowest id.

subisomorphism_rm_permutation <- function(si) {
    # method that identifies the right subisomorphisms to keep
    is_first_min <- function(x)  {  return(x[1] == min(x))  }
    sel <- lapply(si, is_first_min)    # select by list apply
    return(si[unlist(sel)])            # return selected elements from list si
} 


# Function that samples species energies from normal distribution and 
# activation energies from a distribution of the form p(E) = 6/(pi^2*(exp(1/x)-1))

rplancklike <- function(N) {
    r <- c()    

    for(i in 1:N) {
        x <- runif(min=0, max=100,n=1)
        y <- (6/pi^2)/(x^3*(exp(1/x)-1))

        while(y < runif(min=0, max=1, n=1)) {
            x <- runif(min=0, max=100,n=1)
            y <- (6/pi^2)/(x^3*(exp(1/x)-1))
        }
        
        r <- c(r, x)
    }

    return(r)
}
