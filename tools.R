# author: jakob fischer (jakob@automorph.info)
# description: 
# little helper functions that implement general functions that are not sufficiently
# available through the language

sourced_tools <- T

library(pracma)   # for gcd function
library(igraph)

# 'acc' was a self implemented version with the functionality of cumsum before
# TODO If all code using "acc" has been changed to "cumsum", remove

acc <- cumsum


# Removes all flags that indicate that files have been sourced into the current
# session by removing variables starting with "sourced_"

clear_sourced_flags <- function() {
    for(i in ls()[grepl("sourced_", ls())])
        rm(i)
}


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


# After drawing a histogram with the rainbow(20) palette 
# this function returns a vector associating the same colors
# used in the histogramms to the values given.
#
# 'values' vector of values
# 'hi' object returned from the hist-call

get_color_by_hist <- function(hi, values) {
    d <- character()

    for(i in values) {
        if(is.na(i)) {
                d <- c(d, "#FFFFFFFF")
	    } else {
                lisa <- max(which(i > hi$breaks))
                d <- c(d, rainbow(20)[lisa])
	    }
    }

    return(d)
}


# This function calculates the greatest common divisor for all elements in 
# a vector. For this the gcd function of the package 'pracma' is used.

vec_gcd <- function(x) {
    a <- x[1]
    x <- x[-1]
    x <- x[x != 0]

    if(length(x) > 1) 
        for(i in 2:length(x))
            a <- gcd(a,x[i])             
           
    return(a)
}


# Helper-function that scales all rows in a matrix by different scalar values
# given by a vector.

scale_mat_rows <- function(mat, vec) {
    if(nrow(mat) != length(vec)) {
        cat("ERROR: scale_mat_rows: dimension mismatch!\n")
        return(mat)
    }    

    return(mat * matrix(rep(vec, ncol(mat)), ncol=ncol(mat)))
}


# Converts an igraph graph to an adjacency matrix...

graph_to_amatrix <- function(g) {
    N <- length(V(g))
    M <- length(E(g))

    x <- matrix(0, ncol=N, nrow=N)
    for(i in 1:M) {
        a <- ends(g, i)[1,1]
        b <- ends(g, i)[1,2]

        x[a,b] <- x[a,b] + 1
    }  

    return(x)
}
