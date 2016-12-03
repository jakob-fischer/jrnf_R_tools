# author: jakob fischer (jakob@automorph.info)
# description: 
# little helper functions that implement general functions that are not sufficiently
# available through the language.

sourced_tools <- T

library(pracma)   # for gcd function
library(igraph)   # graph library

# 'acc' was a self implemented version with the functionality of cumsum before
# TODO If all code using "acc" has been changed to "cumsum", remove

acc <- cumsum


# Sum function. All NA values are removed. If all values of the vector are NA
# the function returns NA instead of 0 (standard behaviour)
msum <- function(x) {
    if(all(is.na(x)))
        return(NA)

    return(sum(x, na.rm=T))
}


# Multiplicates the matrix 'm' with the vector 'v' 
# If a zero entry in the matrix is multiplied with a NA entry in the vector a 
# value of 0 is taken in the corresponding part of the sum 
mmul <- function(m, v) {
    if(ncol(m) != length(v)) {
        cat("mmul: dimension mismatch!\n")
        return(0)
    }

    res <- rep(0, nrow(m))

    for(i in 1:nrow(m))
        res[i] <- sum((m[i,]*v)[m[i,] != 0])
    
    return(res)
}


# Simple helper function to normalize a vector...

mt_norm <- function(x) {
    return(x/sum(x))
}


# inverse which - function that "inverts" the effect of the which function

invwhich <- function(x, N)  {  a <- rep(F, N);  a[x] <- T;  return(a)  }


# Removes all flags that indicate that files have been sourced into the current
# session by removing variables starting with "sourced_"

clear_sourced_flags <- function() {
    objs <- ls(pos =".GlobalEnv")
    rm(list= objs[grepl("^sourced_", objs)], pos=".GlobalEnv")
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
# activation energies from a distribution of the form p(E) = 6/(pi^2*(exp(1/x)-1)).
# The distribution is sampled in the area [0, d_max]. Each sample is found by
# sampling a x uniformly between [0,d_max] and comparing the value of the distribution 
# at x with a uniformly sampled value y between [0,1] and keeping it if y < P(x).

rplancklike <- function(N, d_max=100) {
    r <- c()    

    for(i in 1:N) {
        x <- runif(min=0, max=d_max,n=1)
        y <- (6/pi^2)/(x^3*(exp(1/x)-1))

        while(y < runif(min=0, max=1, n=1)) {
            x <- runif(min=0, max=d_max,n=1)
            y <- (6/pi^2)/(x^3*(exp(1/x)-1))
        }
        
        r <- c(r, x) # add found sample
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

    for(i in x)
        a <- gcd(a,i)             
           
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
    a <- get.adjacency(g)
    return(matrix(as.numeric(a), ncol=ncol(a)))
}


# Helper function that is given a (weighted) adjacency matrix <M_adj> and a number 
# of modules <N_mod> (has to be divisor of species number) and then calculates a
# modularity score.

get_amatrix_mod <- function(M_adj, N_mod) { 
    mod_size <- ncol(M_adj)/N_mod

    # Accumulate weight inside modules
    s <- 0.0
    for(i in 1:N_mod) {
        x <- (i-1)*mod_size+1:mod_size
        s <- s + sum(M_adj[x,x])/(mod_size**2)
    }

    return(s/(sum(M_adj)/(nrow(M_adj)**2)))
}


# Find a modular reordering (for given module number) by trying a big amount
# of random exchanges and testing if it improves modularity.
#
# parameters:
# <M_adj> - Adjacency matrix (now + colum number must match and be a multiple of <N_mod>)
# <N_mod> - Number of modules

find_modular_reordering <- function(M_adj, N_mod) {
    N <- ncol(M_adj)                     # number of species
    o <- 1:N                             # current reordering
    mod <- get_amatrix_mod(M_adj, N_mod) # modularity

    # sample big number of pairs to exchange
    r1 <- sample(N, N**2*10, replace=T) 
    r2 <- sample(N, N**2*10, replace=T)
   
    for(i in 1:(N**2*10*4)) {
      # calculate modified reordering <o_t>
      o_t <- o
      o_t[r2[i]] <- o[r1[i]]
      o_t[r1[i]] <- o[r2[i]]  

      # recalculate modularity for new reordering
      mod_t <- get_amatrix_mod(M_adj[o_t,o_t], N_mod)
      
      # keep if modularity is improved
      if(mod_t > mod) {
          mod <- mod_t
          o <- o_t
          cat("-> ", mod, "\n")
      }
    }

    return(o)
}


# Transfers between multidimensional (vector) coordinates and one dimensional
# The function creates a list / object with two methods to transform in both
# ways. Parameter <v> is a integer vector with every element containing the
# size of one dimension. Both spaces (linear, multidimensional) are 1-indexed.

b_mdim_transf <- function(v) {
    to_linear <- function(x) {
        if(length(v) == 1)
            return(x[1])
        else
            return(x[1] + v[1]*((b_mdim_transf(v[-1])$to_linear(x[-1]))-1))
    }

    from_linear <- function(x) {
        if(length(v) == 1)
            return(x)
        else
            return(c((x-1) %% v[1] + 1,
                     b_mdim_transf(v[-1])$from_linear((x-1)%/% v[1]+1)))
    }

    return(list(to_linear=to_linear, from_linear=from_linear))
}


# Function for generating matrices in which each element is it's respective row-id or column-id.
# This is useful if one want's to transform from selected parts of a matix to their indices).
rowid_m <- function(nrow, ncol)  { return(matrix(rep(1:nrow, ncol), ncol=ncol) )  }
colid_m <- function(nrow, ncol)  { return(t(matrix(rep(1:ncol, nrow), ncol=nrow)) )  }

# Given a id (<x>) and matrix dimension (<nrow>, <ncol>) the functions return the 
# row and column indices necessary to access the <x>th element in the matrixs linear
# representation.
mat_lin_to_rowid <- function(x, nrow, ncol)  {  return(as.vector(rowid_m(nrow, ncol))[x])  }
mat_lin_to_colid <- function(x, nrow, ncol)  {  return(as.vector(colid_m(nrow, ncol))[x])  }

# Given a boolean matrix <x> the functions return the row and column indices (vectors)
# to access all those matrix elements that are logical TRUE.
mat_sel_to_rowid <- function(x) {  return(as.vector(rowid_m(nrow=nrow(x), ncol=ncol(x)))[as.vector(x)])  }
mat_sel_to_colid <- function(x) {  return(as.vector(colid_m(nrow=nrow(x), ncol=ncol(x)))[as.vector(x)])  }
