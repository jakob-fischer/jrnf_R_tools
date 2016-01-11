# author: jakob fischer (jakob@automorph.info)
# description: 
# File contains code that is necessary for pathway analysis code (pathway_analysis.R)
# but generic enough to be usefull outside of the specific algorithm used there.

sourced_pathway_analysis_tools <- T

if(!exists("sourced_tools"))
    source("tools.R") 

if(!exists("sourced_cycles"))
    source("cycles.R")

if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")


# Function writes elementary mode <em> of network <net>. Reactions that have 
# only one educt and no products or only one product and no educts are ignored 
# (they are probably) pseudo-reactions for the output and not considered when 
# printing the net reaction. 

pa_write_em <- function(net, em, discard_s=T) {
    N <- jrnf_calculate_stoich_mat(net)   
    x <- (apply(N != 0, 2, sum) == 1)     # remove reactions with just one species involved
    if(discard_s)
        em[x] <- 0    

    # list all (remaining) reactions with coefficients 
    for(i in which(em != 0)) {
        cat(em[i], "X : ")
        cat(jrnf_reaction_to_string(net, i), "\n") 
    }

    cat("============================================================\n")
    # now print the net reaction    
    x <- N %*% em
    ed <- which(x < 0)    

    if(length(ed) > 0)
    for(i in 1:length(ed)) {
        if(-x[ed[i]] != 1)
            cat(-x[ed[i]], " ", sep="")

        cat(net[[1]]$name[ed[i]])

        if(i < length(ed))
            cat(" + ")
    }

    cat(" => ")

    prod <- which(x > 0)    

    if(length(prod) > 0)
    for(i in 1:length(prod)) {
        if(x[prod[i]] != 1)
            cat(x[prod[i]], " ", sep="")

        cat(net[[1]]$name[prod[i]])

        if(i < length(prod))
            cat(" + ")
    }

    cat("\n")    
}


# This function writes a list of <N> most relevant elementary modes /
# pathways to the file with filename <filename>.
#
# TODO: cleanup? Right place for this function?   

pa_write_em_set <- function(filename, net, em, exp_f, rates, N=100) {
    if(nrow(em) == 0)
        return()

    N <- min(N, nrow(em))
        

    em <- em[1:N,]
    exp_f <- exp_f[1:N]
    rates <- rates[1:N]

    x <- c()

    for(i in 1:nrow(em)) {
        x <- c(x, capture.output(cat(i, ": coefficient is", rates[i], "   explained fraction", exp_f[i])))
        x <- c(x, capture.output(cat("============================================================")))
        x <- c(x, capture.output(pa_write_em(net, em[i,])))
        x <- c(x, capture.output(cat("\n")))
    }

    writeLines(x, filename)
}


# The helper function takes a matrix (complete list of pathways) as argument
# and identifies those that are actually unique. 

pa_find_elementary_pw <- function(mat) {
    keep_pw <- rep(T, nrow(mat))         
    c_nzero_row <- apply(mat != 0, 1, sum)
        
    for(i in 1:nrow(mat)) {
        for(j in 1:nrow(mat)) 
            if(all(mat[i,] != 0 | mat[j,] == 0) & 
                   c_nzero_row[i] != c_nzero_row[j])
                    keep_pw[i] <- F
        }

    return(keep_pw)
}


# Function extends the pathway representation to one that does not only contain
# the (positive integer) coefficients of all the reactions but also integer 
# coefficients for the change of all species concentration by applying the pathway.
# TODO choose a better name?

pa_extend_pathways <- function(pw, N) {
    if(is.list(pw)) # if <pw> is not a matrix of pathways but a list of pathways
                    # $M and coefficients $coef...
        return(list(M=pa_extend_pathways(pw$M, N), coef=pw$coef))

    if(is.vector(pw))
        return(c(N %*% pw, pw))

    return(t(apply(dec_1$M, 1, function(x)  {  return(c(N %*% x, x))}  )))
}


# Given a reaction network and reaction rates this function adds inflow and outflow
# reactions that would maintain steady state. Function returns a list with new network
# and new rates.

pa_extend_net <- function(net, rates) {
    # calculate rates that balance growth / decrease of concentrations
    cdif_r <- jrnf_calculate_concentration_change(net, rates)

    # add reactions to balance growth / decrease
    for(i in 1:length(cdif_r)) 
        # if species' concentration increases one pseudoreaction has to be included to remove it ("X -> ")
        if(cdif_r[i] > 0) {   
	    net[[2]] <- rbind(net[[2]], data.frame(reversible=factor(c(FALSE)), 
                              c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                              activation=as.numeric(c(0)),educts=I(list(i)), educts_mul=I(list(1)),
                              products=I(list(c())), products_mul=I(list(c()))))
        # if species' concentration decreases one pseudoreaction is included to add it ("-> X ")
        } else { 
	    net[[2]] <- rbind(net[[2]], data.frame(reversible=factor(c(FALSE)), 
                              c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                              activation=as.numeric(c(0)),educts=I(list(c())), educts_mul=I(list(c())),
                              products=I(list(i)), products_mul=I(list(1))))
        }

    # include inflow / outflow rates
    rates_ext <- c(rates, abs(cdif_r))
    return(list(net, rates_ext))
}


# Function checks the reachability of individual species from a set of paths.
# The species that can not be produced and those who can not be consumed are 
# printed out...

pa_check_reachability <- function(N, path_M) { 
    x <- N %*% t(path_M)

    cat("inflow missing: ")
    for(i in 1:nrow(N))
        if(length(which(x[i,] > 0)) == 0)
            cat(i, " ")
    cat("\n")
    

    cat("outflow missing: ")
    for(i in 1:nrow(N))
        if(length(which(x[i,] < 0)) == 0)
            cat(i, " ")
    cat("\n")
}


# Helper function that takes a matrix of pathways and their coefficients,
# identifies duplicated pathways and then removes multiple occurences of
# pathways and calculate new accumulated rates for the remaining pathways.
# Function returns a logical vector identifying pathways to keep and a vector
# of new coefficients.

pa_rm_duprows_accum <- function(mat, coef) {
    keep_pw <- !duplicated(matrix(mat,nrow=nrow(mat)))
        
    for(i in which(!keep_pw)) {
        k <- 1
        while(!all(mat[k,] == mat[i,]))
            k <- k + 1
            coef[k] <- coef[k] + coef[i]
        }

    return(list(keep=keep_pw, coef=coef[keep_pw]))
} 



# Function calculates the reaction rates given a set of pathways and their
# associated rates.

pa_reconstruct_rates <- function(path_M, path_rates) {
    return(as.vector(t(path_M) %*% path_rates))
}


# Helper function for pa_subpath_decompostion
# The function checks in which rows of the matrix <path_M> the same elements
# are non-zero than in the vector <v> and returns a boolean vector of length
# nrow(path_M)

pa_is_row_contained <- function(v, path_M) {
    if(!is.matrix(path_M))
        return(c(F))

    if(nrow(path_M) == 0)
        return(c(F))
  
    return(apply(path_M, 1, function(x)  {  return(!any(x != 0 & v == 0))  }))
}

