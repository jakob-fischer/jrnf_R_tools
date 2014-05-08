# author: jakob fischer (jakob@automorph.info)
# date: 11. December 2013
# description: 
# testing-tool for algorithm which determines pathways / elementary flux modes
# of reaction networks considering a flux vector v. All functions of this 
# file start with prefix "pa_"

library(pracma)   # for gcd function
source("cycles.R")


# Funktion to write reaction <i> of network <net> to standard ouptut. If 
# stoichiometric matrix is not given as third parameter it is calculated. 
# TODO: Calculating the entire matrix just for one column is overkill!

pa_write_rea <- function(net, i, N=c()) {
    if(!is.matrix(N))
        N <- jrnf_calculate_stoich_mat(net)

    t <- N[,i]

    for(i in which(t < 0))
        cat(-t[i], " ", net[[1]]$name[i], "  ")

    cat("=>  ")
    
    for(i in which(t > 0))
        cat(t[i], " ", net[[1]]$name[i], "  ")
}


# Function writes elementary mode <em> of network <net>. Reactions that have 
# only one educt and no products or only one product and no educts are ignored 
# (they are probably) pseudo-reactions for the output and not considered when 
# printing the net reaction. 

pa_write_em <- function(net, em) {
    N <- jrnf_calculate_stoich_mat(net)   
    x <- (apply(N != 0, 2, sum) == 1)     # remove reactions with just one species involved
    em[x] <- 0    

    # subfunction for printing reaction <i> of network <net>
    pa_write_rea <- function(i) {  pa_write_rea(net, i, N)  }

    # list all (remaining) reactions with coefficients 
    for(i in which(em != 0)) {
        cat(em[i], "X : ")
        pa_write_rea(i);      
        cat("\n")
    }

    cat("============================================================\n")
    # now print the net reaction    
    x <- N %*% em
    
    for(i in which(x < 0))
        cat(-x[i], " ", net[[1]]$name[i], "  ")

    cat("=>  ")
    for(i in which(x > 0))
        cat(x[i], " ", net[[1]]$name[i], "  ")

    cat("\n")    
}


# This function calculates the greatest common divisor for all elements in 
# a vector. For this the gcd function of the package 'pracma' is used.

vec_gcd <- function(x) {
    a <- x[1]

    if(length(x) > 1) 
        for(i in 2:length(x))
            a <- gcd(a,x[i])             
           
    return(a)
}


# Function checks the reachability of individual species from a set of paths.
# The species that can not be produced and those who can not be consumed are 
# printed out...

check_reachability <- function(N, path_M) { 
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


# Function calculates the reaction rates given a set of pathways and their
# associated rates.

reconstruct_rates <- function(path_M, path_rates) {
    return(as.vector(t(path_M) %*% path_rates))
}


# Helper function for do_subpath_decompostion

is_row_contained <- function(v, path_M) {
    if(!is.matrix(path_M))
        return(c(F))

    if(nrow(path_M) == 0)
        return(c(F))
  
    return(apply(path_M, 1, function(x)  {  return(!any(x != 0 & v == 0))  }))
}


do_subpath_decomposition <- function(N_orig, path_orig) {
    # calculate the reaction set for which decomposition is done
    sel_rea <- which(path_orig != 0)
    # only species that are created and consumed by above reaction set are 
    # taken as branching species

    N <- matrix(N_orig[,sel_rea], nrow(N_orig), length(sel_rea))
    #branch_sp <- which(apply(N, 1, max) > 0 & apply(N, 1, min) < 0)
    branch_sp <- which(N %*% path_orig[sel_rea] == 0)
        

    if(length(branch_sp) == 0)
        return(path_orig)

    N <- matrix(N[branch_sp,], length(branch_sp))   
    path_M <- matrix(0, length(sel_rea), length(sel_rea))
    for(i in 1:length(sel_rea))
        path_M[i,i] <- 1

    for(i in 1:nrow(N)) {
        # net creator
        x <- N %*% t(path_M)

        cr <- which(x[i,] > 0)
        cr_m <- x[i,cr]
       
        # net consumer
        con <- which(x[i,] < 0)
        con_m <- -x[i,con]

        path_M_new <- matrix(path_M[x[i,] == 0,], ncol=ncol(path_M))

        if(length(cr) != 0 && length(con) != 0)
        for(src in cr) {
            K_src <- abs((N %*% path_M[src,])[i])              # How does the pathway src change i
            for(dest in con) {
                K_dst <- abs((N %*% path_M[dest,])[i])        # How does the pathway dest change i 
                np <- K_dst*path_M[src,] +                          # new pathway
                      K_src*path_M[dest,]
                np <- np/vec_gcd(np)
                hull <- (path_M[dest,] != 0 | path_M[src,] != 0) 

                if(!any(is_row_contained(hull, path_M_new)))
                    path_M_new <- rbind(path_M_new,matrix(np, 1, length(sel_rea)))                
            }
        }
        
        path_M <- path_M_new   
    }

    # transform back to space with all reactions in it!
    path_M_ret <- matrix(0, nrow(path_M), ncol(N_orig))
    path_M_ret[,sel_rea] <- path_M

    return(path_M_ret)
}




check_pathway_present <- function(path_M, path_rates, path_M_el) {
    id <- rep(0, nrow(path_M_el))
    present <- rep(F, nrow(path_M_el))

    if(length(path_rates) != 0 & is.matrix(path_M))
        for(i in 1:nrow(path_M_el)) {
            if(nrow(path_M) != 0)
                for(j in 1:nrow(path_M)) {
                    if(all(path_M[j,] == path_M_el[i,])) {
                        present[i] <- T
                        id[i] <- j
                    }            
                }
        }

    rate <- rep(0, nrow(path_M_el))
    rate[present] <- path_rates[id[present]]
    complexity <- apply(path_M_el != 0, 1, sum)

    return(data.frame(id=id, present=present, rate=rate, complexity=complexity))
}



#
#
#

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

    # include 
    rates_ext <- c(rates, abs(cdif_r))
    return(list(net, rates_ext))
}




pathway_analysis_prototype <- function(net, rates, fexp=0.001, f2=1e-40) {
    # rates of the discarded pathways
    rates_dropped <- rep(0, length(rates))
    
    # Stoichiometric matrix of consumption (N_in), production (N_out) and total (N)
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)
    N <- N_out - N_in

    # rates of pathways + matrix for pathways (one row = one pathway) 
    # Initially every reaction equals one pathway with the reaction's rate being the rate of the pathway
    path_rates <- rates
    path_M <- matrix(0, nrow(net[[2]]), nrow(net[[2]]))
    for(i in 1:nrow(net[[2]]))
        path_M[i,i] <- 1

    # reactions / pathways with rate of zero have to be removed
    path_rates <- path_rates[rates != 0]
    path_M <- path_M[rates != 0,]

    # maybe take generation rate for ordering?  or degree of substrate graph?
    rate_gen <- as.vector(N_in %*% rates)
    g <- jrnf_to_directed_network(net)

    # All species are taken as branching species (higher production rate first)
    #for(i in order(degree(g), decreasing=F)) { 
    for(i in order(rate_gen, decreasing=F)) {
        cat("branching at species ", net[[1]]$name[i], "\n")

        # calculate maximal relative deviation of reaction rates
        mm <- max(abs((reconstruct_rates(path_M, path_rates)-rates)/rates), na.rm=T)
        # which element has this maximal relative deviation?
        w_mm <- which(mm == abs((reconstruct_rates(path_M, path_rates)-rates)/rates))[1]
        # print check:  maximal relative deviation;  reaction id; absolute 
        cat("MCHECK: ",  mm, "  - ", (reconstruct_rates(path_M, path_rates)-rates)[w_mm] , "  - ",  w_mm ,  "\n")

        # net creator
        x <- N %*% t(path_M)
        cr <- which(x[i,] > 0)
        cr_m <- x[i,cr]
        rate_gen <- sum(cr_m*path_rates[cr])
        cr_f <- (cr_m*path_rates[cr])/rate_gen
       
        # net consumer
        #check_reachability(N, path_M)

        con <- which(x[i,] < 0)
        con_m <- -x[i,con]
        con_f <- con_m*path_rates[con]/rate_gen

        z <- cr_f %*% t(con_f)
        n <- z > fexp
        #cat(sum(z), "\n")

        cat("length(cr)=", length(cr), "  length(con)=", length(con), "\n")

        if(length(cr) == 0 & length(con) == 0) {
            cat("species ", i, " not reachable, omitting!\n")
        } else if(length(cr) == 0 | length(con) == 0) {
            cat("Warning, species ", i, " is only partially reachable!\n")
            if(length(cr) != 0)
                cat("unexplained is: ", sum(path_rates[cr]), "\n")
            else
                cat("unexplained is: ", sum(path_rates[con]), "\n")

            cat("Removing pathways that end/start here!\n")
            path_M <- path_M[x[i,] == 0,]
            path_rates <- path_rates[x[i,] == 0]
        } else {
            for(src in 1:length(cr)) {
                a <- z[src,]/sum(z[src,])
                if(length(which(a > fexp | a > 1/(length(a)+1))) == 0) {
                    cat("HAH!\n")
                    return(list(cr, cr_m, cr_f, con, con_m, con_f, path_rates))
                }

                n[src, a > fexp | a > 1/(length(a)+1)] <- T
            }

            for(dest in 1:length(con)) {
                b <- z[,dest]/sum(z[,dest])

                if(length(which(b > fexp | b > 1/(length(b)+1))) == 0) {
                    cat("HUH!\n")
                    return(z)
                }

                n[b > fexp | b > 1/(length(b)+1),dest] <- T
            }


            # Start filling new matrix with pathways that don't net produce or consume i
            #cat("pre rates: ", path_rates, "\n\n")

            path_M_new <- path_M[x[i,] == 0,]
            path_rates_new <- path_rates[x[i,] == 0]
            #cat("int rates: ", path_rates_new, "\n\n")

            # Add contribution of dropped / deleted pathways to vector rates_dropped
            for(src in 1:length(cr)) {
                K_src <- cr_m[src]             # How does the pathway src change i
                cat(".")
                for(dest in 1:length(con)) {
                    #cat("MCHECK: ",  max(abs(reconstruct_rates(path_M_new, path_rates_new)-rates)/rates, na.rm=T), "\n")

                    K_dst <- con_m[dest]       # How does the pathway dest change i 
                    nr <- z[src,dest]*rate_gen/(K_dst*K_src)                    # new rate
                    np <- K_dst*path_M[cr[src],] +                          # new pathway
                          K_src*path_M[con[dest],]  

                    


                    if(n[src,dest]) {                                 # Add combined pathway
                        if((N %*% np)[i] != 0) {
                            cat("ERROR: pathway combination not zero!\n")
                            cat("src=", src, "  dest=", dest, "   cr[src]=", cr[src],   
                                "  con[dest]=", con[dest], "\n")
                            cat("K_src=", K_src, "\n")
                            cat("K_dst=", K_dst, "\n")
                            cat("np=", np, "\n\n")
                            return(N %*% np)
                        }

                        # First calculate pathway decomposition
                        path_M_dec <- do_subpath_decomposition(N, np)
                        #path_M_dec <- matrix(np, ncol=length(np))
 
                        # 
                        p_sort <- check_pathway_present(path_M_new, path_rates_new, path_M_dec)
                        o <- order(-p_sort$rate, p_sort$complexity)
                        p_sort <- p_sort[o,]
                        path_M_dec <- matrix(path_M_dec[o,], ncol=ncol(path_M_dec))

                        s <- np*nr
                        for(k in 1:nrow(p_sort)) {
                            y <- s/path_M_dec[k,]
        
                            if(length(which(y>0)) > 0) {
                                rt <- min(y[path_M_dec[k,]>0])
                                if(rt < 0)
                                    rt <- 0
                                   

                                s <- s - path_M_dec[k,]*rt

                                if(is.na(rt)) { 
                                    cat("NA rate found!\n")
                                    return(list(path_M_dec[k,], s))
                                }
       
                                if(p_sort$present[k]) {
                                    path_rates_new[p_sort$id[k]] <- path_rates_new[p_sort$id[k]] + rt
                                } else {
                                    if(rt > f2) {
                                        path_rates_new[length(path_rates_new)+1] <- rt
                                        path_M_new <- rbind(path_M_new, matrix(path_M_dec[k,], 1, length(path_M_dec[k,])))
                                    } else {
                                        rates_dropped <- rates_dropped + rt*path_M_dec[k,]
                                    }
                                }
                            } 
                        }
                                     
                    } else {                   # Drop pathway combination
                        rates_dropped <- rates_dropped + nr*np
                    }
                }
            }

            # swap / report / ready for next step  
            path_rates <- path_rates_new
            cat("sum(is.na(path_rates))=", sum(is.na(path_rates)), "\n")
            path_M <- path_M_new
            cat("Having ", nrow(path_M), " pathways!\n")
            cat("max dropped: ", max((rates_dropped/rates)[rates != 0]), "\n")
        }
    }

    return(list(path_M, path_rates))
} 



# Trying an alternative version of the pathway analysis tool (calculating elementary modes 
# while including 
#
#

pathway_analysis_alternative <- function(net, rates, f2=1e-40, comp=50) {
    # rates of the discarded pathways
    rates_dropped <- rep(0, length(rates))
    
    # Stoichiometric matrix of consumption (N_in), production (N_out) and total (N)
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)
    N <- N_out - N_in

    # rates of pathways + matrix for pathways (one row = one pathway) 
    # Initially every reaction equals one pathway with the reaction's rate being the rate of the pathway
    path_rates <- rates
    path_M <- matrix(0, nrow(net[[2]]), nrow(net[[2]]))
    for(i in 1:nrow(net[[2]]))
        path_M[i,i] <- 1

    # reactions / pathways with rate of zero have to be removed
    path_rates <- path_rates[rates != 0]
    path_M <- path_M[rates != 0,]

    # maybe take generation rate for ordering?  or degree of substrate graph?
    rate_gen <- as.vector(N_in %*% rates)
    g <- jrnf_to_directed_network(net)

    # All species are taken as branching species (higher production rate first)
    #for(i in order(degree(g), decreasing=F)) { 
    for(i in order(rate_gen, decreasing=F)) {
        cat("branching at species ", net[[1]]$name[i], "\n")

        # net creator
        x <- N %*% t(path_M)
        cr <- which(x[i,] > 0)
        cr_m <- x[i,cr]
        rate_gen <- sum(cr_m*path_rates[cr])
        #cr_f <- (cr_m*path_rates[cr])/rate_gen    # fraction generated by pathways
       
        # net consumer
        con <- which(x[i,] < 0)
        con_m <- -x[i,con]
        rate_con <- sum(con_m*path_rates[con])
        #con_f <- con_m*path_rates[con]/rate_con   # fraction consumed by pathways

        #z <- cr_f %*% t(con_f)    # "transition" probability
        #n <- z > fexp             # which transitions are relevant

        cat("length(cr)=", length(cr), "  length(con)=", length(con), "\n")

        if(length(cr) == 0 & length(con) == 0) {
            cat("species ", i, " not reachable, omitting!\n")
        } else if(length(cr) == 0 | length(con) == 0) {
            cat("Warning, species ", i, " is only partially reachable!\n")
            if(length(cr) != 0)
                cat("unexplained is: ", sum(path_rates[cr]), "\n")
            else
                cat("unexplained is: ", sum(path_rates[con]), "\n")

            cat("Removing pathways that end/start here!\n")
            path_M <- path_M[x[i,] == 0,]
            path_rates <- path_rates[x[i,] == 0]
        } else {
            # The matrix of "relevant transitions" z has to have at least one "true" 
            # element in every row and column. This is ensured and checked here
            #for(src in 1:length(cr)) {
            #    a <- z[src,]/sum(z[src,])
            #    n[src, a > fexp | a > 1/(length(a)+1)] <- T
            #}

            #for(dest in 1:length(con)) {
            #    b <- z[,dest]/sum(z[,dest])
            #    n[b > fexp | b > 1/(length(b)+1),dest] <- T
            #}

            # Start filling new matrix with pathways that don't net produce or consume i
            path_M_new <- path_M[x[i,] == 0,]
            path_rates_new <- path_rates[x[i,] == 0]
 
            # Add contribution of dropped / deleted pathways to vector rates_dropped
            for(src in 1:length(cr)) {
                K_src <- cr_m[src]             # How does the pathway src change i
                cat(".")
                for(dest in 1:length(con)) {
                    K_dst <- con_m[dest]       # How does the pathway dest change i 
                    np <- K_dst*path_M[cr[src],] +                          # new pathway
                          K_src*path_M[con[dest],]  


                    {                                 # Add combined pathway
                        if((N %*% np)[i] != 0) {
                            cat("ERROR: pathway combination not zero!\n")
                            return()
                        }

                        # First calculate pathway decomposition
                        path_M_dec <- do_subpath_decomposition(N, np)
                        
                        # returns boolean array of EMs already present (p_sort$present)
                        p_sort <- check_pathway_present(path_M_new, path_rates_new, path_M_dec)

                        for(k in 1:nrow(p_sort)) {
                            if(!p_sort$present[k]) {

                                y <- rates/path_M_dec[k,]
        
                                if(length(which(y>0)) > 0) {
                                    rt <- min(y[path_M_dec[k,]>0])
                                    
                                    if(rt < 0)
                                        rt <- 0
                                                                        
                                     if(is.na(rt)) { 
                                         cat("NA rate found!\n")
                                         return()
                                     }
       
                                     if(rt > f2 & p_sort$complexity[k] < comp) {
                                         path_rates_new[length(path_rates_new)+1] <- rt
                                         path_M_new <- rbind(path_M_new, matrix(path_M_dec[k,], 1, length(path_M_dec[k,])))
                                     } 
                                }
                            }
                        }
                                     
                    } 
                }
            }

            # swap / report / ready for next step  
            path_rates <- path_rates_new
            cat("sum(is.na(path_rates))=", sum(is.na(path_rates)), "\n")
            path_M <- path_M_new
            cat("Having ", nrow(path_M), " pathways!\n")
            cat("max dropped: ", max((rates_dropped/rates)[rates != 0]), "\n")
        }
    }

    return(list(path_M, path_rates))
} 




# Trying an alternative version of the pathway analysis tool (calculating elementary modes 
# while including 
#
#

pathway_analysis_alternative2 <- function(net, rates, f2=1e-40, comp=50) {
    # rates of the discarded pathways
    rates_dropped <- rep(0, length(rates))
    
    # Stoichiometric matrix of consumption (N_in), production (N_out) and total (N)
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)
    N <- N_out - N_in

    # rates of pathways + matrix for pathways (one row = one pathway) 
    # Initially every reaction equals one pathway with the reaction's rate being the rate of the pathway
    path_rates <- rates
    path_M <- matrix(0, nrow(net[[2]]), nrow(net[[2]]))
    for(i in 1:nrow(net[[2]]))
        path_M[i,i] <- 1

    # reactions / pathways with rate of zero have to be removed
    path_rates <- path_rates[rates != 0]
    path_M <- path_M[rates != 0,]

    # maybe take generation rate for ordering?  or degree of substrate graph?
    rate_gen <- as.vector(N_in %*% rates)
    g <- jrnf_to_directed_network(net)

    # All species are taken as branching species (higher production rate first)
    #for(i in order(degree(g), decreasing=F)) { 
    for(i in order(rate_gen, decreasing=F)) {
        cat("branching at species ", net[[1]]$name[i], "\n")

        # net creator
        x <- N %*% t(path_M)
        cr <- which(x[i,] > 0)
        cr_m <- x[i,cr]
        rate_gen <- sum(cr_m*path_rates[cr])
        #cr_f <- (cr_m*path_rates[cr])/rate_gen    # fraction generated by pathways
       
        # net consumer
        con <- which(x[i,] < 0)
        con_m <- -x[i,con]
        rate_con <- sum(con_m*path_rates[con])
        #con_f <- con_m*path_rates[con]/rate_con   # fraction consumed by pathways

        #z <- cr_f %*% t(con_f)    # "transition" probability
        #n <- z > fexp             # which transitions are relevant

        cat("length(cr)=", length(cr), "  length(con)=", length(con), "\n")

        if(length(cr) == 0 & length(con) == 0) {
            cat("species ", i, " not reachable, omitting!\n")
        } else if(length(cr) == 0 | length(con) == 0) {
            cat("Warning, species ", i, " is only partially reachable!\n")
            if(length(cr) != 0)
                cat("unexplained is: ", sum(path_rates[cr]), "\n")
            else
                cat("unexplained is: ", sum(path_rates[con]), "\n")

            cat("Removing pathways that end/start here!\n")
            path_M <- path_M[x[i,] == 0,]
            path_rates <- path_rates[x[i,] == 0]
        } else {
            # The matrix of "relevant transitions" z has to have at least one "true" 
            # element in every row and column. This is ensured and checked here
            #for(src in 1:length(cr)) {
            #    a <- z[src,]/sum(z[src,])
            #    n[src, a > fexp | a > 1/(length(a)+1)] <- T
            #}

            #for(dest in 1:length(con)) {
            #    b <- z[,dest]/sum(z[,dest])
            #    n[b > fexp | b > 1/(length(b)+1),dest] <- T
            #}

            # Start filling new matrix with pathways that don't net produce or consume i
            path_M_new <- path_M[x[i,] == 0,]
            path_rates_new <- path_rates[x[i,] == 0]
 
            # Add contribution of dropped / deleted pathways to vector rates_dropped
            for(src in 1:length(cr)) {
                K_src <- cr_m[src]             # How does the pathway src change i
                cat(".")
                for(dest in 1:length(con)) {
                    K_dst <- con_m[dest]       # How does the pathway dest change i 
                    np <- K_dst*path_M[cr[src],] +                          # new pathway
                          K_src*path_M[con[dest],]  


                    {                                 # Add combined pathway
                        if((N %*% np)[i] != 0) {
                            cat("ERROR: pathway combination not zero!\n")
                            return()
                        }

                        # First calculate pathway decomposition
                        path_M_dec <- do_subpath_decomposition(N, np)
                        
                        # returns boolean array of EMs already present (p_sort$present)
                        p_sort <- check_pathway_present(path_M_new, path_rates_new, path_M_dec)

                        for(k in 1:nrow(p_sort)) {
                            if(!p_sort$present[k]) {

                                y <- rates/path_M_dec[k,]
        
                                if(length(which(y>0)) > 0) {
                                    rt <- min(y[path_M_dec[k,]>0])
                                    
                                    if(rt < 0)
                                        rt <- 0
                                                                        
                                     if(is.na(rt)) { 
                                         cat("NA rate found!\n")
                                         return()
                                     }
       
                                     new_comp <- sum((N %*% path_M_dec[k,]) != 0) 

                                     if(rt > f2 & new_comp <= comp) {
                                         path_rates_new[length(path_rates_new)+1] <- rt
                                         path_M_new <- rbind(path_M_new, matrix(path_M_dec[k,], 1, length(path_M_dec[k,])))
                                     } 
                                }
                            }
                        }
                                     
                    } 
                }
            }

            # swap / report / ready for next step  
            path_rates <- path_rates_new
            cat("sum(is.na(path_rates))=", sum(is.na(path_rates)), "\n")
            path_M <- path_M_new
            cat("Having ", nrow(path_M), " pathways!\n")
            cat("max dropped: ", max((rates_dropped/rates)[rates != 0]), "\n")
        }
    }

    return(list(path_M, path_rates))
} 


# Trying an alternative version of the pathway analysis tool (calculating elementary modes 
# while including 
#
#

pathway_analysis_alternative3 <- function(net, rates, f2=1e-40, C1=6, C2=10) {
    # rates of the discarded pathways
    rates_dropped <- rep(0, length(rates))
    
    # Stoichiometric matrix of consumption (N_in), production (N_out) and total (N)
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)
    N <- N_out - N_in

    # rates of pathways + matrix for pathways (one row = one pathway) 
    # Initially every reaction equals one pathway with the reaction's rate being the rate of the pathway
    path_rates <- rates
    path_M <- matrix(0, nrow(net[[2]]), nrow(net[[2]]))
    for(i in 1:nrow(net[[2]]))
        path_M[i,i] <- 1

    # reactions / pathways with rate of zero have to be removed
    path_rates <- path_rates[rates != 0]
    path_M <- path_M[rates != 0,]

    # maybe take generation rate for ordering?  or degree of substrate graph?
    rate_gen <- as.vector(N_in %*% rates)
    g <- jrnf_to_directed_network(net)

    # All species are taken as branching species (higher production rate first)
    for(i in order(degree(g), decreasing=F)) { 
    #for(i in order(rate_gen, decreasing=F)) {
        cat("branching at species ", net[[1]]$name[i], "\n")

        # net creator
        x <- N %*% t(path_M)
        cr <- which(x[i,] > 0)
        cr_m <- x[i,cr]
        rate_gen <- sum(cr_m*path_rates[cr])
        #cr_f <- (cr_m*path_rates[cr])/rate_gen    # fraction generated by pathways
       
        # net consumer
        con <- which(x[i,] < 0)
        con_m <- -x[i,con]
        rate_con <- sum(con_m*path_rates[con])
        #con_f <- con_m*path_rates[con]/rate_con   # fraction consumed by pathways

        #z <- cr_f %*% t(con_f)    # "transition" probability
        #n <- z > fexp             # which transitions are relevant

        cat("length(cr)=", length(cr), "  length(con)=", length(con), "\n")

        if(length(cr) == 0 & length(con) == 0) {
            cat("species ", i, " not reachable, omitting!\n")
        } else if(length(cr) == 0 | length(con) == 0) {
            cat("Warning, species ", i, " is only partially reachable!\n")
            if(length(cr) != 0)
                cat("unexplained is: ", sum(path_rates[cr]), "\n")
            else
                cat("unexplained is: ", sum(path_rates[con]), "\n")

            cat("Removing pathways that end/start here!\n")
            path_M <- path_M[x[i,] == 0,]
            path_rates <- path_rates[x[i,] == 0]
        } else {
            # The matrix of "relevant transitions" z has to have at least one "true" 
            # element in every row and column. This is ensured and checked here
            #for(src in 1:length(cr)) {
            #    a <- z[src,]/sum(z[src,])
            #    n[src, a > fexp | a > 1/(length(a)+1)] <- T
            #}

            #for(dest in 1:length(con)) {
            #    b <- z[,dest]/sum(z[,dest])
            #    n[b > fexp | b > 1/(length(b)+1),dest] <- T
            #}


            # Start filling new matrix with pathways that don't net produce or consume i
            path_M_new <- path_M[x[i,] == 0,]
            path_rates_new <- path_rates[x[i,] == 0]
 
            # Add contribution of dropped / deleted pathways to vector rates_dropped
            for(src in 1:length(cr)) {
                K_src <- cr_m[src]             # How does the pathway src change i
                cat(".")
                for(dest in 1:length(con)) {
                    K_dst <- con_m[dest]       # How does the pathway dest change i 
                    np <- K_dst*path_M[cr[src],] +                          # new pathway
                          K_src*path_M[con[dest],]  


                    {                                 # Add combined pathway
                        if((N %*% np)[i] != 0) {
                            cat("ERROR: pathway combination not zero!\n")
                            return()
                        }

                        # First calculate pathway decomposition
                        path_M_dec <- do_subpath_decomposition(N, np)
                        
                        # returns boolean array of EMs already present (p_sort$present)
                        p_sort <- check_pathway_present(path_M_new, path_rates_new, path_M_dec)

                        for(k in 1:nrow(p_sort)) {
                            if(!p_sort$present[k]) {

                                y <- rates/path_M_dec[k,]
        
                                if(length(which(y>0)) > 0) {
                                    rt <- min(y[path_M_dec[k,]>0])
                                    
                                    if(rt < 0)
                                        rt <- 0
                                                                        
                                     if(is.na(rt)) { 
                                         cat("NA rate found!\n")
                                         return()
                                     }
       
                                     new_comp <- sum((N %*% path_M_dec[k,]) != 0) 

                                     if(rt > f2 & new_comp <= C1 & p_sort$complexity[k] < C2) {
                                         path_rates_new[length(path_rates_new)+1] <- rt
                                         path_M_new <- rbind(path_M_new, matrix(path_M_dec[k,], 1, length(path_M_dec[k,])))
                                     } 
                                }
                            }
                        }
                                     
                    } 
                }
            }

            # swap / report / ready for next step  
            path_rates <- path_rates_new
            cat("sum(is.na(path_rates))=", sum(is.na(path_rates)), "\n")
            path_M <- path_M_new
            cat("Having ", nrow(path_M), " pathways!\n")
            cat("max dropped: ", max((rates_dropped/rates)[rates != 0]), "\n")
        }
    }

    return(list(path_M, path_rates))
} 




pathway_analysis_prototype2 <- function(net, rates, fexp=0.001, f2=1e-10, C1=6, C2=10) {
    # rates of the discarded pathways
    rates_dropped <- rep(0, length(rates))
    
    # Stoichiometric matrix of consumption (N_in), production (N_out) and total (N)
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)
    N <- N_out - N_in

    # rates of pathways + matrix for pathways (one row = one pathway) 
    # Initially every reaction equals one pathway with the reaction's rate being the rate of the pathway
    path_rates <- rates
    path_M <- matrix(0, nrow(net[[2]]), nrow(net[[2]]))
    for(i in 1:nrow(net[[2]]))
        path_M[i,i] <- 1

    # reactions / pathways with rate of zero have to be removed
    path_rates <- path_rates[rates != 0]
    path_M <- path_M[rates != 0,]

    # maybe take generation rate for ordering?  or degree of substrate graph?
    rate_gen <- as.vector(N_in %*% rates)
    g <- jrnf_to_directed_network(net)

    # All species are taken as branching species (higher production rate first)
    #for(i in order(degree(g), decreasing=F)) { 
    for(i in order(rate_gen, decreasing=F)) {
        cat("branching at species ", net[[1]]$name[i], "\n")

        # calculate maximal relative deviation of reaction rates
        mm <- max(abs((reconstruct_rates(path_M, path_rates)-rates)/rates), na.rm=T)
        # which element has this maximal relative deviation?
        w_mm <- which(mm == abs((reconstruct_rates(path_M, path_rates)-rates)/rates))[1]
        # print check:  maximal relative deviation;  reaction id; absolute 
        cat("MCHECK: ",  mm, "  - ", (reconstruct_rates(path_M, path_rates)-rates)[w_mm] , "  - ",  w_mm ,  "\n")

        # net creator
        x <- N %*% t(path_M)
        cr <- which(x[i,] > 0)
        cr_m <- x[i,cr]
        rate_gen <- sum(cr_m*path_rates[cr])
        cr_f <- (cr_m*path_rates[cr])/rate_gen
       
        # net consumer
        #check_reachability(N, path_M)

        con <- which(x[i,] < 0)
        con_m <- -x[i,con]
        con_f <- con_m*path_rates[con]/rate_gen

        z <- cr_f %*% t(con_f)
        n <- z > fexp
        #cat(sum(z), "\n")

        cat("length(cr)=", length(cr), "  length(con)=", length(con), "\n")

        if(length(cr) == 0 & length(con) == 0) {
            cat("species ", i, " not reachable, omitting!\n")
        } else if(length(cr) == 0 | length(con) == 0) {
            cat("Warning, species ", i, " is only partially reachable!\n")
            if(length(cr) != 0)
                cat("unexplained is: ", sum(path_rates[cr]), "\n")
            else
                cat("unexplained is: ", sum(path_rates[con]), "\n")

            cat("Removing pathways that end/start here!\n")
            path_M <- path_M[x[i,] == 0,]
            path_rates <- path_rates[x[i,] == 0]
        } else {
            for(src in 1:length(cr)) {
                a <- z[src,]/sum(z[src,])
                if(length(which(a > fexp | a > 1/(length(a)+1))) == 0) {
                    cat("HAH!\n")
                    return(list(cr, cr_m, cr_f, con, con_m, con_f, path_rates))
                }

                n[src, a > fexp | a > 1/(length(a)+1)] <- T
            }

            for(dest in 1:length(con)) {
                b <- z[,dest]/sum(z[,dest])

                if(length(which(b > fexp | b > 1/(length(b)+1))) == 0) {
                    cat("HUH!\n")
                    return(z)
                }

                n[b > fexp | b > 1/(length(b)+1),dest] <- T
            }

            # Start filling new matrix with pathways that don't net produce or consume i
            #cat("pre rates: ", path_rates, "\n\n")

            path_M_new <- path_M[x[i,] == 0,]
            path_rates_new <- path_rates[x[i,] == 0]
            #cat("int rates: ", path_rates_new, "\n\n")

            # Add contribution of dropped / deleted pathways to vector rates_dropped
            for(src in 1:length(cr)) {
                K_src <- cr_m[src]             # How does the pathway src change i
                cat(".")
                for(dest in 1:length(con)) {
                    #cat("MCHECK: ",  max(abs(reconstruct_rates(path_M_new, path_rates_new)-rates)/rates, na.rm=T), "\n")

                    K_dst <- con_m[dest]       # How does the pathway dest change i 
                    nr <- z[src,dest]*rate_gen/(K_dst*K_src)                    # new rate
                    np <- K_dst*path_M[cr[src],] +                          # new pathway
                          K_src*path_M[con[dest],]  

                    


                    if(n[src,dest]) {                                 # Add combined pathway
                        if((N %*% np)[i] != 0) {
                            cat("ERROR: pathway combination not zero!\n")
                            cat("src=", src, "  dest=", dest, "   cr[src]=", cr[src],   
                                "  con[dest]=", con[dest], "\n")
                            cat("K_src=", K_src, "\n")
                            cat("K_dst=", K_dst, "\n")
                            cat("np=", np, "\n\n")
                            return(N %*% np)
                        }

                        # First calculate pathway decomposition
                        path_M_dec <- do_subpath_decomposition(N, np)
                        #path_M_dec <- matrix(np, ncol=length(np))
 
                        # 
                        p_sort <- check_pathway_present(path_M_new, path_rates_new, path_M_dec)
                        o <- order(-p_sort$rate, p_sort$complexity)
                        p_sort <- p_sort[o,]
                        path_M_dec <- matrix(path_M_dec[o,], ncol=ncol(path_M_dec))

                        s <- np*nr
                        for(k in 1:nrow(p_sort)) {
                            y <- s/path_M_dec[k,]
        
                            if(length(which(y>0)) > 0) {
                                rt <- min(y[path_M_dec[k,]>0])
                                if(rt < 0)
                                    rt <- 0
                                   

                                s <- s - path_M_dec[k,]*rt

                                if(is.na(rt)) { 
                                    cat("NA rate found!\n")
                                    return(list(path_M_dec[k,], s))
                                }
       
                                if(p_sort$present[k]) {
                                    path_rates_new[p_sort$id[k]] <- path_rates_new[p_sort$id[k]] + rt
                                } else {
                                     new_comp <- sum((N %*% path_M_dec[k,]) != 0) 
                                     if(rt > f2 & new_comp <= C1 & p_sort$complexity[k] < C2) {
                                        path_rates_new[length(path_rates_new)+1] <- rt
                                        path_M_new <- rbind(path_M_new, matrix(path_M_dec[k,], 1, length(path_M_dec[k,])))
                                    } else {
                                        rates_dropped <- rates_dropped + rt*path_M_dec[k,]
                                    }
                                }
                            } 
                        }
                                     
                    } else {                   # Drop pathway combination
                        rates_dropped <- rates_dropped + nr*np
                    }
                }
            }

            # swap / report / ready for next step  
            path_rates <- path_rates_new
            cat("sum(is.na(path_rates))=", sum(is.na(path_rates)), "\n")
            path_M <- path_M_new
            cat("Having ", nrow(path_M), " pathways!\n")
            cat("max dropped: ", max((rates_dropped/rates)[rates != 0]), "\n")
        }
    }

    return(list(path_M, path_rates))
} 





em_rm_multiple_check_elementary <- function(net, ems) {
    ems <- ems[,1:nrow(net[[2]])]
    em_new <- ems[c(),]
    N <- jrnf_calculate_stoich_mat(net)

    for(i in 1:nrow(ems)) {
        cat(".")
        path_M_dec <- do_subpath_decomposition(N, ems[i,])
        if(nrow(path_M_dec) != 1)
            cat("WARNING: Mode number ", i, " was not elementary!\n")

        for(j in 1:nrow(path_M_dec)) {
            cpp <- check_pathway_present(em_new, rep(1, nrow(em_new)), matrix(path_M_dec[j,],1))
            if(!any(cpp$present))  # adding
                em_new <- rbind(em_new, path_M_dec[j,])
        }
    }
    cat("\n\n")

    return(em_new)
}



em_sort_by_rate <- function(ems, v) {
    ems <- ems[,1:length(v)]

    al <- function(x) {
        y <- v/x
        if(length(which(y>0)) > 0) 
            return(min(y[x>0]))
        else
            return(0)
    }

    a <- apply(ems, 1, al)

    return(ems[order(a,decreasing=T),])
}



em_develop_rates_max <- function(ems, v) {
    ems <- matrix(ems[,1:length(v)], ncol=length(v))
    em_new <- ems[c(),]
    coeff <- c()
    non_zero <- v != 0

    finished <- !is.matrix(ems) | length(v) == 0

    al <- function(x) {
        y <- v/x
        if(length(which(y>0)) > 0) 
            return(min(y[x>0]))
        else
            return(0)
    }

    while(!finished) {
        a <- apply(ems, 1, al)
        m_id <- order(a, decreasing=T)[1]

        #y <- v/ems[m_id,]
        m_rate <- a[m_id]
        #lim_rea <- which(y[ems[m_id,]>0] == m_rate)[1]

        cat("coefficient = ", m_rate, "\n")

        if(m_rate != 0) {
            coeff <- c(coeff, m_rate)
            v <- v - ems[m_id,]*m_rate
            v[v<0] <- 0
            #v[lim_rea] <- 0

            em_new <- rbind(em_new, ems[m_id,])   # add elementary mode to new elementary modes
            ems <- ems[-m_id,]                    # remove it from old list
            #ems <- ems[ems[,lim_rea] == 0,]       # also remove all elementary modes that have the limiting (now zero) reaction
        }

       finished <- !is.matrix(ems) | length(v) == 0 | m_rate == 0
    }

    # Return new list / matrix of elementary modes  +  coefficients  + 
    return(list(em_new, coeff, v))  
}






# This function orders the elementary modes and generates a data frame containing information
# on how much the reordered elementary modes explain the rate vector 'v'. This is done es well
# cumulative (only those parts not explained by previous modes) as well as individual. It is 
# also interesting what is the average fraction and the minimum fraction of reaction rates. 
# Here 'individual' fraction means the maximum fraction of this and all previous elementary modes.
# min_R ( max_{1...i} ( f ) )
#
# For every elementary mode also it's complexity is put into the data frame. 

pa_iterate_rates_max <- function(em, v, order="max", net) {
    # cut columns of data frame or elements of rate vector 'v'
    if(length(v) < ncol(em))
        em <- em[,1:length(v)]
    else
        v <- v[1:ncol(em)]

    N <- jrnf_calculate_stoich_mat(net)      #
    # inflow / outflow reaction are those with exactly one coefficient in a stoichiometric matrix column 
    pseudo_r <- apply(N != 0, 2, sum) == 1   

    em_id <- c()       # index of the elementary mode in matrix 'em'
    coeff <- c()       # coefficient of em for initial v 
    coeff_acc <- c()   # coefficient of em for v - parts explained by previous ems
    exp_f <- c()       # fraction of rate explained with current em
    exp_f_acc <- c()   # fraction of rate explained with current and previous ems
    C1 <- c()          # complexity of current em (sum of all coeffiecients)
    C2 <- c()          # complexity of current em (number of non zero coefficients)
    C3 <- c()          # number of pseudoreactions (inflow or outflow)
    min_f <- c()       # minimum fraction of reaction explained by individual ems
    min_f_acc <- c()   # minimum fraction of reaction explained by this+previous ems 

    exp_v <- rep(0, length(v))   # The fraction of reaction i explained by the best of the previous EMs

    v_ <- v                       # 'v_' is the rate vector minus parts explained by previous modes

    # if elementary modes 'em' is empty or no matrix just return empty data frame
    if(!is.matrix(em) | nrow(em) == 0) 
        return(data.frame(em_id=numeric(), coeff=numeric(), coeff_acc=numeric(), exp_f=numeric(), 
               exp_f_acc=numeric(), C1=numeric(), C2=numeric(), min_f=numeric(), min_f_acc=numeric()))

    # al-function gets one elementary mode and returns it's max coefficient (for v_)
    al <- function(x) {
        y <- v_/x
        if(length(which(y>0)) > 0) 
            return(min(y[x>0]))
        else
            return(0)
    }

    # adds an element to the results. you can give an value for the coefficient co and the
    # accumulated coeffiecient co_a (expansion) ; else NA should be given and the function
    # will calculate them
    add_element <- function(i, co, co_a) {
        if(is.na(co))            # the (maximal) coefficent of this em   
            co <- a_init[i]

        if(is.na(co_a))          # the coefficient of the em considering the previous em's have been "substracted"
            co_a <- al(em[i,])

        em_id <<- c(em_id, i)               # build vector of em's ids
        coeff <<- c(coeff, co)              
        coeff_acc <<- c(coeff_acc, co_a)
        C1 <<- c(C1, sum(em[i,]))           # sum of all reaction's coefficients in the em 
        C2 <<- c(C2, sum(em[i,] != 0))      # number of reactions with non-zero coefficients
        C3 <<- c(C3, sum((em[i,] != 0) & pseudo_r))

        v_ <<- v_ - co_a*em[i,]             # subtracting maximal possible part of this em from v_
        v_[v_<0] <<- 0                      # v_ has to stay non-negative (may become zero because of numerical problems)

        exp_f <<- c(exp_f, sum(co*em[i,])/sum(v))       # (maximum) fraction of activity explained by this reaction 
        exp_f_acc <<- c(exp_f_acc, 1-sum(v_)/sum(v))    # fraction of activity explained when expanding up to this reaction

        tmp <- (co*em[i,])/v                            # (maximum) fraction of flow this reaction explains (for every reaction)
        tmp[v == 0] <- 1                                 
        exp_v <<- pmax(exp_v, tmp, na.rm=T)             # Fraction the worst reaction is explained by its best em      
 
        min_f <<- c(min_f, min(exp_v[v!=0], na.rm=T))                   # build vector
        min_f_acc <<- c(min_f_acc, min(((v-v_)/v)[v != 0], na.rm=T))    # build vector; fraction of the worst explained reaction (from the expansion)
    }


    a_init <- apply(em, 1, al)


    # different modes in how the elementary modes are ordered for developing v / v_    
    # "max": after each step the maximum mode best for the next step is selected
    if(order == "max") {
        prev <- T          # will be assigned the coefficient of the previous em 

        for(i in 1:nrow(em)) {
            if(prev) {
                a <- apply(em, 1, al)
                m_id <- order(a, decreasing=T)
                add_element(m_id[1], NA, a[m_id[1]])
                cat(".")

                if(a[m_id[1]] == 0) {
                    prev <- F
                    m_id <- m_id[-1]
                }
            } else {
                add_element(m_id[1], NA, a[m_id[1]])
                cat(".")
                m_id <- m_id[-1]
            }


        }

    # "initial": initially the elementary pathways are ordered with decreasing coefficient
    } else if(order == "initial") {
        a <- apply(em, 1, al)
        m_id <- order(a, decreasing=T)

        for(i in m_id) { 
            add_element(i, a[i], NA);
            cat(".")
        }

    # -: pathways are ordered as given
    } else {
        for(i in 1:nrow(em)) {
            add_element(i, NA, NA)
            cat(".")
        }
    }

    # return data frame with results
    return(data.frame(em_id=em_id, coeff=coeff, coeff_acc=coeff_acc, exp_f=exp_f, 
                      exp_f_acc=exp_f_acc, C1=C1, C2=C2, C3=C3, min_f=min_f, min_f_acc=min_f_acc))  
}
