library(pracma)   # for gcd function


# This function calculates the greatest common divisor for all elements in 
# a vector. For this the gcd function of the package 'pracma' is used.

vec_gcd <- function(x) {
    a <- x[1]

    if(length(x) > 1) 
        for(i in 2:length(x))
            a <- gcd(a,x[i])             
           
    return(a)
}


# Calculates the 'in' stoichiometric matrix of the jrnf-network net. 
# 'in' means every column contains the left side of one reaction. 
# All values are positive.

jrnf_calculate_stoich_mat_in <- function(net) {
    no_sp <- nrow(net[[1]])
    no_re <- nrow(net[[2]])
    N <- matrix(0, no_sp, no_re)
    
    for(i in 1:no_re) {
        for(j in 1:length(net[[2]]$educts[[i]])) {
            e <- net[[2]]$educts[[i]][j]
            N[e, i] <- N[e, i] + net[[2]]$educts_mul[[i]][j]
        }
    }

    return(N)
}


# Calculates the 'out' stoichiometric matrix of the jrnf-network net. 
# 'out' means every column contains the right side of one reaction. 
# All values are positive.

jrnf_calculate_stoich_mat_out <- function(net) {
    no_sp <- nrow(net[[1]])
    no_re <- nrow(net[[2]])
    N <- matrix(0, no_sp, no_re)
    
    for(i in 1:no_re) {
        for(j in 1:length(net[[2]]$products[[i]])) {
            p <- net[[2]]$products[[i]][j]
            N[p, i] <- N[p, i] + net[[2]]$products_mul[[i]][j]
        }
    }

    return(N)
}


# Calculates the stoichiometric matrix for a jrnf-network
# (matrix contains information on net change of species with reactions)

jrnf_calculate_stoich_mat <- function(net) {
    return(jrnf_calculate_stoich_mat_out(net)-jrnf_calculate_stoich_mat_in(net))
}


# Checks matrix for pathway_analysis_prototype function. Every row and every
# column has to have at least one TRUE element.

check_mat <- function(m) {
    c <- T

    for(i in 1:nrow(m)) 
        if(all(!m[i,]))
            c <- F
        
    for(i in 1:ncol(m))
        if(all(!m[,i]))
            c <- F

    if(!c) 
        cat("ERROR: matrix check failed!\n")

    return(!c)
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
#

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


    #return(list(N, path_M))



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
    complexity <- apply(path_M_el, 1, sum)

    return(data.frame(id=id, present=present, rate=rate, complexity=complexity))
}



                

#
#N_test <- t(matrix(c(1, -1, 0, -1, -1, 0, 0, 0, -1, 0,
#                    0, 0, 0, 1, 1, 0, -1, 0, 0, 0,
#                    0, 1, 1, -1, 0, 0, 0, 0, 0, 0,
#                    0, 0, -1, 1, 0, -1, 0, 0, 0, 0,
#                    0, 0, 0, 0, 0, 1, 1, -1, 0, 1,
#                    0, 0, 0, 0, 0, 0, 0, 0, 1, -1), 10))

#N_test2 <- t(matrix(c(1, -1, 0, 0, 0, 0,
#                      0, 1, -1, 1, -1, 0,
#                      0, 0, 0, 0, 1, -1, 
#                      0, 0, 1, -1, 0, 0), 6))



pathway_analysis_extend <- function(net, rates) {
    # calculate rates that balance growth / decrease of concentrations
    cdif_r <- jrnf_calculate_concentration_change(net, rates)

    # add reactions to balance growth / decrease
    for(i in 1:length(cdif_r)) 
        if(cdif_r[i] > 0) {   
	    net[[2]] <- rbind(net[[2]], data.frame(reversible=factor(c(FALSE)), 
                              c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                              activation=as.numeric(c(0)),educts=I(list(i)), educts_mul=I(list(1)),
                              products=I(list(c())), products_mul=I(list(c()))))
        } else { 
	    net[[2]] <- rbind(net[[2]], data.frame(reversible=factor(c(FALSE)), 
                              c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                              activation=as.numeric(c(0)),educts=I(list(c())), educts_mul=I(list(c())),
                              products=I(list(i)), products_mul=I(list(1))))
        }

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

    # How fast are the species produced / consumed
    #rate_gen <- N_in %*% rates      

    g <- jrnf_to_directed_network(net)

    # maybe take generation rate for ordering?
    rate_gen <- as.vector(N_in %*% rates)

    # All species are taken as branching species (higher production rate first)
    #for(i in order(degree(g), decreasing=T)) { 
    for(i in order(rate_gen, decreasing=T)) {
        cat("branching at species ", net[[1]]$name[i], "\n")
        #cat("CHECK: ", (reconstruct_rates(path_M, path_rates)-rates)/rates, "\n")
        mm <- max(abs((reconstruct_rates(path_M, path_rates)-rates)/rates), na.rm=T)
        w_mm <- which(mm == abs((reconstruct_rates(path_M, path_rates)-rates)/rates))[1]

        cat("MCHECK: ",  mm, "  - ", (reconstruct_rates(path_M, path_rates)-rates)[w_mm] , "  - ",  w_mm ,  "\n")

        #cat("PATHS: ", path_rates, "\n")

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

            if(check_mat(n))  # Every row and every column has to have at least 1 true element
                return()    

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
    ems <- ems[,1:length(v)]
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





em_develop_rates_min <- function(ems, v) {
    ems <- ems[,1:length(v)]
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
        b <- order(a, decreasing=F)
        m_id <- b[a[b] != 0][1]

        #y <- v/ems[m_id,]
        m_rate <- a[m_id]
        #lim_rea <- which(y[ems[m_id,]>0] == m_rate)[1]

        cat("coefficient = ", m_rate, "\n")

        if(!is.na(m_rate) & m_rate != 0) {
            coeff <- c(coeff, m_rate)
            v <- v - ems[m_id,]*m_rate
            v[v<0] <- 0
            #v[lim_rea] <- 0

            em_new <- rbind(em_new, ems[m_id,])   # add elementary mode to new elementary modes
            ems <- ems[-m_id,]                    # remove it from old list
            #ems <- ems[ems[,lim_rea] == 0,]       # also remove all elementary modes that have the limiting (now zero) reaction
        }

       finished <- !is.matrix(ems) | length(v) == 0 | m_rate == 0 | is.na(m_rate)
    }

    # Return new list / matrix of elementary modes  +  coefficients  + 
    return(list(em_new, coeff, v))  
}

