library(pracma)   # for gcd function

vec_gcd <- function(x) {
    a <- x[1]

    if(length(x) > 1) 
        for(i in 2:length(x))
            a <- gcd(a,x[i])             
           
    return(a)
}



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


jrnf_calculate_stoich_mat <- function(net) {
    return(jrnf_calculate_stoich_mat_out(net)-jrnf_calculate_stoich_mat_in(net))
}



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


reconstruct_rates <- function(path_M, path_rates) {
    return(as.vector(t(path_M) %*% path_rates))
}


is_row_contained <- function(v, path_M) {
    if(!is.matrix(path_M))
        return(c(F))

    a <- c()

    for(i in 1:nrow(path_M))
        a <- c(a,!any((path_M[i,] != 0) & (v == 0)))
   
    return(a)
}


do_subpath_decomposition <- function(N_orig, path_M_orig, path_orig) {
    # calculate the reaction set for which decomposition is done
    sel_rea <- which(path_orig != 0)
    # only species that are created and consumed by above reaction set are 
    # taken as branching species

    N <- N_orig[,sel_rea]
    branch_sp <- which(apply(N, 1, max) > 0 & apply(N, 1, min) < 0)
    N <- N[branch_sp,]    
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

        path_M_new <- path_M[x[i,] == 0,]

        for(src in 1:length(cr)) {
            K_src <- abs((N %*% path_M[cr[src],])[i])              # How does the pathway src change i
            for(dest in 1:length(con)) {
                K_dst <- abs((N %*% path_M[con[dest],])[i])        # How does the pathway dest change i 
                np <- K_dst*path_M[cr[src],] +                          # new pathway
                      K_src*path_M[con[dest],]
                np <- np/vec_gcd(np)
                hull <- (path_M[con[dest],] != 0 | path_M[cr[src],] != 0) 
                
                if(!any(is_row_contained(hull, path_M_new)))
                    path_M_new <- rbind(path_M_new,matrix(np, 1, length(sel_rea)))                
            }
        }
        
        path_M <- path_M_new       
    }

    # transform to space with all reactions in it!
    path_M_ret <- matrix(0, nrow(path_M), ncol(N_orig))
    path_M_ret[,sel_rea] <- path_M

    return(path_M_ret)
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




pathway_analysis_prototype <- function(net, rates, fexp=0.001, f2=1e-50) {
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

    cat("-> ", nrow(net[[2]]), "\n")

    # rates of the reactions + rates of discarded pathways
    rates <- c(rates, abs(cdif_r))
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
    rate_gen <- N_in %*% rates        # should be same anyway?
    rate_des <- N_out %*% rates

    # All species are taken as branching species (higher production rate first)
    for(i in order(N_out %*% rates,decreasing=T)) { 
        cat("branching at species ", i, "\n")

        # net creator
        x <- N %*% t(path_M)
        cr <- which(x[i,] > 0)
        cr_m <- x[i,cr]
        cr_f <- (cr_m*path_rates[cr])/rate_gen[i]
       
        # net consumer
        check_reachability(N, path_M)

        con <- which(x[i,] < 0)
        con_m <- -x[i,con]
        con_f <- con_m*path_rates[con]/rate_des[i]

        z <- cr_f %*% t(con_f)
        n <- z > fexp

        cat("length(cr)=", length(cr), "  length(con)=", length(con), "\n")

        if(length(cr) == 0 | length(con) == 0) {
            cat("species ", i, " not reachable, omitting!\n")
        } else {
            for(src in 1:length(cr)) {
                a <- z[src,]/sum(z[src,])
                if(length(which(a > fexp)) == 0) {
                    cat("HAH!\n")
                    return(list(cr, cr_m, cr_f, path_rates))
                }

                n[src, a > fexp] <- T
            }

            for(dest in 1:length(con)) {
                b <- z[,dest]/sum(z[,dest])

                if(length(which(b > fexp)) == 0) {
                    cat("HUH!\n")
                    return(z)
                }

                n[b > fexp,dest] <- T
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
                K_src <- abs((N %*% path_M[cr[src],])[i])              # How does the pathway src change i
                for(dest in 1:length(con)) {
                    K_dst <- abs((N %*% path_M[con[dest],])[i])        # How does the pathway dest change i 
                    nr <- z[src,dest]/(K_dst*K_src)    # new rate
                    if(n[src,dest] && nr > f2) {                        # Add combined pathway
                        
                        np <- K_dst*path_M[cr[src],] +                          # new pathway
                              K_src*path_M[con[dest],]  
                    
                        if((N %*% np)[i] != 0) {
                            cat("ERROR: pathway combination not zero!\n")
                            cat("src=", src, "  dest=", dest, "   cr[src]=", cr[src],   
                                "  con[dest]=", con[dest], "\n")
                            cat("K_src=", K_src, "\n")
                            cat("K_dst=", K_dst, "\n")
                            cat("np=", np, "\n\n")
                            return(N %*% np)
                        }

                        #cat("z= ", z[src,dest], "  K_dst=", K_dst, "  K_src=", K_src, "\n")
                     


                        # TO INCLUDE: pathway decomposition


                        path_rates_new[length(path_rates_new)+1] <- nr
                        path_M_new <- rbind(path_M_new, np)
                                     
                    } else {                   # Drop pathway combination
                        rates_dropped <- rates_dropped + z[src,dest]*(path_M[src,]+path_M[dest,])
                    }
                }
            }

            # swap / report / ready for next step  
            path_rates <- path_rates_new
            path_M <- path_M_new
            cat("Having ", nrow(path_M), " pathways!\n")
            cat("max dropped: ", max((rates_dropped/rates)[rates != 0]), "\n")
        }
    }

    return(list(path_M, path_rates))
} 





