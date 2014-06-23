# author: jakob fischer (jakob@automorph.info)
# date: 11. December 2013
# description: 
# testing-tool for algorithm which determines pathways / elementary flux modes
# of reaction networks considering a flux vector v. All functions of this 
# file start with prefix "pa_"

library(pracma)   # for gcd function
source("cycles.R")
source("jrnf_network.R")


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

    # list all (remaining) reactions with coefficients 
    for(i in which(em != 0)) {
        cat(em[i], "X : ")
        pa_write_rea(net, i, N);      
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


# Given a reaction network <N_orig> and a subset of active reactions <path_orig>
# this function calculates all elementary flux modes in the subset of reactions.
# This function follows the original efm-algorithm (TODO: source) and can only be
# used for small networks. It is used as a helper for the function which is usable
# for bigger networks.

pa_subpath_decomposition <- function(N_orig, path_orig, sel=c()) {
    # calculate the reaction set for which decomposition is done
    sel_rea <- which(path_orig != 0)

    # only species that are created and consumed by above reaction set are 
    # taken as branching species

    N <- matrix(N_orig[,sel_rea], nrow(N_orig), length(sel_rea))
    #branch_sp <- which(N %*% path_orig[sel_rea] == 0)
    if(length(sel) == 0)
        branch_sp <- which(apply(N, 1, max) > 0 & apply(N, 1, min) < 0)
    else
        branch_sp <- which(apply(N, 1, max) > 0 & apply(N, 1, min) < 0 & sel)      
 
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
                K_dst <- abs((N %*% path_M[dest,])[i])         # How does the pathway dest change i 
                np <- K_dst*path_M[src,] +                          # new pathway
                      K_src*path_M[dest,]
                np <- np/vec_gcd(np)
                hull <- (path_M[dest,] != 0 | path_M[src,] != 0) 

                if(!any(pa_is_row_contained(hull, path_M_new)))
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


# The function checks for all pathways in <path_M_el> if they are already in
# <path_M>.

pa_check_pathway_present <- function(path_M, path_rates, path_M_el) {
    id <- rep(0, nrow(path_M_el))
    present <- rep(F, nrow(path_M_el))

    # First check which new pathways (in <path_M_el>) are actually present
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

    # transfer rate and calculate complexit (number of reactions that are non-zero)
    rate <- rep(0, nrow(path_M_el))
    rate[present] <- path_rates[id[present]]
    complexity <- apply(path_M_el != 0, 1, sum)

    # return data frame with id of pathway in <path_M>, logical value if it is there at 
    # all and its rate there as well as complexity
    return(data.frame(id=id, present=present, rate=rate, complexity=complexity))
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


# TODO find a more intelligent method of deciding which pathways to drop
# 
# parameters:
# net  -  the network that is analysed
# rates  -  the reaction rates of the network

pa_analysis <- function(net, rates, fexp=0.1, f2=1e-20, pmin=0.01, dir=F) {
    # flag those rates that the reduction condition is applied to
    flag_red <- rates > f2

    # rates of the discarded pathways
    rates_dropped <- rep(0, length(rates))
    # flag showing which species have been used for branching
    sp_br_flag <- rep(F, nrow(net[[1]]))
    
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
    count <- 1

    # All species are taken as branching species (higher production rate first)
    #for(i in order(degree(g), decreasing=F)) {
    for(i in order(rate_gen, decreasing=dir)) {
        cat("branching at species ", net[[1]]$name[i], "\n")
        sp_br_flag[i] <- T

        # net creator
        x <- N %*% t(path_M)
        cr <- which(x[i,] > 0)
        cr_m <- x[i,cr]
        rate_gen <- sum(cr_m*path_rates[cr])
        cr_f <- (cr_m*path_rates[cr])/rate_gen
       
        # net consumer
        con <- which(x[i,] < 0)
        con_m <- -x[i,con]
        con_f <- con_m*path_rates[con]/rate_gen

        z <- cr_f %*% t(con_f)
        n <- z > fexp

        cat("length(cr)=", length(cr), "  length(con)=", length(con), "\n")

        # Species is not part of any pathway (easy as pie)!
        if(length(cr) == 0 & length(con) == 0) {
            cat("species ", i, " not reachable, omitting!\n")
        # If species is only in pathways as outflow or inflow:
        # remove pathways and print a warning!
        } else if(length(cr) == 0 | length(con) == 0) {
            cat("Warning, species ", i, " is only partially reachable!\n")
            if(length(cr) != 0)
                cat("unexplained is: ", sum(path_rates[cr]), "\n")
            else
                cat("unexplained is: ", sum(path_rates[con]), "\n")

            cat("Removing pathways that end/start here!\n")
            path_M <- path_M[x[i,] == 0,]
            path_rates <- path_rates[x[i,] == 0]

        # Species is produced and consumed. 
        } else {
            # Ensure that in each row and each column there is at least one TRUE
            # (implies that each pathway is at least connected to one other pathway)
            for(src in 1:length(cr)) {
                a <- z[src,]/sum(z[src,])
                if(length(which(a > fexp | a > 1/(length(a)+1))) == 0) {
                    #cat("src=", src, "\n")
                    #cat("z[src,]=", z[src,], "\n")
                    #cat("a=", a, "\n")
                    cat("HAH!\n")
                    #return(list(cr, cr_m, cr_f, con, con_m, con_f, path_rates))
                }

                n[src, a > fexp | a > 1/(length(a)+1)] <- T
            }

            for(dest in 1:length(con)) {
                b <- z[,dest]/sum(z[,dest])

                if(length(which(b > fexp | b > 1/(length(b)+1))) == 0) {
                    #cat("z=", z, "\n\n")
                    #cat("dim(z)=", dim(z), "\n")
                    #cat("z[,dest]=", z[,dest], "\n")
                    #cat("dest=", dest,"\n")
                    #cat("b=", b, "\n")
                    cat("HUH!\n")
                    #return(z)
                }

                n[b > fexp | b > 1/(length(b)+1),dest] <- T
            }

            # Start filling new matrix with pathways that don't net produce or consume i
            path_M_new <- path_M[x[i,] == 0,]
            path_rates_new <- path_rates[x[i,] == 0]
 
            # Add contribution of dropped / deleted pathways to vector rates_dropped
            for(src in 1:length(cr)) {
                K_src <- cr_m[src]             # How does the pathway src change i
                cat(".")
                for(dest in 1:length(con)) {
                    K_dst <- con_m[dest]       # How does the pathway dest change i 
                    nr <- z[src,dest]*rate_gen/(K_dst*K_src)                    # new rate
                    np <- K_dst*path_M[cr[src],] +                          # new pathway
                          K_src*path_M[con[dest],]  

                    # Calculate the minimal contribution of this (coupled) pathway to the total rate
                    # If it is above <pmin> the two pathway are coupled even if the matrix <n> for it
                    # is false.

                    x <- abs((np*nr/rates)[np != 0 & flag_red])
                    if(length(x) == 0)
                        x <- abs((np*nr/rates)[np != 0 & !flag_red])                        

                    mm <- max(x)
                    #cat("mm=", mm, "\n")

                    if(n[src,dest] | mm > pmin) {                                 # Add combined pathway
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
                        path_M_dec <- pa_subpath_decomposition(N, np, sp_br_flag)
                        #path_M_dec <- matrix(np, ncol=length(np))

                        if(nrow(path_M_dec) == 0) {
                            cat("ERROR after subpath decomposition (empty)!")
                            return();
                        }
 
                        # 
                        p_sort <- pa_check_pathway_present(path_M_new, path_rates_new, path_M_dec)
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
                                     # The pathway has to be added
                                     x <- abs((path_M_dec[k,]*rt/rates)[path_M_dec[k,] != 0 & flag_red])
                                     if(length(x) == 0)
                                         x <- abs((np*nr/rates)[np != 0 & !flag_red])     

                                     mm <- max(x)

                                     new_comp <- sum((N %*% path_M_dec[k,]) != 0) 
#                                     if(rt > f2 & new_comp <= C1 & p_sort$complexity[k] < C2 | mm > pmin) {

                                     if(mm > pmin) {
                                        path_rates_new[length(path_rates_new)+1] <- rt
                                        path_M_new <- rbind(path_M_new, matrix(path_M_dec[k,], 1, length(path_M_dec[k,])))
                                    } else {
                                        rates_dropped <- rates_dropped + rt*path_M_dec[k,]
                                    }
                                }
                            } 
                        }
                                     
                    } else {  # Drop pathway combination
                        rates_dropped <- rates_dropped + nr*np
                    }
                }
            }

            # swap / report / ready for next step  
            path_rates <- path_rates_new
            path_M <- path_M_new
        }

        cat("Having ", nrow(path_M), " pathways! (", count, "/",nrow(N_in),")\n")

        # calculate maximal relative deviation of reaction rates
        rec_ra <- pa_reconstruct_rates(path_M, path_rates)
        mm_a <- max(abs(rec_ra-rates)/rates, na.rm=T)
        mm_b <- max((abs(rec_ra-rates)/rates)[flag_red], na.rm=T)
        mm_c <- max(abs(rec_ra-rates), na.rm=T)
        mm_d <- sum(abs(rec_ra-rates))/sum(rates)

        # print check:  maximal relative deviation;  reaction id; absolute        
        cat("check: ",  mm_a, "  - ", mm_b , "  - ",  mm_c , "  - ", mm_d,  "\n")
        count <- count + 1
    }

    return(list(path_M, path_rates, mm_b, mm_d))
} 


# The function checks in a matrix of elementary modes if they really are elementary
# (by decomposing it with pa_subpath_decomposition) and removes multiple occurences.

pa_rm_multiple_check_elementary <- function(net, ems) {
    # cut these reactions of elementary modes that are not in reaction network
    ems <- ems[,1:nrow(net[[2]])]

    # new matrix / ems that are checked are put here
    em_new <- ems[c(),]
    N <- jrnf_calculate_stoich_mat(net)     # stoichiometric matrix

    for(i in 1:nrow(ems)) {
        cat(".")
        path_M_dec <- pa_subpath_decomposition(N, ems[i,])

        # Show warning if decomposition into elementary modes gives more than em back.
        # This doesn't mean that ems[i,] was not an elementary mode. Could be the same
        # set of reactions is involved in different elementary modes (different coefficients)!
        if(nrow(path_M_dec) != 1)
            cat("WARNING: Mode number ", i, " may not have been elementary!\n")

        # If not already in em_new add the pathway
        for(j in 1:nrow(path_M_dec)) {
            cpp <- pa_check_pathway_present(em_new, rep(1, nrow(em_new)), matrix(path_M_dec[j,],1))
            if(!any(cpp$present))  # adding
                em_new <- rbind(em_new, path_M_dec[j,])
        }
    }
    cat("\n\n")

    return(em_new)
}


# This function orders the elementary modes and generates a data frame containing information
# on how much the reordered elementary modes explain the rate vector 'v'. This is done es well
# cumulative (only those parts not explained by previous modes) as well as individual. It is 
# also interesting what is the average fraction and the minimum fraction of reaction rates. 
# Here 'individual' fraction means the maximum fraction of this and all previous elementary modes.
# min_R ( max_{1...i} ( f ) )
#
# For every elementary mode also it's complexity is put into the data frame. 

pa_iterate_rates_max <- function(em, v, order="max", net, coef=c()) {
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


    if(length(coef) == 0)
        a_init <- apply(em, 1, al)
    else
        a_init <- coef

    # different modes in how the elementary modes are ordered for developing v / v_    
    # "max": after each step the maximum mode best for the next step is selected
    if(order == "max_step") {
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
    } else if(order == "max_init") {
        a <- apply(em, 1, al)
        m_id <- order(a, decreasing=T)

        for(i in m_id) { 
            add_element(i, a[i], NA);
            cat(".")
        }

    # -: pathways are ordered as given / if coef is given it is used as accumulative coefficients
    } else {
        for(i in 1:nrow(em)) {
            add_element(i, NA, a_init[i])
            cat(".")
        }
    }

    # return data frame with results
    return(data.frame(em_id=em_id, coeff=coeff, coeff_acc=coeff_acc, exp_f=exp_f, 
                      exp_f_acc=exp_f_acc, C1=C1, C2=C2, C3=C3, min_f=min_f, min_f_acc=min_f_acc))  
}



# Function derives properties of all elementary modes in a matrix and returns
# them in a data frame. The function returns them 
#
#

pa_em_derive <- function(em_matrix, net, c_max=15) {
    # generating empty vectors
    C <- matrix(0, ncol=c_max, nrow=nrow(em_matrix))     # Number of cycles of different length
    C_s <- matrix(0, ncol=c_max, nrow=nrow(em_matrix))   # Counting cycles by subgraph isomorphism (nodes!)
    Deg <- rep(0, nrow(em_matrix))           # mean degree of all associated species
    Deg_int <- rep(0, nrow(em_matrix))       # mean degree of all associated species ignoring other reactions
    Deg_max <- rep(0, nrow(em_matrix))       # max degree of all associated species
    Deg_max_int <- rep(0, nrow(em_matrix))   # max degree of all associated species ignoring other reactions
    Sp_no <- rep(0, nrow(em_matrix))         # number of species taking part in elementary mode
    Re <- rep(0, nrow(em_matrix))            # number of reactions taking part in elementary mode
    Re_s <- rep(0, nrow(em_matrix))          # number of reactions (counting each only once)
    Ex <- rep(0, nrow(em_matrix))            # Number of exchanges with environment
    Ex_s <- rep(0, nrow(em_matrix))          # Number of exchanges (counting each species only once)
    In <- rep(0, nrow(em_matrix))            # Number of input from environment
    In_s <- rep(0, nrow(em_matrix))          # Number of input (counting each species only once)
    Out <- rep(0, nrow(em_matrix))           # Number of output to environment
    Out_s <- rep(0, nrow(em_matrix))         # Number of output to environment (each species only once) 


    # calculating degree of <net> independently of em's
    net_deg <- as.vector(degree(jrnf_to_undirected_network(net)))
    N <- jrnf_calculate_stoich_mat(net)

    # iterating all pathways
    for(i in 1:nrow(em_matrix)) {
        cat(".")

        rev <- c()                   # is the reaction reverted in the specific elementary mode
        em_abs <- abs(em_matrix[i,]) 
        sel <- c()                   # which reactions are in the elementary mode (id's of ones with 
        sel_s <- c()                 # higher coefficients are put there multiple times)

        for(j in which(em_abs != 0)) {
            sel_s <- c(sel_s, j)

            for(k in 1:em_abs[j]) {
                sel <- c(sel, j)
                rev <- c(rev, em_matrix[i,j] < 0)
            }           
        }
            
        # build temporary networks
        net_tmp <- list(net[[1]], net[[2]][sel,])
        net_tmp_s <- list(net[[1]], net[[2]][sel_s,])
        N_tmp <- jrnf_calculate_stoich_mat(net_tmp)  
        N_tmp_in <- jrnf_calculate_stoich_mat_in(net_tmp)  
        N_tmp_out <- jrnf_calculate_stoich_mat_out(net_tmp)  

        # generate stoichometric matrix with all 1 species reactions (inflow / outflow) removed
        N_int <- N
        for(j in which(apply(N_int != 0, 2, sum) == 1)) 
            N_int[,j] <- 0
  
        em_bilance <- N_int %*% em_matrix[i,]        
        

        # find cycles
        g_tmp <- jrnf_to_directed_network_d(net_tmp, rev)
        for(k in 1:c_max) {
            C[i,k] <- get_n_cycles_directed(g_tmp,k)[[1]]
            C_s[i,k] <- get_n_cycles_directed_A(g_tmp,k)[[1]]
        }
            
        # calculate various other complexity measures
        species_con <- which(apply(N_tmp_in != 0, 1, any) | apply(N_tmp_out != 0, 1, any))
        degree_global <- net_deg[ species_con ]
        degree_internal <- as.vector(degree( jrnf_to_undirected_network(net_tmp) ))
       

        Deg[i] <- mean( degree_global )
        Deg_int[i] <- mean( degree_internal )
        Deg_max[i] <- max( degree_global )
        Deg_max_int[i] <- max( degree_internal )
        Sp_no[i] <- length( species_con )

        Re[i] <- nrow( net_tmp[[2]] )
        Re_s[i] <- nrow( net_tmp_s[[2]] )
        Ex[i] <- sum(abs(em_bilance))
        Ex_s[i] <- length(which(em_bilance != 0))
        In[i] <- sum(abs(pmin(em_bilance,0)))
        In_s[i] <- length(which(abs(pmin(em_bilance,0)) != 0))
        Out[i] <- sum(abs(pmax(em_bilance,0)))
        Out_s[i] <- length(which(abs(pmax(em_bilance,0)) != 0))
    }    

    # sum of cycles is calculated
    C_sum <- apply(C, 1, sum)
    C_s_sum <- apply(C_s, 1, sum)


    # assemble data frame and return it
    return(data.frame(C_sum=C_sum, C_s_sum=C_s_sum, C_1=C[,1], C_2=C[,2], C_3=C[,3], C_1_s=C_s[,1], 
                      C_2_s=C_s[,2], C_3_s=C_s[,3], Deg=Deg, Deg_int=Deg_int, Deg_max=Deg_max, 
                      Deg_max_int=Deg_max_int, Sp_no=Sp_no, Re=Re, Re_s=Re_s, Ex=Ex, Ex_s=Ex_s,
                      In=In, In_s=In_s, Out=Out, Out_s=Out_s))
}



# Uses an existing expansion of a network <net>'s steady state flow <v> and 
# calculates how good it is absolute (using all elementary modes) and cummulative
# (using elementary modes 1 to <n>). The expansion / elementary modes have to be
# ordered in advance. If coefficients of em's are negative this means that the 
# reactions of the network are reversed before starting evaluating the expansion.
# For this it is checked if there is no sign missmatch for all the pathways which 
# have <em_rates> > 0.

pa_ana_expansion <- function(em_matrix, em_rates, net, v, em_df=c()) {
    # First check the sign condition if rates of elementary modes are nonzero all coefficients of 
    # the respective elementary modes have to have the same sign as the reactions rate (v)
    

    # If sign conditions are fulfilled reactions, reaction rates and coefficients of 
    # elementary modes are reversed for negative rate reactions.
    
    
    v_sum <- sum(v) 
    v_ <- rep(0, length(v))                # exansion for ems 1 to <n>

    coeff <- em_rates                      # coefficient of the em
    exp_f <- rep(0, length(em_rates))      # fraction of all rates explained by current em            
    exp_f_acc <- rep(0, length(em_rates))  # fractions of all rates explained by ems 1 to <n>
    err_f_50p <- rep(1, length(em_rates))  # fraction of reaction which is explained by less than 50% by ems 1 to <n>
    err_rmax_50p <- rep(0, length(em_rates)) # highest rate of reaction which is explained by less than 50% by ems 1 to <n> (argmax)
    err_rates <- list()                    # list of data frames. Each data frame contains absolute and relative error by the expansion for each reaction  

    #
    for(i in 1:length(em_rates)) {
        cat(".")
        dv <- em_rates[i]*em_matrix[i,]
        v_ <- v_ + dv
        exp_f[i] <- sum(dv)/v_sum      
        exp_f_acc[i] <- sum(v_)/v_sum
        # total error 
        err_abs <- v_ - v
        err_rel <- abs((v_-v)/v)
        err_rel[v == 0] <- 0

        err_f_50p[i] <- length(which(err_rel > 0.5)) / length(err_rel)
        err_rmax_50p[i] <- max(c(v[err_rel > 0.5],0))

        err_rates[[i]] <- data.frame(v=v, err_abs=err_abs, err_rel=err_rel, v_acc=v_, dv=dv) 
    }

    return(data.frame(exp_f=exp_f, exp_f_acc=exp_f_acc,                                                   
                      err_f_50p=err_f_50p, err_rmax_50p=err_rmax_50p, 
                      err_rates=I(err_rates)))
}


