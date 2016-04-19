# author: jakob fischer (jakob@automorph.info)
# description: 
# testing-tool for algorithm which determines pathways / elementary flux modes
# of reaction networks and a given flux (stea vector v. All functions of this 
# file start with prefix "pa_"

sourced_pathway_analysis <- T

if(!exists("sourced_tools"))
    source("tools.R") 

if(!exists("sourced_cycles"))
    source("cycles.R")

if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")

if(!exists("sourced_composition_analysis"))
    source("composition_analysis.R")

if(!exists("sourced_pathway_analysis_tools"))
    source("pathway_analysis_tools.R")

if(!exists("sourced_pathway_analysis_eval"))
    source("pathway_analysis_eval.R")


# Given a reaction network <N_orig> and a subset of active reactions <path_orig>
# this function calculates all elementary flux modes in the subset of reactions.
# This function follows the original efm-algorithm (TODO: source) and can only be
# used for small sub-networks. It is used as a helper for the function which is usable
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






# This is the older pathway decomposition function that decomposes the
# steady state in one function call and uses "pa_subpath_decomposition"
# as a helper.
# The newer one is implemented in pa_initialize and pa_step and uses 
# "pa_decompose" for decomposition. A "one-function call" can be made 
# through pa_analysis
# 
# parameters:
# net  -  the network that is analysed
# rates  -  the reaction rates of the network
# TODO This function might be improved by finding a more intelligent criterion 
#      on which pathways to drop. This is not a simple implementation issue
#      and would require some research and testing.

pa_analysis_legacy <- function(net, rates, fexp=0.1, pmin=0.01, do_decomposition=T) {
    # Flag those rates that the reduction condition is applied to, even if this is
    # actually all reactions this makes still sense as some rates could
    # be zero and you want to exclude them.
    flag_red <- rates > 0

    # Accounting for the reaction rates of the discarded pathways
    rates_dropped <- rep(0, length(rates))
    # Array containing flags showing which species have been used for branching
    sp_br_flag <- rep(F, nrow(net[[1]]))
    
    # Stoichiometric matrix of consumption (N_in), production (N_out) and total (N)
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)
    N <- N_out - N_in

    # Rates of pathways + matrix for pathways (one row = one pathway) 
    # Initially every reaction equals one pathway with the pathways initial rate 
    # being the rate of the reaction.
    path_rates <- rates
    path_M <- matrix(0, nrow(net[[2]]), nrow(net[[2]]))
    for(i in 1:nrow(net[[2]]))
        path_M[i,i] <- 1

    # reactions / pathways with rate of zero have to be removed
    path_rates <- path_rates[rates != 0]
    path_M <- path_M[rates != 0,]

    # maybe take generation rate for ordering?  or degree of substrate graph?
    turnover <- as.vector(N_in %*% rates)
    count <- 1

    # All species are taken as branching species (lower turnover rate first)
    for(i in order(turnover, decreasing=F)) {
        cat("branching at species ", net[[1]]$name[i], "\n")
        sp_br_flag[i] <- T

        # net creator
        x <- N %*% t(path_M)
        cr <- which(x[i,] > 0)
        cr_m <- x[i,cr]
        turnover <- sum(cr_m*path_rates[cr])
        cr_f <- (cr_m*path_rates[cr])/turnover 
       
        # net consumer
        con <- which(x[i,] < 0)
        con_m <- -x[i,con]
        con_f <- con_m*path_rates[con]/turnover

        z <- cr_f %*% t(con_f)
        n <- z > fexp

        #cat("length(cr)=", length(cr), "  length(con)=", length(con), "\n")

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
            path_M <- matrix(path_M[x[i,] == 0,], ncol=ncol(path_M))
            path_rates <- path_rates[x[i,] == 0]

        # Species is produced and consumed. 
        } else {
            # Ensure that in each row and each column there is at least one TRUE
            # (implies that each pathway is at least connected to one other pathway)
            for(src in 1:length(cr)) {
                a <- z[src,]/sum(z[src,])
                if(length(which(a > fexp | a > 1/(length(a)+1))) == 0) 
                    cat("ERROR: assortion violated (A)!\n")                 

                n[src, a > fexp | a > 1/(length(a)+1)] <- T
            }

            for(dest in 1:length(con)) {
                b <- z[,dest]/sum(z[,dest])

                if(length(which(b > fexp | b > 1/(length(b)+1))) == 0) 
                    cat("ERROR: assortion violated (A)!\n")

                n[b > fexp | b > 1/(length(b)+1),dest] <- T
            }

            # Start filling new matrix with pathways that don't net produce or consume i
            path_M_new <- matrix(path_M[x[i,] == 0,], ncol=ncol(path_M))
            path_rates_new <- path_rates[x[i,] == 0]
                                       
            # Add contribution of dropped / deleted pathways to vector rates_dropped
            for(src in 1:length(cr)) {
                K_src <- cr_m[src]             # How does the pathway src change i
                cat(".")
                for(dest in 1:length(con)) {
                    K_dst <- con_m[dest]       # How does the pathway dest change i 
                    nr <- z[src,dest]*turnover/(K_dst*K_src)                # new rate
                    np <- K_dst*path_M[cr[src],] +                          # new pathway
                          K_src*path_M[con[dest],]  

                    # Calculate the minimal contribution of this (coupled) pathway to the total rate
                    # If it is above <pmin> the two pathway are coupled even if the matrix <n> for it
                    # is false.

                    x <- abs((np*nr/rates)[np != 0 & flag_red])
                    if(length(x) == 0)
                        x <- abs((np*nr/rates)[np != 0 & !flag_red])                        

                    mm <- max(x)

                    if(n[src,dest] | mm > pmin) {                                 # Add combined pathway
                        if((N %*% np)[i] != 0) 
                            cat("ERROR: assortion violated (C)!\n")

                        # TEST if speedup if check if in list BEFORE doing subpath decomposition (TODO CHECK / REMOVE (NOTE))
                        # before kasting 220 with 
                        # system.time(x <- pa_analysis(net_kasting_220_rr_ext, v_kasting_220_rr_ext, 0.1, 1e-5))
                        # took 21406.46 seconds and resulted in 16247 pathways
                        #if(any(apply(path_M_new == np, 1, all)) || !do_decomposition)
                        {
                            # First calculate pathway decomposition
                            if(!do_decomposition)
                                path_M_dec <- matrix(np, nrow=1)    # don't do subpath stuff (SLOW)
                            else
                                path_M_dec <- pa_subpath_decomposition(N, np, sp_br_flag)
                        }

                        if(nrow(path_M_dec) == 0) 
                            cat("ERROR: assortion violated (D)!\n")
 
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

                                if(is.na(rt))  
                                    cat("ERROR: assortion violated (E)!\n")
       
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

            # swap (ready for next step ) 
            path_rates <- path_rates_new
            path_M <- path_M_new
        }

        cat("Having ", nrow(path_M), " pathways! (", count, "/",nrow(N_in),")\n")

        # calculate maximal relative deviation of reaction rates

        rec_ra <- pa_reconstruct_rates(path_M, path_rates)
        mm_a <- max(abs(rec_ra-rates)/rates, na.rm=T)
        mm_c <- max(abs(rec_ra-rates), na.rm=T)
        mm_d <- sum(abs(rec_ra-rates))/sum(rates)

        # print check:  maximal relative deviation;  reaction id; absolute        
        cat("check: ",  mm_a, "  - ",  mm_c , "  - ", mm_d,  "\n")
        count <- count + 1
    }

    return(list(path_M, path_rates, mm_a, mm_d))
} 






# --------------------------------------------------------------------------------
# CODE SPECIFIC TO NEW PATHWAY ANALYSIS CODE BELOW
# --------------------------------------------------------------------------------
#
# source("~/tmp/test.R", chdir=T)
# PA_new_10_T <- pa_analysis(net_rr_ext, rates_rr_ext, 0.5, 1e-3, T)     - 12
# PA_new_10_F <- pa_analysis(net_rr_ext, rates_rr_ext, 0.5, 1e-3, F)     - 12
# PA_new_20_T <- pa_analysis(net_rr_ext, rates_rr_ext, 0.5, 0, T)        - 183
# PA_new_20_F <- pa_analysis(net_rr_ext, rates_rr_ext, 0.5, 0, F)        - 183
# PA_classic_10_T <- pa_analysis_legacy(net_rr_ext, rates_rr_ext, 0.5, 1e-3, T)   - 2592
# PA_classic_10_F <- pa_analysis_legacy(net_rr_ext, rates_rr_ext, 0.5, 1e-3, F)   - 2838
# PA_new_01_T <- pa_analysis(net_rr_ext, rates_rr_ext, 0, 1e-1, T)       - 5
# PA_new_01_F <- pa_analysis(net_rr_ext, rates_rr_ext, 0, 1e-1, F)       - 5
# PA_classic_01_F <- pa_analysis_legacy(net_rr_ext, rates_rr_ext, 0, 1e-1, F)     - 399
# PA_classic_01_T <- pa_analysis_legacy(net_rr_ext, rates_rr_ext, 0, 1e-1, T)     - 391




pa_backlog_apply <- function(bl, x) {
    # if there is no backlog yet there is nothing to apply
    if(is.null(bl))
        return(x)

    M <- x$M
    coef <- x$coef

    #print(M)
    #print(bl$key)

    # first find backlog keys and pathway entries that overlap
    bl_ol <- duplicated(rbind(M, bl$key))[-(1:nrow(M))]    
    M_ol <- duplicated(rbind(bl$key, M))[-(1:nrow(bl$key))]

    M_ <- matrix(M[!M_ol,], ncol=ncol(M))   # Nothing to do for these
    coef_ <- coef[!M_ol]

    M <- matrix(M[M_ol,], ncol=ncol(M))
    coef <- coef[M_ol]    

    # now find for each element in M the (id of the) exact matching key
    sel <- rep(1, nrow(M))
    
    if(length(sel) > 0) {
        for(i in 1:length(sel)) 
            while(!all(bl$key[which(bl_ol)[sel[i]],] == M[i,]))
                sel[i] <- sel[i] + 1

        x0 <- scale_mat_rows(matrix(bl$coef[which(bl_ol)[sel],], ncol=ncol(bl$coef)), coef)
        x1 <- apply(x0, 2, sum) 
        x2 <- apply(matrix(bl$sel[which(bl_ol)[sel],], ncol=ncol(bl$sel)), 2, any)
    
        coef__ <- x1[x2]
        M__ <- bl$M[x2,]

        coef_ <- c(coef_, coef__)
        M_ <- rbind(M_, M__)
    } 

    y <- pa_rm_duprows_accum(M_, coef_)


    return(list(M=M_[y$keep,], coef=y$coef))
}


pa_backlog_add <- function(bl, key, x) {
    if(is.null(bl)) {
        # If backlist is empty, generating a new one is easy

        # Only thing to consider: key is part of it's own representation?
        tmp <- which(duplicated(rbind(matrix(key, nrow=1), x$M))[-1])

        if(length(tmp) == 1) {
            sc <- 1/(1-x$coef[tmp])
            x$M <- x$M[-tmp,]
            x$coef <- x$coef[-tmp]*sc
        }

        bl <- list(key=matrix(key, nrow=1),
                   M=x$M, coef=matrix(x$coef, nrow=1), 
                   sel=matrix(rep(T, length(x$coef)), nrow=1))
    } else {
        # First check that the new key isn't already in the backlist
        if(any(duplicated(rbind(bl$key, key)))) {
            cat("WARNING: pa_backlog_add - added key multiple times!\n")
            return(bl)
        }

        # Then apply the backlist to the representation x
        x <- pa_backlog_apply(bl, x)

        # Key might be used in representation (remove and rescale)
        tmp <- which(duplicated(rbind(key, x$M))[-1])
         
        if(length(tmp) == 1) {
            sc <- 1/(1-x$coef[tmp])
            x$M <- x$M[-tmp,]
            x$coef <- x$coef[-tmp]*sc
        }

        # Identify pathways in x that are already in bl$M
        sel <- duplicated(rbind(bl$M, x$M))[-(1:nrow(bl$M))]
        n_sel <- sum(sel)
        n_nsel <- sum(!sel)
        id_sel <- rep(1, n_sel)
        if(n_sel > 0)
            for(i in 1:n_sel) 
                while(!all(x$M[which(sel)[i],] == bl$M[id_sel[i],]))
                    id_sel[i] <- id_sel[i] + 1

        # ADD key, ADD empty columns in bl$sel, bl$coef and ADD pathways to bl$M
        bl$key <- rbind(bl$key, matrix(key, nrow=1))

        if(n_nsel > 0) {
            bl$M <- rbind(bl$M, x$M[!sel,])
            bl$sel <- cbind(bl$sel, matrix(F, ncol=n_nsel, nrow=nrow(bl$sel)))
            bl$coef <- cbind(bl$coef, matrix(0, ncol=n_nsel, nrow=nrow(bl$sel)))
        }

        # Construct new row and add to bl$sel and bl$coef 
        sel_ <- c(rep(F,ncol(bl$sel)-n_nsel), rep(T, n_nsel))
        if(n_sel > 0)
            sel_[id_sel] <- T
        coef_ <- c(rep(0, ncol(bl$coef)-n_nsel), x$coef[!sel]) 
        if(n_sel > 0)
            coef_[id_sel] <- x$coef[sel]

        bl$sel <- rbind(bl$sel, matrix(sel_, nrow=1))
        bl$coef <- rbind(bl$coef, matrix(coef_, nrow=1))


        # If the newly added key is used as representation of the other keys update
        # these representations and remove the key from M (+ delete the columns in sel and coef)
        tmp <- which(duplicated(rbind(matrix(key, nrow=1), bl$M))[-1])
        
         
        if(length(tmp) == 1) {
            # first use the last row's coefficients (just calculated) to insert
            # in all the representations of the other keys
                
            for(i in which(bl$sel[,tmp])) {
                bl$coef[i,] <- bl$coef[i,] + bl$coef[i,tmp]*coef_   
                bl$sel[i,] <- bl$sel[i,] | sel_
            }

            bl$M <- matrix(bl$M[-tmp,], ncol=ncol(bl$M))
            bl$sel <- bl$sel[,-tmp]
            bl$coef <- bl$coef[,-tmp]
        }

        return(bl)
    }
}


# TODO This function need a way to handle numerical problems
# (Especially when used with fuzzy input) no pathway but steady state vector
# these can really slow down the calculation...
#

pa_decompose_plain <- function(N, pw_init , branch_sp, do_backlog=T, cutoff=0) {
    rates <- pw_init
    turnover <- 0.5*as.vector(abs(N) %*% rates) 
    rea_col <- -(1:nrow(N))
    M <- cbind(t(N), diag(ncol(N)))

    #cat("CALL DECOMPOSE PLAIN\n")
    for(m in 1:length(branch_sp)) {
        i <- branch_sp[m]
        #cat("i= ", i, "\n")

        sel_keep <- which(M[,i] == 0)
        sel_produce <- which(M[,i] > 0)
        sel_consume <- which(M[,i] < 0) 


        f <- function(x) {
            x <- x - 1 
            k <- sel_produce[(x %% length(sel_produce)) + 1]
            l <- sel_consume[as.integer(x / length(sel_produce)) + 1] 

            cre <- M[k,i]
            cre_f <- cre*rates[k]/turnover[i]           
            con <- -M[l,i]   
            con_f_x_turnover <- con*rates[l]  

            pw <- con*M[k,] + cre*M[l,]
            gcd_l <- vec_gcd(pw[nrow(N)+(1:ncol(N))])
            pw <- pw / gcd_l 
            new_rate <- cre_f*con_f_x_turnover*(gcd_l/(cre*con))
            if(is.nan(new_rate)) # if new rate is nan (mostly because of turnover == 0 set it to 0
                new_rate <- 0 

            return(matrix(c(pw, new_rate), nrow=1))
        }

        if(length(sel_produce)*length(sel_consume) != 0)
            M_ <- rbind(cbind(matrix(M[sel_keep,], ncol=ncol(M)), matrix(rates[sel_keep], ncol=1)),
                        do.call("rbind", lapply(1:(length(sel_produce)*length(sel_consume)), f)) )
        else
            M_ <- cbind(matrix(M[sel_keep,], ncol=ncol(M)), matrix(rates[sel_keep], ncol=1))

        rates_ <- as.vector(M_[,ncol(M_)])
        M_ <- matrix(M_[,-ncol(M_)], nrow=nrow(M_))

        # Before thinking about not elementary pathways first remove all duplicates
        # (nontrivial because their coefficients / rates have to be added up) 
        x <- pa_rm_duprows_accum(matrix(M_[,rea_col], nrow=nrow(M_)), rates_)
        M_ <- matrix(M_[x$keep,], ncol=ncol(M_))        
        rates_ <- x$coef

        if(cutoff != 0) {
            M_ <- matrix(M_[rates_ > cutoff,], ncol=ncol(M_)) 
            rates_ <- rates_[rates_ > cutoff]         
        }

        backlog <- c()
        repeat{
            keep_pw <- pa_find_elementary_pw(matrix(M_[,rea_col], nrow=nrow(M_)))
            
            if(all(keep_pw) | !do_backlog){
                break
            }

            #cat("having ", length(keep_pw), " pathways!\n")
            #cat(sum(!keep_pw), " have to be dropped.\n")

            for(k in which(!keep_pw)) {
                # recursive call just keep the backlog

                y <- pa_decompose_plain(N, M_[k,rea_col], branch_sp[1:m], F, cutoff)
                
                key_ext <- pa_extend_pathway_representation(M_[k,rea_col],N)
                backlog <- pa_backlog_add(backlog, key_ext, list(M=y$M, coef=y$coef))   
            }

            x <- pa_backlog_apply(backlog, list(M=M_, coef=rates_))
            M_ <- x$M
            rates_ <- x$coef
            
        }

        # And replace old by new (rates + pathways)
        M <- M_
        rates <- rates_

        # Turnover has to be updated to remove part that is cycling in pathways!
        # (This is due to the fact that turnover is used to calculate branching
        #  probability and already closed cycles don't contribute any more.)

        for(i in 1:nrow(N))     # for each species
            turnover[i] <- 0.5*sum(abs(M[,i])*rates)
    }


    return(list(M=M, coef=rates))
}




# BLA 
# TODO: documentation

pa_decompose <- function(N_orig, path_orig, do_backlog=T, branch_all=F, cutoff=0) {
    # calculate the reaction set for which decomposition is done
    sel_rea <- which(path_orig != 0)

    # only species that are created and consumed by above reaction set are 
    # taken as branching species
    N <- matrix(N_orig[,sel_rea], nrow(N_orig), length(sel_rea))
    if(!branch_all)
        branch_sp <- which(apply(N, 1, max) > 0 & apply(N, 1, min) < 0 & N %*% path_orig[sel_rea] == 0)
    else
        branch_sp <- which(apply(N, 1, max) > 0 & apply(N, 1, min) < 0)

    if(length(branch_sp) == 0)
        return(list(M=matrix(path_orig, nrow=1), coef=c(1)))

    # reduce N matrix even further
    N <- matrix(N[branch_sp,], nrow=length(branch_sp))   
    turnover <- pa_calculate_turnover(scale_mat_rows(t(N), path_orig[sel_rea]))

    #if(branch_all)
    #    cutoff <- max(abs(N %*% path_orig[sel_rea])) / (1e4)     
    #else 
    #    cutoff <- 0

    # because there might be reactions like "2A -> 2B" and the subsequent subfunctions
    # only operate onto the effective change of concentration (but pathways coefficients
    # should be integers)
    s <- apply(N, 2, vec_gcd)
    N <- t(scale_mat_rows(t(N), 1/s))

    x <- pa_decompose_plain(N, path_orig[sel_rea]*s, order(turnover), do_backlog, cutoff)

    x$M <- x$M[,-(1:nrow(N))]

    for(i in 1:nrow(x$M)) {
        sk <- x$M[i,] != 0
        p <- prod(s[sk])
        x$M[i,sk] <- x$M[i,sk]*p/s[sk]
        x$coef[i] <- x$coef[i]/p
        sh <- vec_gcd(x$M[i,sk])
        if(sh!=1) {
            x$coef[i] <- x$coef[i]*sh
            x$M[i,sk] <- x$M[i,sk] / sh
        }   
    }
 

    # transform back to space with all reactions in it!
    M_ret <- matrix(0, nrow=nrow(x$M), ncol=ncol(N_orig))
    M_ret[,sel_rea] <- x$M

    return(list(M=M_ret, coef=x$coef))
}



pa_calculate_turnover <- function(M, N=c()) {
    if(!is.null(N))
        M <- matrix(M[,1:nrow(N)], nrow=nrow(M))

    return(pmax(apply(pmax(M, 0), 2, sum),
                apply(pmax(-M, 0), 2, sum)))
}


pa_initialize <- function(net, rates, co_branch=0, co_exp_rea=0, co_exp_turnover=0, prep=F, decompose=T) {
    # if preparation is done here: remove reverse reactions and add pseudoreactions
    if(prep) {
        x <- jrnf_remove_reverse_pairs(net, rates)
        net_rr <- x[[1]]
        rates_rr <- x[[2]]
    
        x <- pa_extend_net(net_rr, rates_rr)
        net <- x[[1]]
        rates <- x[[2]]
    }

    N <- jrnf_calculate_stoich_mat(net)
    
    M_init <- cbind(t(N), diag(ncol(N)))
    rates_active <- rates
    rates_dumped <- rep(0, length(rates))

    turnover <- pa_calculate_turnover(scale_mat_rows(t(N),rates))

    turnover_active <- turnover
    turnover_eliminated <- rep(0, length(turnover))

    coefficients <- rates 
    sp_planned <- order(turnover)
    sp_done <- c()
    
    active_f <- rep(T, nrow(M_init))
     
    # The pathway explains a fraction of the steady state of each reaction 
    # exp_rates is the maximum of these values over all reactions 
    exp_rates <- rep(1, nrow(M_init))
    # throughput is is the (again maximum) fraction of the throughput of all
    # species that is explained by this pathway. This includes generated or
    # consumed species.
    calc_tp <- function(x) {  
        sp <- x[1:nrow(N)]
        r <- x[length(x)]
        return(max(c(0, (abs(sp*r)/turnover)[sp != 0] ), na.rm=T))  
    }
    exp_throughput <- apply(cbind(M_init, matrix(rates,ncol=1)), 1, calc_tp)
    exp_interact <- exp_throughput

    return(list(pathways=list(M=M_init, coefficients=coefficients,  
                              active_f=active_f, exp_rates=exp_rates, 
                              exp_turnover=exp_throughput, exp_interact=exp_interact),
                state=list(rates_active=rates_active, rates_dumped=rates_dumped,
                           turnover_active=turnover_active, sp_done=sp_done, sp_planned=sp_planned,
                           turnover_eliminated=turnover_eliminated), 
                parameters=list(co_branch=co_branch, net=net, N=N, 
                                co_exp_rea=co_exp_rea, co_exp_turnover=co_exp_turnover, 
                                rates=rates, turnover=turnover, do_decompose=decompose)))
}


pa_limit_sp_planned <- function(x, bound=1e-20) {
    sel <- as.logical(abs(x$parameters$N %*% x$parameters$rates) < bound)
    x$state$sp_planned <- x$state$sp_planned[x$state$sp_planned %in% which(sel)]

    # TODO REMOVE NEXT LINE / just for testing
    x$state$sp_planned <- sort(x$state$sp_planned, decreasing=T)

    return(x)
}


pa_is_done <- function(obj) {
    if(length(obj$state$sp_planned) == 0)
        return(T)
    else 
        return(F)
}




pa_step <- function(obj, i=c()) {
    gc()

    if(is.null(i))
        if(is.null(obj$state$sp_planned) || size(obj$state$sp_planned) == 0) {
            cat("ERROR: pa_step called with no species to collapse!\n")
            return(obj)  
        } else {
            i <- obj$state$sp_planned[1]
        } 

    # separate pathways in ones that are kept, one that produce species i and 
    # those that consume it
    sel_keep <- which(obj$pathways$M[,i] == 0 | !obj$pathways$active_f)
    sel_produce <- which(obj$pathways$M[,i] > 0 & obj$pathways$active_f)
    sel_consume <- which(obj$pathways$M[,i] < 0 & obj$pathways$active_f)

    # temporary values (easier access)
    N <- obj$parameters$N
    net <- obj$parameters$net

    # Build a combined matrix of M, rate (coefficient), active flag,  
    M <- obj$pathways$M   
    rates <- obj$pathways$coefficient
    turnover <- obj$state$turnover_active   
    M_ext <- cbind(M, 
                   matrix(obj$pathways$coefficient, ncol=1),
                   matrix(obj$pathways$exp_rates, ncol=1),
                   matrix(obj$pathways$exp_turnover, ncol=1),
                   matrix(obj$pathways$exp_interact, ncol=1),
                   matrix(obj$pathways$active_f, ncol=1),
                   matrix(1, ncol=1, nrow=nrow(M)))

    # APPLY calculations
    # first combine all of sel_produce with all of sel_consume
    # TODO: MEANINGFUL NAME
    f <- function(x) {
        x <- x - 1 
        k <- sel_produce[(x %% length(sel_produce)) + 1]
        l <- sel_consume[as.integer(x / length(sel_produce)) + 1] 

        cre <- M[k,i]
        cre_f <- cre*rates[k]/turnover[i]           
        con <- -M[l,i]   
        con_f_x_turnover <- con*rates[l]  

        pw <- con*M[k,] + cre*M[l,]
        gcd_l <- vec_gcd(pw[nrow(N)+(1:ncol(N))])
        pw <- pw / gcd_l 
        new_rate <- cre_f*con_f_x_turnover*(gcd_l/(cre*con))
        if(is.nan(new_rate)) # if new rate is nan (mostly because of turnover == 0 set it to 0
            new_rate <- 0 

        branch_p <- cre_f*con_f_x_turnover/turnover[i]
        active <- T
        # TODO include creation / consumtion fraction::: branch_p > obj$parameters$co_branch
        #cat("bp=", branch_p, "\n")

        return(matrix(c(pw, new_rate, NA, NA, NA, active, 0), nrow=1))
    }

    # Calculate explained rates, turnover and interaction for all
    # rows in which at least one of those is NA
    # TODO GIVE MEANINGFUL NAME
    g <- function(x) {
        exp_rates <- x[nrow(N)+ncol(N)+2]
        exp_turnover <- x[nrow(N)+ncol(N)+3]
        exp_interact <- x[nrow(N)+ncol(N)+4]
        
        if(any(is.na(exp_rates),
               is.na(exp_turnover),
               is.na(exp_interact)) | T) {
            sp_c <- abs(x[1:nrow(N)])
            re_c <- x[nrow(N)+1:ncol(N)]
            rate <- x[nrow(N)+ncol(N)+1]
            active <- x[nrow(N)+ncol(N)+5]
            y <- pa_calculate_turnover(scale_mat_rows(t(N),re_c))

            exp_rates <- max(c(0, (re_c*rate/obj$parameters$rates)[re_c != 0]) ,na.rm=T)
            exp_turnover <- max(c(0, (abs(sp_c*rate)/obj$state$turnover_active)[sp_c != 0] ), na.rm=T)
            exp_interact <- max(c(0, (abs(y*rate)/obj$parameter$turnover)[y != 0] ),na.rm=T)
            
            x[nrow(N)+ncol(N)+2] <- exp_rates 
            x[nrow(N)+ncol(N)+3] <- exp_turnover
            x[nrow(N)+ncol(N)+4] <- exp_interact

            if(active | is.na(active)) {
                x[nrow(N)+ncol(N)+5] <- exp_rates > obj$parameters$co_exp_rea 
                if(exp_rates <= obj$parameters$co_exp_rea)
                    cat("dumped: exp_rates=", exp_rates, "\n")
            } 
        } 
        return(x)
    }

    # Decompose pathway to elementary pathways
    # TODO GIVE MEANINGFUL NAME
    h <- function(x) {
        # inactive (don't decompose)
        if(x[nrow(N)+ncol(N)+5] < 0.5)
            return(matrix(x, ncol=length(x)))

        # also don't decompose if it has been checked already
        if(x[nrow(N)+ncol(N)+6] > 0.5)
            return(matrix(x, ncol=length(x)))

        active <- x[nrow(N)+ncol(N)+5]
        coef <- x[nrow(N)+ncol(N)+1]

        rt <- x[nrow(N)+1:ncol(N)]
        tmp <- pa_decompose(N, rt)

        if(nrow(tmp$M) <= 1)         
            return(matrix(x, ncol=length(x)))

        return(cbind(t(N %*% t(tmp$M)), 
                     tmp$M,
                     matrix(tmp$coef*coef, ncol=1),
                     matrix(NA, ncol=3, nrow=nrow(tmp$M)),
                     matrix(active, ncol=1, nrow=nrow(tmp$M)),
                     matrix(1, ncol=1, nrow=nrow(tmp$M))))
    }

    cat("working ", obj$parameters$net[[1]]$name[i], "(", length(obj$state$sp_done)+1, "/") 
    cat(length(obj$state$sp_planned)+length(obj$state$sp_done), ") - combining ", length(sel_produce))
    cat("x", length(sel_consume), "=",length(sel_produce)*length(sel_consume),  "pathways!\n")
    cat("leaving ", sum(M_ext[sel_keep,ncol(N)+nrow(N)+5]>0.5), " active pathways unchanged - ", sum(M_ext[sel_keep,ncol(N)+nrow(N)+5]<0.5), " inactive.\n")

    if(length(sel_produce)*length(sel_consume) != 0) {
        M_ext <- rbind(matrix(M_ext[sel_keep,], ncol=ncol(M_ext)),
                       do.call("rbind", lapply(1:(length(sel_produce)*length(sel_consume)), f)) )
    } else if(length(sel_produce) != 0 || length(sel_consume) != 0) {
        M_ext[sel_produce,ncol(N)+nrow(N)+5] <- 0
        M_ext[sel_consume,ncol(N)+nrow(N)+5] <- 0       

        cat("=> Species is not consumed and produced: ", length(sel_produce) + length(sel_consume),  
            "pathways are dumped!\n")
    } 
        


    rm_duprows <- function(M) {
        if(nrow(M) < 1)
            return(M)

        # first split the matrix (M) in a part that doesn't contain duplicates and one that doesn't
        uniq_sec <- nrow(N)+1:ncol(N)
        dups <- duplicated(M[,uniq_sec]) | duplicated(M[nrow(M):1,uniq_sec])[nrow(M):1]
        M_ <- matrix(M[dups,], ncol=ncol(M))
        M <- matrix(M[!dups,], ncol=ncol(M))

        # No duplicates, just return M
        if(nrow(M_) == 0)
            return(M)

        order_matrix <- function(x) {  return(do.call(order, lapply(1:ncol(x), function(i) x[, i])))  }
 
        # Order matrix and build vector with indices of the rows that are non duplicates.
        # These are now essentially the ones that are kept. The ones inbetween are duplicates
        # and their coefficients are just summed up accordingly.
        M_ <- M_[order_matrix(M_[,uniq_sec]),]
        dd <-which(!duplicated(M_[,uniq_sec]))

        # function sums up all the coefficients and combines all parameters for the
        # x-th unique pathway / row in M_
        k <- function(x) {
            first <- dd[x]
            if(x != length(dd))
                last <- dd[x+1]-1
            else
                last <- nrow(M_)

            pw <- M_[first,1:(ncol(N)+nrow(N))]
            new_rate <- sum(M_[first:last, nrow(N)+ncol(N)+1])
            active <- any(M_[first:last, nrow(N)+ncol(N)+5] > 0.5, na.rm=T)
            elementary <- any(M_[first:last, nrow(N)+ncol(N)+6] > 0.5, na.rm=T)

            return(matrix(c(pw, new_rate, NA, NA, NA, as.numeric(active), as.numeric(elementary)), nrow=1))
        }
        #

        if(length(dd) > 0)
            M_ <- do.call("rbind", lapply(1:length(dd), k))  

        return(rbind(M, M_))
    }

    # Take care of duplicates for the first time
    M_ext <- rm_duprows(M_ext)
    M_ext <- t(apply(M_ext, 1, g))
    
    a_F <- sum(M_ext[,nrow(N)+ncol(N)+5] > 0.5)
    t_F <- nrow(M_ext)
    cat("having t_F=", t_F, " with a_F=", a_F, " active!\n")

    if(obj$parameters$do_decompose) {
        # Expand to elementary pathways (those that are still active)
        M_ext <- apply(M_ext, 1, h) 

        # If one call of "h" returns 
        if(is.list(M_ext))                    
            M_ext <- do.call("rbind", M_ext)
        else 
            M_ext <- t(M_ext)
    }

    a_F <- sum(M_ext[,nrow(N)+ncol(N)+5] > 0.5)
    t_F <- nrow(M_ext)
    cat("having t_F=", t_F, " with a_F=", a_F, " active!\n")

    # Remove duplicates for the second time
    M_ext <- rm_duprows(M_ext)
    M_ext <- t(apply(M_ext, 1, g))

    a_F <- sum(M_ext[,nrow(N)+ncol(N)+5] > 0.5)
    t_F <- nrow(M_ext)
    cat("having t_F=", t_F, " with a_F=", a_F, " active!\n")

    M_ext <- matrix(M_ext[M_ext[,nrow(N)+ncol(N)+5] > 0.5,], ncol=ncol(M_ext)) # only keep active ones

    cat("having t_F=", t_F, " with a_F=", a_F, " active!\n")
    
    # copy back to "obj"
    # pathway specific
    obj$pathways$M <- matrix(M_ext[,1:(ncol(N)+nrow(N))], nrow=nrow(M_ext))
    obj$pathways$coefficients <- M_ext[,ncol(N)+nrow(N)+1]
    obj$pathways$exp_rates <- M_ext[,ncol(N)+nrow(N)+2]
    obj$pathways$exp_turnover <- M_ext[,ncol(N)+nrow(N)+3]
    obj$pathways$exp_interact <- M_ext[,ncol(N)+nrow(N)+4]
    obj$pathways$active_f <- as.logical(M_ext[,ncol(N)+nrow(N)+5])
    # Apply additional constraint on active (exp_rates, exp_turnover, exp_interact)
    obj$pathways$active_f <- obj$pathways$active_f & obj$pathways$exp_rates > obj$parameters$co_exp_rea

    # state specific
    obj$state$sp_done <- c(obj$state$sp_done, obj$state$sp_planned[1])
    obj$state$sp_planned <- obj$state$sp_planned[-1]


    tmp <- scale_mat_rows(obj$pathways$M, obj$pathways$coefficients)
    obj$state$rates_active <- apply(matrix(tmp[obj$pathways$active_f, nrow(N)+1:ncol(N)],ncol=ncol(N)), 2, sum)
    obj$state$rates_dumped <- apply(matrix(tmp[!obj$pathways$active_f, nrow(N)+1:ncol(N)],ncol=ncol(N)), 2, sum)
    obj$state$turnover_active <- pa_calculate_turnover( tmp[, 1:nrow(N)] ) 
    obj$state$turnover_eliminated <- abs(obj$parameters$turnover - obj$state$turnover_active)

    return(obj)
}




pa_analysis <- function(net, rates, fexp=0.1, pmin=0.01, do_decomposition=T) {
    
    xx <- pa_initialize(net, rates, pmin, fexp, 0, F, do_decomposition)

    while(!pa_is_done(xx))
        xx <- pa_step(xx)

    return(list(xx$pathways$M[xx$pathways$active_f,], xx$pathways$coefficients[xx$pathways$active_f], xx))
} 






