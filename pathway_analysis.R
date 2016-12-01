# author: jakob fischer (jakob@automorph.info)
# description: 
# Tool with algorithm which determines pathways / elementary flux modes of
# reaction networks. The main functionality for calculating reaction
# pathways is in "pa_initialize" and "pa_step". It can also be called at once
# through "pa_analysis".

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


# Function calculates the turnover for all species. The turnover is defined as
# the total amount of a species that is consumed or produced by the species 
# through all reactions. If total prouction differs from total consumption the
# maximum is returned.
#
# parameters:
# <N>  - Stoichiometric matrix
# <v>  - (steady state) rate vector

pa_calculate_turnover <- function(N, v) {
    N_p <- matrix(pmax(N, 0), ncol=ncol(N))
    N_n <- matrix(pmax(-N, 0), ncol=ncol(N))

    return(as.numeric(pmax(abs(N_p %*% v), abs(N_n %*% v))))
}


#
#
# parameters:
# <N>      - Stoichiometric matrix
# <M>      - Elementary modes / pathway (matrix)
# <alpha>  - Pathways' coefficients

pa_calculate_interaction <- function(N, M, alpha) {
    NM <- N %*% t(M)
    NM_p <- matrix(pmax(NM, 0), ncol=ncol(NM)) 
    NM_n <- matrix(pmax(-NM, 0), ncol=ncol(NM))

    return(as.numeric(pmax(abs(NM_p %*% alpha), abs(NM_n %*% alpha))))
}


# TODO comment + test
#
# parameters:

pa_decompose_plain <- function(N, pw_init , branch_sp) {
    rea_col <- -(1:nrow(N))
    M <- cbind(t(N), diag(ncol(N)))

    for(m in 1:length(branch_sp)) {
        i <- branch_sp[m]

        sel_keep <- which(M[,i] == 0)
        sel_produce <- which(M[,i] > 0)
        sel_consume <- which(M[,i] < 0) 


        f <- function(x) {
            x <- x - 1 
            k <- sel_produce[(x %% length(sel_produce)) + 1]
            l <- sel_consume[as.integer(x / length(sel_produce)) + 1] 

            cre <- M[k,i]         
            con <- -M[l,i]    
            pw <- con*M[k,] + cre*M[l,]
            gcd_l <- vec_gcd(pw[nrow(N)+(1:ncol(N))])
            pw <- pw / gcd_l 
            return(matrix(pw, nrow=1))
        }

        if(length(sel_produce)*length(sel_consume) != 0)
            M_ <- rbind(matrix(M[sel_keep,], ncol=ncol(M)),
                        do.call("rbind", lapply(1:(length(sel_produce)*length(sel_consume)), f)) )
        else
            M_ <- matrix(M[sel_keep,], ncol=ncol(M))

        # Before thinking about not elementary pathways first remove all duplicates
        # (nontrivial because their coefficients / rates have to be added up) 
        x_keep <- !duplicated(matrix(M_,nrow=nrow(M_)))
        M_ <- matrix(M_[x_keep,], ncol=ncol(M_))        
        
        keep_pw <- pa_find_elementary_pw(matrix(M_[,rea_col], nrow=nrow(M_)))

        # And replace old by new (rates + pathways)
        M <- matrix(M_[keep_pw], ncol=ncol(M_))
    }

    return(M)
}




# BLA 
# TODO: documentation + test

pa_decompose <- function(N_orig, path_orig, branch_all=F, rnd_o=F) {
    # calculate the reaction set for which decomposition is done
    sel_rea <- which(path_orig != 0)

    # Only species that are created and consumed by above reaction set are taken 
    # as branching species. 
    N <- matrix(N_orig[,sel_rea], nrow(N_orig), length(sel_rea))
    if(!branch_all)
        branch_sp <- which(apply(N, 1, max) > 0 & apply(N, 1, min) < 0 & N %*% path_orig[sel_rea] == 0)
    else
        branch_sp <- which(apply(N, 1, max) > 0 & apply(N, 1, min) < 0)

    if(length(branch_sp) == 0)
        return(list(M=matrix(path_orig, nrow=1), coef=c(1)))

    # reduce N matrix even further
    N <- matrix(N[branch_sp,], nrow=length(branch_sp))   
    turnover <- pa_calculate_turnover(N, path_orig[sel_rea])

    # because there might be reactions like "2A -> 2B" and the subsequent subfunctions
    # only operate onto the effective change of concentration (but pathways coefficients
    # should be integers)
    s <- apply(abs(N), 2, vec_gcd)
    N <- t(scale_mat_rows(t(N), 1/s))

    # Branch in random order or ordered with increasing turnover?
    ord <- 0
    if(rnd_o)
        ord <- sample(length(branch_sp))
    else 
        ord <- order(turnover, decreasing=T)

    x_M <- pa_decompose_plain(N, path_orig[sel_rea]*s, ord)

    x_M <- matrix(x_M[,-(1:nrow(N))], nrow=nrow(x_M))

    for(i in 1:nrow(x_M)) {
        sk <- x_M[i,] != 0
        p <- prod(s[sk])
        x_M[i,sk] <- x_M[i,sk]*p/s[sk]
        sh <- vec_gcd(x_M[i,sk])
        if(sh!=1) 
            x_M[i,sk] <- x_M[i,sk] / sh
    }
 

    # transform back to space with all reactions in it!
    M_ret <- matrix(0, nrow=nrow(x_M), ncol=ncol(N_orig))
    M_ret[,sel_rea] <- x_M

    return(M_ret)
}


#
#
# parameters:
# <net>             - the reaction network to decompose  
# <rates>           - steady state rate of the network
# <co_branch>       - cutoff parameter for branching (which pw "to combine")
# <co_exp_rea>      - 
# <co_exp_turnover> -
# <prep>            - should network be prepared by remove reverse pairs and 
#                     add pseudo in- / outflow reactions
# <decompose>       - 
# <dump>            -
# <decreasing>      -

pa_initialize <- function(net, rates, co_branch=0, co_exp_rea=0, co_exp_turnover=0, prep=F, decompose=T, dump=T, decreasing=T) {
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

    turnover <- pa_calculate_turnover(N,rates)

    coefficients <- rates 
    sp_planned <- order(turnover, decreasing=decreasing)
    sp_planned <- sp_planned[0 != turnover[sp_planned]]
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
    exp_turnover <- apply(cbind(M_init, matrix(rates,ncol=1)), 1, calc_tp)
    exp_interact <- exp_turnover

    # explanation / definition of pa-structure:
    # pathways (part): state containing pathways or data associated with them
    #    M            : Matrix containing the pathways (each row one pathway)
    #                 : first a integer for each species, followed by an positive
    #                 : integer for each reaction
    #    coefficients : Floating point coefficient (positive) for each pathway
    #    active_f     : Boolean flag if pathway is still active (if do_dump is
    #                 : false pathways that are dropped are kept but set inactive)
    #    exp_rates    : Maximul fraction of a reaction's rate explained by all the
    #                 : pathways (with current set of coefficients)
    #    exp_turnover : Maximum fraction of turnover (creation, consumption,
    #                   creation+consumption - absolut, not effective) explained
    #                   by this pathway (with current set of coefficients)
    #    exp_interact : Maximum fraction of interaction (effective creation or
    #                   consumption) explained by this pathway
    #
    # state (part): state not directly associated with pathways (
    #               but changing with every step of algorithm)
    #    rates_active       : Current expansion of steady state with all active
    #                       : pathways and their coefficients
    #    turnover_active    : Turnover (see above) of current expansion
    #    interaction_active : Interaction (see above) of current expansion
    #    sp_done            : Vector of species already used for branching
    #    sp_planned         : Vector of species that are scheduled for branching 
    # 
    # parameters (part): parameters of algorithm that do not change after init
    #    co_branch       : cutoff parameter for branching probability
    #    co_exp_rea      : cutoff parameter for explained rates
    #    co_exp_turnover : cutoff parameter for explained turnover / interaction
    #    rates           : steady state rates for which expansion is to be calculated
    #    turnover        : turnover of all species (equivalent to interaction at initial state)
    #    do_decompose    : Flag indicates if decomposition (to elementary modes) is done
    #    do_dump         : Flag indicating if dropped (inactive) pathways are dumped (deleted)

    return(list(pathways=list(M=M_init, coefficients=coefficients,  
                              active_f=active_f, exp_rates=exp_rates, 
                              exp_turnover=exp_turnover, exp_interact=exp_interact),
                state=list(rates_active=rates,
                           turnover_active=turnover, interaction_active=turnover,
                           sp_done=sp_done, sp_planned=sp_planned), 
                parameters=list(co_branch=co_branch, net=net, N=N, 
                                co_exp_rea=co_exp_rea, co_exp_turnover=co_exp_turnover, 
                                rates=rates, turnover=turnover, do_decompose=decompose, 
                                do_dump=dump)))
}


# Limits the species that are planned for a pathway analysis object. For this 
# the given rate vector is used to find out which species are not in steady 
# state or have inflow or outflow outside of the reactions described in the network.
# Species are excluded if their change rate (assuming the 'rates' known are a
# steady state) is above 'bound'.

pa_limit_sp_planned <- function(obj, bound=1e-20) {
    sel <- as.logical(abs(obj$parameters$N %*% obj$parameters$rates) < bound)
    obj$state$sp_planned <- obj$state$sp_planned[obj$state$sp_planned %in% which(sel)]

    return(obj)
}


# Function returns true if no further species for branching are available 

pa_is_done <- function(obj) {
    return(length(obj$state$sp_planned) == 0)
}


# Make one step for pathway analysis on pa-object <obj>. Step means that
# pathways (net-) producing species <i> are combined with those that consume
# it. On this way with minor contribution to steady state rate are dropped.
# If no <i> parameter is driven the function takes the next planned species
# of the pa-object.
#
# TODO exp_interact

pa_step <- function(obj, i=c()) {
    gc()  # call garbage collector first because function needs a lot of memory

    if(is.null(i))
        if(is.null(obj$state$sp_planned) || size(obj$state$sp_planned) == 0) {
            cat("ERROR: pa_step called with no species to collapse!\n")
            return(obj)  
        } else {
            i <- obj$state$sp_planned[1]
        } 

    # separate pathways in ones that are kept, those that produce species i and 
    # those that consume it
    sel_keep <- which(obj$pathways$M[,i] == 0 | !obj$pathways$active_f)
    sel_produce <- which(obj$pathways$M[,i] > 0 & obj$pathways$active_f)
    sel_consume <- which(obj$pathways$M[,i] < 0 & obj$pathways$active_f)

    # temporary values (easier access)
    N <- obj$parameters$N
    net <- obj$parameters$net

    # Build a combined matrix of M, rate (coefficient), explained rate, explained 
    # turnover,explained interaction and active flag,  
    M <- obj$pathways$M   
    rates <- obj$pathways$coefficient
    turnover <- obj$state$turnover_active   
    M_ext <- cbind(M, 
                   matrix(obj$pathways$coefficient, ncol=1),
                   matrix(obj$pathways$exp_rates, ncol=1),
                   matrix(obj$pathways$exp_turnover, ncol=1),
                   matrix(obj$pathways$exp_interact, ncol=1),
                   matrix(obj$pathways$active_f, ncol=1),
                   matrix(NA, ncol=1, nrow=nrow(M)))

    # APPLY calculations
    # first combine all of sel_produce with all of sel_consume
    gen_comb_pw <- function(x) {
        k <- sel_produce[mat_lin_to_rowid(x, length(sel_produce), length(sel_consume))]
        l <- sel_consume[mat_lin_to_colid(x, length(sel_produce), length(sel_consume))] 

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
        return(matrix(c(pw, new_rate, NA, NA, NA, 1, branch_p), nrow=1))
    }

     # Calculate explained rates, turnover and interaction for all
    # rows in which at least one of those is NA
    calculate_scores <- function(x) {
        exp_rates <- x[nrow(N)+ncol(N)+2]
        exp_turnover <- x[nrow(N)+ncol(N)+3]
        exp_interact <- x[nrow(N)+ncol(N)+4]
        active <- x[nrow(N)+ncol(N)+5]
        
        # calculate explained fractions if one of them is NA
        if(any(is.na(exp_rates),
               is.na(exp_turnover),
               is.na(exp_interact))) {
            sp_c <- abs(x[1:nrow(N)])
            re_c <- x[nrow(N)+1:ncol(N)]
            rate <- x[nrow(N)+ncol(N)+1]
            y <- pa_calculate_turnover(N,re_c)
            z <- abs(x[1:nrow(N)])

            exp_rates <- max(c(0, (re_c*rate/obj$parameters$rates)[re_c != 0]) ,na.rm=T)
            exp_turnover <- max(c(0, (y*rate/obj$parameters$turnover)[y != 0] ), na.rm=T)
            exp_interact <- max(c(0, (z*rate/obj$parameters$turnover)[sp_c != 0] ),na.rm=T)
            
            x[nrow(N)+ncol(N)+2] <- exp_rates 
            x[nrow(N)+ncol(N)+3] <- exp_turnover
            x[nrow(N)+ncol(N)+4] <- exp_interact
        } 

        # "first run" - while branching probability is available only dropping
        #               in combination with it is legal
        if(!is.na(x[nrow(N)+ncol(N)+6]) && ( active > 0.5 || is.na(active))) {
            x[nrow(N)+ncol(N)+5] <- ( x[nrow(N)+ncol(N)+6] > obj$parameters$co_branch ||
                                      exp_rates > obj$parameters$co_exp_rea || 
                                      exp_turnover > obj$parameters$co_exp_turnover )
        #
        # 
        } else if(active > 0.5 || is.na(active)) {
            x[nrow(N)+ncol(N)+5] <- ( exp_rates > obj$parameters$co_exp_rea || 
                                      exp_turnover > obj$parameters$co_exp_turnover )
        }
        return(x)
    }

    # Decompose pathway to elementary pathways
    decompose_to_elementary <- function(x) {
        # inactive (don't decompose)
        if(x[nrow(N)+ncol(N)+5] < 0.5)
            return(matrix(x, ncol=length(x)))

        # also don't decompose if it has been checked already
        if(is.na(x[nrow(N)+ncol(N)+6]))
            return(matrix(x, ncol=length(x)))

        rt <- x[nrow(N)+1:ncol(N)]
        tmp <- pa_decompose(N, rt, rnd_o=T)

        if(nrow(tmp) <= 1) {        
            x[nrow(N)+ncol(N)+6] <- NA
            return(matrix(x, ncol=length(x)))
        }

        return(cbind(t(N %*% t(tmp)), 
                     tmp,
                     matrix(NA, ncol=4, nrow=nrow(tmp)),
                     matrix(1, ncol=1, nrow=nrow(tmp)),
                     matrix(NA, ncol=1, nrow=nrow(tmp))))
    }


    rm_duprows <- function(M) {
        if(nrow(M) <= 1)
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
            if(all(is.na(M_[first:last, nrow(N)+ncol(N)+5])))
                active <- NA
            elementary <- sum(M_[first:last, nrow(N)+ncol(N)+6])

            return(matrix(c(pw, new_rate, NA, NA, NA, as.numeric(active), as.numeric(elementary)), nrow=1))
        }
        #

        if(length(dd) > 0)
            M_ <- do.call("rbind", lapply(1:length(dd), k))  

        return(rbind(M, M_))
    }


    # simple diagnostic function prints state of intermediate pathways (in X)
    print_diag <- function(X) {
        a_F <- sum(X[,nrow(N)+ncol(N)+5] > 0.5)
        t_F <- nrow(X)
        cat("having t_F=", t_F, " with a_F=", a_F, " active!\n")
    }



    # Parameters that should be considered: 
    #  co_branch, co_exp_rea, co_exp_turnover
    #  do_decompose, do_dump


    # Quite extensive (diagnostics) TODO simplify when function is tested
    cat("working ", obj$parameters$net[[1]]$name[i], "(", length(obj$state$sp_done)+1, "/") 
    cat(length(obj$state$sp_planned)+length(obj$state$sp_done), ") - combining ", length(sel_produce))
    cat("x", length(sel_consume), "=",length(sel_produce)*length(sel_consume),  "pathways!\n")
    cat("leaving ", sum(M_ext[sel_keep,ncol(N)+nrow(N)+5]>0.5), " active pathways unchanged - ", sum(M_ext[sel_keep,ncol(N)+nrow(N)+5]<0.5), " inactive.\n")


    # Function is fundamentally different depending on whether there are actually 
    # pathways to combine or not. If not (especially if the species is only created
    # or only consumed these are simply dumped or dropped depending do_dump
    # 
    if(length(sel_produce)*length(sel_consume) != 0) {
        # calculating branching probability (of fraction that is still exchanged)

        cat("A")
        # Create list of new (combined pathways)
        if(obj$parameters$co_branch > 0) {
            v_in <- abs(M_ext[sel_produce,i])*rates[sel_produce] 
            v_out <- abs(M_ext[sel_consume,i])*rates[sel_consume] 
            v_in <- v_in / sum(v_in)
            v_out <- v_out / sum(v_out)

            o_in <- order(v_in, decreasing=T)
            s_in <- which(cumsum(v_in[o_in]) > (1-obj$parameters$co_branch))[1]
            o_out <- order(v_out, decreasing=T)
            s_out <- which(cumsum(v_out[o_out]) > (1-obj$parameters$co_branch))[1]

            mb <- matrix(F, nrow=length(v_in), ncol=length(v_out))
            mb[,o_out[1:s_out]] <- T
            mb[o_in[1:s_in],] <- T

            cat("_", length(which(as.vector(mb))), "_")
            M_new <- do.call("rbind", lapply(which(as.vector(mb)), gen_comb_pw)) 
        } else {
            M_new <- do.call("rbind", lapply(1:(length(sel_produce)*length(sel_consume)), gen_comb_pw)) 
        }
  
        M_ext <- matrix(M_ext[sel_keep,], ncol=ncol(M_ext))

        cat("B")
        # Remove duplicates of new (combined) pathways.
        M_new <- rm_duprows(M_new)

        # Decompose new pathways and add to existing ones (unchanged ones)
        if(obj$parameters$do_decompose) {

            cat("E") 
            list_result <- lapply(split(M_new,seq(nrow(M_new))),decompose_to_elementary)
            M_new <- do.call(rbind,list_result)
        }

        cat("F-")
        M_ext <- rbind(M_ext, M_new)
        cat(nrow(M_ext))

        cat("-G-")
        # Make unique again
        M_ext <- rm_duprows(M_ext)
        cat(nrow(M_ext))

        cat("-H")
        # Recalculate rates of decomposition from skratch
        act_id <- as.logical(M_ext[,ncol(N)+nrow(N)+5])
        coef_new <- pa_calc_coefficients(M_ext[act_id,nrow(N)+1:ncol(N)], obj$parameters$rates, con_fb=F)
        M_ext[act_id, nrow(N)+ncol(N)+1] <- coef_new$coef
        M_ext[act_id, nrow(N)+ncol(N)+2:5] <- NA

        cat("I")
        # Drop pathways with too low rates / exp fractions
        M_ext <- t(apply(M_ext, 1, calculate_scores))
      
        cat("J\n")
  
    } else if(length(sel_produce) != 0 || length(sel_consume) != 0) {
        # If there is only an input or an output to the branching species these
        # pathways have to be droped (as there is nothing to combine them with).
        M_ext[sel_produce,ncol(N)+nrow(N)+5] <- 0
        M_ext[sel_consume,ncol(N)+nrow(N)+5] <- 0       

        cat("=> Species is not consumed and produced: ", length(sel_produce) + length(sel_consume),  
            "pathways are dropped!\n")
    } 


    if(obj$parameters$do_dump)
        M_ext <- matrix(M_ext[M_ext[,ncol(N)+nrow(N)+5] > 0.5,], ncol=ncol(M_ext)) # only keep active ones

    cat("having", nrow(M_ext), "pathways after dump!\n")
    
    # copy back to "obj"
    # pathway specific
    obj$pathways$M <- matrix(M_ext[,1:(ncol(N)+nrow(N))], nrow=nrow(M_ext))
    obj$pathways$coefficients <- M_ext[,ncol(N)+nrow(N)+1]
    obj$pathways$exp_rates <- M_ext[,ncol(N)+nrow(N)+2]
    obj$pathways$exp_turnover <- M_ext[,ncol(N)+nrow(N)+3]
    obj$pathways$exp_interact <- M_ext[,ncol(N)+nrow(N)+4]
    obj$pathways$active_f <- as.logical(M_ext[,ncol(N)+nrow(N)+5])

    # state specific 
    obj$state$sp_done <- c(obj$state$sp_done, i)
    obj$state$sp_planned <- obj$state$sp_planned[obj$state$sp_planned != i]

    # recalculate state parameters
    active <- obj$pathways$active_f
    obj$state$rates_active <- as.vector(t(obj$pathways$M[active,nrow(N)+1:ncol(N)]) %*% obj$pathways$coefficients[active])

    obj$state$turnover_active <- pa_calculate_turnover( N, obj$state$rates_active)
    obj$state$interaction_active <- pa_calculate_interaction( N, obj$pathways$M[active,nrow(N)+1:ncol(N)], 
                                                              obj$pathways$coefficients[active] )

    return(obj)
}


#
#
# parameters:
# <net>              - Reaction network
# <rates>            - reaction network's steady state vector
# <fexp>             -
# <pmin>             -
# <do_decomposition> - Can be used to switch off intermediate decomposition to
#                      elementary pathways.
# <decreasing>       - Order intermediate species with decreasing turnover?

pa_analysis <- function(net, rates, fexp=1e-2, pmin=1e-3, do_decomposition=T, decreasing=T) {
    # Initialize - Don't use explained turnover for cutoff (not implemented as of April 16)
    xx <- pa_initialize(net, rates, pmin, fexp, fexp, prep=F, decompose=do_decomposition, decreasing=decreasing)

    # Do step till nothing more to do
    while(!pa_is_done(xx))
        xx <- pa_step(xx)

    # Return list with active pathways as first element and coefficients as second 
    return(list(xx$pathways$M[xx$pathways$active_f,-(1:nrow(net[[1]]))], xx$pathways$coefficients[xx$pathways$active_f], xx))
} 

