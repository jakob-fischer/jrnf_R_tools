# author: jakob fischer (jakob@automorph.info)
# description: 
# This file contains code that allows deeper analysis after reaction networks have 
# been calculated and a decomposition (coefficients) of the steady state has been
# found.

sourced_pathway_analysis_eval <- T

if(!exists("sourced_tools"))
    source("tools.R") 

if(!exists("sourced_cycles"))
    source("cycles.R")

if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")


# Function derives properties of all elementary modes in a matrix and returns
# them in a data frame. The function returns them 
#
# parameters:
# <em_matrix>  - Elementary modes
# <net>        - Network object
# <c_max>      - Maximum size of cycles that are calculated.
# <ignore_hv>  - Don't count "hv" species as exchange species (Ex, In, Out)
# 
# TODO At the moment his function assumes to get a network without pseudo-exchange
#      species and a set of elementary modes that fit the number of reaction. For
#      future releases it should be extended to cut the network to pathway size and
#      also drop all exchange pseudoreactions automaticly.      

pa_em_derive <- function(em_matrix, net, c_max=4, ignore_hv=F) {
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
    Hv <- rep(0, nrow(em_matrix))            # Number of hv species consumed (generated?) by this pathway
    Hv_s <- rep(0, nrow(em_matrix))          # Hv input present? (1 iff, 0 else)

    # these species are not counted for Ex, In, Out!
    ex_exclude_f <- c()
    ex_hv <- c(which(net[[1]]$name == "hv'" | net[[1]]$name == "hv"))

    if(ignore_hv) 
        ex_exclude_f <- ex_hv

    # calculating degree of <net> independently of em's
    net_deg <- as.vector(degree(jrnf_to_undirected_network(net)))
    N <- jrnf_calculate_stoich_mat(net)

    # iterating all pathways
    for(i in 1:nrow(em_matrix)) {
        cat(".")

        rev <- em_matrix[i,] < 0       # is the reaction reverted in the specific elementary mode
        em_abs <- abs(em_matrix[i,]) 
        sel <- which(em_abs != 0)       # which reactions are in the elementary mode (id's of ones with 
                                        # higher coefficients are not considered

        # build temporary networks
        net_tmp <- list(net[[1]], net[[2]][sel,])
        N_tmp <- jrnf_calculate_stoich_mat(net_tmp)  
        N_tmp_in <- jrnf_calculate_stoich_mat_in(net_tmp)  
        N_tmp_out <- jrnf_calculate_stoich_mat_out(net_tmp)  

        # generate stoichometric matrix with all 1 species reactions (inflow / outflow) removed
        N_int <- N
        for(j in which(apply(N_int != 0, 2, sum) == 1)) 
            N_int[,j] <- 0
  
        em_bilance <- N_int %*% em_matrix[i,]   
        Hv[i] = sum(abs(em_bilance[ex_hv]))
        em_bilance[ex_exclude_f] <- 0

        # find cycles
        g_tmp <- jrnf_to_directed_network_d(net_tmp, rev[sel])
        for(k in 1:c_max) {
            C[i,k] <- get_n_cycles_directed(g_tmp,k)[[1]]
            C_s[i,k] <- get_n_cycles_directed_V(g_tmp,k)[[1]]
        }
            
        # Calculate various other complexity measures. "species_con" are those 
        # species that are in- or output species of at least one reaction in the
        # pathway.
        species_con <- which(apply(N_tmp_in != 0, 1, any) | apply(N_tmp_out != 0, 1, any))
        degree_global <- net_deg[ species_con ]
        degree_internal <- as.vector(degree( jrnf_to_undirected_network(net_tmp) ))[species_con]

        Deg[i] <- mean( degree_global )
        Deg_int[i] <- mean( degree_internal )
        Deg_max[i] <- max( degree_global )
        Deg_max_int[i] <- max( degree_internal )
        Sp_no[i] <- length( species_con )

        Re[i] <- sum(em_abs) 
        Re_s[i] <- sum(em_abs != 0)
        Ex[i] <- sum(abs(em_bilance))
        Ex_s[i] <- sum(em_bilance != 0)
        In[i] <- sum(abs(pmin(em_bilance,0)))
        In_s[i] <- sum(abs(pmin(em_bilance,0)) != 0)
        Out[i] <- sum(abs(pmax(em_bilance,0)))
        Out_s[i] <- sum(abs(pmax(em_bilance,0)) != 0)
    }    

    # sum of cycles is calculated
    C_sum <- apply(C, 1, sum)
    C_s_sum <- apply(C_s, 1, sum)
    Hv_s = sign(Hv)


    # assemble data frame and return it
    return(data.frame(C_sum=C_sum, C_s_sum=C_s_sum, C_1=C[,1], C_2=C[,2], C_3=C[,3], C_1_s=C_s[,1], 
                      C_2_s=C_s[,2], C_3_s=C_s[,3], Deg=Deg, Deg_int=Deg_int, Deg_max=Deg_max, 
                      Deg_max_int=Deg_max_int, Sp_no=Sp_no, Re=Re, Re_s=Re_s, Ex=Ex, Ex_s=Ex_s,
                      In=In, In_s=In_s, Out=Out, Out_s=Out_s, Hv=Hv, Hv_s=Hv_s))
}


# Uses an existing expansion of a network <net>'s steady state flow <v> and 
# calculates how good it is absolute (using all elementary modes) and cummulative
# (using elementary modes 1 to <n>). The expansion / elementary modes have to be
# ordered in advance. If coefficients of reactions inside em's are negative tese 
# reactions of the network are reversed before starting evaluating the expansion.
# For this it is checked if there is no sign missmatch for all the pathways which 
# have <em_coef> > 0.
#
# parameters:
# <em_matrix>  - Elementary modes
# <em_coef>    - Coefficients of the elementary modes
# <net>        - Reaction network
# <v>          - Steady state rate vector
# <em_df>      - Data frame with extended pathway information (from "pa_em_derive")

pa_ana_expansion <- function(em_matrix, em_coef, net, v, em_df=c()) {
    # First check the sign condition if rates of elementary modes are nonzero all coefficients of 
    # the respective elementary modes have to have the same sign as the reactions rate (v)
    
    em_matrix_red <- em_matrix[em_coef != 0, ]

    ma <- apply(em_matrix_red, 2, max)
    mi <- apply(em_matrix_red, 2, min)

    if(!all(-1 != sign(ma)*sign(mi))) {
        cat("Reaction's signs in elementary modes did not match!\n") 
        return()
    }

    if(!all(-1 != sign(ma)*sign(v))) {
        cat("Reaction's signs in elementary modes did not match with reaction rates'!\n") 
        return()
    }


    # If sign conditions are fulfilled reactions, reaction rates and coefficients of 
    # elementary modes are reversed for negative rate reactions.
    for(i in 1:ncol(em_matrix)) 
        if(v[i] < 0) 
            em_matrix[,i] <- -em_matrix[,i]
    
    v <- abs(v)


    # strictly speaking only non pseudo reactions are part of the network!
    N <- jrnf_calculate_stoich_mat(net)
    rea_id_np <- (1 != apply(N != 0, 2, sum))

    v_sum <- sum(v[rea_id_np]) 
    v_ <- rep(0, length(v))                # expansion for ems 1 to <n>

    coeff <- em_coef                      # coefficient of the em
    exp_f <- rep(0, length(em_coef))      # fraction of all rates explained by current em            
    exp_f_acc <- rep(0, length(em_coef))  # fractions of all rates explained by ems 1 to <n>
    err_f_50p <- rep(1, length(em_coef))  # fraction of reaction which is explained by less than 50% by ems 1 to <n>
    err_rmax_50p <- rep(0, length(em_coef)) # highest rate of reaction which is explained by less than 50% by ems 1 to <n> (argmax)
    err_rates <- list()                    # list of data frames. Each data frame contains absolute and relative error by the expansion to this point for each reaction 

    # If <em_df> is given the information on the used elementary modes is accumulated with the expansion of the steady state
    C_sum_acc <- rep(0, length(em_coef))                     # Number of cycles of different length
    C_s_sum_acc <- rep(0, length(em_coef))                   # Counting cycles by subgraph isomorphism (nodes!)
    Deg_acc <- rep(0, length(em_coef))                       # mean degree of all associated species
    Deg_int_acc <- rep(0, length(em_coef))                   # mean degree of all associated species ignoring other reactions
    Deg_max_acc <- rep(0, length(em_coef))                   # max degree of all associated species
    Deg_max_int_acc <- rep(0, length(em_coef))               # max degree of all associated species ignoring other reactions
    Sp_no_acc <- rep(0, length(em_coef))                     # number of species taking part in elementary mode
    Re_acc <- rep(0, length(em_coef))            # number of reactions taking part in elementary mode
    Re_s_acc <- rep(0, length(em_coef))          # number of reactions (counting each only once)
    Ex_acc <- rep(0, length(em_coef))            # Number of exchanges with environment
    Ex_s_acc <- rep(0, length(em_coef))          # Number of exchanges (counting each species only once)
    In_acc <- rep(0, length(em_coef))            # Number of input from environment
    In_s_acc <- rep(0, length(em_coef))          # Number of input (counting each species only once)
    Out_acc <- rep(0, length(em_coef))           # Number of output to environment
    Out_s_acc <- rep(0, length(em_coef))         # Number of output to environment (each species only once) 

    # Accumulate data for all pathways in the expansion
    for(i in 1:length(em_coef)) {
        cat(".")
        dv <- em_coef[i]*em_matrix[i,]
        v_ <- v_ + dv
        exp_f[i] <- sum(dv[rea_id_np])/v_sum      
        exp_f_acc[i] <- sum(v_[rea_id_np])/v_sum
        # total error 
        err_abs <- v_ - v
        err_rel <- abs((v_-v)/v)
        err_rel[v == 0] <- 0

        err_f_50p[i] <- length(which(err_rel > 0.5)) / length(err_rel)
        err_rmax_50p[i] <- max(c(v[err_rel > 0.5],0))

        err_rates[[i]] <- data.frame(v=v, err_abs=err_abs, err_rel=err_rel, v_acc=v_, dv=dv) 
 
        # if em derived data (em_df) is available their accumulated quantities are derived
        C_sum_acc[i] <- sum(em_df$C_sum[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        C_s_sum_acc[i] <- sum(em_df$C_s_sum[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        Deg_acc[i] <- sum(em_df$Deg[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        Deg_int_acc[i] <- sum(em_df$Deg_int[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        Deg_max_acc[i] <- sum(em_df$Deg_max[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        Deg_max_int_acc[i] <- sum(em_df$Deg_max_int[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        Sp_no_acc[i] <- sum(em_df$Sp_no[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        Re_acc[i] <- sum(em_df$Re[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        Re_s_acc[i] <- sum(em_df$Re_s[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        Ex_acc[i] <- sum(em_df$Ex[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        Ex_s_acc[i] <- sum(em_df$Ex_s[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        In_acc[i] <- sum(em_df$In[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        In_s_acc[i] <- sum(em_df$In_s[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        Out_acc[i] <- sum(em_df$Out[1:i]*exp_f[1:i])/sum(exp_f[1:i])
        Out_s_acc[i] <- sum(em_df$Out_s[1:i]*exp_f[1:i])/sum(exp_f[1:i])

    }

    if(length(em_df) == 0)
        return(data.frame(exp_f=exp_f, exp_f_acc=exp_f_acc,                                                   
                          err_f_50p=err_f_50p, err_rmax_50p=err_rmax_50p, 
                          err_rates=I(err_rates)))
    else
        return(data.frame(exp_f=exp_f, exp_f_acc=exp_f_acc,                                                   
                          err_f_50p=err_f_50p, err_rmax_50p=err_rmax_50p, C_sum_acc=C_sum_acc, C_s_sum_acc=C_s_sum_acc,
                          Deg_acc=Deg_acc, Deg_int_acc=Deg_int_acc, Deg_max_acc=Deg_max_acc, 
                          Deg_max_int_acc=Deg_max_int_acc, Sp_no_acc=Sp_no_acc, Re_acc=Re_acc,
                          Re_s_acc=Re_s_acc, Ex_acc=Ex_acc, Ex_s_acc=Ex_s_acc, In_acc=In_acc,
                          In_s_acc=In_s_acc, Out_acc=Out_acc, Out_s_acc=Out_s_acc, 
                          err_rates=I(err_rates)))
}


# Build matrices that contain the information (for each pathway) how each species
# is occuring in it / influenced by it.

pa_build_species_incidence <- function(em_matrix, net_N) {
    # calculate stoichiometric matrix (if necessary)
    N <- jrnf_calculate_stoich_mat(net_N)

    # change of species by each elementary mode is one column in matrix x1
    x1 <- t(N%*%t(em_matrix)) 
    # the matrix x2 contains the cummulative absolute change (if one unit of
    # a chemical species is created by one reaction and removed by another 
    x2 <- t(abs(N)%*%t(em_matrix))

    # return in a list
    # first: relative exchange with environment (with sign)
    # second: absolute change of all reactions added (non-negative)
    # third: absolute change that is not explained by exchange with environment
    return(list(change=x1,prod_cons=x2,internal=(x2-abs(x1))))
}


# Build matrices that contain the information (for each pathway) how each 
# elementary component is occuring in it / influenced by it.
#
# parameters:
# <em_matrix> - Reaction pathways / elementary modes
# <net_N>     - Reaction network / stoichiometric matrix
# <comp_mat>  - Composition matrix
# <si>        - Species incidence (pa_build_species_incidence) for same 
#               network and same pathways
 
pa_build_elements_incidence <- function(em_matrix, net_N, comp_mat=c(), si=c()) {
    if(is.matrix(net_N) && is.null(comp_mat)) {
        cat("Error in pa_build_elements_incidence: ")
        cat("Have to have either network object or composition matrix!")
        return(FALSE)
    }

    N <- jrnf_calculate_stoich_mat(net_N)

    # If no composition matrix is given ("comp_mat") it is calculated from the
    # network object.
    if(is.null(comp_mat))
        comp_mat <- build_composition_matrix(net_N)

    # If species incidence is not given it is calculated 
    if(is.null(si))
        si <- pa_build_species_incidence(em_matrix, N);  

    x1 <- t(t(comp_mat)%*%t(abs(si[[1]])))
    x2 <- t(t(comp_mat)%*%t(si[[2]]))

    return(list(change=x1,          # Number of parts of the elementary component
                                    # that is contained in all species that are 
                                    # effectively produced or consumed 
                prod_cons=x2,       # Accumulative turnover of all elementary
                                    # components
                internal=(x2-x1)))  # The part of turnover that is not explained  
}                                   # by in- or outflow


# Calculate of how many cycles of a given size each species is part of in each 
# pathway. For calculating cycles each reaction is only taken once (even if it 
# has a bigger coefficient in the pathway).
#
# parameters:
# <em_matrix> - Elementary modes / reaction pathways
# <net>       - Reaction network
# <c_max>     - Maximal cycle size

pa_cycle_incidence_sp <- function(em_matrix, net, c_max=5) {
    # generate empty matrices
    incidence_cycles_sp <- list()
    for(i in 1:c_max)
        incidence_cycles_sp[[i]] <- matrix(0, ncol=nrow(net[[1]]), nrow=nrow(em_matrix))


    # iterate all pathways
    for(i in 1:nrow(em_matrix)) {
        cat(".")

        rev <- em_matrix[i,] < 0     # is the reaction reverted in the specific elementary mode
        sel <- em_matrix[i,] != 0    # Which reactions are in the elementary mode? 
                                     # Each reaction is just taken once.
            
        # build temporary networks
        net_tmp <- list(net[[1]], net[[2]][sel,])
        N_tmp <- jrnf_calculate_stoich_mat(net_tmp)  
        N_tmp_in <- jrnf_calculate_stoich_mat_in(net_tmp)  
        N_tmp_out <- jrnf_calculate_stoich_mat_out(net_tmp)  

        # find cycles
        g_tmp <- jrnf_to_directed_network_d(net_tmp, rev[sel])
        for(k in 1:c_max) {
            x <- get_n_cycles_directed(g_tmp,k,F)
            # save / collect cycle number (for each species)
            incidence_cycles_sp[[k]][i,] <- x[[2]] 
        }
    }    

    # assemble data frame and return it
    return(incidence_cycles_sp)
}


# Calculate of how many cycles of a given size each elementary component is part 
# of in each pathway. An elementary component being part of a cycle means being
# presen in EVERY species that is forming the cycle! For calculating cycles each 
# reaction is only taken once (even if it has a bigger coefficient in the pathway).
#
# parameters:
# <em_matrix> - Elementary modes / reaction pathways
# <net>       - Reaction network
# <c_max>     - Maximal cycle size

pa_cycle_incidence_el <- function(em_matrix, net, c_max=5) {
    # generating empty vectors
    incidence_cycles_el <- list()
    for(i in 1:c_max)
        incidence_cycles_el[[i]] <- matrix(0, ncol=length(get_element_names()), nrow=nrow(em_matrix))

    comp_m <- build_composition_matrix(net)

    # iterating all pathways
    for(i in 1:nrow(em_matrix)) {
        cat(".")

        rev <- em_matrix[i,] < 0     # Which reactions are reversed?
        sel <- em_matrix[i,] != 0    # Which reactions are in the elementary mode. All

            
        # build temporary networks
        net_tmp <- list(net[[1]], net[[2]][sel,])
        N_tmp <- jrnf_calculate_stoich_mat(net_tmp)  
        N_tmp_in <- jrnf_calculate_stoich_mat_in(net_tmp)  
        N_tmp_out <- jrnf_calculate_stoich_mat_out(net_tmp)  

        # find cycles
        g_tmp <- jrnf_to_directed_network_d(net_tmp, rev[sel])
        for(k in 1:c_max) {
            x <- get_n_cycles_directed(g_tmp,k,T)
            cycles <- x[[3]]
    
            if(nrow(cycles) != 0)
                for(j in 1:nrow(cycles)) {
                    u <- which(cycles[j,] != 0)
                    if(length(u) != 0) {
                        # increment for those elementary components that are 
                        # part of ALL species of the cycle
                        t <- as.numeric(apply(matrix(comp_m[u,] != 0,nrow=k), 2, all))
                        incidence_cycles_el[[k]][i,] <- incidence_cycles_el[[k]][i,] + t
                    }
                 }  
        }
    }    

    # assemble data frame and return it
    return(incidence_cycles_el)
}


# For a cycle incidence obtained from "pa_cycle_incidence_sp" or 
# "pa_cycle_incidence_el" this function calculates an cycling 
# score for each species / element if given an explained fraction
# (<exp_f>) that weights the reaction pathways.

pa_calculate_sp_el_cycling <- function(ci_sp_el_l, exp_f) {
    n_cycle <- list()
    accum <- rep(0, ncol(ci_sp_el_l[[1]]))

    for(X in ci_sp_el_l) {
        y <- apply(scale_mat_rows(X, exp_f), 2, sum)
        accum <- accum + y
        n_cycle[[length(n_cycle)+1]] <- y
    }

    # return in two elemen list
    return(list(n_cycle=n_cycle, accum=accum))
}


