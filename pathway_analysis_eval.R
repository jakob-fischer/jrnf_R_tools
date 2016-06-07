# author: jakob fischer (jakob@automorph.info)
# description: 

sourced_pathway_analysis_eval <- T

if(!exists("sourced_tools"))
    source("tools.R") 

if(!exists("sourced_cycles"))
    source("cycles.R")

if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")


# Function Transforms a number of pathways in a reaction network into a latex table 
# that then is written to file. 
#
# net      - network that is transformed to table 
# pw       - matrix describing the pathways (1 pathway = 1 row)           
#          - 
# filename - latex table is written to this file
# style=1  - standard style (others not available yet)
# add_info - 1 (numeric) : means that counting is added as additional information
#          - c()         : nothing is added
#          -             : Data frame which has the same number of rows as pathway
#                        : matrix has rows

pa_pathways_to_ltable <- function(filename, pw, net, add_info=1, style=1, longtable=F) {
    con <- file(filename, "w")
    N <- jrnf_calculate_stoich_mat(net)

    # Standard argument for add_info (==1) is interpreted as having to replace it
    # by a data frame that adds a number for each pathway as single string. 
    if(is.numeric(add_info))
        add_info <- pa_h_add_count_column(nrow(pw))

    # fill up spaces (adds spaces to string to make the table code look smoother)
    fus <- function(x, s=25) {
        if(nchar(x) >= s)
            return(x)
        else 
            return(paste(c(x, rep(" ", s-nchar(x))), collapse=""))
    }

    # insert function writes string to file...
    i <- function(...)  {  writeLines(paste(list(...), collapse=""), con)  }

    # Returns string for the main layout of table 
    # <multiplicity> & <educts> & => & products ...
    get_layout_head <- function() {  
        b <- "l r c l"

        if(is.data.frame(add_info))
            b <- paste(c(b, rep(" c", ncol(add_info)), " "), collapse="")
            
        return(b)  
    }
 
    # draws (write to file) the header
    draw_header <- function()  {
        if(is.data.frame(add_info) && 
           (ncol(add_info) != 1 || names(add_info)[[1]] != "EMPTY")) {
            x <- " & & & "
            for(k in 1:ncol(add_info)) {
                n <- names(add_info)[[k]]
                if(n == "EMPTY") 
                    n <- ""
                x <- paste(x, "& ", n, " ", sep="")
            }

            i(paste(x, "\\\\", sep=""))
            i("\\hline")
        }
    }

    # builds (return as string) the specific reaction with specific multiplicity
    # (if multiplicity is negative it means direction is reversed)
    build_rea <- function(rea, mul) {
        educts_s <- jrnf_educts_to_string(net, rea, jrnf_species_name_to_latex_rep, "\\,+\\,")
        products_s <- jrnf_products_to_string(net, rea, jrnf_species_name_to_latex_rep, "\\,+\\,") 
        x <- ""

        if(mul > 1)
            x <- paste(as.character(mul), "X & $", fus(educts_s), "$ & $\\rightarrow $ &  $",
                       fus(products_s), " $ ", sep="")        
        else if (mul == 1) 
            x <- paste("  & $", fus(educts_s), "$ & $\\rightarrow $ &  $",
                       fus(products_s), " $ ", sep="")        
        else if(mul < -1)
            x <- paste(as.character(-mul), "X & $", fus(products_s), "$ & $\\rightarrow $ &  $",
                       fus(educts_s), " $ ", sep="")        
        else if (mul == -1) 
            x <- paste("  & $", fus(products_s), "$ & $\\rightarrow $ &  $",
                       fus(educts_s), " $ ", sep="")        
        return(x)
    }

    # builds (return as string) the additional columns / information for pathway
    # 'l'. If first is numberic, the data (multirow) is written, else empty columns. 
    build_info <- function(l, first=F) {
        x <- ""

        if(is.data.frame(add_info)) {
            for(k in 1:ncol(add_info)) {
                if(is.numeric(first))
                    x <- paste(x, "& ", fus(paste("\\multirow{", as.character(first) ,
                                                  "}{*}{", as.character(add_info[l,k]),
                                                  "}", sep="")), 
                               " ", sep="")
                else
                    x <- paste(x, "& ", fus(" "), " ", sep="")
            }
        }
 
        return(x)
    }

    #
    draw_pw <- function(j) {
        # first list all the reactions in pathway
        x <- which(pw[j,] != 0)
        first <- length(x)+1
        accum_s <- N %*% pw[j,]

        for(k in x) {
            i(paste(build_rea(k, pw[j,k]), build_info(j, first), "\\\\", sep=""))
            first <- F
        }

        # now draw effective reaction
        lhs <- jrnf_side_to_string(which(accum_s < 0),
                                   abs(accum_s[accum_s < 0]),
                                   net, jrnf_species_name_to_latex_rep, "\\,+\\,")
        
        rhs <- jrnf_side_to_string(which(accum_s > 0),
                                   accum_s[accum_s > 0],
                                   net, jrnf_species_name_to_latex_rep, "\\,+\\,")
        
        i("\\cline{1-4}")
        i(paste("net: & $", lhs, " $ & $\\rightarrow $ &  $", rhs, " $ ", build_info(j, first), "\\\\", sep="")) 
        i("\\hline")
    }


    con <- file(filename, "w")
    i("% created by pa_pathways_to_ltable!")
    if(!longtable) {
        i("% please don't forget to include \\usepackage{multirow}.")
        i("\\begin{table}[h]")
        i("%\\caption[table title]{long table caption}")
        i("%\\label{tab:mytabletable}")
        i("\\begin{tabular}{ ", get_layout_head(), " }")
        i("\\hline")
        draw_header()
        for(k in 1:nrow(pw)) 
            draw_pw(k)
        i("\\end{tabular}")
        i("\\end{table}")
    } else {
        i("% (longtable version - make sure longtable package is active)")
        i("% please don't forget to include \\usepackage{multirow}.")
        i("\\begin{center}")
        i("\\begin{longtable}{ ", get_layout_head(), " }")
        i("\\caption[table title]{long table caption} \\label{tab:mylabel} \\\\")
        i("\\hline")
        draw_header()
        i("\\hline")
        i("\\endfirsthead")
        i("\\caption[]{(continued)}\\\\")
        i("\\hline")
        draw_header()
        i("\\hline")
        i("\\endhead")
        i("\\hline")
        i(paste("\\multicolumn{", as.character(3+ncol(add_info)) ,
                "}{c}{(continued on next page)}"))
        i("\\endfoot")
        i("\\hline")
        i("\\endlastfoot")
        for(k in 1:nrow(pw)) 
            draw_pw(k)
        i("\\end{longtable}")
        i("\\end{center}")
    }
    close(con)
}


# This function writes a list of <N> most relevant elementary modes /
# pathways to the file with filename <filename>.
#
# TODO: cleanup? Right place for this function?   

pa_write_em_set <- function(filename, net, em, exp_f, rates, N=100) {
    # Don't have to do anything if no em's present
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
        x <- c(x, capture.output(pa_print_em(net, em[i,])))
        x <- c(x, capture.output(cat("\n")))
    }

    writeLines(x, filename)
}


# TODO 
# Function decomposes all pathways to be elementary and distributes the rates 
# accordingly. For this the same method as in pa_analysis is used 
# (first decompose the pathway, order all elementary pathways with increasing 
#  complexity and decreasing current coefficients, and then develop the 
#  (non-elementary) pathway with the elementary pathways.)
# 
# Maybe this function can late even be called as a subfunction through pa_analysis. 
# (For this the possibility to specify which pathways to check should be given, 
#  also a flag to exclude species from the subpath composition would be necessary.)

pa_make_elementary <- function(net, ems, coef) {
    # cut these reactions of elementary modes that are not in reaction network
    ems <- ems[,1:nrow(net[[2]])]

    # new matrix / ems that have been checked are put here
    ems_new <- ems[c(),]
    coef_new <- c()
    N <- jrnf_calculate_stoich_mat(net)     # stoichiometric matrix

    # The information on which of the original pathways are elementary is written in this vector:
    is_el <- rep(F, nrow(ems))

    for(i in 1:nrow(ems)) {
        cat(".")

        subpath_dec <- pa_subpath_decomposition(N, ems[i,])

        if(nrow(subpath_dec) == 0) 
           cat("ERROR: assortion violated (G)!\n")
 
        p_sort <- pa_check_pathway_present(ems_new, coef_new, subpath_dec)
        o <- order(p_sort$complexity, -p_sort$rate)  # low complexity is more important than high rate
        p_sort <- p_sort[o,]
        subpath_dec <- matrix(subpath_dec[o,], ncol=ncol(subpath_dec))

        s <- ems[i,]*coefficients[i]

        if(nrow(p_sort) == 1) {
            is_el[i] = T


        } 
                        
        for(k in 1:nrow(p_sort)) {
            y <- s/path_M_dec[k,]
        
            if(length(which(y>0)) > 0) {
                rt <- max(min(y[path_M_dec[k,]>0]), 0)
                s <- s - path_M_dec[k,]*rt

                if(is.na(rt))  
                    cat("ERROR: assortion violated (E)!\n")
       
                if(p_sort$present[k]) {
                    path_rates_new[p_sort$id[k]] <- path_rates_new[p_sort$id[k]] + rt
                } else {
                    path_rates_new[length(path_rates_new)+1] <- rt
                    path_M_new <- rbind(path_M_new, matrix(path_M_dec[k,], 1, length(path_M_dec[k,])))
                }
            } 
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
            C_s[i,k] <- get_n_cycles_directed_V(g_tmp,k)[[1]]
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
    
    em_matrix_red <- em_matrix[em_rates != 0, ]

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
    v_ <- rep(0, length(v))                # exansion for ems 1 to <n>

    coeff <- em_rates                      # coefficient of the em
    exp_f <- rep(0, length(em_rates))      # fraction of all rates explained by current em            
    exp_f_acc <- rep(0, length(em_rates))  # fractions of all rates explained by ems 1 to <n>
    err_f_50p <- rep(1, length(em_rates))  # fraction of reaction which is explained by less than 50% by ems 1 to <n>
    err_rmax_50p <- rep(0, length(em_rates)) # highest rate of reaction which is explained by less than 50% by ems 1 to <n> (argmax)
    err_rates <- list()                    # list of data frames. Each data frame contains absolute and relative error by the expansion to this point for each reaction 

    # If <em_df> is given the information on the used elementary modes is accumulated with the expansion of the steady state
    C_sum_acc <- rep(0, length(em_rates))                     # Number of cycles of different length
    C_s_sum_acc <- rep(0, length(em_rates))                   # Counting cycles by subgraph isomorphism (nodes!)
    Deg_acc <- rep(0, length(em_rates))                       # mean degree of all associated species
    Deg_int_acc <- rep(0, length(em_rates))                   # mean degree of all associated species ignoring other reactions
    Deg_max_acc <- rep(0, length(em_rates))                   # max degree of all associated species
    Deg_max_int_acc <- rep(0, length(em_rates))               # max degree of all associated species ignoring other reactions
    Sp_no_acc <- rep(0, length(em_rates))                     # number of species taking part in elementary mode
    Re_acc <- rep(0, length(em_rates))            # number of reactions taking part in elementary mode
    Re_s_acc <- rep(0, length(em_rates))          # number of reactions (counting each only once)
    Ex_acc <- rep(0, length(em_rates))            # Number of exchanges with environment
    Ex_s_acc <- rep(0, length(em_rates))          # Number of exchanges (counting each species only once)
    In_acc <- rep(0, length(em_rates))            # Number of input from environment
    In_s_acc <- rep(0, length(em_rates))          # Number of input (counting each species only once)
    Out_acc <- rep(0, length(em_rates))           # Number of output to environment
    Out_s_acc <- rep(0, length(em_rates))         # Number of output to environment (each species only once) 

    #
    for(i in 1:length(em_rates)) {
        cat(".")
        dv <- em_rates[i]*em_matrix[i,]
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
        #
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



#
# Build matrices that contain the information (for each pathway) which species
# are 
#
#

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
    return(list(x1,x2,(x2-abs(x1))))
}

 
pa_build_elementary_components_incidence <- function(em_matrix, net_N, comp_mat=c(), si=c()) {
    if(is.matrix(net_N) && is.null(comp_mat)) {
        cat("Error in pa_build_elementary_components_incidence_1:")
        cat("Have to have either network object or composition matrix!")
        return(FALSE)
    }

    N <- jrnf_calculate_stoich_mat(net_N)

    # if no composition matrix is given ("comp_mat") it is calculated from
    if(is.null(comp_mat))
        comp_mat <- build_composition_matrix(net_N)

    # if species incidence is not given it is calculated 
    if(is.null(si))
        si <- pa_build_species_incidence(em_matrix, N);  

    # 
    #
    x1 <- t(t(comp_mat)%*%t(abs(si[[1]])))
    x2 <- t(t(comp_mat)%*%t(si[[2]]))

    return(list(x1,x2,(x2-x1)))    
}



pa_cycle_incidence_sp <- function(em_matrix, net, c_max=5) {
    # generating empty vectors
    incidence_cycles_sp <- list()
    for(i in 1:c_max)
        incidence_cycles_sp[[i]] <- matrix(0, ncol=nrow(net[[1]]), nrow=nrow(em_matrix))


    # iterating all pathways
    for(i in 1:nrow(em_matrix)) {
        cat(".")

        rev <- c()                   # is the reaction reverted in the specific elementary mode
        em_abs <- abs(em_matrix[i,]) 
        sel <- c()                   # which reactions are in the elementary mode (id's of ones with 
                                     # higher coefficients are put there multiple times)

        for(j in which(em_abs != 0)) {
            for(k in 1:em_abs[j]) {
                sel <- c(sel, j)
                rev <- c(rev, em_matrix[i,j] < 0)
            }           
        }
            
        # build temporary networks
        net_tmp <- list(net[[1]], net[[2]][sel,])
        N_tmp <- jrnf_calculate_stoich_mat(net_tmp)  
        N_tmp_in <- jrnf_calculate_stoich_mat_in(net_tmp)  
        N_tmp_out <- jrnf_calculate_stoich_mat_out(net_tmp)  

        # find cycles
        g_tmp <- jrnf_to_directed_network_d(net_tmp, rev)
        for(k in 1:c_max) {
            x <- get_n_cycles_directed_B_fast(g_tmp,k,F)
            incidence_cycles_sp[[k]][i,] <- x[[2]] 
        }
    }    

    # assemble data frame and return it
    return(incidence_cycles_sp)
}


pa_cycle_incidence_el <- function(em_matrix, net, c_max=5) {
    # generating empty vectors
    incidence_cycles_el <- list()
    for(i in 1:c_max)
        incidence_cycles_el[[i]] <- matrix(0, ncol=length(get_element_names()), nrow=nrow(em_matrix))

    cm <- build_composition_matrix(net)

    # iterating all pathways
    for(i in 1:nrow(em_matrix)) {
        cat(".")

        rev <- c()                   # is the reaction reverted in the specific elementary mode
        em_abs <- abs(em_matrix[i,]) 
        sel <- c()                   # which reactions are in the elementary mode (id's of ones with 
                                     # higher coefficients are put there multiple times)

        for(j in which(em_abs != 0)) {
            for(k in 1:em_abs[j]) {
                sel <- c(sel, j)
                rev <- c(rev, em_matrix[i,j] < 0)
            }           
        }
            
        # build temporary networks
        net_tmp <- list(net[[1]], net[[2]][sel,])
        N_tmp <- jrnf_calculate_stoich_mat(net_tmp)  
        N_tmp_in <- jrnf_calculate_stoich_mat_in(net_tmp)  
        N_tmp_out <- jrnf_calculate_stoich_mat_out(net_tmp)  

        # find cycles
        g_tmp <- jrnf_to_directed_network_d(net_tmp, rev)
        for(k in 1:c_max) {
            x <- get_n_cycles_directed_B_fast(g_tmp,k,T)
            cycles <- x[[3]]
    
            if(nrow(cycles) != 0)
                for(j in 1:nrow(cycles)) {
                    u <- which(cycles[j,] != 0)
                    if(length(u) != 0) {
                        t <- as.numeric(apply(matrix(cm[u,] != 0,nrow=k), 2, all))
                        incidence_cycles_el[[k]][i,] <- incidence_cycles_el[[k]][i,] + t
                    }
                 }

             
        }
    }    

    # assemble data frame and return it
    return(incidence_cycles_el)
}





# Try to build an heuristic. The main assumption is that all reactions have same thermodynamic
# disequilibrium. First all elementary modes without hv-inflow are removed and their thermodynamic
# efficiency is set to zero. Also those that have only hv-inflow (and no other exchange with boundary)
# are removed.
#

pa_calc_hv_efficiency <- function(em_matrix, em_rates, net) {
    # First identify all exchange (pseudo) reaction
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)
    N <- N_out - N_in
    hv_id <- which(net[[1]]$name == "hv")

    pseudo_rea <- (1 == apply(N != 0, 2, sum) & 1 >= apply(N_in != 0, 2, sum) & 1 >= apply(N_out != 0, 2, sum))
    hv_rea <- (N[hv_id,] != 0) & pseudo_rea

    # 
    hv_count <- em_matrix[,hv_rea]
    ex_count <- apply(em_matrix[,pseudo_rea], 1, sum) - hv_count
    np_count <- apply(em_matrix[,!pseudo_rea], 1, sum)

    efficiency <- ex_count/(ex_count+np_count)
    eff_zero <- which(hv_count != 0 & ex_count == 0 | hv_count == 0 & ex_count != 0)   # efficiency zero ems
    efficiency[eff_zero] <- 0

    return(efficiency)
}


# 
pa_calc_species_cycling <- function(em_matrix, em_rates, net, ci_sp_l) {
    accum <- rep(0, nrow(net[[1]]))
    net <- (jrnf_calculate_stoich_mat(net))
    net_a <- abs(jrnf_calculate_stoich_mat(net))

    # first calculate one matrix containing information whether the species are part 
    # of a cycle in all the path
    c_f <- ci_sp_l[[2]] != 0
    ci_sp_l <- ci_sp_l[-(1:2)]

    for(x in ci_sp_l) {
        c_f <- c_f | x != 0
    }

    # now iterate all the pathways
    for(i in 1:nrow(em_matrix)) {
        a <- as.vector(net %*% em_matrix[i,])
        b <- as.vector(net_a %*% em_matrix[i,])
        sel <- b != 0 & a == 0 & c_f[i,]
 
        accum <- accum + 0.5*b*sel*em_rates[i]
    }

    return(accum)
}

