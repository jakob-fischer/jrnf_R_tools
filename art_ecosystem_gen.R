# author: jakob fischer (jakob@automorph.info)
# description: 

sourced_art_ecosystem_gen <- T

if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")

# Function creates artificial ecosystem with <N> species... 
#


# Function draws energies for artificial ecosystem
# Activation energies are drawn from "rplancklike" distribution and standard chemical
# potentials are drawn from a gauß distribution. Except the chemical potential for "hv" 
# which is either 50 or the maximum of all other chemical potentials + 5...
# If no <flat_energies> (all zero) are drawn additionally the constraint that 
# photoreactions increase the energy is preserved.

jrnf_ae_draw_energies <- function(net, flat_energies=F, limit_AE=F) {
    lAE <- NA
    if(limit_AE) lAE <- 3
    if(is.numeric(limit_AE)) lAE <- limit_AE

    net[[1]]$energy <- rnorm(nrow(net[[1]]))    
    net[[2]]$activation <- pmin(rplancklike(nrow(net[[2]])), lAE, na.rm=T)

    if(flat_energies) {
        net[[1]]$energy <- rep(0, nrow(net[[1]]))    
        net[[2]]$activation <- rep(1, nrow(net[[2]]))    
    }
  
    if(any(net[[1]]$name == "hv")) {
        if(!flat_energies) {   # If having hv-species and no flat energies have
                               # to maintain property of photochemical reactions
                               # only increasing energy levels
            hv_id <- which(net[[1]]$name == "hv")
            N <- jrnf_calculate_stoich_mat(net)

            hv_lhs <- which(N[hv_id,] < 0)
            hv_rhs <- which(N[hv_id,] > 0)
            N[,hv_id] <- 0

            is_v <- function() {
                energy <- net[[1]]$energy
                energy[hv_id] <- 0
                d_E <- as.numeric(t(N) %*% matrix(energy, ncol=1))
                return(all(d_E[hv_lhs] > 0) && all(d_E[hv_rhs] < 0))
            }
 
            while(!is_v()) 
                net[[1]]$energy <- rnorm(nrow(net[[1]])) 
        }

        net[[1]]$energy[net[[1]]$name == "hv"] <- max(50, max(net[[1]]$energy)+5)
    }

    return(net)
}


# Helper to create names consistent with elementary constituents
hcae_create_name <- function(comp, c_names, ex_names, empty_name="X", prefix="") {
    name <- prefix

    # first build name
    if(sum(comp) == 0) 
        name <- paste(name, empty_name, sep="")
    else 
        for(i in 1:length(comp))
            if(comp[i] != 0) 
                name <- paste(name, c_names[i], comp[i], sep="")    


    # second step - unify name (by appending "_2", "_3", ...
    if(any(name == ex_names)) {    # name already used
        i <- 2
        while(any(paste(name, "_", i, sep="") == ex_names))
            i <- i+1

        name <- paste(name, "_", i, sep="")
    }

    return(name)   # TODO implement ;)
}

# Checks if a reaction is possible from elementary constituents
hcae_check_rea_constituents <- function(rea, comp, N) {
    comp <- comp$composition
    # Returns composition of <i>th participant in reaction
    get_c <- function(i) {
        if(rea[i] == 1)
            return(rep(0, ncol(comp)))
        else
            return(comp[rea[i]-1,])
    }

    # All components have to match (empty species rea[x] == 1) is mapped to
    # zero of all components by the get_c-helper function
    return(all(get_c(1)+get_c(2) == get_c(3)+get_c(4)))
}


# This function draws an elementary composition for <N> species that are build from
# the number of elementary components in comp$components. Parameters is the maximal
# allowed number of duplicates (<max_tot>), the maximum number of duplicates for a single
# species (<max_single>), maximal stoichiometric coefficient (<max_stoich>) and the flag 
# controlling massless species creation (<allow_massless>). The parameter <lambda> 
# corresponds to this parameter in the poisson distribution from which stoichiometric
# coefficents are sampled.

hcae_extend_elem_composition <- function(N, comp, max_tot=c(), max_single=c(), 
                                         lambda=0.1, allow_massless=F, max_stoich=10, sp_n_pre="") {
    comp_no <- ncol(comp$composition)  # number of components derived by existing composition matrix
    N_new <- N + nrow(comp$composition)

    if(comp_no == 1) {                 # if only one component is allowed and max_tot, max_single
        if(is.null(max_tot))           # are unset they are set to total number of species
            max_tot <- N_new

        if(is.null(max_single))
            max_single <- N_new
    } else {
        if(is.null(max_tot))           # else (more than one component) the standard maximum is 
            max_tot <- (N_new)/3       # two duplicates per composition and one third in total

        if(is.null(max_single))
            max_single <- 2
    }

    composition <- matrix(0, nrow=N_new, ncol=comp_no)
    energy <- c(comp$energy, rnorm(N))         
    component_names <- c("C", "N", "O", "H", "P")
    component_names <- component_names[1:comp_no]

    v <- rep(max_stoich+1, comp_no)
    trans <- b_mdim_transf(v) 

    calc_prop <- function(x)  {
        p_x <- function(y)  {  return((lambda ** y)/factorial(y) * exp(-lambda))  }
        return( prod( sapply(trans$from_linear(x), p_x) ) )   
    } 

    prob <- sapply(1:((max_stoich+1)**comp_no), calc_prop)
    if(!allow_massless)
        prob[trans$to_linear(rep(1, comp_no))] <- 0   
    cnt <- rep(0, length(prob))
    cnt_t <- 0

    # Following loop is run for all existing and newly created compositions. 
    # The reason it is rerun for existing ones is that the duplicates have 
    # to be counted in 'cnt' and 'cnt_t'.

    for(i in 1:N_new) { 
        if(i > (N_new-N)) {
            x <- sample(1:length(prob), 1, F, prob/sum(prob))   
            composition[i,] <- trans$from_linear(x)-1
        } else {  
            composition[i,] <- comp$composition[i,]
            x <- trans$to_linear(composition[i,]+1)
        }

        cnt[x] <- cnt[x] + 1
        if(cnt[x] > 1) {
            cnt_t <- cnt_t + 1
             
            if(cnt_t >= max_tot)
                prob[cnt != 0] <- 0
        }
        
        if(cnt[x] >= max_single)
            prob[x] <- 0
    }

    # now (re)calculate names
    name <- comp$name
    for(i in (nrow(comp$composition)+1):N_new)
        name <- c(name, hcae_create_name(composition[i,], component_names,
                                         name, prefix=sp_n_pre))   
 
    return(list(composition=composition, name=name, energy=energy))
}


# wrapper function / drawing an elementary composition is a special case of 
# adding new species to an existing composition. All parameters are explained in
# the corresponding function (above).

hcae_draw_elem_composition <- function(N, comp_no, max_tot=c(), max_single=c(), 
                                         lambda=0.1, allow_massless=F, max_stoich=10) {
    comp <- list(composition=matrix(0, nrow=0, ncol=comp_no), name=c(), energy=c())

    return(hcae_extend_elem_composition(N, comp, max_tot, max_single, lambda, allow_massless, max_stoich))
}


# If reaction is valid (by composition) returns max mass transfer
hcae_calc_rea_transfer <- function(rea, comp, N) {
    comp <- comp$composition
    get_c <- function(i) {
        if(rea[i] == 1)
            return(rep(0, ncol(comp)))
        else
            return(comp[rea[i],])
    }

    # the components that are hv are set to zero / empty species first
    for(i in 1:4)
        if(rea[i] == N+2) 
            rea[i] <- 1
   
    if(rea[1] != 1)
        return(min(sum(abs(get_c(1) - get_c(3))), 
                   sum(abs(get_c(1) - get_c(4)))))

    if(rea[2] != 1)
        return(min(sum(abs(get_c(2) - get_c(3))), 
                   sum(abs(get_c(2) - get_c(4)))))

    return(0)    
}

# Check further conditions on reactions 
hcae_check_rea_conditions <- function(rea, N, comp) {
    hv_id <- which(comp$name == "hv") + 1
    empty_id <- 1

    empty_ed  <- (rea[1] == empty_id) + (rea[2] == empty_id)
    empty_pro  <- (rea[3] == empty_id) + (rea[4] == empty_id)

    # Avoid double reactions (through reverse reaction)
    if(rea[2] > rea[1] || rea[4] > rea[3] || rea[3] > rea[1])
        return(F) 

    hv_ed <- (rea[1] == hv_id) + (rea[2] == hv_id)
    hv_pro <- (rea[3] == hv_id) + (rea[4] == hv_id)

    # One side is completely empty
    if(empty_ed == 2 || empty_pro == 2)
        return(F)

    # Reaction doesn't do anything  
    if(rea[1] == rea[3] && rea[2] == rea[4] || rea[1] == rea[4] && rea[2] == rea[3])
        return(F)


    if(length(hv_id) != 1)  # If there are no photoreactions: can finish here
        return(T)


    # hv is on both sides  OR  one side only has hv
    if(hv_ed+hv_pro > 1 || hv_ed+empty_ed > 1 || hv_pro+empty_pro > 1)
        return(F)


    # Reaction has the form "hv + A -> A"
    # (complicated - permutations because hv_id can be bigger or smaller than id(A))
    if(rea[1] == hv_id && rea[2] == rea[3] && rea[4] == empty_id ||
       rea[2] == hv_id && rea[1] == rea[3] && rea[4] == empty_id ||
       rea[3] == hv_id && rea[4] == rea[1] && rea[2] == empty_id ||
       rea[4] == hv_id && rea[3] == rea[1] && rea[2] == empty_id)
        return(F)

    # Exclude drawing "hv + A -> B" AND "hv + B -> A"
    # 
    # Question: should reactions like "hv + A -> B + C" be 
    # forbidden if mu_A > mu_B + mu_C 
    energy <- c(0, comp$energy)
    energy[hv_id] <- 0
    delta_energy <- sum(energy[c(rea[3], rea[4])], -energy[c(rea[1], rea[2])])

    # iff no energy difference "linear" photochemical reactions are only allowed
    # from lower to higher species id's - photochemical reactions with three 
    # components are always allowed
    if(delta_energy == 0) {
        if(rea[1] == hv_id && rea[4] == empty_id && rea[2] >= rea[3] ||
           rea[2] == hv_id && rea[4] == empty_id && rea[1] >= rea[3] ||
           rea[3] == hv_id && rea[2] == empty_id && rea[4] >= rea[1] ||
           rea[4] == hv_id && rea[2] == empty_id && rea[3] >= rea[1])
            return(F)
    # if energy decreases in forward direction the reaction is discarded if the
    # photon occurs on the LHS (independently of the number of components)
    } else if (delta_energy < 0) { 
        if(rea[1] == hv_id || rea[2] == hv_id)
            return(F)
    # if energy increases in forward direction the reaction is discarded if the
    # photon occurs on the RHS (independently of the number of components)
    } else if (delta_energy > 0) {
        if(rea[3] == hv_id || rea[4] == hv_id)
            return(F)
    }

    return(T)
}



# Helper function that is given a (weighted) adjacency matrix and a number of modules
# (has to be divisor of species number) and then tries to maximize the modularity
# (weight of edges inside the modules) by randomly reordering the species ids.


gmr_get_inner_density <- function(M_adj, N_mod) {
    s <- 0
    N <- ncol(M_adj)    
    mod_size <- N/N_mod

    for(i in 1:N_mod) {
        x <- (i-1)*mod_size+1:mod_size
        s <- s + sum(M_adj[x,x])/(mod_size**2)
    }

    return(s/N_mod)
}


jrnf_get_modular_reordering <- function(M_adj, N_mod) {
    #
    #



    N <- ncol(M_adj)                    # number of species
    o <- 1:N                            # current reordering
    density <- sum(M_adj)/(N**2)        # density



    r1 <- sample(N, N**2*10, replace=T) 
    r2 <- sample(N, N**2*10, replace=T)
   
    for(i in 1:(N**2*10)) {
      o_t <- o
      o_t[r2[i]] <- o[r1[i]]
      o_t[r1[i]] <- o[r2[i]]  

      density_t <- gmr_get_inner_density(M_adj[o_t,o_t], N_mod)
      
      if(density_t > density) {
          density <- density_t
          o <- o_t
          cat("-> ", density, "\n")
      }
    }

    return(o)
}


jrnf_analyze_ecosystem_constituents <- function(names) {
    if(is.list(names))
        names <- names[[1]]$name

    if("hv" == names[length(names)])
        names <- names[-length(names)]

    constituents <- c()

    # First remove tailing "_<number>" and extract components names
    for(i in 1:length(names)) {
        names[i] <- strsplit(names[i], "_", fixed=TRUE)[[1]][1]
        constituents <- c(constituents, strsplit(names[i], "[0-9]+")[[1]])
    }

    constituents <- sort(unique(constituents))

    m <- matrix(0, ncol=length(constituents), nrow=length(names))


    # Now extract component's names again / including multiplicity
    for(i in 1:length(names)) {
        con <- strsplit(names[i], "[0-9]+")[[1]]
        mul <- as.numeric(strsplit(names[i], "[A-Za-z]+")[[1]][-1])

        if(length(con) != length(mul))
            cat("error: length of con and mul have to match (jrnf_analyze_ecosystem_constituents)!\n")

        for(j in 1:length(con)) {
            sel <- which(constituents == con[j])
            m[i,sel] <- mul[j]
        }
    }

    return(list(constituents, m))
} 


# Function creates an artificial ecosystem randomly 
# <N> - number of species
# <M> - total number of reactions
# <no_2fold> - number of linear reactions
# <no_hv> - number of photochemical reactions
# <cat_as_lin> - flag that indicates whether catalytical reactions (A + C -> B + C)
#                are treated identical to linear reactions in terms of duplicates 
#                and in terms of the reactions of the different types that are drawn
# <type_spec_dup> - duplicates are treated specific to the three types and not 
#                   accross them
# <comp> - number of elementary components the set of species has 
#          (elementary composition is generated by hcae_draw_elem_composition
#          using standard parameters.) Alternatively one can give a composition
#          created by hcae_draw_elem_composition directly. Note that this method
#          won't add a hv-pseudospecies then but it has to be included.
# <allow_direct_backflow> - Makes it possible to forbid direct backflow if it would
#                           be otherways possible. Direct backflow means a combination
#                           of reactions like: hv + A -> B;  B -> A
# <rm_dup> - Remove duplicates. Only allow one reaction for each column in the 
#            stoichiometric matrix. Influenced by <type_spec_dup> and <cat_as_lin>
#
# TODO Maybe function can be reengineered #
#      Sort different reactions in multiple equivalency groups (same effective change) 
#      and types (linear, photoreaction, rest). 
#
#      At the moment <N> does not include the hv-pseudospecies. If a composition is given 
#      (<comp>) <N> is set minus one of the number of species in this composition.

jrnf_ae_create <- function(N, M, no_2fold, no_hv, comp=c(), cat_as_lin=F, 
                                             type_spec_dup=T, allow_direct_backflow=T, rm_dup=T) { 
    hv_name <- "hv"

     if(is.numeric(comp)) {
        # draw compostition (el. constituents) and energy of species
        cat("Drawing elementary composition")
        comp <- hcae_draw_elem_composition(N, comp)
        cat(".\n")

        # Add species for photons / energy source
        comp$composition <- rbind(comp$composition, matrix(0, nrow=1, ncol=ncol(comp$composition)))
        comp$energy <- c(comp$energy, max(max(comp$energy)+5, 50))
        comp$name <- c(comp$name, hv_name) 
    } else 
        N <- nrow(comp$composition)-1

    hv_id <- which(comp$name == hv_name)
    trans <- b_mdim_transf(rep(N+2, 4))   

    composition <- comp$composition
    energy <- comp$energy
    name <- comp$name

    
    is_rea_possible <- function(i) {
        x <- trans$from_linear(i)

        if(!hcae_check_rea_constituents(x, comp, N) ||
           !hcae_check_rea_conditions(x, N, comp))
            return(F)
            
        return(T)
    }

    is_hv_rea <- function(i) {
        x <- trans$from_linear(i)
        return(any(x-1 == hv_id)) 
    }

    # Check if the reaction with id <i> is linear (definition depends on whether 
    # catalytic reactions are considered linear)
    is_lin_rea <- function(i) {
        x <- trans$from_linear(i)
        if(!cat_as_lin)   
            return(x[2] == 1 && x[4] == 1) 
        else
            return(x[1] == x[3] || x[1] == x[4] || x[2] == x[3] || x[2] == x[4])
    }

    # Subfunction transforms a linear reaction id <v> into an effective change
    # of concentration. Change is normalized to first nonzero component being
    # greater zero. If direct backflow is forbidden change of hv is set to zero
    # (so backflow is a duplicate of "forward flow") 
    fval_to_stoichcol <- function(v) {
        v <- trans$from_linear(v)
        x <- matrix(0, ncol=N+1, nrow=1)
        if(v[1] != 1) x[v[1]-1] <- x[v[1]-1] - 1
        if(v[2] != 1) x[v[2]-1] <- x[v[2]-1] - 1
        if(v[3] != 1) x[v[3]-1] <- x[v[3]-1] + 1
        if(v[4] != 1) x[v[4]-1] <- x[v[4]-1] + 1
        # unique representation - first nonzero entry has to be positive
        m <- which(x[1,] != 0)[1]

        if(x[1,m] < 0)
            x <- -x       

        if(!allow_direct_backflow)
            x[1,hv_id] <- 0

        return(x)
    }

    # Sample one element from x
    s <- function(x) {
        return(x[sample(length(x))])
    }

    r <- rev

    # Find the duplicates (in terms of effective change of species) reactions
    dup <- function(xx) {
        return(duplicated(t(sapply(xx, fval_to_stoichcol))))
    }


    reactions <- which(sapply(1:((N+2)**4), is_rea_possible))
    cat("found", length(reactions), "valid reactions!\n")
    reactions <- s(reactions)


    if(!type_spec_dup && rm_dup) {
        reactions <- s(reactions)
        x <- dup(reactions)
        reactions <- reactions[!x]
        cat("after removing (effective) doubles, there are ", length(reactions), "valid reactions!\n")
    }

    sel <- which(sapply(reactions, is_hv_rea))
    if(length(sel) == 0) {
        rea_hv <- c()
    } else {
        rea_hv <- reactions[sel]
        reactions <- reactions[-sel]
    }

    sel <- sapply(reactions, is_lin_rea)
    if(length(sel) == 0) {
        rea_lin <- c()
        rea_nonl <- reactions
    } else {
        rea_lin <- reactions[sel]
        rea_nonl <- reactions[!sel]
    }

    cat(length(rea_hv), "photochemical,", length(rea_lin), "linear (no loops) and", 
        length(rea_nonl), "nonlinear reactions!\n")

    if(type_spec_dup && rm_dup) {
        # first remove duplicates inside of the different types
        rea_hv <- s(rea_hv)
        y <- dup(rea_hv)
        rea_hv <- rea_hv[!y]

        rea_lin <- s(rea_lin) 
        y <- dup(rea_lin)
        rea_lin <- rea_lin[!y]

        rea_nonl <- s(rea_nonl)
        y <- dup(rea_nonl)
        rea_nonl <- rea_nonl[!y]

        if(!allow_direct_backflow) {
            # combine linear and hv reactions and shuffle
            #rea_lin_hv <- s(c(rea_lin, rea_hv))
            rea_lin_hv <- s(c(rea_lin, rea_hv))
            # identify duplicates to remove (x) and to keep (x_)
            x <- dup(rea_lin_hv)
            x_ <- r(dup(r(rea_lin_hv)))
            # 
            dp_rm <- rea_lin_hv[x]   # duplicates that are removed
            dp_keep <- rea_lin_hv[x_]  # duplicates that are kept

            rea_lin <- rea_lin[! rea_lin %in% dp_rm]
            rea_hv_keep <- rea_hv[rea_hv %in% dp_keep]
            rea_hv <- rea_hv[! (rea_hv %in% dp_rm | rea_hv %in% dp_keep)]

            # similar as above / combine nonlinear with (undecided!) hv
            rea_nonl_hv <- c(rea_hv_keep, s(c(rea_hv, rea_nonl)))
            #cat("rea_nonl_hv=", rea_nonl_hv, "\n")
            dp_rm_ <- rea_nonl_hv[dup(rea_nonl_hv)]

            rea_nonl <- rea_nonl[! rea_nonl %in% dp_rm_]
            rea_hv <- c(rea_hv_keep, rea_hv[! rea_hv %in% dp_rm_])
        }
    
        cat("After removing type specific duplicates there are", length(rea_hv), 
            "photochemical,", length(rea_lin), "linear (no loops) and", 
            length(rea_nonl), "nonlinear reactions!\n")
    }

    if(length(rea_lin) < no_2fold)
        no_2fold <- length(rea_lin)

    if(length(rea_hv) < no_hv)
        no_hv <- length(rea_hv)

     rea_lin <- rea_lin[sample(length(rea_lin), no_2fold)]
     rea_hv <- rea_hv[sample(length(rea_hv), no_hv)]
     if(M-no_hv-no_2fold > 0 && M-no_hv-no_2fold < length(rea_nonl))
         rea_nonl <- rea_nonl[sample(length(rea_nonl), M-no_hv-no_2fold)]


    possible_reas_to_jrnf <- function(s_rea) {
        M_ <- length(s_rea)

        # create network object and return it...
        e <- list()
        em <- list()
        p <- list()
        pm <- list()

        for(i in 1:length(s_rea)) {
            e_ <- c()
            em_ <- c()
            p_ <- c()
            pm_ <- c()

            x <- trans$from_linear(s_rea[i])
            if(x[1] != 1) {
                e_ <- c(e_, x[1]-1)
                em_ <- c(em_, 1)
            }

            if(x[2] != 1) {
                e_ <- c(e_, x[2]-1)
                em_ <- c(em_, 1)
            }

            if(x[3] != 1) {
                p_ <- c(p_, x[3]-1)
                pm_ <- c(pm_, 1)
            }

            if(x[4] != 1) {
                p_ <- c(p_, x[4]-1)
                pm_ <- c(pm_, 1)
            }
            e[[i]] <- e_
            em[[i]] <- em_
            p[[i]] <- p_
            pm[[i]] <- pm_
        }

        

        species <- data.frame(type=as.integer(rep(0, N+1)), name=as.character(name), 
                              energy=as.numeric(energy),
                              constant=as.logical(rep(F, N+1)),
                              stringsAsFactors=FALSE)

        species$constant[hv_id] <- T

        reactions <- data.frame(reversible=as.logical(rep(T, M_)),
                                c=as.numeric(rep(0, M_)), 
                                k=as.numeric(rep(0, M_)),
                                k_b=as.numeric(rep(0,M_)), 
                                activation=as.numeric(rplancklike(M_)), 
                                educts=I(e), educts_mul=I(em),
                                products=I(p), products_mul=I(pm))

        return(list(species, reactions))
    }
   
    net <- possible_reas_to_jrnf(c(rea_hv, rea_lin, rea_nonl))
    net$composition <- comp$composition
    return(net)
}


# legace name  // TODO remove when not needed any more
jrnf_create_artificial_ecosystem <- jrnf_ae_create


#
# NEW CODE FOR ARTIFICIAL ECOSYSTEM EVOLUTION
#


# Function creates an anorganic core. In principle this corresponds to creating
# an artificial ecosystem with jrnf_ae_create.
#
#

jrnf_ae_create_anorganic_core <- function(N, M, no_2fold, no_hv, comp_no,
                                          oi_N=c(), o_N=c(), o_M=c(), 
                                          o_no_2fold=c(), o_no_hv=c()) {
    net <- jrnf_ae_create(N, M, no_2fold, no_hv, comp_no) 
    net$para <- list(an = list(N=N, M=M, no_2fold=no_2fold, no_hv=no_hv, comp_no=comp_no),
                     org = list(i_N=oi_N, N=o_N, M=o_M, no_2fold=o_no_2fold, 
                                no_hv=o_no_hv, next_id=1))
    net$assoc <- list(sp=rep(0,nrow(net[[1]])), re=rep(0,nrow(net[[2]])))

    return(net)   
}



jrnf_ae_add_organism <- function(net, 
                                 oi_N=c(), o_N=c(), o_M=c(), 
                                 o_no_2fold=c(), o_no_hv=c()) {
    prefix_table <- c("X_", "Y_", "Z_", "U_", "V_", "W_", "K_", "L_", "M_", "I_") 

    if(is.null(net$para$org$i_N)) # if parameters are not available yet - update
        net$para$org <- list(i_N=oi_N, N=o_N, M=o_M, no_2fold=o_no_2fold, 
                             no_hv=o_no_hv, next_id=net$para$org$next_id)

    next_id <- net$para$org$next_id; N <- nrow(net[[1]]); M <- nrow(net[[2]])

    # select net$para$org$i_N interacting species
    i <- sample(which(net$assoc$sp == 0), net$para$org$i_N)
    # construct composition object for current net - and extend
    comp <- list(composition=net$composition, name=net[[1]]$name, energy=net[[1]]$energy)
    cc <- list(composition=matrix(comp$composition[i,], ncol=ncol(comp$composition)), name=comp$name[i], energy=comp$energy[i])
    c <- hcae_extend_elem_composition(net$para$org$N, cc, 
                                      sp_n_pre=prefix_table[net$para$org$next_id])

    cat("extended composition (private part) contains ", nrow(c$composition) - nrow(cc$composition), "\n")

    # create subnet for organism
    net_ <- jrnf_ae_create(net$para$org$i_N+net$para$org$N, net$para$org$M, net$para$org$no_2fold, 
                           net$para$org$no_hv, c)

    # FOR DIAGNOSTICS TODO REMOVE
    #net$net_old <- net
    #net$net_new <- net_
    #net$c_new <- c

    # join new organism subnet with existing network
    net <- jrnf_merge_net(net, net_)

    # update $composition as well as $assoc
    net$composition <- rbind(comp$composition, matrix(c$composition[-(1:length(i)),], ncol=ncol(c$composition)))
    net$assoc$sp <- c(net$assoc$sp, rep(next_id,nrow(net[[1]])-N))
    net$assoc$re <- c(net$assoc$re, rep(next_id,nrow(net[[2]])-M))

    # update next free organism id
    x <- net$assoc$sp
    net$para$org$next_id <- min(which(!(1:(max(x)+1) %in% x)))

    return(net)
}


jrnf_ae_remove_organism <- function(net, id) {
    # first remove all associated reactions
    keep_re <- net$assoc$re != id
    net[[2]] <- net[[2]][keep_re,]
    net$assoc$re <- net$assoc$re[keep_re]  

    # now one can simply remove all associated species
    keep_sp <- net$assoc$sp != id
    net$assoc$sp <- net$assoc$sp[keep_sp]
    net$composition <- matrix(net$composition[keep_sp,], ncol=ncol(net$composition))
    net <- jrnf_subnet(net, keep_sp)

    # update next free organism id
    x <- net$assoc$sp
    net$para$org$next_id <- min(which(!(1:(max(x)+1) %in% x)))
    return(net)
}


jrnf_ae_replace_organism <- function(net, id, 
                                     oi_N=c(), o_N=c(), o_M=c(), 
                                     o_no_2fold=c(), o_no_hv=c()) {
    net <- jrnf_ae_remove_organism(net, id)

    return(jrnf_ae_add_organism(net, 
                                oi_N, o_N, o_M, o_no_2fold, o_no_hv))
}







#- Function that adds an additional organism to an network out of environment + organisms
#  Parameters are: number of interacting species, number of additional species, number of reactions, number of linear reactions, number of photochemical ones



#- Function for removing an organism



