# author: jakob fischer (jakob@automorph.info)
# description: 
# Generation of artificial ecosystems. AE are reaction networks with an 
# elementary composition of species that are driven thermodynamically by 
# photochemical reactions. The second part of this module allows to 
# generate networks that consist of an "anorganic" part and various 
# organisms that can be replaced separately. This allows to evolve these
# systems (-> "simulation_builder_evol.R"). 
#    The general philosophy for artifical ecosystem generation here is to
# first generate an elementary composition. Afterwards ALL legal reactions
# with two or less educts and two or less products are determined and the
# network is creating by drawing a subset. This limits the size of the 
# networks that can be created to a size of around 50 species.

sourced_art_ecosystem_gen <- T

if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")


# Function draws energies for artificial ecosystems. Activation energies are 
# drawn from "rplancklike" distribution and standard chemical potentials are 
# drawn from a gauß distribution. Except the chemical potential for "hv" 
# which is either 50 or the maximum of all other chemical potentials + 5...
#
# parameters:
# <net>           - Reaction network
# <flat_energies> - Simple energies drawn? "0" for chemical pot. "1" for activation.
# <limit_AE>      - Set an limit for activation energy. Either give "F" for "no limit",
#                   "T" for limit "3" or a numeric value of the limit.

jrnf_ae_draw_energies <- function(net, flat_energies=F, limit_AE=F) {
    lAE <- NA
    if(limit_AE) lAE <- 3
    if(is.numeric(limit_AE)) lAE <- limit_AE

    # draw energies
    if(flat_energies) {
        net[[1]]$energy <- rep(0, nrow(net[[1]]))    
        net[[2]]$activation <- rep(1, nrow(net[[2]]))    
    } else {
        net[[1]]$energy <- rnorm(nrow(net[[1]]))    
        net[[2]]$activation <- rplancklike(nrow(net[[2]]), lAE)
    }
  
    # energy of hv is 50 + at least 5 higher than all other species
    if(any(net[[1]]$name == "hv"))
        net[[1]]$energy[net[[1]]$name == "hv"] <- max(50, max(net[[1]]$energy)+5)

    return(net)
}

# HELPER for jrnf_ae_extend_el_comp
# Function creates names consistent with elementary constituents.
#
# parameters:
# <comp>       - Composition vector (how much of each elementary component)
# <el_names>   - Names of elementary components (same length than <comp>)
# <ex_names>   - Vector of existing names
# <empty_name> - Name for empty species (NO elementary components)
# <prefix>     - Prefix that is put in front of the name

hcae_create_name <- function(comp, el_names, ex_names, empty_name="X", prefix="") {
    name <- prefix

    # first build name
    if(sum(comp) == 0) 
        name <- paste(name, empty_name, sep="")
    else 
        for(i in 1:length(comp))
            if(comp[i] != 0) 
                name <- paste(name, el_names[i], comp[i], sep="")    


    # second step - unify name (by appending "_2", "_3", ...
    if(any(name == ex_names)) {    # name already used
        i <- 2
        while(any(paste(name, "_", i, sep="") == ex_names))
            i <- i+1

        name <- paste(name, "_", i, sep="")
    }

    return(name)  
}


# This function draws an elementary composition for <N> species that are build from
# the number of elementary components in <comp$components>. 
#
# parameters:
# <N>              - Number of species
# <comp>           - Composition object (data.frame)
# <max_tot>        - Maximal allowed number of duplicates
# <max_singles>    - Maximum number of duplicates for a single species
# <lambda>         - Parameter in the poisson distribution from which 
#                    stoichiometric coefficents are sampled
# <allow_massless> - Are massless species allowed
# <max_stoich>     - Maximal stoichiometric coefficient
# <sp_n_pre>       - Prefix for species names (see "hcae_create_name")

jrnf_ae_extend_el_comp <- function(N, comp, max_tot=c(), max_single=c(), 
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

jrnf_ae_draw_el_comp <- function(N, comp_no, max_tot=c(), max_single=c(), 
                                         lambda=0.1, allow_massless=F, max_stoich=10) {
    comp <- list(composition=matrix(0, nrow=0, ncol=comp_no), name=c(), energy=c())

    return(jrnf_ae_extend_el_comp(N, comp, max_tot, max_single, lambda, allow_massless, max_stoich))
}


# Function reconstructs a composition object from the network <net> similar to 
# the one returned from "jrnf_ae_draw_el_comp". Additionally a fourth element 
# "elements" is added that contains the names of the elementary components.	

jrnf_ae_reconstruct_comp <- function(net) {
    names <- net[[1]]$name

    constituents <- c()

    # First remove tailing "_<number>" and extract components names
    for(i in 1:length(names)) {
        names[i] <- strsplit(names[i], "_", fixed=TRUE)[[1]][1]
        constituents <- c(constituents, strsplit(names[i], "[0-9]+")[[1]])
    }

    constituents <- sort(unique(constituents))
    constituents <- constituents[constituents != "hv"]

    m <- matrix(0, ncol=length(constituents), nrow=length(names))


    # Now extract component's names again / including multiplicity
    for(i in 1:length(names)) 
        if(names[i] != "hv") {
            con <- strsplit(names[i], "[0-9]+")[[1]]
            mul <- as.numeric(strsplit(names[i], "[A-Za-z]+")[[1]][-1])

            if(length(con) != length(mul))
                cat("error: length of con and mul have to match (jrnf_analyze_ecosystem_constituents)!\n")

            for(j in 1:length(con)) {
                sel <- which(constituents == con[j])
                m[i,sel] <- mul[j]
            }
        }

    return(list(composition=m, name=net[[1]]$name, energy=net[[1]]$energy, elements=constituents))
} 


# HELPER for jrnf_ae_create
# Checks if a reaction is possible in terms of elementary constituents.
#
# parameters:
# <rea>    - Reaction (integer representation)
# <comp>   - Components data frame

hcae_check_rea_constituents <- function(rea, comp) {
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


# HELPER TO jrnf_ae_create
# Check further conditions on reactions besides elementary constituents.
#
# parameters:
# <rea>    - Reaction (integer representation)
# <comp>   - Components data frame

hcae_check_rea_conditions <- function(rea, comp) {
    hv_id <- which(comp$name == "hv") + 1
    empty_id <- 1

    # how many of the eucts and products are empty
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

    return(T)
}


# HELPER for jrnf_ae_create
# This function receives a list of reactions representated as a vector of
# integers <s_rea> and transforms it into an jrnf-network.
#
# parameters:
# <s_rea>  - List of reactions (vector of integers)
# <comp>   - Components data frame
# <AE_max> - Maximum activation energy 

hcae_reas_to_jrnf <- function(s_rea, comp, AE_max) {
        # rebuild transformation
        trans <- b_mdim_transf(rep(nrow(comp$composition)+1, 4))  
        M_ <- length(s_rea)   # shortcut

        # create lists for educts, educts multipliers, products, 
        # products multipliers
        e <- list()
        em <- list()
        p <- list()
        pm <- list()

        for(i in 1:M_) {
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

        # Now build network object. Fields for reaction constants are left zero.
        species <- data.frame(type=as.integer(rep(0, nrow(comp$composition))), 
                              name=as.character(comp$name), 
                              energy=as.numeric(comp$energy),
                              constant=as.logical(rep(F, nrow(comp$composition))),
                              stringsAsFactors=FALSE)

        species$constant[comp$name == "hv"] <- T

        reactions <- data.frame(reversible=as.logical(rep(T, M_)),
                                c=as.numeric(rep(0, M_)), 
                                k=as.numeric(rep(0, M_)),
                                k_b=as.numeric(rep(0,M_)), 
                                activation=as.numeric(rplancklike(M_, AE_max)), 
                                educts=I(e), educts_mul=I(em),
                                products=I(p), products_mul=I(pm))

        return(list(species, reactions))
    }




# Function creates an artificial ecosystem randomly 
# <N> - number of species
# <M> - total number of reactions
# <no_2fold> - number of linear reactions
# <no_hv> - number of photochemical reactions
# <comp> - number of elementary components the set of species has 
#          (elementary composition is generated by jrnf_ae_draw_el_comp
#          using standard parameters.) Alternatively one can give a composition
#          created by jrnf_ae_draw_el_comp directly. Note that this method
#          won't add a hv-pseudospecies then but it has to be included.
# <cat_as_lin> - flag that indicates whether catalytical reactions (A + C -> B + C)
#                are treated identical to linear reactions in terms of duplicates 
#                and in terms of the reactions of the different types that are drawn
# <type_spec_dup> - duplicates are treated specific to the three types and not 
#                   accross them
# <rm_dup> - Remove duplicates. Only allow one reaction for each column in the 
#            stoichiometric matrix. Influenced by <type_spec_dup> and <cat_as_lin>
#
# TODO At the moment <N> does not include the hv-pseudospecies. If a composition is given 
#      (<comp>) <N> is set minus one of the number of species in this composition.

jrnf_ae_create <- function(N, M, no_2fold, no_hv, comp=c(), cat_as_lin=F, 
                           type_spec_dup=F, rm_dup=T, AE_max=3) { 
    hv_name <- "hv"

    # If no cmposition is given draw compostition (el. constituents) and energy of species
     if(is.numeric(comp)) {    
        cat("Drawing elementary composition")
        comp <- jrnf_ae_draw_el_comp(N, comp)
        cat(".\n")

        # Add species for photons / energy source
        comp$composition <- rbind(comp$composition, 
                                  matrix(0, nrow=1, ncol=ncol(comp$composition)))
        comp$energy <- c(comp$energy, max(max(comp$energy)+5, 50))
        comp$name <- c(comp$name, hv_name) 
    } else 
        # This assumes
        N <- nrow(comp$composition)-1

    hv_id <- which(comp$name == hv_name)  # shortcut
    # Create transformation object from integer to vector (each element 
    # one reactant) representation.
    trans <- b_mdim_transf(rep(nrow(comp$composition)+1, 4))   


    # SUBFUNCTIONS (access previously defined values)

    # Check if the reaction is valid
    is_rea_possible <- function(i) {
        x <- trans$from_linear(i)
        return(hcae_check_rea_constituents(x, comp) &&
               hcae_check_rea_conditions(x, comp))
    }

    # Check <i> is photoreaction
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
    # greater zero.  
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

        # Change of hv is set to zero (so backflow is a duplicate of 
        # "forward flow" and "hv + A -> B" and "hv + B -> A" won't cooccur).
         x[1,hv_id] <- 0

        return(x)
    }

    # Randomly permutate x
    s <- function(x) {
        return(x[sample(length(x))])
    }

    # Find the duplicates (in terms of effective change of species) reactions
    dup <- function(xx) {
        return(duplicated(t(sapply(xx, fval_to_stoichcol))))
    }

    # SUBFUNCTIONS END

   
    # Check for all possible reaction if they are legal.
    reactions <- which(sapply(1:((nrow(comp$composition)+1)**4), is_rea_possible))
    cat("found", length(reactions), "valid reactions!\n")
    # Permutate legal reactions
    reactions <- s(reactions)

    if(!type_spec_dup && rm_dup) {
        reactions <- s(reactions)
        x <- dup(reactions)
        reactions <- reactions[!x]
        cat("after removing (effective) doubles, there are ", length(reactions), "valid reactions!\n")
    }

    # Separate photoreactions, linear and nonlinear reactions
    sel <- sapply(reactions, is_hv_rea)
    rea_hv <- reactions[sel]
    reactions <- reactions[!sel]

    sel <- sapply(reactions, is_lin_rea)
    rea_lin <- reactions[sel]
    rea_nonl <- reactions[!sel]

    cat(length(rea_hv), "photochemical,", length(rea_lin), "linear (no loops) and", 
        length(rea_nonl), "nonlinear reactions!\n")

    if(type_spec_dup && rm_dup) {
        # Remove duplicates inside of the different types
        rea_hv <- s(rea_hv)
        y <- dup(rea_hv)
        rea_hv <- rea_hv[!y]

        rea_lin <- s(rea_lin) 
        y <- dup(rea_lin)
        rea_lin <- rea_lin[!y]

        rea_nonl <- s(rea_nonl)
        y <- dup(rea_nonl)
        rea_nonl <- rea_nonl[!y]
    
        cat("After removing type specific duplicates there are", length(rea_hv), 
            "photochemical,", length(rea_lin), "linear (no loops) and", 
            length(rea_nonl), "nonlinear reactions!\n")
    }

    # Ensure number of linear reactions and photoreactions not to exceed the
    # number of available (legal) reactions
    if(length(rea_lin) < no_2fold)
        no_2fold <- length(rea_lin)

    if(length(rea_hv) < no_hv)
        no_hv <- length(rea_hv)

    # Sample reactions for all three different types
    rea_lin <- rea_lin[sample(length(rea_lin), no_2fold)]
    rea_hv <- rea_hv[sample(length(rea_hv), no_hv)]
    if(M-no_hv-no_2fold > 0 && M-no_hv-no_2fold < length(rea_nonl))
        rea_nonl <- rea_nonl[sample(length(rea_nonl), M-no_hv-no_2fold)]

    # This should not really happen, but ensure consistency....
    if(M <= no_hv+no_2fold) 
        rea_nonl <- c()
   
    # transform chosen reaction into jrnf network format, add composition and
    # return the compiled object
    net <- hcae_reas_to_jrnf(c(rea_hv, rea_lin, rea_nonl), comp, AE_max)
    net$composition <- comp$composition
    return(net)
}


#
# CODE FOR GENERATING EVOLVABLE ARTIFICIAL ECOSYSTEMS
#

# Function creates an anorganic core. In principle this corresponds to creating
# an artificial ecosystem with jrnf_ae_create. But also parameters for creation
# of organisms can be appended to the network. And arrays that indicate to which
# organism (anorganic network) species and reactions belong to are added. In 
# principle parameter match with "jrnf_ae_create". For those not present here 
# "jrnf_ae_create" is called with standard values.
#
# parameters:
# <N>           - number of anorganic species 
# <M>           - total number of (anorganic) reactions (if possible)
# <no_2fold>    - number of (anorganic) reactions (if possible)
# <no_hv>       - number of (anorganic) photoreactions (if possible)
# <comp_no>     - number of elemenetary components
# <oi_N>        - number of (anorganic) species that will be choosen as interacting
#                 species for an organism (can take part in organisms reactions)
# <o_N>         - number of organisms internal species 
# <o_M>         - total number of reactions (of one organism)
# <o_no_2fold>  - number of linear reactions (of one organism)
# <o_no_hv>     - number of photoreactions (if possible / org. has to have "hv" 
#                 as interacting sp)
# <AE_max>      - maximum activation energy 

jrnf_ae_create_anorganic_core <- function(N, M, no_2fold, no_hv, comp_no,
                                          oi_N=c(), o_N=c(), o_M=c(), 
                                          o_no_2fold=c(), o_no_hv=c(), AE_max=3) {
    net <- jrnf_ae_create(N, M, no_2fold, no_hv, comp_no, AE_max=AE_max) 
    net$para <- list(an = list(N=N, M=M, no_2fold=no_2fold, no_hv=no_hv, comp_no=comp_no),
                     org = list(i_N=oi_N, N=o_N, M=o_M, no_2fold=o_no_2fold, 
                                no_hv=o_no_hv, next_id=1),
                     AE_max=AE_max)
    net$assoc <- list(sp=rep(0,nrow(net[[1]])), re=rep(0,nrow(net[[2]])))

    return(net)   
}


# Function that adds an additional organism to an network out of environment + organisms
# Parameters are taken from the network object ("net$para$org") if not given directly. 
# See "jrnf_ae_create_anorganic_core" for a description of the parameter. 

jrnf_ae_add_organism <- function(net, 
                                 oi_N=c(), o_N=c(), o_M=c(), 
                                 o_no_2fold=c(), o_no_hv=c()) {
    prefix_table <- c("X_", "Y_", "Z_", "U_", "V_", "W_", "K_", "L_", "M_", "I_", "J_", 
                      "x_", "y_", "z_", "u_", "v_", "w_", "k_", "l_", "m_", "i_", "j_") 

    if(is.null(net$para$org$i_N)) # if parameters are not available yet - update
        net$para$org <- list(i_N=oi_N, N=o_N, M=o_M, no_2fold=o_no_2fold, 
                             no_hv=o_no_hv, next_id=net$para$org$next_id)

    next_id <- net$para$org$next_id; N <- nrow(net[[1]]); M <- nrow(net[[2]])

    # select net$para$org$i_N interacting species
    i <- sample(which(net$assoc$sp == 0), net$para$org$i_N)
    # construct composition object for current net - and extend
    comp <- list(composition=net$composition, name=net[[1]]$name, energy=net[[1]]$energy)
    cc <- list(composition=matrix(comp$composition[i,], ncol=ncol(comp$composition)), name=comp$name[i], energy=comp$energy[i])
    c <- jrnf_ae_extend_el_comp(net$para$org$N, cc, 
                                      sp_n_pre=prefix_table[net$para$org$next_id])

    cat("extended composition (private part) contains ", nrow(c$composition) - nrow(cc$composition), "\n")

    # create subnet for organism
    net_ <- jrnf_ae_create(net$para$org$i_N+net$para$org$N, net$para$org$M, net$para$org$no_2fold, 
                           net$para$org$no_hv, c, AE_max=net$para$AE_max)

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


# Function for removing the organism with id <id>. 

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


# Simply replace an organism by removing it first and then adding a new one
# that is randomly generated. If no further parameters besides <id> are given
# the ones present under "net$para$org" are used. For a description of the 
# meaning of the parameters se above ("jrnf_ae_create_anorganic_core").

jrnf_ae_replace_organism <- function(net, id, 
                                     oi_N=c(), o_N=c(), o_M=c(), 
                                     o_no_2fold=c(), o_no_hv=c()) {
    net <- jrnf_ae_remove_organism(net, id)

    return(jrnf_ae_add_organism(net, 
                                oi_N, o_N, o_M, o_no_2fold, o_no_hv))
}














