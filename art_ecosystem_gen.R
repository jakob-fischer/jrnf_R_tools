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

jrnf_ae_draw_energies <- function(net, flat_energies=F) {
    net[[1]]$energy <- rnorm(nrow(net[[1]]))    
    net[[2]]$activation <- rplancklike(nrow(net[[2]]))

    if(flat_energies) {
        net[[1]]$energy <- rep(0, nrow(net[[1]]))    
        net[[2]]$activation <- rep(1, nrow(net[[2]]))    
    }
  
    if(any(net[[1]]$name == "hv")) 
        net[[1]]$energy[net[[1]]$name == "hv"] <- max(50, max(net[[1]]$energy)+5)

    return(net)
}


# Helper to create names consistent with elementary constituents
hcae_create_name <- function(comp, c_names, ex_names, empty_name="X") {
    name <- ""

    # first build name
    if(sum(comp) == 0) 
        name <- empty_name
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
    #cat("hcrc rea=", rea, "  comp=", dim(comp), "  N=", N, "\n")
    get_c <- function(i) {
        if(rea[i] == 1)
            return(rep(0, ncol(comp)))
        else
            return(comp[rea[i]-1,])
    }

    # the components that are hv are set to zero / empty species first
    for(i in 1:4)
        if(rea[i] == N+2) 
            rea[i] <- 1

    return(all(get_c(1)+get_c(2) == get_c(3)+get_c(4)))
}


# TODO document

hcae_draw_elem_composition <- function(N, comp_no, max_tot=N/3, max_single=2, 
                                         lambda=0.1, allow_massless=F, max_stoich=10) {
    composition <- matrix(0, nrow=N, ncol=comp_no)
    energy <- rnorm(N)         
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

    for(i in 1:nrow(composition)) { 
        x <- sample(1:length(prob), 1, F, prob/sum(prob))   

        composition[i,] <- trans$from_linear(x)-1

        cnt[x] <- cnt[x] + 1
        if(cnt[x] > 1) {
            cnt_t <- cnt_t + 1
             
            if(cnt_t >= max_tot)
                prob[cnt != 0] <- 0
        }
        
        if(cnt[x] >= max_single)
            prob[x] <- 0
    }

    # now calculate names
    name <- c()
    for(i in 1:N)
        name <- c(name, hcae_create_name(composition[i,], component_names,
                                         name, empty_name))   
 
    return(list(composition=composition, name=name, energy=energy))
}


# If reaction is valid (by composition) returns max mass transfer
hcae_calc_rea_transfer <- function(rea, comp, N) {
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
hcae_check_rea_conditions <- function(rea, N) {
    hv_ed <- sum((rea[1] == N+2) + (rea[2] == N+2))
    hv_pro <- sum((rea[3] == N+2) + (rea[4] == N+2))
    empty_ed  <- sum((rea[1] == 1) + (rea[2] == 1))
    empty_pro  <- sum((rea[3] == 1) + (rea[4] == 1))

    # hv is on both sides
    if(hv_ed+hv_pro > 1)
        return(F)

    # One side of reaction has only hv
    if(hv_ed+empty_ed > 1 || hv_pro+empty_pro > 1)
        return(F)

    # One side is completely empty
    if(empty_ed == 2 || empty_pro == 2)
        return(F)

    # Reaction doesn't do anything  (TODO second part of condition excluded by next check?)
    if(rea[1] == rea[3] && rea[2] == rea[4] || rea[1] == rea[4] && rea[2] == rea[3])
        return(F)

    # Also avoid double reactions (through reverse reaction)
    if(rea[2] > rea[1] || rea[4] > rea[3] || rea[3] > rea[1])
        return(F) 

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


#
#
#
# TODO Function needs much cleanup #
#      Especially the way the elementary composition is drawn needs to be made easier (parameter of distribution!)

jrnf_create_artificial_ecosystem <- function(N, M, no_2fold, no_hv, comp_no, mod_no=0, mod_f=1, 
                                             no_reordering=T, enforce_transfer=T) { 
    hv_name <- "hv"
    trans <- b_mdim_transf(rep(N+2, 4))    

    # draw compostition (el. constituents) and energy of species
    cat("Drawing elementary composition")
    x <- hcae_draw_elem_composition(N, comp_no)
    cat(".\n")

    composition <- x$composition
    energy <- x$energy
    name <- x$name

    # Add species for photons / energy source
    energy <- c(energy, max(max(energy)+5, 50))
    name <- c(name, hv_name) 

    # TODO check parameters   (N/mod_no has to be a natural number)!
    # module of each species / first entry is (NO SPECIES and last hv)
    #sp_mod_id <- c(mod_no+1, floor((1:N-1)/(N/mod_no)), mod_no+2)    
    #is_mod= mod_no!=0 && mod_no != 1

    is_rea_possible <- function(i) {
        x <- trans$from_linear(i)

        if(!hcae_check_rea_constituents(x, composition, N) ||
           !hcae_check_rea_conditions(x, N))
            return(F)
            
        return(T)
    }

    reactions <- which(sapply(1:((N+2)**4), is_rea_possible))

    cat("found", length(reactions), "valid reactions!\n") 


    possible_reas_to_jrnf <- function(s_rea, s) {
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

            x <- trans$from_linear()
            if(a != 1) {
                e_ <- c(e_, a-1)
                em_ <- c(em_, 1)
            }

            if(b != 1) {
                e_ <- c(e_, b-1)
                em_ <- c(em_, 1)
            }

            if(c != 1) {
                p_ <- c(p_, c-1)
                pm_ <- c(pm_, 1)
            }

            if(d != 1) {
                p_ <- c(p_, d-1)
                pm_ <- c(pm_, 1)
            }
            e[[i]] <- e_
            em[[i]] <- em_
            p[[i]] <- p_
            pm[[i]] <- pm_
        }

        species <- data.frame(type=as.integer(rep(0, N+1)), name=as.character(name), 
                              energy=as.numeric(energy),
                              constant=as.logical(c(rep(F, N), T)),
                              stringsAsFactors=FALSE)

        reactions <- data.frame(reversible=as.logical(rep(T, M_)),
                                c=as.numeric(rep(0, M_)), 
                                k=as.numeric(rep(0, M_)),
                                k_b=as.numeric(rep(0,M_)), 
                                activation=as.numeric(rplancklike(M_)), 
                                educts=I(e), educts_mul=I(em),
                                products=I(p), products_mul=I(pm))

        return(list(species, reactions))
    }
}
