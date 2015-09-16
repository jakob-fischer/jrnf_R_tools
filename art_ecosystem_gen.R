

source("jrnf_network.R")

# Function creates artificial ecosystem with <N> species... 
#



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
        if(rea[i] == 0)
            return(rep(0, ncol(comp)))
        else
            return(comp[rea[i],])
    }

    # the components that are hv are set to zero / empty species first
    for(i in 1:4)
        if(rea[i] == N+1) 
            rea[i] <- 0

    return(all(get_c(1)+get_c(2) == get_c(3)+get_c(4)))
}


# TODO DOCUMENTATION

hcae_draw_elem_composition <- function(N, comp_no, max_dup=c(), tries=1000) {
    composition <- matrix(0, nrow=N, ncol=comp_no)
    energy <- rnorm(N)       
    component_names <- c("C", "N", "O", "H", "P")
    component_names <- component_names[1:comp_no]
    if(is.null(max_dup))
        max_dup = N
        
    draw <- function(i) {
        if(comp_no == 1)
            composition[i,] <<- rpois(comp_no, 0.5)+1
        else
            composition[i,] <<- rpois(comp_no, 0.5)+sample(c(0,1), comp_no, T, c(comp_no-1,1))

        if(sum(composition[i,]) == 0)
            draw(i)
    }

    get_dup <- function() {
        return(which(duplicated(composition)))
    }

    for(i in 1:N)
        draw(i)

    x <- get_dup()
    while(length(x) > max_dup) {
        for(i in x[-(1:max_dup)])
            draw(i)
        x <- get_dup()
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
    #cat("hcrc rea=", rea, "  comp=", dim(comp), "  N=", N, "\n")
    get_c <- function(i) {
        if(rea[i] == 0)
            return(rep(0, ncol(comp)))
        else
            return(comp[rea[i],])
    }

    # the components that are hv are set to zero / empty species first
    for(i in 1:4)
        if(rea[i] == N+1) 
            rea[i] <- 0
   
    if(rea[1] != 0)
        return(min(sum(abs(get_c(1) - get_c(3))), 
                   sum(abs(get_c(1) - get_c(4)))))

    if(rea[2] != 0)
        return(min(sum(abs(get_c(2) - get_c(3))), 
                   sum(abs(get_c(2) - get_c(4)))))

    return(0)    
}

# Check further conditions on reactions 
hcae_check_rea_conditions <- function(rea, N) {
    hv_ed <- sum((rea[1] == N+1) + (rea[2] == N+1))
    hv_pro <- sum((rea[3] == N+1) + (rea[4] == N+1))
    empty_ed  <- sum((rea[1] == 0) + (rea[2] == 0))
    empty_pro  <- sum((rea[3] == 0) + (rea[4] == 0))

    # hv is on both sides
    if(hv_ed+hv_pro > 1)
        return(F)

    # One side of reaction has only hv
    if(hv_ed+empty_ed > 1 || hv_pro+empty_pro > 1)
        return(F)

    # Reaction doesn't do anything
    if(rea[1] == rea[3] && rea[2] == rea[4] || rea[1] == rea[4] && rea[2] == rea[3])
        return(F)

    # If one educt / product is empty it has to be the second one (avoid double occurence of probable reactions)
    #if(rea[1] == 0 && rea[2] != 0 || rea[3] == 0 && rea[4] != 0)
    #    return(F)

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
                                             no_reordering=F, enforce_transfer=T) { 

    # TODO check parameters   (N/mod_no has to be a natural number)!
    sp_mod_id <- floor((1:N-1)/mod_no)    # module of each species
    mod= mod_no!=0

    eval_possible_reas <- function(s) { 
        #cat("call to eval_possible_reas with s=", s, "\n")
        p <- rep(1, s^4) 
        has_hv <- rep(F, s^4)
        rea_no <- rep(0, s^4)
        transfer <- rep(0, s^4)
        elementary_valid <- rep(T, s^4)
        mod_structure <- rep(F, s^4)
        double_species <- rep(F, s^4)
        self_loop <- rep(F, s^4)

        # TODO: speed this up with apply?
        for(i in 1:length(p)) {
            j <- i-1
            a <- j%%s
            b <- ((j-a)/s)%%s
            c <- ((j-a-b*s)/s^2)%%s
            d <- ((j-a-b*s-c*s^2)/s^3)%%s 

            #cat("a=", a, " b=", b, " c=", c, " d=", d, "\n")

            # Increase probability of reactions "in" modular structure
            if(mod & a != 0 && c != 0 && a != N+1 && c != N+1 && sp_mod_id[a] == sp_mod_id[c] ||
               mod & b != 0 && d != 0 && b != N+1 && d != N+1 &&  sp_mod_id[b] == sp_mod_id[d])
                mod_structure[i] <- T

            # Decrease probability of autocatalytic self loops
            if(mod & a != 0 && b != 0 && c != 0 && d != 0 && 
               (a == c || a == d || b == c || b == d))
                self_loop[i] <- T
                
            # Check if reaction works by elementary constituents
            if(!hcae_check_rea_constituents(c(a,b,c,d), composition, N) ||
               !hcae_check_rea_conditions(c(a,b,c,d), N))
                elementary_valid[i] <- F
                

            # If reaction is valid, calculate transfer
            if(elementary_valid[i] != 0)
                transfer[i] <- hcae_calc_rea_transfer(c(a,b,c,d), composition, N)

            if(a == N+1 || b == N+1 || c == N+1 || d == N+1)
                has_hv[i] <- T

            if(a == b || c == d)
                double_species[i] <- T

            # CAREFUL: reactants occuring on both sides are not counted (have to be subtracted)
            rea_no[i] <- sum(c(a != 0 , b != 0, c != 0, d != 0)) - 
                         2*sum(c(a == c & a != 0, a == d & a != 0, b == c & b != 0, b == d & b != 0))

            if(mod_structure[i])
                p[i] <- p[i]/mod_f

            if(self_loop[i])
                p[i] <- p[i]/N

            if(double_species[i])
                p[i] <- p[i]/N

            if(!elementary_valid[i])
                p[i] <- 0

            if(enforce_transfer && !has_hv[i] && transfer[i] == 0 && rea_no[i] != 2)
                p[i] <- 0
        }

        p[rea_no == 0] <- 0
        p[rea_no == 1] <- 0

        return(list(p=p, has_hv=has_hv, rea_no=rea_no, el_transfer=transfer, 
                    el_valid=elementary_valid, mod_structure=mod_structure,
                    double_species=double_species, self_loop=self_loop))
    }


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

            j <- s_rea[i]-1
            a <- j%%s
            b <- ((j-a)/s)%%s
            c <- ((j-a-b*s)/s^2)%%s
            d <- ((j-a-b*s-c*s^2)/s^3)%%s

            if(a != 0) {
                e_ <- c(e_, a)
                em_ <- c(em_, 1)
            }

            if(b != 0) {
                e_ <- c(e_, b)
                em_ <- c(em_, 1)
            }

            if(c != 0) {
                p_ <- c(p_, c)
                pm_ <- c(pm_, 1)
            }

            if(d != 0) {
                p_ <- c(p_, d)
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

    hv_name <- "hv"
    
    # draw compostition (el. constituents) and energy of species
    cat("Drawing elementary composition")
    if(comp_no == 1)
        x <- hcae_draw_elem_composition(N, comp_no, as.integer(N*0.8))
    else
        x <- hcae_draw_elem_composition(N, comp_no, as.integer(N/4))

    composition <- x$composition
    energy <- x$energy
    name <- x$name

    # Add species for photons / energy source
    energy <- c(energy, max(max(energy)+5, 50))
    name <- c(name, hv_name) 
 
    # Now investigate all possible reactions up to 2x2 and build one 
    # vector containing information whether the reaction is possible
    x <- eval_possible_reas(N+2)
    possible_reas <- x[[1]]
    has_hv <- x[[2]]
    rea_no <- x[[3]]
    transfer <- x[[4]]

    if(mod_no != 0 && !no_reordering) {
        cat("DOING MODULARIZATION...\n")
        cat("having ", sum(possible_reas != 0), " possible reactions out of ", (N+2)^4, "!\n")
        cat(" ", sum(has_hv & possible_reas != 0), " of them are photoreactions.\n") 
        cat(" ", sum(rea_no[possible_reas != 0] == 2), " have 2 reactants, ", sum(rea_no[possible_reas != 0] == 3), 
            " have 3 reactants and ", sum(rea_no[possible_reas != 0] == 4), " have 4!\n")

        net <- possible_reas_to_jrnf(which(possible_reas != 0), N+2) 
        g <- jrnf_to_undirected_network(net)
        mat <- jrnf_graph_to_amatrix(g)
        mat <- mat[1:N,1:N]      

        o <- jrnf_get_modular_reordering(mat, mod_no)

        cat("REORDERING o=", o, "\n")  
        cat("nrow(composition) = ", nrow(composition), "\n")
        composition <- matrix(composition[o,], nrow=length(o))
        name <- c(name[o], name[length(name)])   

        # reevaluating
        x <- eval_possible_reas(N+2)
        possible_reas <- x[[1]]
        has_hv <- x[[2]]
        rea_no <- x[[3]]   
        transfer <- x[[4]]  
    }



    cat("having ", sum(possible_reas != 0), " possible reactions out of ", (N+2)^4, "!\n")
    cat(" ", sum(has_hv & possible_reas != 0), " of them are photoreactions.\n") 
    cat(" ", sum(rea_no[possible_reas != 0] == 2), " have 2 reactants, ", sum(rea_no[possible_reas != 0] == 3), 
        " have 3 reactants and ", sum(rea_no[possible_reas != 0] == 4), " have 4!\n")

    # now draw M reactions     # allow to select which fraction should be 1x1 reactions or how many photoreactions should be choosen

    {   # Now sample; we are sampling hv-containing reactions
        l_hv_reas <- which(possible_reas != 0 & has_hv)
        if(length(l_hv_reas) > no_hv)
            s_rea_hv <- sample(l_hv_reas, no_hv, F, possible_reas[l_hv_reas]/sum(possible_reas[l_hv_reas]))
        else {
            s_rea_hv <- l_hv_reas
            no_hv <- length(l_hv_reas)   # update no_hv so other types can take more
        }

        l_reas_2 <- which(possible_reas & rea_no == 2 & !has_hv)
        if(length(l_reas_2) > no_2fold)
            s_rea_2 <- sample(l_reas_2, no_2fold, F, possible_reas[l_reas_2]/sum(possible_reas[l_reas_2]))
        else {
            s_rea_2 <- l_reas_2      
            no_2fold <- length(l_reas_2)   # update no_2fold so other types can take more
        }

        l_reas_other <- which(possible_reas & rea_no > 2 & !has_hv)
        if(length(l_reas_other) > M-no_2fold-no_hv)
            s_rea_other <- sample(l_reas_other, M-no_2fold-no_hv, F, possible_reas[l_reas_other]/sum(possible_reas[l_reas_other]))
        else
            s_rea_other <- l_reas_other

        s_rea <- c(s_rea_hv, s_rea_2, s_rea_other)  
    }

    cat(" -> ", sum(has_hv[s_rea]), " photoreactions drawn!\n")
    cat(" -> ", sum(rea_no[s_rea] == 2), " reactions with (effectively) 2 reactants!\n")
    cat(" -> ", sum(rea_no[s_rea] == 3), " reactions with (effectively) 3 reactants!\n")
    cat(" -> ", sum(rea_no[s_rea] == 4), " reactions with (effectively) 4 reactants!\n")
    cat(" -> ", sum(transfer[s_rea] == 0), " reactions without (mass) transfer!\n")
    cat(sum(duplicated(composition)), " species are duplicates in terms of elementary composition!\n")

    net <- possible_reas_to_jrnf(s_rea, N+2)
    
    if(mod_no != 0) {
        g <- jrnf_to_directed_network(net)
        mat <- jrnf_graph_to_amatrix(g)[1:N,1:N]
        f <- gmr_get_inner_density(mat, mod_no)/(sum(mat)/N**2)
        cat("modularity factor is ", f, "\n")
    }

    # If no hv reactions are created on purpose, hv is removed from the
    # species list 
    if(no_hv == 0) 
        net <- list(net[[1]][1:N,], net[[2]])

    return(net)
}
