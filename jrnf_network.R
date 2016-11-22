# author: jakob fischer (jakob@automorph.info)
# description: 
# R-tools for handling jrnf-network data and using parts of igraph
# for network analysis and visualisation.
#
# Internally reaction networks are stored in the following format:
# A list of two data frames
# - The first data frame contains all species with 'type' (integer),
# energy (numeric), 'name' (character) and 'constant' (boolean).
# - The second data frame contains all reactions with information whether
# it is 'reversible' (boolean), the reaction constants 'c', 'k', 'k_b', 
# 'activation' energy (numeric), lists on 'products' and 'educts' and 
# their multiplicity in the reaction ('educts_mul' and 'products_mul')
#
# TODO: For some specific cases additional entries might be added to the list.
#       This is of course not saved in the official file format but only when
#       object is saved ad R-object. For artificial ecosystems for example the
#       species composition is saved in the $composition field as matrix. 
#       Unresolved issue is how methods that transform the network handle this.
#       Many rebuild the network and drop the additional data. This could be 
#       changed but then the metadata might become inconsistent.

sourced_jrnf_network <- T

library(igraph)

if(!exists("sourced_tools"))
    source("tools.R")    # general helper functions

if(!exists("sourced_cycles"))
    source("cycles.R")

if(!exists("sourced_jrnf_network_io"))
    source("jrnf_network_io.R")


# Transform all reactions to extended form where all multiplicities are 1

jrnf_expand <- function(net) {
    # separate species and reactions dataframe (easier to access)
    species <- net[[1]]
    reactions <- net[[2]]    
 
    for(i in 1:nrow(reactions)) { 
        educts <- reactions$educts[[i]]
        products <- reactions$products[[i]]

        if(length(educts) > 0)
        for(ed in 1:length(educts)) {
            # If multiplicity is higher 1 add aditional entries (with multiplicity 1)...
            if(as.integer(reactions$educts_mul[[i]][ed]) > 1) {
                for(k in 2:as.integer(reactions$educts_mul[[i]][ed])) {
                    reactions$educts[[i]] <- c(reactions$educts[[i]], reactions$educts[[i]][ed])
                    reactions$educts_mul[[i]] <- c(reactions$educts_mul[[i]], 1) 
                }    
            } 
            # ...and set multiplicity of original entry to 1.
            reactions$educts_mul[[i]][ed] <- 1
        }            


        if(length(products) > 0)
        for(prod in 1:length(products)) {
            # If multiplicity is higher 1 add aditional entries (with multiplicity 1)...
            if(as.integer(reactions$products_mul[[i]][prod]) > 1) {
                for(k in 2:as.integer(reactions$products_mul[[i]][prod])) {
                    reactions$products[[i]] <- c(reactions$products[[i]], reactions$products[[i]][prod])
                    reactions$products_mul[[i]] <- c(reactions$products_mul[[i]], 1) 
                }
            }
            # ...and set multiplicity of original entry to 1.
            reactions$products_mul[[i]][prod] <- 1
        }

    }

    # reassemble network
    return(list(species, reactions))
}


# Reverse of jrnf_expand. If some reactions are internally formulated
# like "A + A -> B + C + C" they are collapsed to "2 A -> B + 2 C". 

jrnf_collapse <- function(net) {   
    if(nrow(net[[2]]) > 0)
    for(i in 1:nrow(net[[2]])) { 
        educts <- unique(net[[2]]$educts[[i]])
        educts_mul <- net[[2]]$educts_mul[[i]][!duplicated(net[[2]]$educts[[i]])]
        products <- unique(net[[2]]$products[[i]])
        products_mul <- net[[2]]$products_mul[[i]][!duplicated(net[[2]]$products[[i]])]
        
        if(length(educts) > 0) {
            for(ed in which(duplicated(net[[2]]$educts[[i]]))) {
                # find educts and add up
                j <- which(educts == net[[2]]$educts[[i]][ed])
                educts_mul[j] <- educts_mul[j] + net[[2]]$educts_mul[[i]][ed]
            }  

            net[[2]]$educts[[i]] <- educts
            net[[2]]$educts_mul[[i]] <- educts_mul           
        }

        if(length(products) > 0) {
            for(prod in which(duplicated(net[[2]]$products[[i]]))) {
                # find products and add up
                j <- which(products == net[[2]]$products[[i]][prod])
                products_mul[j] <- products_mul[j] + net[[2]]$products_mul[[i]][prod]
            }      

            net[[2]]$products[[i]] <- products
            net[[2]]$products_mul[[i]] <- products_mul          
        }
    }

    return(net)    
}


# The function simplifies the representation of the reactions in the way that
# all second order terms (multiple occurence of same educt or product) are 
# simplified to first order. For example, the equation: "A + 3 B -> 3 A + C" is
# simplified to "A + B -> A + C".
# 
# This only considers the jrnf-list representation. If a reaction is saved in the
# form 1*X + 1*X -> 1*Y it is not simplified! For this use 
# "jrnf_simplify_multiplicity" (below).

jrnf_simplify_multiplicity_rep <- function(net) {
    species <- net[[1]]
    reactions <- net[[2]]    
 
    for(i in 1:nrow(reactions)) { 
        educts <- reactions$educts[[i]]
        products <- reactions$products[[i]]

        if(length(educts) > 0)
        for(ed in 1:length(educts)) 
            reactions$educts_mul[[i]][ed] <- 1  

        if(length(products) > 0)
        for(prod in 1:length(products)) 
            reactions$products_mul[[i]][prod] <- 1
    }

    return(list(species, reactions))
}


# The function simplifies reactions in the way that all second order terms 
# (multiple occurence of same educt or product) are simplified to first order. 
# For example, the equation: "A + 3 B -> 3 A + C" is simplified to 
# "A + B -> A + C".

jrnf_simplify_multiplicity <- function(net) {
    net <- jrnf_collapse(net)
    return(jrnf_simplify_multiplicity_rep(net))
}


# Randomize the reaction network by replacing each species id in each
# reaction with a randomly choosen chemical species (uniformly)
jrnf_randomize_species <- function(net) {
    net <- jrnf_expand(net)

    species <- net[[1]]
    reactions <- net[[2]]    
    N <- nrow(species) 

    for(i in 1:nrow(reactions)) { 
        educts <- reactions$educts[[i]]
        products <- reactions$products[[i]]

        # set ids of all educts of reaction 'i' randomly
        if(length(educts) > 0)
        for(ed in 1:length(educts))
            reactions$educts[[i]][ed] <- sample(1:N, 1)  

        # set ids of all products of reaction 'i' randomly
        if(length(products) > 0)
        for(prod in 1:length(products))  
            reactions$products[[i]][prod] <- sample(1:N, 1)  
    }

    return(list(species, reactions))
}


# Calculates the 'in' stoichiometric matrix of the jrnf-network net. 
# 'in' means every column contains the left side of one reaction. 
# All values are positive.

jrnf_calculate_stoich_mat_in <- function(net) {
    no_sp <- nrow(net[[1]])
    no_re <- nrow(net[[2]])
    N <- matrix(0, no_sp, no_re)
    
    for(i in 1:no_re) {
        l <- length(net[[2]]$educts[[i]])
        if(l != 0)
            for(j in 1:l) {
                e <- net[[2]]$educts[[i]][j]
                N[e, i] <- N[e, i] + net[[2]]$educts_mul[[i]][j]
            }
        }

    return(N)
}


# Calculates the 'out' stoichiometric matrix of the jrnf-network net. 
# 'out' means every column contains the right side of one reaction. 
# All values are positive.

jrnf_calculate_stoich_mat_out <- function(net) {
    no_sp <- nrow(net[[1]])
    no_re <- nrow(net[[2]])
    N <- matrix(0, no_sp, no_re)
    
    for(i in 1:no_re) {
        l <- length(net[[2]]$products[[i]])
        if(l != 0)
            for(j in 1:l) {
                p <- net[[2]]$products[[i]][j]
                N[p, i] <- N[p, i] + net[[2]]$products_mul[[i]][j]
            }
        }

    return(N)
}


# Calculates the stoichiometric matrix for a jrnf-network
# (matrix contains information on net change of species with reactions)
# The method has the additional semantic of, if given a matrix instead of a network 
# object, just returning it. This allows other functions which take network objects 
# or stoichiometric matrices as parameters using this function without having to 
# check if the object is a matrix themselfes.

jrnf_calculate_stoich_mat <- function(net_N) {
    if(is.matrix(net_N))
        return(net_N)

    return(jrnf_calculate_stoich_mat_out(net_N)-jrnf_calculate_stoich_mat_in(net_N))
}


# Calculate the rate of concentration change assuming the reactions are
# happening at the rates

jrnf_calculate_concentration_change <- function(network, rates) {
    N <- jrnf_calculate_stoich_mat(network)
    return(N %*% rates)
}


# Calculates the flow / reaction rates for a given concentration vector
# 

jrnf_calculate_flow <- function(network, concentrations) {
    cf_a <- function(reaction) {
        f_forward <- reaction$k
        f_backward <- reaction$k_b
        energy_dif <- 0 
 

        # Iterate educts
        for(j in 1:length(reaction$educts)) {
            id <- reaction$educts[j]
            mul <- reaction$educts_mul[j]
            energy_dif <- energy_dif - network[[1]]$energy[id]*mul
            f_forward <- f_forward * concentrations[id]^mul
        }

        # Iterate products
        for(j in 1:length(reaction$products)) {
            id <- reaction$products[j]
            mul <- reaction$products_mul[j]
            energy_dif <- energy_dif + network[[1]]$energy[id]*mul
            f_backward <- f_backward * concentrations[id]^mul
	}

        f_effective <- f_forward - f_backward
        energy_f <- energy_dif*f_effective
        if(is.finite(f_forward/f_backward) && f_forward != 0)
            entrop_p <- (f_forward-f_backward)*log(abs(f_forward/f_backward))   
        else
            entrop_p <- 0

        return(data.frame(flow_effective=as.numeric(f_effective), 
                          flow_forward=as.numeric(f_forward), 
                          flow_backward=as.numeric(f_backward), 
                          energy_flow=as.numeric(energy_f), 
                          entropy_prod=as.numeric(entrop_p)))

    }

    f_list <- apply(network[[2]], 1, cf_a)
    return(do.call("rbind", f_list))
}






calculate_flow_dif <- function(x, network, flow_effective) {
    dif <- 0

    for(i in 1:nrow(network[[2]])) {
        # Iterate educts
        for(j in 1:length(network[[2]]$educts[[i]])) {
            id <- network[[2]]$educts[[i]][j]
            if(id == x) {
                mul <- network[[2]]$educts_mul[[i]][j]
                dif <- dif + flow_effective[i]*mul  
                if(length(dif) != 1) {
                    cat("BADERROR\n")
                    cat("i=", i, "j=", j, "id=", id, "\n")
                }
            }
        }

        # Iterate products
        for(j in 1:length(network[[2]]$products[[i]])) {
            id <- network[[2]]$products[[i]][j]
            if(id == x) {
                mul <- network[[2]]$products_mul[[i]][j]
                dif <- dif - flow_effective[i]*mul
                if(length(dif) != 1) {
                    cat("BADERROR\n")
                    cat("i=", i, "j=", j, "id=", id, "\n")
                }
            }     
        }
    }

    return(dif)
}




#
#
#

jrnf_reverse_reactions <- function(net, rev) {
    if(!is.logical(rev))
        rev <- (rev < 0)

    for(i in 1:nrow(net[[2]])) {
        if(rev[i]) {
            # reverse topological part
            e <- net[[2]]$educts[[i]]
            e_m <- net[[2]]$educts_mul[[i]]
            p <- net[[2]]$products[[i]]
            p_m <- net[[2]]$products_mul[[i]]
         
            if(is.null(p)) 
                net[[2]]$educts[i] <- list(NULL)
            else
                net[[2]]$educts[[i]] <- I(p)

            if(is.null(p_m)) 
                net[[2]]$educts_mul[i] <- list(NULL)              
            else
                net[[2]]$educts_mul[[i]] <- I(p_m)

            if(is.null(e))
                net[[2]]$products[i] <- list(NULL)
            else
                net[[2]]$products[[i]] <- I(e)

            if(is.null(e_m))
                net[[2]]$products_mul[i] <- list(NULL)                    
            else
                net[[2]]$products_mul[[i]] <- I(e_m)

            # reverse reaction constants
            tmp <- net[[2]]$k[i]
            net[[2]]$k[i] <- net[[2]]$k_b[i]
            net[[2]]$k_b[i] <- tmp
        }
    }
    
    return(net)
}


# randomizes the reaction direction...

jrnf_randomize_dir <- function(net) {
    N <- nrow(net[[2]])
    x <- runif(N)*2-1
    return(jrnf_reverse_reactions(net, x))
}




jrnf_find_duplicate_reactions <- function(net) {
    if(nrow(net[[2]]) < 2)
        return(list())

    p <- list()
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)
    N <- N_out - N_in

    return(duplicated(N_in) & duplicated(N_out) & duplicated(N))
}


# Function returns a list of 2-element vectors of reactions that are the 
# reverse of each other 
 
jrnf_find_duplicate_reactions <- function(net) {
    if(nrow(net[[2]]) < 2)
        return(list())

    p <- list()
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)
    N <- N_out - N_in

    for(i in 1:(nrow(net[[2]])-1)) {
        for(j in (i+1):(nrow(net[[2]]))) {
            if(all(N_in[,i] == N_out[,j]) && all(N_out[,i] == N_in[,j]))
                p[[length(p)+1]] <- c(i,j)             
        }
    }

    return(p)


}


jrnf_find_reverse_pairs <- function(net) {
    if(nrow(net[[2]]) < 2)
        return(list())

    p <- list()
    N <- jrnf_calculate_stoich_mat(net)

    for(i in 1:(nrow(net[[2]])-1)) {
        for(j in (i+1):(nrow(net[[2]]))) {
            if(all(N[,i] == -N[,j]))
                p[[length(p)+1]] <- c(i,j)             
        }
    }

    return(p)
}


#
 
jrnf_remove_reverse_pairs <- function(net, rates) {
    sel <- rep(T, nrow(net[[2]]))        # flag of those reactions which are selected

    rp <- jrnf_find_reverse_pairs(net)   # list of pairs that are reverse of each other

    for(x in rp) {
        if(rates[x[1]] > rates[x[2]]) {
            sel[x[2]] <- F
            rates[x[1]] <- rates[x[1]] - rates[x[2]]
        } else {
            sel[x[1]] <- F
            rates[x[2]] <- rates[x[2]] - rates[x[1]]
        }
    }

    return(list(list(net[[1]], net[[2]][sel,]), rates[sel]))
}


associate_reaction_to_species <- function(network, sp_prop) {
    t <- list()
    for(i in 1:nrow(network[[1]]))
        t[[i]] <- as.numeric(c())

    for(i in 1:nrow(network[[2]])) {
        # Iterate educts
        for(j in 1:length(network[[2]]$educts[[i]])) {
            sp <- network[[2]]$educts[[i]][j]
            t[[sp]] <- c(t[[sp]], i)
        }

        # Iterate products
        for(j in 1:length(network[[2]]$products[[i]])) {
            sp <- network[[2]]$products[[i]][j]
            t[[sp]] <- c(t[[sp]], i)
	}
    }

    res <- c()

    for(x in t) {
        m <- sp_prop[x]
        if(length(m) == 0)
            res <- c(res, 0)
        else
            res <- c(res, mean(m))
    }
    
    return(res)   
}



jrnf_get_educts_mul <- function(net) {
    e <- c()

    for(i in 1:nrow(net[[2]])) {
        e <- c(e, sum(net[[2]]$educts_mul[[i]]))
    }
    
    return(e)
}

jrnf_get_products_mul <- function(net) {
    p <- c()

    for(i in 1:nrow(net[[2]])) {
        p <- c(p, sum(net[[2]]$products_mul[[i]]))
    }
    
    return(p)
}



# This function transforms a jrnf-network to a directed network by connecting
# two species / nodes if there is a reaction with the first as educt and the
# second as product

jrnf_to_directed_network <- function(jrnf_data, rnd=F) {
    jrnf_species <- jrnf_data[[1]]
    jrnf_reactions <- jrnf_data[[2]]    
    g <- graph.empty()
    g <- add.vertices(g, nrow(jrnf_species), 
                          name=as.vector(jrnf_species$name))   
 
    for(i in 1:nrow(jrnf_reactions)) { 
        educts <- jrnf_reactions$educts[[i]]
        products <- jrnf_reactions$products[[i]]
        educts_m <- jrnf_reactions$educts_mul[[i]]
        products_m <- jrnf_reactions$products_mul[[i]]

        if(length(educts) != 0 && length(products) != 0)    
        for(ed in 1:length(educts)) {
            for(prod in 1:length(products)) {
                if(rnd) {
                    b <- runif(1) > 0.5
                    if(b)
                        g <- add.edges(g, c(educts[ed], products[prod]))
                    else 
                        g <- add.edges(g, c(products[prod], educts[ed]))
                } else {
                    X <- as.integer(educts_m[ed])*as.integer(products_m[prod])
                    if(X != 0)
                    for(u in 1:(X))
                        g <- add.edges(g, c(educts[ed], products[prod]))
                }
            }
        } 
    }

    V(g)$label <- V(g)$name

    return(g)    
}




#
#
#

jrnf_to_directed_network_d <- function(jrnf_data, direction) {
    jrnf_species <- jrnf_data[[1]]
    jrnf_reactions <- jrnf_data[[2]]    
    g <- graph.empty()
    g <- add.vertices(g, nrow(jrnf_species), 
                          name=as.vector(jrnf_species$name))   

    for(i in 1:nrow(jrnf_reactions)) {
        for(ed in jrnf_reactions$educts[[i]]) {
            for(prod in jrnf_reactions$products[[i]]) {
                if(direction[i] > 0)
                    g <- add.edges(g, c(ed, prod))
                else
                    g <- add.edges(g, c(prod, ed))
            }
        } 
    }

    V(g)$label <- V(g)$name

    return(g)    
}









# TODO comment

jrnf_to_undirected_network <- function(jrnf_network) {
    return(as.undirected(jrnf_to_directed_network(jrnf_network), mode="each"))
}



# Transforms a jrnf-network (pair of data frames) to an undirected 
# igraph network where for every pair of species in every reaction on
# node is created.

jrnf_to_dense_network <- function(jrnf_data) {
    jrnf_species <- jrnf_data[[1]]
    jrnf_reactions <- jrnf_data[[2]]    
    g <- graph.empty()
    g <- add.vertices(g, nrow(jrnf_species), 
                          name=as.vector(jrnf_species$name))   

    for(i in 1:nrow(jrnf_reactions)) {
        for(ed in jrnf_reactions$educts[[i]]) {
            for(prod in jrnf_reactions$products[[i]]) {
                g <- add.edges(g, c(ed, prod))
            }
        } 
    }

    for(i in 1:nrow(jrnf_reactions)) {
        for(ed1 in 1:length(jrnf_reactions$educts[[i]])) {
            for(ed2 in 1:length(jrnf_reactions$educts[[i]])) {
                if(ed1 != ed2) {
                    g <- add.edges(g, c(jrnf_reactions$educts[[i]][ed1], 
                                        jrnf_reactions$educts[[i]][ed2]))
		}
            }
        } 
    }

    for(i in 1:nrow(jrnf_reactions)) {
        for(prod1 in 1:length(jrnf_reactions$products[[i]])) {
            for(prod2 in 1:length(jrnf_reactions$products[[i]])) {
                if(prod1 != prod2) {
                    g <- add.edges(g, c(jrnf_reactions$products[[i]][prod1], 
                                        jrnf_reactions$products[[i]][prod2]))
		}
            }
        } 
    }

    return(g)    
}


# Takes a jrnf-network and a boolean vector indicating which species to keep and 
# calculates the reduced network. The parameter 'rm_reaction' indicates whether 
# to remove every reaction that contains a removed species ('r') or if just the species
# should be removed from all reactions ('s'). If rm_reaction is neither 'r' or 's' reactions
# are removed if they have no educts or products but had them before.
# The 'list_changes' tells the function to write which species are removed.

jrnf_subnet <- function(jrnf_network, keep_flag, rm_reaction="r", list_changes=FALSE) {
    new_to_old <- which(keep_flag)
    old_to_new <- rep(0, nrow(jrnf_network[[1]]))
    for(i in 1:length(new_to_old)) old_to_new[new_to_old[i]]=i 

    # First subset the species and copy the reactions (remove and reindex later)
    new_species <- jrnf_network[[1]][keep_flag,]
    new_reactions <- jrnf_network[[2]]

    # If wanted the species that are deleted are described
    if(list_changes) {
        cat("Removing species: ")
        cat(as.vector((jrnf_network[[1]])$name)[!keep_flag], "\n") 
    }

    # Submethod which computes if the specific reaction is being kept
    calculate_keep <- function(reaction) {
        e <- as.vector(reaction$educts)
        em <- as.vector(reaction$educts_mul)
        p <- as.vector(reaction$products)
        pm <- as.vector(reaction$products_mul)

        e_f <- keep_flag[e]   # Logical vectors, which reaction are kept
        p_f <- keep_flag[p]   

        if(rm_reaction=="r") {
            if(sum(e_f) == length(e) && 
               sum(p_f) == length(p))
                return(TRUE)
            else
                return(FALSE)
        } else if(rm_reaction == "s") {
            if(sum(e_f) == 0 && sum(p_f) == 0)
                return(FALSE)
            else
                return(TRUE)
        } else {
            if(sum(e_f) == 0 && length(e) != 0 ||
               sum(p_f) == 0 && length(p) != 0)
                return(FALSE)
            else
                return(TRUE)
        }

        return(TRUE)
    }

    
    # Create a logic vector on which reactions to delete
    keep_r <- as.vector(apply(new_reactions, MARGIN=1, calculate_keep))
    new_reactions <- as.data.frame(new_reactions[keep_r,])

    if(list_changes) {
        cat("Removing reactions: ")
        cat(as.vector(which(!(keep_r))), "\n")
    }

    for(j in 1:nrow(new_reactions)) {


        e <- unlist(new_reactions[j,]$educts)
        em <- unlist(new_reactions[j,]$educts_mul)
        p <- unlist(new_reactions[j,]$products)
        pm <- unlist(new_reactions[j,]$products_mul)

        e_f <- keep_flag[e]
        p_f <- keep_flag[p]

        new_reactions$educts[j] <- list(old_to_new[e[e_f]])
        new_reactions$educts_mul[j] <- list(em[e_f])
        new_reactions$products[j] <- list(old_to_new[p[p_f]])
        new_reactions$products_mul[j] <- list(pm[p_f])
    }

    jrnf_network[[1]] <- new_species
    jrnf_network[[2]] <- new_reactions
    return(jrnf_network)
}



#
# Subnet function that removes all reactions indicated. After that it additionally 
# removes all species that are not occuring in any reaction...
#
#

jrnf_subnet_r <- function(jrnf_network, keep_flag) {
    net <- list(jrnf_network[[1]], jrnf_network[[2]][keep_flag,])
    degs <- degree(jrnf_to_undirected_network(net))
    net <- jrnf_subnet(net, degs != 0)   
    return(net)
}


# For a jrnf-reaction network object. This function calculates the associated
# substrate graph / network (using jrnf_to_directed_network), calculates
# the biggest strongly connected subgraph and returns a boolean vector indicating
# which nodes are in this subgraph.

jrnf_get_s_con_subnet <- function(jrnf_network) {
    directed_substrate_net <- jrnf_to_directed_network(jrnf_network)
    strong_clusters <- clusters(directed_substrate_net, mode="strong")
    keep_flag <- (strong_clusters$membership == which.max(strong_clusters$csize))

    return(keep_flag)
}


# TODO comment 
#

jrnf_get_w_con_subnet <- function(jrnf_network) {
    directed_substrate_net <- jrnf_to_directed_network(jrnf_network)
    strong_clusters <- clusters(directed_substrate_net, mode="weak")
    keep_flag <- (strong_clusters$membership == which.max(strong_clusters$csize))

    return(keep_flag)
}


# Simplifies a atmospheric reaction network in a common form. 
# N2, hv, CH4 and CO2 are kept even if they are not produced / consumed inside of the network.
# For hv a inflow reaction is added. The function returns the reduced subnetwork.
#
# If all reactions containing one species not belonging to the first strongly connected 
# subnet are removed this may lead again to a not strongly connected subset. If recursive
# is true this leads to applying the reduction until this does not happen. 

jrnf_simplify_AC_RN <- function(jrnf_network, recursive=TRUE, inflow=c("hv", "O2"), outflow=c("CO2"), 
                                keep=c("N2", "hv", "CH4", "CO2"), remove_M=TRUE) {
    id_external <- c()

    # Handling inflow species
    for(i in inflow) {
        x <- which(jrnf_network[[1]]$name == i)

        # adding to external list if necessary
        if(length(x) > 0 && length(which(id_external == x)) == 0)
            id_external <- c(id_external, x)

        # adding inflow reaction to network        
	    jrnf_network[[2]] <- rbind(jrnf_network[[2]], data.frame(reversible=factor(c(FALSE)), 
                                   c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                                   activation=as.numeric(c(0)),educts=I(list(c())), educts_mul=I(list(c())),
                                   products=I(list(x)), products_mul=I(list(1))))
    }

    # Handling outflow species
    for(i in outflow) {
        x <- which(jrnf_network[[1]]$name == i)

        # adding to external list if necessary
        if(length(x) > 0 && length(which(id_external == x)) == 0)
            id_external <- c(id_external, x)

        # adding outflow reaction to network
	    jrnf_network[[2]] <- rbind(jrnf_network[[2]], data.frame(reversible=factor(c(FALSE)), 
                                   c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                                   activation=as.numeric(c(0)),educts=I(list(c(x))), educts_mul=I(list(c(1))),
                                   products=I(list(c())), products_mul=I(list(c()))))
    }

    # Species that are kept (without adding inflow or outflow reaction)
    for(i in keep) {
        x <- which(jrnf_network[[1]]$name == i)
        if(length(x) > 0 && length(which(id_external == x)) == 0)
            id_external <- c(id_external, x)
    }

     
    kp <- jrnf_get_s_con_subnet(jrnf_network)
    kp[id_external] <- TRUE
    jrnf_sn <- jrnf_subnet(jrnf_network, kp, rm_reaction="r", list_changes=FALSE)

    # 
    if(sum(kp) < length(kp) && recursive)
        return(jrnf_simplify_AC_RN(jrnf_sn, TRUE, c(), c(), c(inflow, outflow, keep), remove_M))
    else
        if(remove_M)
            return(jrnf_subnet(jrnf_network, jrnf_network[[1]]$name != "M", rm_reaction="s", list_changes=FALSE))
        else
            return(jrnf_sn)
}



# new version of function above (should be faster through using apply...)
# TODO CHECK!

jrnf_calc_reaction_r <- function(network, kB_T=1) {
    calc_kkb <- function(x) {
        e <- unlist(x$educts)
        e_m <- unlist(x$educts_mul)
        p <- unlist(x$products)
        p_m <- unlist(x$products_mul)

        E_e <- sum(network[[1]]$energy[e]*e_m)
        E_p <- sum(network[[1]]$energy[p]*p_m)
        E_a <- max(E_e, E_p) + x$activation
        k <- exp(-(E_a-E_e)/kB_T)
  
        if(x$reversible)
            k_b <- exp(-(E_a-E_p)/kB_T)
        else
            k_b <- 0   

        return(data.frame(k=as.numeric(k), k_b=as.numeric(k_b)))
   }

   sp_list <- apply(network[[2]], 1, calc_kkb)

   sp <- do.call("rbind", sp_list)

   network[[2]]$k <- sp$k
   network[[2]]$k_b <- sp$k_b

   return(network)
}


# TODO document
#

jrnf_Ea_rel_to_abs <- function(network, Ea_rel) {
    calc_abs <- function(x) {
        e <- unlist(x$educts)
        e_m <- unlist(x$educts_mul)
        p <- unlist(x$products)
        p_m <- unlist(x$products_mul)

        E_e <- sum(network[[1]]$energy[e]*e_m)
        E_p <- sum(network[[1]]$energy[p]*p_m)
        E_a <- max(E_e, E_p) + x$activation_rel

        return(E_a)
   }

   network[[2]]$activation_rel <- Ea_rel
   sp <- as.numeric(apply(network[[2]], 1, calc_abs))
   return(sp)
}


jrnf_Ea_abs_to_rel <- function(network, Ea_abs) {
    calc_abs <- function(x) {
        e <- unlist(x$educts)
        e_m <- unlist(x$educts_mul)
        p <- unlist(x$products)
        p_m <- unlist(x$products_mul)

        E_e <- sum(network[[1]]$energy[e]*e_m)
        E_p <- sum(network[[1]]$energy[p]*p_m)
        E_a <- x$activation - max(E_e, E_p)

        return(E_a)
   }

   network[[2]]$activation <- Ea_abs
   sp <- as.numeric(apply(network[[2]], 1, calc_abs))
   return(sp)
}


#
#

jrnf_side_to_string <- function(r, r_m, net, f=function(x) {  x  }, sum_c=" + ") {
    a <- ""

    if(!is.null(r) && length(r) != 0)
        for(i in 1:length(r)) {
            if(r_m[i] != 1)
                a <- paste(a, as.numeric(r_m[i]), " ", sep="")

            a <- paste(a, f(net[[1]]$name[r[i]]), sep="")

            if(i < length(r))
                a <- paste(a, sum_c, sep="")
        }

    return(a)
}


#
#

jrnf_educts_to_string <- function(net, re_id, f=function(x) {  x  }, sum_c=" + ") {
    return(jrnf_side_to_string(net[[2]]$educts[re_id][[1]],
                               net[[2]]$educts_mul[re_id][[1]],
                               net, f, sum_c))
}


#
#

jrnf_products_to_string <- function(net, re_id, f=function(x) {  x  }, sum_c=" + ") {
    return(jrnf_side_to_string(net[[2]]$products[re_id][[1]],
                               net[[2]]$products_mul[re_id][[1]],
                               net, f, sum_c))
}


# Function creates a string representation of the reaction 're_id' from the 
# network 'net'

jrnf_reaction_to_string <- function(net, re_id, f=function(x) {  x  }, sum_c=" + ") {
    educts <- jrnf_educts_to_string(net, re_id, f, sum_c)
    products <- jrnf_products_to_string(net, re_id, f, sum_c)

    return(paste(educts, " => ", products, sep=""))
}



# Create a initial file for a certain jrnf network which can 
# then be used as a starting point for solving an ode. In default every 
# concentration is initialized with random values drawn from a gaussian
# distribution, added to one. 
# Additional 'bc_id' can be given a vector of species id names or a vector
# of (1-indexed) ids. If this is done, then 'bc_v' should be set to a 
# numeric vector of the same size indicating the respective concentrations.
#
# setBCE0  -  setting the boundary condition species energy to zero

jrnf_create_initial <- function(jrnf_network, init_file, network_file=NA, bc_id=NA, bc_v=NA, kB_T=1, setBCE0=TRUE) {
    jrnf_species <- jrnf_network[[1]]
    jrnf_reactions <- jrnf_network[[2]]    

    # ensure bc_id (if given) is numeric
    if(!is.na(bc_id) && !is.numeric(bc_id)) 
        for(i in 1:length(bc_id)) 
            bc_id[i] <- which(jrnf_species$name == bc_id[i])

    df <- data.frame(time=as.numeric(0),msd=as.numeric(0))
    df[as.vector(jrnf_species$name)] <- abs(rnorm(length(jrnf_species$name), mean=mean(bc_v), sd=sqrt(mean((bc_v- mean(bc_v))^2))))

    if(!is.na(bc_id) && !is.na(bc_v)) {
        if(length(bc_id) != length(bc_v)) {
            cat("Error in jrnf_create_initial; bc_id,bc_v length missmatch!")
            return() 
        }

        df[1,bc_id+2] <- bc_v 
    } 

    # and write 
    write.csv(df, init_file, row.names=FALSE)

    if(is.na(network_file)) 
        return()

    # Writing a network with the boundary species set constant
    if(!is.na(network_file) && !is.na(bc_id))
        jrnf_network[[1]]$constant[bc_id] <- TRUE;

    if(!is.na(network_file) && !is.na(bc_id) && setBCE0) {
        jrnf_network[[1]]$energy[bc_id] <- 0
        jrnf_network <- jrnf_calc_reaction_r(jrnf_network, kB_T)
    }

    write_jrnf(network_file, jrnf_network)
}


# This function calculates information on the network topology using igraph and 
# writes it to two files. 
# - One file (pfile) contains information specific to pairs of nodes
# shortest path, number of different shortest paths, ...
# - The second file (nfile) contains information specific to nodes, like degree,
# betweenness and whether the node belongs to the biggest cluster of the network.
# This two outputs can be enabled separately. The topological analysis is done
# on the directed graph one gets from 'jrnf_to_directed_network'
# TODO implement

jrnf_create_pnn_file <- function(jrnf_network, pfile=NA, nfile=NA) {
    N <- nrow(jrnf_network[[1]])
    g <- jrnf_to_undirected_network(jrnf_network)
    g_s <- simplify(g, remove.multiple=TRUE, remove.loops=FALSE)

    sps_g <- shortest.paths(g)

    if(!is.na(pfile)) {
        df_1 <- data.frame(from=as.numeric(rep(0, N*N)), to=as.numeric(rep(0, N*N)), 
                           shortest_path=as.numeric(rep(0, N*N)), sp_multiplicity=as.numeric(rep(0, N*N)), 
                           sp_multiplicity_s=as.numeric(rep(0, N*N)),
                           stringsAsFactors=FALSE)
        df_1$from <- floor((1:(N*N)-1)/N) + 1
        df_1$to <- (1:(N*N)-1) %% N + 1
            
        df_1$shortest_path <- as.vector(sps_g)
        #for(i in 1:(N*N)) {
        #    df_1$shortest_path[i] <- sps_g[floor((i-1)/N) + 1, (i-1) %% N + 1]
        #}
        df_1$sp_multiplicity <- 0
        df_1$sp_multiplicity_s <- 0


        for(k in 1:N) {
            if(k %% 10 == 1)
                cat(".")
            sp <- get.all.shortest.paths(g, from=k, mode="out")$res
            sp_s <- get.all.shortest.paths(g_s, from=k, mode="out")$res                     

            sp_x <- c()
            for(x in sp)
                sp_x <- c(sp_x, x[length(x)])

            adf <- as.data.frame(table(sp_x))
            df_1$sp_multiplicity[(k-1)*N+as.numeric(as.vector(adf$sp_x))] <- adf$Freq


            sp_s_x <- c()
            for(x in sp_s)
                sp_s_x <- c(sp_s_x, x[length(x)])

            adf <- as.data.frame(table(sp_s_x))
            df_1$sp_multiplicity_s[(k-1)*N+as.numeric(as.vector(adf$sp_s_x))] <- adf$Freq
        }
        cat("\n")

        write.csv(df_1, pfile, row.names=FALSE)
    }


    if(!is.na(nfile)) {
        df_2 <- data.frame(node=as.numeric(1:N), deg_in=as.numeric(rep(NA,N)), deg_out=as.numeric(rep(NA,N)), deg_all=as.numeric(rep(NA,N)),
                           deg_s_in=as.numeric(rep(NA,N)), deg_s_out=as.numeric(rep(NA,N)), deg_s_all=as.numeric(rep(NA,N)), betweenness=as.numeric(rep(NA,N)), betweenness_s=as.numeric(rep(NA,N)), main=as.logical(rep(NA,N)), main_w=as.logical(rep(NA,N)), stringsAsFactors=FALSE) 

        df_2$deg_in <- as.numeric(degree(g, mode="in"))
        df_2$deg_out <- as.numeric(degree(g, mode="out"))
        df_2$deg_all <- as.numeric(degree(g, mode="all"))
        df_2$deg_s_in <- as.numeric(degree(g_s, mode="in"))
        df_2$deg_s_out <- as.numeric(degree(g_s, mode="out"))
        df_2$deg_s_all <- as.numeric(degree(g_s, mode="all"))
        df_2$betweenness <- as.numeric(betweenness(g))
        df_2$betweenness_s <- as.numeric(betweenness(g_s))

        strong_clusters <- clusters(g, mode="strong")
        weak_clusters <- clusters(g, mode="weak")
        df_2$main <- (strong_clusters$membership == which.max(strong_clusters$csize))
        df_2$main_w <- (weak_clusters$membership == which.max(weak_clusters$csize))

        write.csv(df_2, nfile, row.names=FALSE)
    }

}



# a faster variant of pfile generation. Only calculates thes pairs that are given in 
# b_list structure...
# TODO

jrnf_create_pfile_bignet <- function(jrnf_network, b_list, pfile, calc_sp_mul=TRUE) {
    L <- length(b_list)
    g <- jrnf_to_undirected_network(jrnf_network)
    g_s <- simplify(g, remove.multiple=TRUE, remove.loops=FALSE)

    df_1 <- data.frame(from=as.numeric(rep(0, L)), to=as.numeric(rep(0, L)), 
                       shortest_path=as.numeric(rep(0, L)), sp_multiplicity=as.numeric(rep(0, L)), 
                       sp_multiplicity_s=as.numeric(rep(0, L)),
                       stringsAsFactors=FALSE)

    for(x in 1:length(b_list)) {
        from <- b_list[[x]][1]
        to <- b_list[[x]][2]
        df_1$from[x] <- from
        df_1$to[x] <- to
            
        cat(".")

        df_1$shortest_path[x] <- as.numeric(shortest.paths(g, from, to))
        df_1$sp_multiplicity[x] <- length(get.all.shortest.paths(g, from=from, to=to, mode="all")$res)        
        df_1$sp_multiplicity_s[x] <- length(get.all.shortest.paths(g_s, from=from, to=to, mode="all")$res)
    }        
    cat("\n")

    write.csv(df_1, pfile, row.names=FALSE)
}






#
# TODO implement + comment
# copy parameters from original reactions, even if it doesnt make sense

jrnf_copy_linearize <- function(infile, outfile, C) {
    net_in <- jrnf_read(infile)

    s_in <- net_in[[1]]
    r_in <- net_in[[2]]

    flag_11 <- logical()
    flag_22 <- logical()

    for(i in 1:nrow(r_in)) {
        flag_11 <- c(flag_11, sum(unlist(r_in$educts_mul[i])) == 1 & sum(unlist(r_in$products_mul[i])) == 1)
        flag_22 <- c(flag_22, sum(unlist(r_in$educts_mul[i])) == 2 & sum(unlist(r_in$products_mul[i])) == 2)
    }
    
    if(sum(!xor(flag_11, flag_22)) != 0) {
        cat("jrnf_copy_linearize - found non 1-1 / 2-2 reaction. Aborting!\n")
        return()
    }

    if(sum(flag_22) < C) {
        cat("Less 2-2 reactions than desired C!\n")
        return()
    }

    flag_keep <- logical(nrow(r_in)) 
    flag_keep[sample(which(flag_22), C)] <- TRUE

    r_out <- rbind(r_in[flag_keep,], r_in[flag_11,])

    # Transform and add linear version of substrate graph representation
    for(i in which(!flag_keep & !flag_11)) {
        educts <- unlist(r_in$educts[i])
        products <- unlist(r_in$products[i])
        ff <- educts[1]
        fs <- educts[length(educts)]
        sf <- products[1]
        ss <- products[length(products)]

        r1 <- data.frame(reversible=factor(TRUE), c=as.numeric(0), k=as.numeric(0),
                         k_b=as.numeric(0), activation=as.numeric(0), 
                         educts=as.numeric(ff), educts_mul=as.numeric(1),
                         products=as.numeric(sf), products_mul=as.numeric(1))

        r2 <- data.frame(reversible=factor(TRUE), c=as.numeric(0), k=as.numeric(0),
                         k_b=as.numeric(0), activation=as.numeric(0), 
                         educts=as.numeric(ff), educts_mul=as.numeric(1),
                         products=as.numeric(ss), products_mul=as.numeric(1))

        r3 <- data.frame(reversible=factor(TRUE), c=as.numeric(0), k=as.numeric(0),
                         k_b=as.numeric(0), activation=as.numeric(0), 
                         educts=as.numeric(fs), educts_mul=as.numeric(1),
                         products=as.numeric(sf), products_mul=as.numeric(1))

        r4 <- data.frame(reversible=factor(TRUE), c=as.numeric(0), k=as.numeric(0),
                         k_b=as.numeric(0), activation=as.numeric(0), 
                         educts=as.numeric(fs), educts_mul=as.numeric(1),
                         products=as.numeric(ss), products_mul=as.numeric(1))

        r_out <- rbind(r_out, r1)
        r_out <- rbind(r_out, r2)
        r_out <- rbind(r_out, r3)
        r_out <- rbind(r_out, r4)

    } 

    #return(list(s_in, r_out))

    net_out <-   
    jrnf_write(outfile, list(s_in, r_out))
}


# Function samples energies (activation energies and chemical standard potential)
# for the reaction network 'net' and calculates the assoziated rate constants.
# The function allows tu set all energies to zero instead of using the standard 
# distributions (zero = TRUE). Also if only 1-1 and 2-2 reactions are present
# the parameter v can be used to virtually increase the density of the system
# and increase the rates of the 2-2 reactions.

jrnf_sample_energies <- function(net, kB_T=1, v=1, zero=FALSE) {
    N <- nrow(net[[1]])
    M <- nrow(net[[2]])

    rea_is_1on1 <- (as.numeric(lapply(net[[2]]$educts, FUN=length)) == 1) &&
                   (as.numeric(lapply(net[[2]]$educts_mul, FUN=length)) == 1) &&
                   (as.numeric(lapply(net[[2]]$products, FUN=length)) == 1) &&
                   (as.numeric(lapply(net[[2]]$products_mul, FUN=length)) == 1) 


    # Draw numbers for species energy and activation energy   
    if(zero) {
        net[[1]]$energy <- rep(0,N)
        net[[2]]$activation <- rep(0,M)
    } else {
        net[[1]]$energy <- rnorm(N)
        net[[2]]$activation <- rplancklike(M)+log(v)*as.numeric(!rea_is_1on1)
    }

    return (jrnf_calc_reaction_r(net, kB_T))
}


#
#
#

jrnf_build_linear_stability_analyzer <- function(net, beta=1) {
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)
    N <- N_out-N_in
    k_f <- net[[2]]$k
    k_b <- net[[2]]$k_b
    mu0 <- net[[1]]$energy
    Ea <- pmax((mu0 %*% N_in)[1,], (mu0 %*% N_out)[1,]) + net[[2]]$activation         
    e_bs_in <- exp(beta*(mu0 %*% N_in)[1,]);
    e_bs_out <- exp(beta*(mu0 %*% N_out)[1,]);
    e_m_bEa <- exp(-beta*Ea)


    calculate <- function(x) {
        M <- matrix(0, ncol=nrow(N), nrow=nrow(N))
        for(i in 1:nrow(N)) 
            if(!net[[1]]$constant[i]) 
                for(l in 1:nrow(N)) {
                    for(k in 1:ncol(N)) {
                        p1 <- prod(x[-k]**N_in[-k,k])
                        p2 <- prod(x[-k]**N_in[-k,k])
                        
                        M[i,l] <- M[i,l] + N[i,k]*e_m_bEa[k]*(e_bs_in[k]*N_in[l,k]*(x[l]**(N_in[l,k]-1))*p1 - 
                                                              e_bs_out[k]*N_out[l,k]*(x[l]**(N_out[l,k]-1))*p2)
                    }

                    if(is.na(M[i,l]))
                        M[i,l] <- 0
                }
        return(M)
    }

    return(list(calculate=calculate, N_in, N_out, N, mu0, Ea,e_bs_in, e_bs_out, e_m_bEa))
}


# Function for merging two networks. Implementation works by adding the second
# network to the first one. Duplicates of species are recognized (by name) 
# and removed. Possible duplicates of reactions are not removed!

jrnf_merge_net <- function(net1, net2) {
    d <- duplicated(c(net1[[1]]$name, net2[[1]]$name))[-(1:nrow(net1[[1]]))]
    # create new species data
    net1[[1]] <- rbind(net1[[1]], net2[[1]][!d,])

    # for each species of the second network calculate the according id in the new network 
    new_id <- rep(0, nrow(net2[[1]]))

    for(i in 1:nrow(net2[[1]])) 
        new_id[i] <- which(net1[[1]]$name == net2[[1]]$name[i])

    # Transform indices for "added" reactions
    new_reactions <- net2[[2]]

    for(j in 1:nrow(new_reactions)) {
        e <- unlist(new_reactions[j,]$educts)
        p <- unlist(new_reactions[j,]$products)

        new_reactions$educts[j] <- list(new_id[e])
        new_reactions$products[j] <- list(new_id[p])
    }

    net1[[2]] <- rbind(net1[[2]], new_reactions)
    return(net1)
}



# Legacy references / can be removed if no longer needed

jrnf_simplify_X <- jrnf_simplify_multiplicity
jrnf_randomize <- jrnf_randomize_species
calculate_flow <- jrnf_calculate_flow
jrnf_graph_to_amatrix <- graph_to_amatrix

