# author: jakob fischer (jakob@automorph.info)
# date: 23. January 2013
# description: 
# R-tools for handling jrnf-network data and using parts of igraph
# for network analysis and visualisation

library(igraph)


# Loads a jrnf-file and returns a list of two data frames
# - The first data frame contains all species with 'type' (integer),
# energy (numeric), 'name' (character) and 'constant' (boolean).
# - The second data frame contains all reactions with information whether
# it is 'reversible' (boolean), the reaction constants 'c', 'k', 'k_b', 
# 'activation' energy (numeric), lists on 'products' and 'educts' and 
# their multiplicity in the reaction ('educts_mul' and 'products_mul')

load_jrnf <- function(filename) {
    # Open file and verify the right header is there
    con <- file(filename, 'r');
    if(readLines(con, 1) != "jrnf0003") {
        cat("File ", filename, "is not in valid jrnf format!")
        return(0)
    }    

    # Declaring empty data frames for species and reactions
    species <- data.frame(type=factor(integer()),
                 name=character(), 
                 energy=integer(),
                 constant=factor(logical()),
                 stringsAsFactors=FALSE)
    reactions <- data.frame(reversible=factor(logical()),
                            c=numeric(), k=numeric(),
                            k_b=numeric(), activation=numeric(), 
                            educts=list(), educts_mul=list(),
                            products=list(), products_mul=list())

    # Read line using ' ' to separate and put parts in 'last_line'-vector
    last_line <- scan(textConnection(readLines(con,1)), sep=" ", quiet="TRUE")
    if(length(last_line) != 2) {
        cat("Error at reading header of ", filename, "!")
        return(0)
    }  

    no_species <- last_line[1]
    no_reactions <- last_line[2]

    # Input species data
    for(i in 1:no_species) {
        last_line <- strsplit(readLines(con,1), " ", fixed = TRUE)[[1]]
        if(length(last_line) != 4) {
            cat("Error at reading species", i, "!")
            return(last_line)
        }

	species <- rbind(species, data.frame(type=factor(as.integer(last_line[1])), name=last_line[2], 
                                             energy=as.numeric(last_line[3]), 
                                             constant=as.logical(as.integer(last_line[4]))))
    }

    # Input data on reactions
    for(i in 1:no_reactions) {
        # read entire line and separate into last_line
        last_line <- strsplit(readLines(con,1), " ", fixed = TRUE)[[1]]
        if(length(last_line) < 7) {
            cat("Error at reading reaction", i, "!")
            return(0)
        }

        e_number <- as.numeric(last_line[6])
        e <- c()
        em <- c()
        p_number <- as.numeric(last_line[7])
        p <- c()
        pm <- c()

        # Fill educts ... 
        for(j in 1:e_number) {
            e <- c(e, as.numeric(last_line[7+j*2-1])+1)
            em <- c(em, as.numeric(last_line[7+j*2])) 
	}

        # ...and product vectors
        for(j in 1:p_number) {
            p <- c(p, as.numeric(last_line[7+e_number*2+j*2-1])+1)
            pm <- c(pm, as.numeric(last_line[7+e_number*2+j*2])) 
	}

        # Add it all togother to reactions data frame
	reactions <- rbind(reactions, data.frame(reversible=factor(as.logical(as.integer(last_line[1]))), 
                                                 c=as.numeric(last_line[2]), k=as.numeric(last_line[3]),
                                                 k_b=as.numeric(last_line[4]), activation=as.numeric(last_line[5]),
                                                 educts=I(list(e)), educts_mul=I(list(em)),
                                                 products=I(list(p)), products_mul=I(list(pm))))
    }

    # close connection and concatenate it to jrnf-object
    close(con)
    return(list(species, reactions))
}


# Writes reaction network data frames to file

write_jrnf <- function(filename, data) {
    # open file and write header
    con <- file(filename, 'w');
    species <- data[[1]]
    reactions <- data[[2]]
    writeLines("jrnf0003", con=con)
    writeLines(paste(as.character(nrow(species)), "\ ", 
               as.character(nrow(reactions)), sep=""), con=con)

    # write species
    for(i in 1:nrow(species)) {
        if(species$constant[i])  # constant?
            sc <- "\ 1"
        else
            sc <- "\ 0"

        # species type, name and energy
        writeLines(paste(as.character(as.vector(species$type)[i]),
                         "\ ", as.character(species$name[i]),"\ ",
                         as.character(species$energy[i]), sc, sep=""), con=con)
    }

    # write reactions
    for(i in 1:nrow(reactions)) {
        # Rename (for convenience)
        educts <- unlist(reactions$educts[[i]])
        educts_mul <- unlist(reactions$educts_mul[[i]])
        products <- unlist(reactions$products[[i]])
        products_mul <- unlist(reactions$products_mul[[i]])

        # reversible?
        if(as.vector(reactions$reversible)[i])
            line <- "1\ "
        else 
            line <- "0\ "

        # Reaction constants and number of educts and products
        line <- paste(line, as.character(as.numeric(reactions$c[i])), " ",
                      as.character(as.numeric(reactions$k[i])), " ",
        	      as.character(as.numeric(reactions$k_b[i])), " ",
        	      as.character(as.numeric(reactions$activation[i])), 
                      " ", as.character(length(educts)), 
                      " ",as.character(length(products)), sep="")

        # Add educts and products
        if(length(educts) != 0)
	    for(j in 1:length(educts)) 
                line <- paste(line, " ", as.character(educts[j]-1), 
                              " ", as.character(educts_mul[j]), sep="")

        if(length(products) != 0)
            for(j in 1:length(products)) 
                line <- paste(line, " ", as.character(products[j]-1), 
                              " ", as.character(products_mul[j]), sep="")
	
        # write reaction to file
	writeLines(line, con=con)

    }
    
    # close and good bye
    close(con)
}


# Calculates the koefficients k and k_b for the networks from
# energies of educts and products and activation energy of the
# reaction.

calculate_rev_koef <- function (network, kT=1.0) {

    for(i in 1:nrow(network[[2]])) {
        e_educts <- 0.0
        e_products <- 0.0
        e_activation <- network[[2]]$activation[i]  

        # Iterate educts
        for(j in 1:length(network[[2]]$educts[[i]])) {
            id <- network[[2]]$educts[[i]][j]
            mul <- network[[2]]$educts_mul[[i]][j]
            e_educts <- e_educts + network[[1]]$energy[id]*mul
        }

        # Iterate products
        for(j in 1:length(network[[2]]$products[[i]])) {
            id <- network[[2]]$products[[i]][j]
            mul <- network[[2]]$products_mul[[i]][j]
            e_products <- e_products + network[[1]]$energy[id]*mul
	}

        network[[2]]$k[i] <- exp((e_educts-e_activation)/kT)
        network[[2]]$k_b[i] <- exp((e_products-e_activation)/kT)
    }

    return(network)
}


# Calculates the flow / reaction rates for a given concentration vector
# TODO: calculate energy dif, check

calculate_flow <- function(network, concentrations) {
    result <- data.frame(flow_effective=numeric(), flow_forward=numeric(), flow_backward=numeric(),
                         energy_flow=numeric(), entropy_prod=numeric())

    for(i in 1:nrow(network[[2]])) {
        f_forward <- network[[2]]$k[i]
        f_backward <- network[[2]]$k_b[i]
        energy_dif <- 0 
 

        # Iterate educts
        for(j in 1:length(network[[2]]$educts[[i]])) {
            id <- network[[2]]$educts[[i]][j]
            mul <- network[[2]]$educts_mul[[i]][j]
            energy_dif <- energy_dif - network[[1]]$energy[id]*mul
            f_forward <- f_forward * concentrations[id]^mul
        }

        # Iterate products
        for(j in 1:length(network[[2]]$products[[i]])) {
            id <- network[[2]]$products[[i]][j]
            mul <- network[[2]]$products_mul[[i]][j]
            energy_dif <- energy_dif + network[[1]]$energy[id]*mul
            f_backward <- f_backward * concentrations[id]^mul
	}

        f_effective <- f_forward - f_backward
        energy_f <- energy_dif*f_effective
        entrop_p <- (f_forward-f_backward)*log(f_forward/f_backward)

        #cat(f_effective, ", ", f_backward, ", ", f_forward, "\n")
        result <- rbind(result, 
                        data.frame(flow_effective=as.numeric(f_effective), 
                                   flow_forward=as.numeric(f_forward), 
                                   flow_backward=as.numeric(f_backward), 
                                   energy_flow=as.numeric(energy_f), 
                                   entropy_prod=as.numeric(entrop_p)))
    }

    return(result)
}


# This function associates a propertie that is given for reactions to
# the species for a jrnf-network. This is done by interpolating on
# all reactions that are using up or producing the specific species

associate_reaction_to_species <- function(network, sp_prop) {
    test <- (1:nrow(network[[1]]))*0

    for(i in 1:nrow(network[[2]])) {
        vec <- c()
        mul <- c()

        # Iterate educts
        for(j in 1:length(network[[2]]$educts[[i]])) {
            vec <- c(vec, network[[2]]$educts[[i]][j])
            mul <- c(mul, network[[2]]$educts_mul[[i]][j])
        }

        # Iterate products
        for(j in 1:length(network[[2]]$products[[i]])) {
            vec <- c(vec, network[[2]]$products[[i]][j])
            mul <- c(mul, network[[2]]$products_mul[[i]][j])
	}

        factor <- 1.0/sum(mul)
        for(j in 1:length(vec)) {
            test[vec[j]] <- test[vec[j]] + sp_prop[i]*factor*mul[j]
        }
    }

    return(test)   
}


# This function transforms a jrnf-network to a directed network by connecting
# two species / nodes if there is a reaction with the first as educt and the
# second as product

jrnf_to_directed_network <- function(jrnf_data) {
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

    V(g)$label <- V(g)$name

    return(g)    
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

    return(list(new_species, new_reactions))
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


# Simplifies a atmospheric reaction network in a common form. 
# N2, hv, CH4 and CO2 are kept even if they are not produced / consumed inside of the network.
# The function returns the reduced subnetwork.
# TODO: Add parameter that allows creation of addition reactions for some species (hv, CH4)
# CAUTION: the pseudo species "M" cannot be removed as the other by deleting all reactions containing
# it without changing the topologic structure signigicantly
# TODO: Removing reactions related with removed species may lead to the result again being not
# strongly connected?

jrnf_simplify_AC_RN <- function(jrnf_network) {
    id_hv <- which(jrnf_network[[1]]$name == "hv")
    id_M <- which(jrnf_network[[1]]$name == "M")
    id_CH4 <- which(jrnf_network[[1]]$name == "CH4")
    id_CO2 <- which(jrnf_network[[1]]$name == "CO2")
    id_N2 <- which(jrnf_network[[1]]$name == "N2")
    id_O2 <- which(jrnf_network[[1]]$name == "O2")

    id_externals <- c(id_hv, id_CO2, id_N2, id_O2)

    keep <- jrnf_get_s_con_subnet(jrnf_network)
    keep[id_externals] <- TRUE

    return(jrnf_subnet(jrnf_network, keep, rm_reaction="r", list_changes=FALSE))
}



# After drawing a histogram with the rainbow(20) palette 
# this function returns a vector associating the same colors
# used in the histogramms to the values given.
#
# 'values' vector of values
# 'hi' object returned from the hist-call

get_color_by_hist <- function(hi, values) {
    d <- character()

    for(i in values) {
        if(is.na(i)) {
            d <- c(d, "#FFFFFFFF")
	} else {
            lisa <- max(which(i > hi$breaks))
            d <- c(d, rainbow(20)[lisa])
	}
    }

    return(d)
}


# Function takes a vector of chemical species names and a data frame of
# entropy-enthaly information and return a vector associating each species
# a formation enthalpy (kJ/mol) - or NA if no information is available...

get_fenthalpy_species <- function(names, ent) {
    res <- numeric()

    for(i in names) {
        t <- which(i == ent$name)
	    if(length(t) == 1) {
            res <- c(res,ent$dHf[t[1]]) 
	    } else {
            res <- c(res, NA)
	    }
    }

    return(res)
}


# Function takes a vector of chemical species names and a data frame of
# entropy-enthaly information and return a vector associating each species
# a standard entropy (J/(mol K)) - or NA if no information is available...

get_sentropy_species <- function(names, ent) {
    res <- numeric()

    for(i in names) {
        t <- which(i == ent$name)
	if(length(t) == 1) {
            res <- c(res,ent$S0[t[1]]) 
	} else {
            res <- c(res, NA)
	}
    }

    return(res)
}

