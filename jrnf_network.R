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

jrnf_read <- function(filename) {
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


# backward compatibility

load_jrnf <- jrnf_read


# Writes reaction network data frames to file

jrnf_write <- function(filename, data) {
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


# Backward compatibility

write_jrnf <- jrnf_write


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

jrnf_simplify_AC_RN <- function(jrnf_network, recursive=TRUE, inflow=c("hv"), outflow=c(""), 
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


# TODO comment + test
#

jrnf_calc_reaction_r <- function(network, kB_T=1) {
    M <- nrow(network[[2]])

    for(i in 1:M) {
        e <- unlist(network[[2]]$educts[i])
        e_m <- unlist(network[[2]]$educts_mul[i])
        p <- unlist(network[[2]]$products[i])
        p_m <- unlist(network[[2]]$products_mul[i])

        E_e <- sum(network[[1]]$energy[e]*e_m)   # Energy of educts
        E_p <- sum(network[[1]]$energy[p]*p_m)   # Energy of products
        E_a <- max(E_e, E_p) + network[[2]]$activation[i]        # (absolute) activation energy

        network[[2]]$k[i] <- exp(-(E_a-E_e)/kB_T)

        if(network[[2]]$reversible[i]) 
            network[[2]]$k_b[i] <- exp(-(E_a-E_p)/kB_T)
        else 
            network[[2]]$k_b[i] <- 0    
    }

    return(network)
}



# Create a initial file for a certain jrnf network which can 
# then be used as a starting point for solving an ode. In default every 
# concentration is initialized with random values drawn from a gaussian
# distribution, added to one. 
# Additional 'bc_id' can be given a vector of species id names or a vector
# of (1-indexed) ids. If this is done, then 'bc_v' should be set to a 
# numeric vector of the same size indicating the respective concentrations.
#
#

jrnf_create_initial <- function(jrnf_network, init_file, network_file=NA, bc_id=NA, bc_v=NA, kB_T=1, setBCE0=FALSE) {
    jrnf_species <- jrnf_network[[1]]
    jrnf_reactions <- jrnf_network[[2]]    

    # ensure bc_id (if given) is numeric
    if(!is.na(bc_id) && !is.numeric(bc_id)) 
        for(i in 1:length(bc_id)) 
            bc_id[i] <- which(jrnf_species$name == bc_id[i])

    df <- data.frame(time=as.numeric(0))
    df[as.vector(jrnf_species$name)] <- 0

    df[1,] <- abs(1+rnorm(length(jrnf_species$name)+1))
    df[1,1] <- 0  # set time zero again

    if(!is.na(bc_id) && !is.na(bc_v)) {
        if(length(bc_id) != length(bc_v)) {
            cat("Error in jrnf_create_initial; bc_id,bc_v length missmatch!")
            return() 
        }

        df[1,bc_id+1] <- bc_v 
    } 

    # and write 
    write.csv(df, init_file, row.names=FALSE)

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
    g <- jrnf_to_directed_network(jrnf_network)
    g_s <- simplify(g, remove.multiple=TRUE, remove.loops=FALSE)

    sps_g <- shortest.paths(g)


    if(!is.na(pfile)) {
        df_1 <- data.frame(from=as.numeric(rep(NA, N*N)), to=as.numeric(rep(NA, N*N)), 
                           shortest_path=as.numeric(rep(NA, N*N)), sp_multiplicity=as.numeric(rep(0, N*N)), 
                           sp_multiplicity_s=as.numeric(rep(0, N*N)),
                           stringsAsFactors=FALSE)

        for(i in 1:(N*N)) {
            fr <- floor((i-1)/N) + 1
            t <- (i-1) %% N + 1

            df_1$from[i] <- fr
            df_1$to[i] <- t
            df_1$shortest_path[i] <- sps_g[fr, t]
            df_1$sp_multiplicity[i] <- 0
            df_1$sp_multiplicity_s[i] <- 0
        }


        for(k in 1:N) {
            sp <- get.all.shortest.paths(g, from=k, mode="out")$res
            sp_s <- get.all.shortest.paths(g_s, from=k, mode="out")$res                     

            for(x in sp) {
                sel <- which(df_1$from == k & df_1$to == x[length(x)])
                df_1$sp_multiplicity[sel] <- df_1$sp_multiplicity[sel] + 1
            }

            for(x in sp_s) {
                sel <- which(df_1$from == k & df_1$to == x[length(x)])
                df_1$sp_multiplicity_s[sel] <- df_1$sp_multiplicity_s[sel] + 1 
            }
        }

        write.csv(df_1, pfile, row.names=FALSE)
    }


    if(!is.na(nfile)) {
        df_2 <- data.frame(node=as.numeric(1:N), deg_in=as.numeric(rep(NA,N)), deg_out=as.numeric(rep(NA,N)), deg_all=as.numeric(rep(NA,N)),
                           deg_s_in=as.numeric(rep(NA,N)), deg_s_out=as.numeric(rep(NA,N)), deg_s_both=as.numeric(rep(NA,N)), betweenness=as.numeric(rep(NA,N)), betweenness_s=as.numeric(rep(NA,N)), main=as.logical(rep(NA,N)), main_w=as.logical(rep(NA,N)), stringsAsFactors=FALSE) 

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







# Function that samples species energies from normal distribution and 
# activation energies from a distribution of the form p(E) = 6/(pi^2*(exp(1/x)-1))

rplancklike <- function(N) {
    r <- c()    

    for(i in 1:N) {
        x <- runif(min=0, max=100,n=1)
        y <- (6/pi^2)/(x^3*(exp(1/x)-1))

        while(y < runif(min=0, max=1, n=1)) {
            x <- runif(min=0, max=100,n=1)
            y <- (6/pi^2)/(x^3*(exp(1/x)-1))
        }
        
        r <- c(r, x)
    }

    return(r)
}


jrnf_sample_energies <- function(jrnf_network, kB_T=1) {
    N <- nrow(jrnf_network[[1]])
    M <- nrow(jrnf_network[[2]])

    # Draw numbers for species energy and activation energy   
    jrnf_network[[1]]$energy <- rnorm(N)
    jrnf_network[[2]]$activation <- rplancklike(M)

    return (jrnf_calc_reaction_r(jrnf_network, kB_T))
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

