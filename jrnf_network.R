# author: jakob fischer (jakob@automorph.info)
# date: 23. January 2013
# description: 
# R-tools for handling jrnf-network data and using parts of igraph
# for network analysis and visualisation

library(igraph)


source("cycles.R")


# Converts an igraph graph to an adjacency matrix...
jrnf_graph_to_amatrix <- function(g) {
    N <- length(V(g))
    M <- length(E(g))

    x <- matrix(0, ncol=N, nrow=N)
    for(i in 1:M) {
        a <- get.edge(g, i)[1]
        b <- get.edge(g, i)[2]

        x[a,b] <- x[a,b] + 1
    }  

    return(x)
}




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
    lines <- readLines(con)
    close(con)


    if(lines[1] != "jrnf0003") {
        cat("File ", filename, "is not in valid jrnf format!")
        return(0)
    }    

    # Read line using ' ' to separate and put parts in 'last_line'-vector
    last_line <- scan(textConnection(lines[2]), sep=" ", quiet="TRUE")
    if(length(last_line) != 2) {
        cat("Error at reading header of ", filename, "!")
        return(0)
    }  

    no_species <- last_line[1]
    no_reactions <- last_line[2]


    # Declaring apply methods 
    sp_a <- function(sp_string) {
        last_line <- strsplit(sp_string, " ", fixed = TRUE)[[1]]
        if(length(last_line) != 4) 
            cat("Error at reading species!")
        
        return(data.frame(type=as.integer(last_line[1]), name=as.character(last_line[2]), 
                          energy=as.numeric(last_line[3]),
                          constant=as.logical(as.logical(as.integer(last_line[4]))),
                          stringsAsFactors=FALSE))
    }

    re_a <- function(re_string) {
        last_line <- strsplit(re_string, " ", fixed = TRUE)[[1]]
        if(length(last_line) < 7) {
            cat("Error at reading reaction:\n")
            cat(last_line, "\n")
        }

        e_number <- as.numeric(last_line[6])
        e <- as.numeric((rep(0, e_number)))
        em <- as.numeric((rep(0, e_number)))
        p_number <- as.numeric(last_line[7])
        p <- as.numeric((rep(0, p_number)))
        pm <- as.numeric((rep(0, p_number)))

        # Fill educts ... 
        for(j in 1:e_number) {
            e[j] <- as.numeric(last_line[7+j*2-1])+1
            em[j] <- as.numeric(last_line[7+j*2]) 
	}

        # ...and product vectors
        for(j in 1:p_number) {
            p[j] <- as.numeric(last_line[7+e_number*2+j*2-1])+1
            pm[j] <- as.numeric(last_line[7+e_number*2+j*2]) 
	}

        if(is.na(e[1])) {
            e <- c()
            em <- c()
        }

        if(is.na(p[1])) {
            p <- c()
            pm <- c()
        }

        return(data.frame(reversible=as.logical(as.integer(last_line[1])),
                          c=as.numeric(last_line[2]), 
                          k=as.numeric(last_line[3]),
                          k_b=as.numeric(last_line[4]), 
                          activation=as.numeric(last_line[5]), 
                          educts=I(list(e)), educts_mul=I(list(em)),
                          products=I(list(p)), products_mul=I(list(pm))))

    }

    # Input species data
    sp_list <- lapply(lines[3:(2+no_species)], sp_a)
    species <- do.call("rbind", sp_list)

    # Input data on reactions
    re_list <- lapply(lines[(3+no_species):(2+no_species+no_reactions)], re_a)

    reactions <- do.call("rbind", re_list)


    # concatenate it to jrnf-object
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


# write the topologic part of the network structure to rea-format
jrnf_write_to_rea <- function(filename, data) {
    # open file and write header
    con <- file(filename, 'w');
    species <- data[[1]]
    reactions <- data[[2]]
    writeLines("# generated by jrnf_write_to_rea", con=con)

    writeLines("# Number of Components", con=con)
    writeLines( as.character(nrow(species)), con=con)
    
    writeLines("# Components", con=con)
    for(i in 1:nrow(species)) 
        writeLines(as.character(species$name[i]) , con=con)
    

    writeLines("# Number of Reactions", con=con)
    writeLines( as.character(nrow(reactions)), con=con)

    writeLines("# Reactions", con=con)
    for(i in 1:nrow(reactions)) {
        # Rename (for convenience)
        educts <- unlist(reactions$educts[[i]])
        educts_mul <- unlist(reactions$educts_mul[[i]])
        products <- unlist(reactions$products[[i]])
        products_mul <- unlist(reactions$products_mul[[i]])

        line <- ""

        # Add educts and products
        if(length(educts) != 0) {
	    for(j in 1:length(educts)) 
                line <- paste(line, as.character(educts_mul[j]), " ", 
                              as.character(species$name[educts[j]]), " ", sep="")
        }

        line <- paste(line, "->", sep="")

        if(length(products) != 0){
	    for(j in 1:length(products)) 
                line <- paste(line, " ", as.character(products_mul[j]), " ", 
                              as.character(species$name[products[j]]), sep="")
        }
	
        # write reaction to file
	writeLines(line, con=con)
    }

    close(con)
}


# Calculates the 'in' stoichiometric matrix of the jrnf-network net. 
# 'in' means every column contains the left side of one reaction. 
# All values are positive.

jrnf_calculate_stoich_mat_in <- function(net) {
    no_sp <- nrow(net[[1]])
    no_re <- nrow(net[[2]])
    N <- matrix(0, no_sp, no_re)
    
    for(i in 1:no_re) {
        for(j in 1:length(net[[2]]$educts[[i]])) {
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
        for(j in 1:length(net[[2]]$products[[i]])) {
            p <- net[[2]]$products[[i]][j]
            N[p, i] <- N[p, i] + net[[2]]$products_mul[[i]][j]
        }
    }

    return(N)
}


# Calculates the stoichiometric matrix for a jrnf-network
# (matrix contains information on net change of species with reactions)

jrnf_calculate_stoich_mat <- function(net) {
    return(jrnf_calculate_stoich_mat_out(net)-jrnf_calculate_stoich_mat_in(net))
}



# Calculate the rate of concentration change assuming the reactions are
# happening at the rates

jrnf_calculate_concentration_change <- function(network, rates) {
    N <- jrnf_calculate_stoich_mat(network)
    return(N %*% rates)
}


# Calculates the flow / reaction rates for a given concentration vector
# TODO: calculate energy dif, check

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
            entrop_p <- (f_forward-f_backward)*log(abs(f_forward/f_backward))    # TODO check if abs is right
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


# Backward compatibility
calculate_flow <- jrnf_calculate_flow




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
            e <- net[[2]]$educts[[i]]
            e_m <- net[[2]]$educts_mul[[i]]
            p <- net[[2]]$products[[i]]
            p_m <- net[[2]]$products_mul[[i]]

            net[[2]]$educts[[i]] <- p
            net[[2]]$educts_mul[[i]] <- p_m
            net[[2]]$products[[i]] <- e
            net[[2]]$products_mul[[i]] <- e_m
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


# Function returns a list of 2-element vectors of reactions that are the 
# reverse of each other 
 
jrnf_find_reverse_pairs <- function(net) {
    if(nrow(net[[2]]) < 2)
        return(list())

    p <- list()
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)

    for(i in 1:(nrow(net[[2]])-1)) {
        for(j in (i+1):(nrow(net[[2]]))) {
            if(all(N_in[,i] == N_out[,j]) & all(N_in[,j] == N_out[,i]))
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



# This function associates a propertie that is given for reactions to
# the species for a jrnf-network. This is done by interpolating on
# all reactions that are using up or producing the specific species

#associate_reaction_to_species <- function(network, sp_prop) {
#    test <- (1:nrow(network[[1]]))*0
#
#    for(i in 1:nrow(network[[2]])) {
#        vec <- c()
#        mul <- c()
#
#        # Iterate educts
#        for(j in 1:length(network[[2]]$educts[[i]])) {
#            vec <- c(vec, network[[2]]$educts[[i]][j])
#            mul <- c(mul, network[[2]]$educts_mul[[i]][j])
#        }
#
#        # Iterate products
#        for(j in 1:length(network[[2]]$products[[i]])) {
#            vec <- c(vec, network[[2]]$products[[i]][j])
#            mul <- c(mul, network[[2]]$products_mul[[i]][j])
#	}
#
#        factor <- 1.0/sum(mul)
#        for(j in 1:length(vec)) {
#            test[vec[j]] <- test[vec[j]] + sp_prop[i]*factor*mul[j]
#        }
#    }
#
#    return(test)   
#}


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
        for(ed in jrnf_reactions$educts[[i]]) {
            for(prod in jrnf_reactions$products[[i]]) {
                if(rnd) {
                    b <- runif(1) > 0.5
                    if(b)
                        g <- add.edges(g, c(ed, prod))
                    else 
                        g <- add.edges(g, c(prod, ed))
                } else {
                    g <- add.edges(g, c(ed, prod))
                }
            }
        } 
    }

    V(g)$label <- V(g)$name

    return(g)    
}


#
#  commented out because content is same than jrnf_to_directed_network_d?
#

#jrnf_to_directed_network_d_mul <- function(jrnf_data, direction) {
#    jrnf_species <- jrnf_data[[1]]
#    jrnf_reactions <- jrnf_data[[2]]    
#    g <- graph.empty()
#    g <- add.vertices(g, nrow(jrnf_species), 
#                          name=as.vector(jrnf_species$name))   
#
#    for(i in 1:nrow(jrnf_reactions)) {
#        for(ed in jrnf_reactions$educts[[i]]) {
#            for(prod in jrnf_reactions$products[[i]]) {
#                if(direction[i] > 0)
#                    g <- add.edges(g, c(ed, prod))
#                else
#                    g <- add.edges(g, c(prod, ed))
#            }
#        } 
#    }
#
#    V(g)$label <- V(g)$name
#
#    return(g)    
#}



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

    return(list(new_species, new_reactions))
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


# TODO comment + test
#

jrnf_calc_reaction_r_OLD <- function(network, kB_T=1) {
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


#
#
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



# TODO comment + test
#
jrnf_reaction_to_string <- function(jrnf_network, reaction_id) {
    return("TODO")
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

    #cat("HALLO bc_id=", bc_id, " bc_v=", bc_v, "\n")

    # ensure bc_id (if given) is numeric
    if(!is.na(bc_id) && !is.numeric(bc_id)) 
        for(i in 1:length(bc_id)) 
            bc_id[i] <- which(jrnf_species$name == bc_id[i])

    #cat("TAT\n")

    df <- data.frame(time=as.numeric(0),msd=as.numeric(0))
    df[as.vector(jrnf_species$name)] <- abs(rnorm(length(jrnf_species$name), mean=mean(bc_v), sd=sqrt(mean((bc_v- mean(bc_v))^2))))


    #cat("BIB\n")

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


jrnf_sample_energies <- function(jrnf_network, kB_T=1, v=1, zero=FALSE) {
    N <- nrow(jrnf_network[[1]])
    M <- nrow(jrnf_network[[2]])

    rea_is_1on1 <- (as.numeric(lapply(jrnf_network[[2]]$educts, FUN=length)) == 1) &&
                   (as.numeric(lapply(jrnf_network[[2]]$educts_mul, FUN=length)) == 1) &&
                   (as.numeric(lapply(jrnf_network[[2]]$products, FUN=length)) == 1) &&
                   (as.numeric(lapply(jrnf_network[[2]]$products_mul, FUN=length)) == 1) 


    # Draw numbers for species energy and activation energy   
    if(zero) {
        jrnf_network[[1]]$energy <- rep(0,N)
        jrnf_network[[2]]$activation <- rep(0,M)
    } else {
        jrnf_network[[1]]$energy <- rnorm(N)
        jrnf_network[[2]]$activation <- rplancklike(M)+log(v)*as.numeric(!rea_is_1on1)
    }



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



# Function creates artificial ecosystem with <N> species... 
#

# Helper to create names consistent with elementary constituents
hcae_create_name <- function(comp, c_names, ex_names, empty_name) {
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

jrnf_create_artificial_ecosystem <- function(N, M, comp_no, mod_no=0, mod_f=1, no_reordering=F, no_hv=2, f_2fold=1/3) { 
    no_2fold = as.integer( M*f_2fold )

    # TODO check parameters   (N/mod_no has to be a natural number)!
    sp_mod_id <- floor((1:N-1)/mod_no)    # module of each species

    eval_possible_reas <- function(s, mod=F) { 
        #cat("call to eval_possible_reas with s=", s, "\n")
        possible_reas <- rep(1, s^4) 
        has_hv <- rep(F, s^4)
        rea_no <- rep(0, s^4)
        transfer <- rep(0, s^4)

        # TODO: speed this up with apply?
        for(i in 1:length(possible_reas)) {
            j <- i-1
            a <- j%%s
            b <- ((j-a)/s)%%s
            c <- ((j-a-b*s)/s^2)%%s
            d <- ((j-a-b*s-c*s^2)/s^3)%%s 

            #cat("a=", a, " b=", b, " c=", c, " d=", d, "\n")

            # Increase probability of reactions "in" modular structure
            if(mod & a != 0 && c != 0 && a != N+1 && c != N+1 && sp_mod_id[a] == sp_mod_id[c])
                possible_reas[j] <- possible_reas[j]/mod_f

            if(mod & b != 0 && d != 0 && b != N+1 && d != N+1 &&  sp_mod_id[b] == sp_mod_id[d])
                possible_reas[j] <- possible_reas[j]/mod_f

            # Decrease probability of autocatalytic self loops
            if(mod & a != 0 && b != 0 && c != 0 && d != 0 && (a == c || a == d || b == c || b == d))
                possible_reas[j] <- possible_reas[j]/N

            # Check if reaction works by elementary constituents
            if(!hcae_check_rea_constituents(c(a,b,c,d), composition, N) ||
               !hcae_check_rea_conditions(c(a,b,c,d), N))
                possible_reas[i] <- 0

            # If reaction is valid, calculate transfer
            if(possible_reas[i] != 0)
                transfer[i] <- hcae_calc_rea_transfer(c(a,b,c,d), composition, N)

            if(a == N+1 || b == N+1 || c == N+1 || d == N+1)
                has_hv[i] <- T

            # CAREFUL: reactants occuring on both sides are not counted (have to be subtracted)
            rea_no[i] <- sum(c(a != 0 , b != 0, c != 0, d != 0)) - 2*sum(c(a == c & a != 0, a == d & a != 0, b == c & b != 0, b == d & b != 0))
        }

        possible_reas[rea_no == 0] <- 0
        possible_reas[rea_no == 1] <- 0

        return(list(possible_reas, has_hv, rea_no, transfer))
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



    empty_name <- "X"
    hv_name <- "hv"
    component_names <- c("C", "N", "O", "H", "P")
    component_names <- component_names[1:comp_no]
    
    # draw compostition (el. constituents) and energy of species
    composition <- matrix(rpois(N*comp_no, 0.8), ncol=comp_no)

    cat("Drawing elementary composition")
    while(sum(duplicated(composition)) > 2) {
        cat(".")
        composition <- matrix(rpois(N*comp_no, 1.2), ncol=comp_no)
    }
    cat("\n")

    energy <- rnorm(N)   

    while(any(apply(composition, 1, sum) == 0)) {
        for(i in which(apply(composition, 1, sum) == 0)) {
            composition[i,] <- rpois(comp_no, 0.8)
        }
    }

    # find names (derived from constituents)
    name <- c()
    for(i in 1:N)
        name <- c(name, hcae_create_name(composition[i,], component_names,
                                         name, empty_name))
 
    # Add species for photons / energy source
    energy <- c(energy, max(max(energy)+5, 50))
    name <- c(name, hv_name) 
 
    # Now investigate all possible reactions up to 2x2 and build one 
    # vector containing information whether the reaction is possible
    x <- eval_possible_reas(N+2, no_reordering)
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
        composition <- composition[o,]
        name <- c(name[o], name[length(name)])   

        # reevaluating
        x <- eval_possible_reas(N+2, T)
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

    return(net)
}
