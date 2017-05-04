# author: jakob fischer (mail@jakobfischer.eu)
# description: Basic functionality for using elementary composition of chemical
#              species for network. At the moment this only contains loading of
#              a elementary composition database and a function for simply testing
#              if chemical reactions violate mass conservation of elementary 
#              constitutents.

sourced_composition_analysis <- T

if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")

ca_elementary_composition_db <- read.csv("elementary_composition_db.csv", stringsAsFactors=FALSE)


# Returns list of species that are in the elementary composition database

ca_get_element_names <- function() {
    return(unique(ca_elementary_composition_db$element))
}


# Calculates a composition matrix that describes the elementary composition of all
# species in the network (assuming all species are in the database).

ca_build_comp_matrix <- function(net) {
    N <- jrnf_calculate_stoich_mat(net)

    x <- matrix(0, ncol=length(ca_get_element_names()), nrow=nrow(N))

    for(i in 1:nrow(ca_elementary_composition_db)) {
        c <- which(ca_elementary_composition_db$element[i] == ca_get_element_names())
        r <- which(ca_elementary_composition_db$species[i] == net[[1]]$name)
        x[r,c] <- ca_elementary_composition_db$c[i]
    }

    return(x)
}


# Assuming all species are contained in the database this function calculates the 
# effective change of all elementary components by the networks reactions. For
# mass conservation the matrix that is returned (rows correspond to species,
# columns to elementary components) has to be zero everywhere.

check_reactions_composition <- function(net_N, comp_mat=c()) {
    # Can only calculate the composition matrix if the network isn't given as a
    # stoichiometric matrix. Else the composition matrix has to be given.
    if(is.matrix(net_N) && is.null(comp_mat)) {
        cat("Error in jrnf_check_reactions_composition:")
        cat("Have to have either network object or composition matrix!")
        return(FALSE)
    }

    N <- jrnf_calculate_stoich_mat(net_N)

    # if no composition matrix is given ("comp_mat") it is calculated from
    if(is.null(comp_mat))
        comp_mat <- ca_build_comp_matrix(net_N)

    x <- matrix(0, nrow=ncol(N), ncol=ncol(comp_mat))

    for(c in 1:ncol(comp_mat)) {
        for(r in 1:ncol(N)) {
            x[r,c] <- sum(N[,r]*comp_mat[,c])
        }
    }

    return(x)
}
