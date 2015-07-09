source("~/development/jrnf_R_tools/jrnf_network.R", chdir=T)

elementary_composition_db <- read.csv("elementary_composition_db.csv", stringsAsFactors=FALSE)


get_element_names <- function() {
    return(unique(elementary_composition_db$element))
}


build_composition_matrix <- function(net) {
    N <- jrnf_calculate_stoich_mat(net)

    x <- matrix(0, ncol=length(get_element_names()), nrow=nrow(N))

    for(i in 1:nrow(elementary_composition_db)) {
        c <- which(elementary_composition_db$element[i] == get_element_names())
        r <- which(elementary_composition_db$species[i] == net[[1]]$name)
        x[r,c] <- elementary_composition_db$c[i]
    }

    return(x)
}


check_reactions_composition <- function(net_N, comp_mat=c()) {
    if(is.matrix(net_N) && is.null(comp_mat)) {
        cat("Error in jrnf_check_reactions_composition:")
        cat("Have to have either network object or composition matrix!")
        return(FALSE)
    }

    N <- jrnf_calculate_stoich_mat(net_N)

    # if no composition matrix is given ("comp_mat") it is calculated from
    if(is.null(comp_mat))
        comp_mat <- jrnf_build_composition_matrix(net_N)

    x <- matrix(0, nrow=ncol(N), ncol=ncol(comp_mat))

    for(c in 1:ncol(comp_mat)) {
        for(r in 1:ncol(N)) {
            x[r,c] <- sum(N[,r]*comp_mat[,c])
        }
    }



    return(x)
}
