# author: jakob fischer (jakob@automorph.info)
# description: 
# R-tools for handling export of jrnf-networks to calculate elementary modes with 
# metatool and handle import of results. Most of this is redundant as the 
# "pathway_analysis" codeset contains a simple implementation of the elementary
# modes algorithm and the code below is not able to interact with the newer (matlab)
# version of metatool. The code requires metatool to be locatet in the subdirectory 
# "metatool" of the current directory if used.
# see: http://pinguin.biologie.uni-jena.de/bioinformatik/networks/metatool/metatool.html
# This is still here out of legacy reasons, but please consider it as DEPRECATED.

sourced_interface_topana <- T


# helper for 'write_jrnf_to_metatool' (used by 'jwtm_m_reaction_to_string'):
# Transforms a vector of products or educts id and a vector
# of asocciated multiplicities to an string that can be read
# from metatool. 
# 'co' - vector of id's
# 'co_mul' - vector of multiplicites (same size!)
# 'species' - species data frame (species$name contains species names)

jwtm_m_rpart_to_string <- function(co, co_mul, species) {
    bla <- ""

    if(length(co))
        for(i in 1:length(co)) {
            for(j in 1:co_mul[i]) {
                bla <- paste(bla, " ", as.character(species$name[co[i]]), " ", sep="")

                if(i != length(co) || j !=co_mul[i])
                    bla <- paste(bla, " + ", sep="")
            }  
        }

    return(bla)
}


# helper for 'write_jrnf_to_metatool':
# Transforms a reaction row out of the jrnf 'reaction' data-frame to a string 
# readable by metatool. The second argument is the species data-frame 'species'.

jwtm_m_reaction_to_string <- function(reaction, species) {
    bla <- ""

    bla <- paste(bla,
                 jwtm_m_rpart_to_string(reaction$educts[[1]],
                                        reaction$educts_mul[[1]],
                                        species), 
                 sep="")

    bla <- paste(bla, " = ", sep="")

    bla <- paste(bla,
                 jwtm_m_rpart_to_string(reaction$products[[1]],
                                        reaction$products_mul[[1]],
                                        species), 
                 sep="")

    return(bla)
}


# Writes a jrnf-network given by 'network' to file named 'filename'. The vectors
# 'input' and 'output' contain the id's of those species that metatool will consider
# to be input / output. 

jrnf_write_to_metatool <- function(network, filename, input, output) {
    sp <- network[[1]]
    re <- network[[2]]

    con <- file(filename, 'w');

    writeLines("-ENZREV", con=con)
    writeLines("", con=con)
    
    writeLines("-ENZIRREV", con=con)
    ac <- ""

    nxt <- nrow(re) + length(input) + length(output)
    for(i in 1:nxt) 
        ac <- paste(ac, "R", as.character(i), " ", sep="")

    writeLines(ac, con=con)
    writeLines("", con=con)

    writeLines("-METINT", con=con)
    writeLines( paste(as.character(sp$name), collapse=" "), con=con)
    writeLines("", con=con)

    writeLines("-METEXT", con=con)
    writeLines("", con=con)

    writeLines("-CAT", con=con)
    for(i in 1:nrow(re))
        writeLines(paste("R", as.character(i), " : ", jwtm_m_reaction_to_string(re[i,], sp), " ", sep=""), con=con)

    nxt <- nrow(re)+1
    for(i in input) {
         writeLines(paste("R", as.character(nxt), " : = ", sp$name[i], sep=""), con=con)
         nxt <- nxt+1
    }

    for(i in output) {
         writeLines(paste("R", as.character(nxt), " : ", sp$name[i], "   = ", sep=""), con=con)
         nxt <- nxt+1
    }

    close(con)
}


# helper for 'jrnf_write_to_expa' (used by 'jwtm_m_reaction_to_string'):
# Transforms a number of species 'co' including their multiplier 'co_mul' 
# into a string compatible with expa-format. 
#

jwte_e_rpart_to_string <- function(co, co_mul, sgn, species) {
    bla <- ""

    if(length(co))
        for(i in 1:length(co)) {
            for(j in 1:co_mul[i]) {
                if(sgn)
                    bla <- paste(bla, " 1 ", as.character(species$name[co[i]]), " ", sep="")
                else
                    bla <- paste(bla, " -1 ", as.character(species$name[co[i]]), " ", sep="")
            }  
        }

    return(bla)

}


# helper for 'jrnf_write_to_expa':
# Transforms one row of reaction data frame in the reaction string compatible with 
# expa-file format.

jwte_e_reaction_to_string <- function(reaction, species) {
    bla <- ""

    bla <- paste(bla,
                 jwte_e_rpart_to_string(reaction$educts[[1]],
                                        reaction$educts_mul[[1]],
                                        FALSE,
                                        species), 
                 sep="")

    bla <- paste(bla,
                 jwte_e_rpart_to_string(reaction$products[[1]],
                                        reaction$products_mul[[1]],
                                        TRUE,
                                        species), 
                 sep="")

    return(bla)

}


# Writes the stoichiometry / topology of an jrnf reaction network to expa
# format. For a definition of expa format, please look here:
# http://www.ce4csb.org/applications/jexpa/expa.html

jrnf_write_to_expa <- function(network, filename, input, output) {
    sp <- network[[1]]
    re <- network[[2]]

    con <- file(filename, 'w');

    writeLines("(Internal fluxes)", con=con)

    for(i in 1:nrow(re))
        writeLines(paste("R", as.character(i), "  I  ", jwte_e_reaction_to_string(re[i,], sp), " ", sep=""), con=con)

    nxt <- nrow(re)+1

    # write all the input / output reactions
    for(i in input) {
         writeLines(paste("R", as.character(nxt), " I  1 ", sp$name[i], sep=""), con=con)
         nxt <- nxt+1
    }

    for(i in output) {
         writeLines(paste("R", as.character(nxt), " I  -1 ", sp$name[i], sep=""), con=con)
         nxt <- nxt+1
    }

    # no external fluxes here
    writeLines("(External fluxes)\n", con=con)
    close(con)
}

 
# Function reads the output from metatool from file 'filename' and looks for
# the part in the ouput that contains the elementary modes and then extracts
# them as a matrix.

read_elementary_modes <- function(filename) {
    con <- file(filename, "r")
    lines <- readLines(con)
    first_em <- grep("ELEMENTARY MODE", lines)+3
    em_head <- lines[first_em-1]
    em_no <- as.numeric(strsplit(strsplit(em_head, split="r")[[1]][3], split=" x")[[1]][1])


    if(length(first_em) && is.na(first_em)) {
        cat("Did not find elementary mode in input file!")
        close(con)
        return(NA)
    }

    suppressWarnings( last_line <- as.numeric(unlist(strsplit(lines[first_em], split=" "))) )
    first_size <- length(last_line)
    em <- numeric()

    if(first_size == 0) {
        cat("Elementary modes set found empty!")
        close(con)
        return(NA)
    }

    trans <- function(x) {
        return(as.numeric(unlist(strsplit(x, split=" "))))
    }

    x <- unlist(lapply(lines[first_em:(first_em+em_no-1)],trans))

    return(matrix(x, ncol=first_size, byrow=T))
}


# Function invokes metatool to calculate all elementary modes. Parameters are the 
# network 'network' in jrnf-format, and names of species that are input, output or
# both (external species). The function expects metatool to be locatet at
# "metatool/fasttool" (from the current directory).
# TODO: Simplify name-to-id transformation with one function that centralizes error handling!

jrnf_calculate_em_metatool <- function(network, in_names=c(), out_names=c(), ext_names=c()) {
    in_vector <- logical(length=nrow(network[[1]]))
    out_vector <- logical(length=nrow(network[[1]]))
    
    for(i in in_names) {
        li <- which(i == network[[1]]$name)

        if(length(li) > 0)
            in_vector[li] <- TRUE
        else 
            cat("Error in calculate_elementary_modes: could not find ", i) 
    }

    for(i in out_names) {
        li <- which(i == network[[1]]$name)

        if(length(li) > 0)
            out_vector[li] <- TRUE
        else 
            cat("Error in calculate_elementary_modes: could not find ", i) 
    }

    for(i in ext_names) {
        li <- which(i == network[[1]]$name)

        if(length(li) > 0) {
            in_vector[li] <- TRUE
            out_vector[li] <- TRUE
        } else 
            cat("Error in calculate_elementary_modes: could not find ", i) 
    }


    # write network in metatool format to temporary file 'tmp.dat'
    jrnf_write_to_metatool(network, "tmp.dat", which(in_vector), which(out_vector))

    # invoke metatool
    system("./metatool/fasttool tmp.dat tmp_out.dat")

    # read back results and return them
    em <- read_elementary_modes("tmp_out.dat")
    return(em)
}



