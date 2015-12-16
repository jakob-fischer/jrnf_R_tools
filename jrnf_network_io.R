# author: jakob fischer (jakob@automorph.info)
# description: 
# Part of the jrnf R codebase that is responsible for loading and writing networks 
# to the file system (in different file formats). Internally the reaction networks
# are stored in a list of two data frames as described in the beginning of
# 'jrnf_network.R'.

sourced_jrnf_network_io <- T

if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")          


# Loads from jrnf-file (see 'jrnf_description' for file format specification)

jrnf_read <- function(filename) {
    # Open file and read in all lines (and close file again) 
    con <- file(filename, 'r');
    lines <- readLines(con)
    close(con)

    # check file header
    if(lines[1] != "jrnf0003") {
        cat("File ", filename, "is not in valid jrnf format!")
        return(0)
    }    

    # Read line using ' ' to separate and put parts in 'last_line'-vector
    last_line <- scan(textConnection(lines[2]), sep=" ", quiet="TRUE")
    if(length(last_line) != 2) {  # second line should consist of two numbers (see below)
        cat("Error at reading header of ", filename, "!")
        return(0)
    }  

    # species number is first entry in second line and reaction number second (separated by ' ') 
    no_species <- last_line[1]
    no_reactions <- last_line[2]


    # Method that is applied to strings and converts them to one
    # row in the species dataframe.
    sp_a <- function(sp_string) {
        last_line <- strsplit(sp_string, " ", fixed = TRUE)[[1]]  # separate string
        if(length(last_line) != 4)       # all species strings consist of four parts
            cat("Error at reading species!")
        
        return(data.frame(type=as.integer(last_line[1]), name=as.character(last_line[2]), 
                          energy=as.numeric(last_line[3]),
                          constant=as.logical(as.logical(as.integer(last_line[4]))),
                          stringsAsFactors=FALSE))
    }

    # Declaring method that is applied to strings and converts them to one
    # row in the reaction data frame. Similar to 'sp_a', but reaction strings
    # can consist of parts (' '-separated) of any number larger than 6...
    re_a <- function(re_string) {
        last_line <- strsplit(re_string, " ", fixed = TRUE)[[1]] 
        if(length(last_line) < 7) {
            cat("Error at reading reaction:\n")
            cat(last_line, "\n")
        }

        # educt and product number are important to know which parts of the
        # last_line vector to identify with which parts of the reaction.
        e_number <- as.numeric(last_line[6])
        e <- as.numeric((rep(0, e_number)))
        em <- as.numeric((rep(0, e_number)))
        p_number <- as.numeric(last_line[7])
        p <- as.numeric((rep(0, p_number)))
        pm <- as.numeric((rep(0, p_number)))

        # Fill educts ... 
        if(e_number != 0)
        for(j in 1:e_number) {
            e[j] <- as.numeric(last_line[7+j*2-1])+1
            em[j] <- as.numeric(last_line[7+j*2]) 
	}

        # ...and product vectors 
        if(p_number != 0)
        for(j in 1:p_number) {
            p[j] <- as.numeric(last_line[7+e_number*2+j*2-1])+1
            pm[j] <- as.numeric(last_line[7+e_number*2+j*2]) 
	}

        # necessary for reactions that have empty input...
        if(is.na(e[1])) {
            e <- c()
            em <- c()
        }

        # ...or output
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


# Writes reaction network to jrnf-file 
# (see 'jrnf_description' for file format specification)

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


# Loads reaction networks in the format of the krome-chemistry-package.
# See http://kromepackage.org/ for details.

jrnf_read_krome <- function(filename) {
    # Open file and read all lines
    con <- file(filename, 'r');
    lines <- readLines(con)
    close(con)

    # remove lines starting with "#" or "@"
    lines <- lines[!grepl("^(#|@)", lines)]

    # species names are collected here in the order they occur (unique)
    species_names <- c()

    # reaction data frame is filled step by step
    reactions <- data.frame(reversible=logical(), c=numeric(), 
                            k=numeric(), k_b=numeric(), 
                            activation=as.numeric(), 
                            educts=list(), educts_mul=list(),
                            products=list(), products_mul=list())

    for(i in 1:length(lines)) {
        last_line <- strsplit(lines[i], ",")[[1]]

        if(lines[i] == "" || lines[i] == " ") {
        } else {
            if(length(last_line) == 8) {
                e <- last_line[2:4]
                p <- last_line[5:7]
            } else if(length(last_line) == 5) {
                e <- c("hv", last_line[2])
                p <- last_line[3:4]
            } 

            for(l in c(e,p))
                if(l != "" && !(l %in% species_names))
                    species_names <- c(species_names, l)

            e_id <- c()
            p_id <- c()

            for(l in e) 
                if(l != "")
                    e_id <- c(e_id, which(l == species_names))

            for(l in p) 
                if(l != "")
                    p_id <- c(p_id, which(l == species_names))

            # add reaction to reaction dataframe
            reactions <- rbind(reactions,
                               data.frame(reversible=as.logical(FALSE),
                                          c=as.numeric(0), 
                                          k=as.numeric(0),
                                          k_b=as.numeric(0), 
                                          activation=as.numeric(1), 
                                          educts=I(list(e_id)), 
                                          educts_mul=I(list(rep(1, length(e_id)))),
                                          products=I(list(p_id)), 
                                          products_mul=I(list(rep(1, length(p_id))))))
        }
    }

    # species dataframe is created last, because all species names are only available 
    # after processing all reactions. 
    species <- data.frame(type=as.integer(rep(1, length(species_names))), 
                          name=as.character(species_names), 
                          energy=as.numeric(rep(0, length(species_names))),
                          constant=rep(FALSE, length(species_names)),
                          stringsAsFactors=FALSE)

    return(list(species, reactions))
}


# Writes jrnf-network ('data') in a human readable format to file 'filename'.
# For exact description of text representation of reactions look into
# 'jrnf_reaction_to_string'.

jrnf_write_to_text <- function(filename, data) {
    con <- file(filename, 'w');
    if(nrow(data[[2]]) > 0) 
        for(i in 1:nrow(data[[2]]))
            writeLines(jrnf_reaction_to_string(data, i), con=con)
    close(con)
}


# Write the topologic part of the network structure to rea-format.
# "Topologic" here means that all thermodynamic information is not saved 
# (but only stoichiometric information / reaction equations).
# TODO add reference to rea-format

jrnf_write_to_rea <- function(filename, data) {
    # open file and write header
    con <- file(filename, 'w');
    species <- data[[1]]
    reactions <- data[[2]]
    writeLines("# generated by jrnf_write_to_rea", con=con)

    writeLines("# Number of Components", con=con)
    writeLines( as.character(nrow(species)), con=con)
    
    # Write out the name of all species / components
    writeLines("# Components", con=con)
    for(i in 1:nrow(species)) 
        writeLines(as.character(species$name[i]) , con=con)
    
    writeLines("# Number of Reactions", con=con)
    writeLines( as.character(nrow(reactions)), con=con)

    # Write out the information on all the individual reactions
    writeLines("# Reactions", con=con)
    for(i in 1:nrow(reactions)) {
        # Rename (for convenience)
        educts <- unlist(reactions$educts[[i]])
        educts_mul <- unlist(reactions$educts_mul[[i]])
        products <- unlist(reactions$products[[i]])
        products_mul <- unlist(reactions$products_mul[[i]])

        line <- ""

        # Add educts and products, each single entry consists of a pair: one number
        # (multiplicity) and the name of the species.
        if(length(educts) != 0) {
	    for(j in 1:length(educts)) 
                line <- paste(line, as.character(educts_mul[j]), " ", 
                              as.character(species$name[educts[j]]), " ", sep="")
        }

        # educt and products entries are separated by "->"
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


# Function loads network data from rea-file format.
# TODO add link to rea specification
 
jrnf_read_rea <- function(filename) {
    # Open file, read all lines and close file object.
    con <- file(filename, 'r');
    lines <- readLines(con)
    close(con)

    # remove all lines that are comments (contain '#')
    is_comment <- function(x)  { return(grepl("#", x))  }
    x <- unlist(lapply(lines, is_comment))
    lines <- lines[!x]


    no_species <- as.numeric(lines[1])
    no_reactions <- as.numeric(lines[2+no_species])
    
    # No apply needed here to load species dataframe, because the data (strings)
    # for species only contains one entry (specie's name).
    species <- data.frame(type=rep(as.integer(1), no_species), 
                          name=as.character(lines[2:(1+no_species)]), 
                          energy=rep(as.numeric(0), no_species),
                          constant=rep(FALSE, no_species),
                          stringsAsFactors=FALSE)

    # Subroutine uses the already created 'species' data frame to match the
    # species identifier to their id. If not unique, an error messae is printed.
    sp_n_to_id <- function(name) {
        id <- which(species$name == name)[1]
               
        if(length(id) != 1)
            cat("ID = ", id, "\n")
               
        if(is.na(id))
            cat("ID is NA: >", ed[i*2], "<\n")
 
        return(id)
    }

    # Method transforms one reaction string into one row of jrnf-reaction dataframe.
    # For this the method accesses the subroutine above to translate species identifier
    # to species ids.
    re_a <- function(re_string) {
        # Split educts and products part (by '->') and individual entries (by ' ')
        educts_str <- strsplit(re_string, "->")[[1]][1]
        products_str <- strsplit(re_string, "->")[[1]][2]
        ed <- strsplit(educts_str, " ")[[1]]
        pro <- strsplit(products_str, " ")[[1]][-1]

        e <- c()
        em <- c()
        p <- c()
        pm <- c()

       # separate strigs 
       if(length(ed) != 0)
           for(i in 1:(length(ed)/2)) {
               e <- c(e, sp_n_to_id(ed[i*2]))
               em <- c(em, as.numeric(ed[i*2-1]))
           }

       if(length(pro) != 0)
           for(i in 1:(length(pro)/2)) {
               p <- c(p, sp_n_to_id(pro[i*2]))
               pm <- c(pm, as.numeric(pro[i*2-1]))
           } 

        return(data.frame(reversible=as.logical(FALSE),
                          c=as.numeric(0), 
                          k=as.numeric(0),
                          k_b=as.numeric(0), 
                          activation=as.numeric(1), 
                          educts=I(list(e)), educts_mul=I(list(em)),
                          products=I(list(p)), products_mul=I(list(pm))))
    }

    # Input data on reactions
    reactions <- do.call("rbind", lapply(lines[(3+no_species):(2+no_species+no_reactions)], re_a))

    # concatenate it to jrnf-object
    return(list(species, reactions))
}


# Legacy references / can be removed if no longer needed
load_jrnf <- jrnf_read
write_jrnf <- jrnf_write

