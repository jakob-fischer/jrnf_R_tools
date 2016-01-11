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


# Loading database for transforming chemical species names into latex representations 
# (to make more beautiful tables of reaction networks)
species_latex_names <- read.csv("species_latex_names_db.csv", stringsAsFactors=F)


# Function transforms species name into a more beautiful representation for
# rendering in latex. For specific (special names) the database in 
# species_latex_names_db.csv is used. For all other names the species are
# transformed according to the following pattern:
#
# "1O2B" -> "_1 \mathrm{O}_2 \mathrm{B}"

jrnf_species_name_to_latex_rep <- function(sp_n) {
    x <- which(sp_n == species_latex_names$input)
    if(length(x) == 1)
        return(species_latex_names$output[x])

    h <- NA
    # looking for "_<number>"
    if(grepl("_", sp_n)) {
        h <- strsplit(sp_n, "_", fixed=TRUE)[[1]][2]
        sp_n <- strsplit(sp_n, "_", fixed=TRUE)[[1]][1]
    }

    # looking for numeric at the start
    y <- regexp(sp_n, "[0-9]+")
    if(!is.null(y$start) && y$start[1] == 1) {
        sp_n <- substring(sp_n, first=y$end[1]+1)
        h <- y$match[1]
    }
  
    y1 <- regexp(sp_n, "[[:alpha:]]+")$match
    y1 <- paste("\\mathrm{", y1, "}", sep="")
    y2 <- regexp(sp_n, "[[:digit:]]+")$match
    if(!is.null(y2))
        y2_ <- paste("_", y2, " ", sep="")

    if(length(y2_) != length(y1))
        y2_ <- c(y2_, "")
    y <- c(rbind(y1,y2_))
 
    if(!is.na(h))
        y <- c("_", h, " ", y)

    return(paste(y, collapse=""))
}


# Function adds one column to dataframe with name "EMPTY" and the content of 
# strings counting from (1) ...  This additional information is used by the
# functions below to add the number / id of reactions or pathways to the 
# generated latex table.
#
# Parameter (<x>) is normally a data frame, but if one only wants to use the
# index as additional information the numbers of rows of the new data frame
# can be given.

pa_h_add_count_column <- function(x) {
    if(!is.data.frame(x))
        return(data.frame(EMPTY=paste("(", as.character(1:x), ")", sep="")))
    else 
        return(cbind(x,
                     data.frame(EMPTY=as.character(paste("(", as.character(1:nrow(x)), ")", sep="")))))
}


# Function Transforms a reaction network into a latex table that then is written 
# to file. 
#
# net      - network that is transformed to table 
# filename - latex table is written to this file
# style=1  - standard style (others not available yet)
# add_info - 1 (numeric) : means that counting is added as additional information
#          - c()         : nothing is added
#          -             : Data frame which has the same number of rows as the 
#                          network has reactions
# marked   - boolean array or vector of id's (numeric) that marks special reactions
#            (not implemented yet but they will be probably hilighted grey - depending
#             on style)
# sep      - vector of id's of reactions after which a separation is drawn ("\hline")
#            (not implemented yet)

jrnf_network_to_ltable <- function(filename, net, add_info=1, marked=c(), sep=c(), style=1) {
    con <- file(filename, "w")

    # Standard argument for add_info (==1) is interpreted as having to replace it
    # by a data frame that adds a number for each reaction as single string.
    if(is.numeric(add_info))
        add_info <- pa_h_add_count_column(nrow(net[[2]]))

    # fill up spaces (adds spaces to string to make the table code look smoother
    fus <- function(x, s=25) {
        if(nchar(x) >= s)
            return(x)
        else 
            return(paste(c(x, rep(" ", s-nchar(x))), collapse=""))
    }

    # insert function writes string to file...
    i <- function(...)  {  writeLines(paste(list(...), collapse=""), con)  }

    # 
    get_layout_head <- function() {  
        b <- "l c r"

        if(is.data.frame(add_info))
            b <- paste(c(b, rep(" c", ncol(add_info)), " "), collapse="")
            
        return(b)  
    }

    #
    draw_header <- function()  {
        if(is.data.frame(add_info) && 
           (ncol(add_info) != 1 || names(add_info)[[1]] != "EMPTY")) {
            x <- " & & "
            for(k in 1:ncol(add_info)) {
                n <- names(add_info)[[k]]
                if(n == "EMPTY") 
                    n <- ""
                x <- paste(x, "& ", n, " ", sep="")
            }

            i(paste(x, "\\\\", sep=""))
            i("\\hline")
        }
    }

    #
    draw_reaction <- function(j) {
        educts_s <- jrnf_educts_to_string(net, j, jrnf_species_name_to_latex_rep, "\\,+\\,")
        products_s <- jrnf_products_to_string(net, j, jrnf_species_name_to_latex_rep, "\\,+\\,") 

        x <- paste("$", fus(educts_s), "$ & $\\rightarrow $ &  $",
                   fus(products_s), " $ ", sep="")
        
        if(is.data.frame(add_info)) 
            for(k in 1:ncol(add_info)) 
                x <- paste(x, "& ", fus(as.character(add_info[j,k])), " ", sep="")

        i(paste(x, "\\\\", sep=""))
    }

    i("% created by pa_network_to_ltable!")
    i("\\begin{table}[h]")
    i("%\\caption[table title]{long table caption}")
    i("%\\label{tab:mytabletable}")
    i("\\begin{tabular}{ ", get_layout_head(), " }")
    i("\\hline")
    draw_header()
    for(k in 1:nrow(net[[2]])) 
        draw_reaction(k)
    i("\\hline")
    i("\\end{tabular}")
    i("\\end{table}")
    close(con)
}



# Legacy references / can be removed if no longer needed
load_jrnf <- jrnf_read
write_jrnf <- jrnf_write

