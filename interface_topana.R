# author: jakob fischer (jakob@automorph.info)
# date: 27. February 2013
# description: 
# R-tools for handling export of jrnf-networks to 
# calculate elementary modes with metatool and handle
# import of results....


wjtm_m_rpart_to_string <- function(co, co_mul, species) {
    bla <- ""

    for(i in 1:length(co)) {
        for(j in 1:co_mul[i]) {
            bla <- paste(bla, " ", as.character(species$name[co[i]]), " ", sep="")

            if(i != length(co) || j !=co_mul[i])
                bla <- paste(bla, " + ", sep="")
        }  
    }

    return(bla)
}

wjtm_m_reaction_to_string <- function(reaction, species) {
    bla <- ""

    bla <- paste(bla,
                 wjtm_m_rpart_to_string(reaction$educts[[1]],
                                        reaction$educts_mul[[1]],
                                        species), 
                 sep="")

    bla <- paste(bla, " = ", sep="")

    bla <- paste(bla,
                 wjtm_m_rpart_to_string(reaction$products[[1]],
                                        reaction$products_mul[[1]],
                                        species), 
                 sep="")

    return(bla)
}



# TODO

write_jrnf_to_metatool <- function(network, filename, input, output) {
    sp <- network[[1]]
    re <- network[[2]]

    con <- file(filename, 'w');

    writeLines("-ENZREV", con=con)
    writeLines("", con=con)
    
    writeLines("-ENZIRREV", con=con)
    ac <- ""

    nxt <- nrow(sp) + length(input) + length(output)
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
    for(i in 1:nrow(sp)) 
        writeLines(paste("R", as.character(i), " : ", wjtm_m_reaction_to_string(re[i,], sp), " ", sep=""), con=con)

    nxt <- nrow(sp)+1
    for(i in input) {
         writeLines(paste("R", as.character(nxt), " : = ", sp$name[i], sep=""), con=con)
         nxt <- nxt+1
    }

    for(i in 1:output) {
         writeLines(paste("R", as.character(nxt), " :", sp$name[i], "   = ", sep=""), con=con)
         nxt <- nxt+1
    }

    close(con)
}



read_elementary_modes <- function(filename) {
    con <- file(filename, "r")
    lines <- readLines(con)
    first_em <- grep("ELEMENTARY MODE", lines)+3

    if(is.na(first_em)) {
        cat("Did not find elementary mode in input file!")
        return(NA)
    }

    suppressWarnings( last_line <- as.numeric(unlist(strsplit(lines[first_em], split=" "))) )
    first_size <- length(last_line)
    em <- numeric()

    if(first_size == 0) {
        cat("Elementary modes set found empty!")
        return(NA)
    }

    while(length(last_line) != 0 && length(last_line) == first_size) {
        em <- c(em, last_line)
        first_em <- first_em + 1
        tryCatch(last_line <- as.numeric(unlist(strsplit(lines[first_em], split=" "))),
                 warning=function(x) {  last_line<<-c()} )
    }

     

    return(matrix(em, ncol=first_size, byrow=T))
}




calculate_elementary_modes <- function(network, in_names, out_names, ext_names) {
    in_vector <- logical(length=nrow(network[[1]]))
    out_vector <- logical(length=nrow(network[[1]]))
    
    for(i in in_names) {
        li <- which(i == network[[1]]$name)

        if(length(li) > 0)
            in_vector[li] = TRUE
        else 
            cat("Error in calculate_elementary_modes: could not find ", i) 
    }

    for(i in out_names) {
        li <- which(i == network[[1]]$name)

        if(length(li) > 0)
            out_vector[li] = TRUE
        else 
            cat("Error in calculate_elementary_modes: could not find ", i) 
    }


    write_jrnf_to_metatool(network, "tmp.dat", in_vector, out_vector)

    system("./metatool/fasttool tmp.dat tmp_out.dat")

    em <- read_elementary_modes("tmp_out.dat")
}
