# author: jakob fischer (jakob@automorph.info)
# description: 
# File contains code that is necessary for pathway analysis code (pathway_analysis.R)
# but generic enough to be usefull outside of the specific algorithm used there.

sourced_pathway_analysis_tools <- T

if(!exists("sourced_tools"))
    source("tools.R") 

if(!exists("sourced_cycles"))
    source("cycles.R")

if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")


# Function prints elementary mode <em> of network <net> to standard output. If
# the option <discard_s> is active reactions that have only one educt and no 
# products or only one product and no educts are ignored (they are probably) 
# pseudo-reactions for in- / outflow and thus concatenated.

pa_print_em <- function(net, em, discard_s=T) {
    em <- as.numeric(em)
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)   
    N <- N_out - N_in
    x <- ( (apply(N_in != 0, 2, sum) + 
            apply(N_out != 0, 2, sum)) == 1)      # remove reactions with just one species involved
    if(discard_s)
        em[x] <- 0    

    # list all (remaining) reactions with coefficients 
    for(i in which(em != 0)) {
        cat(em[i], "X : ")
        cat(jrnf_reaction_to_string(net, i), "\n") 
    }

    cat("============================================================\n")
    # now print the net reaction    
    x <- N %*% em
    ed <- which(x < 0)    

    if(length(ed) > 0)
    for(i in 1:length(ed)) {
        if(-x[ed[i]] != 1)
            cat(-x[ed[i]], " ", sep="")

        cat(net[[1]]$name[ed[i]])

        if(i < length(ed))
            cat(" + ")
    }

    cat(" => ")

    prod <- which(x > 0)    

    if(length(prod) > 0)
    for(i in 1:length(prod)) {
        if(x[prod[i]] != 1)
            cat(x[prod[i]], " ", sep="")

        cat(net[[1]]$name[prod[i]])

        if(i < length(prod))
            cat(" + ")
    }

    cat("\n")    
}


# This function writes a list of <N> most relevant elementary modes /
# pathways to the file with filename <filename>.
# Additional parameters:
# <em>    - elementary modes in a matrix
# <coef>  - coefficients of elementary modes
# <exp_f> - explained fraction of elementary modes (calculated if not given) 

pa_write_em_set <- function(filename, net, em, coef, exp_f=c(), N=100) {
    # Don't have to do anything if no em's present
    if(nrow(em) == 0)
        return()

    if(is.null(exp_f)) {
       exp_f <- apply(abs(scale_mat_rows(em, coef)), 1, sum)
       exp_f <- exp_f / sum(exp_f)
    }

    N <- min(N, nrow(em))
    x <- c()

    for(i in 1:N) {
        x <- c(x, capture.output(cat(i, ": coefficient is", coef[i], "   explained fraction", exp_f[i])))
        x <- c(x, capture.output(cat("============================================================")))
        x <- c(x, capture.output(pa_print_em(net, em[i,])))
        x <- c(x, capture.output(cat("\n")))
    }

    writeLines(x, filename)
}


# Calculates a set of coefficients for a given set of pathways and a steady 
# state flux vector by "distributing" the steady state flux to as many pathways
# as possible.
#
# Parameters:
# <M>       - Matrix of elementary modes / pathways
# <v>       - Steady state flux vector
# <con_fb>  - Give console feedback on pathway / decomposition quality.

pa_calculate_coef <- function(M, v, con_fb=F) {
    coef <- rep(0, nrow(M))
    v_rem <- v      # REMaining flux vector (not explained by pw + coef till now)
    v_active <- v_rem != 0                       # rates / reactions that are not saturated yet
    pw_active <- !apply(matrix(M[,!v_active] != 0, nrow=nrow(M)), 1, any)   # pathways compatible with non_saturated pathways
    
    while(any(v_active) && any(pw_active)) {
       a <- apply(matrix(M[pw_active,], ncol=ncol(M)), 2, sum)
       b <- (v_rem / a)
       b[a == 0] <- NA

       coef[pw_active] <- coef[pw_active] + min(b, na.rm=T)
       v_active[which.min(b)] <- F
       pw_active <- !apply(matrix(M[,!v_active] != 0, nrow=nrow(M)), 1, any)   # pathways compatible with non_saturated pathways

       # recalculate remaining rate part
       if(any(pw_active)) {
           v_rem <- pmax(0, v - apply(scale_mat_rows(M, coef), 2, sum))
           v_rem[!v_active] <- 0   # counter numerical errors
            
           if(con_fb)
               cat(".")
       }
    }

    v_res <- t(M) %*% coef # apply(scale_mat_rows(M, coef), 2, sum)
    if(con_fb)
        cat("\nPerfomance evaluation:", v_res / v, "\n")
    score_b <- max(abs(1-v_res/v)[v != 0])
    return(list(coef=coef, score_b=score_b))
}


# The helper function takes a matrix (complete list of pathways) as argument
# and identifies those that are actually elementary by checking which ones are 
# "contained" in other pathways from the complete set of pathway. 

pa_find_elementary_pw <- function(mat) {
    keep_pw <- rep(T, nrow(mat))         
    c_nzero_row <- apply(mat != 0, 1, sum)
        
    for(i in 1:nrow(mat)) {
        for(j in 1:nrow(mat)) 
            if(all(mat[i,] != 0 | mat[j,] == 0) & 
                   c_nzero_row[i] != c_nzero_row[j])
                    keep_pw[i] <- F
        }

    return(keep_pw)
}


# Function extends the pathway representation to one that does not only contain
# the (positive integer) coefficients of all the reactions but also integer 
# coefficients for the change of all species concentration by applying the pathway.
# This representation first contains nrow(N) coefficients encoding the change of all
# species through the pathway followed by the ncol(N) coefficients for the reactions.
# Function handles pathways as vectors, matrices (multiple pathways) or lists composed
# of pathway matrices and a vector of their respective coefficients.

pa_extend_pathway_representation <- function(pw, N) {
    if(is.list(pw)) # if <pw> is not a matrix of pathways but a list of pathways
                    # $M and coefficients $coef...
        return(list(M=pa_extend_pathway_representation(pw$M, N), coef=pw$coef))

    if(is.vector(pw))
        return(c(N %*% pw, pw))

    return(t(apply(dec_1$M, 1, function(x)  {  return(c(N %*% x, x))}  )))
}


# Given a reaction network and reaction rates this function adds inflow and outflow
# reactions that would maintain steady state. Function returns a list with new network
# and new rates.
# <bound> is used to determined if there is an actual in- / outflow to a species or if
#         it is actually zero. For species without in- / outflow no pseudoreaction is 
#         added. If bound is boolean (vector) with the length equaling the number of 
#         species it is used to select those species for which an exchange pseudoreaction
#         is added. The actual change is only considered to calculate the new / extended
#         rate vector.
# <unique_dir> is set iff all added (pseudo)-reactions are defined to be outflow reactions
#              and for inflow reactions their rate will be set negative. Else the 
#              reaction direction is choosen the way that the rates are positive.

pa_extend_net <- function(net, rates, bound=0, unique_dir=F) {
    # calculate rates that balance growth / decrease of concentrations
    cdif_r <- jrnf_calculate_concentration_change(net, rates)
    
    if(!is.logical(bound)) 
        bound <- abs(cdif_r) >= bound

    # add reactions to balance growth / decrease
    for(i in which(bound))
        # if species' concentration increases one pseudoreaction has to be included to remove it ("X -> ")
        if(cdif_r[i] > 0 | unique_dir) {   
	    net[[2]] <- rbind(net[[2]], data.frame(reversible=factor(c(FALSE)), 
                              c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                              activation=as.numeric(c(0)),educts=I(list(i)), educts_mul=I(list(1)),
                              products=I(list(c())), products_mul=I(list(c()))))
            rates <- c(rates, cdif_r[i])
        # if species' concentration decreases one pseudoreaction is included to add it ("-> X ")
        } else { 
	    net[[2]] <- rbind(net[[2]], data.frame(reversible=factor(c(FALSE)), 
                              c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                              activation=as.numeric(c(0)),educts=I(list(c())), educts_mul=I(list(c())),
                              products=I(list(i)), products_mul=I(list(1))))
            rates <- c(rates, -cdif_r[i])
        }

    return(list(net=net, rates=rates))
}


# Function checks the reachability of individual species from a set of paths.
# The species that are not produced or consumed by any pathway are printed out.
# This is of interest if one has a list of pathways and wants to 
# "eliminate" a chemical species by connecting all producing with all consuming
# pathways. 

pa_check_reachability <- function(N, path_M) { 
    x <- N %*% t(path_M)

    cat("inflow missing: ")
    for(i in 1:nrow(N))
        if(length(which(x[i,] > 0)) == 0)
            cat(i, " ")
    cat("\n")
    

    cat("outflow missing: ")
    for(i in 1:nrow(N))
        if(length(which(x[i,] < 0)) == 0)
            cat(i, " ")
    cat("\n")
}


# Helper function that takes a matrix of pathways and their coefficients,
# identifies duplicated pathways and then removes multiple occurences of
# pathways and calculate new accumulated rates for the remaining pathways.
# Function returns a logical vector identifying pathways to keep and a vector
# of new coefficients.

pa_rm_duprows_accum <- function(mat, coef) {
    keep_pw <- !duplicated(matrix(mat,nrow=nrow(mat)))
        
    for(i in which(!keep_pw)) {
        k <- 1
        while(!all(mat[k,] == mat[i,]))
            k <- k + 1
            coef[k] <- coef[k] + coef[i]
        }

    return(list(keep=keep_pw, coef=coef[keep_pw]))
} 



# Function calculates the reaction rates given a set of pathways and their
# associated rates.

pa_reconstruct_rates <- function(path_M, path_rates) {
    return(as.vector(t(path_M) %*% path_rates))
}


# Helper function for pa_subpath_decompostion
# The function checks in which rows of the matrix <path_M> the same elements
# are non-zero as in the vector <v> and returns a boolean vector of length
# nrow(path_M)

pa_is_row_contained <- function(v, path_M) {
    if(!is.matrix(path_M))
        return(c(F))

    if(nrow(path_M) == 0)
        return(c(F))
  
    return(apply(path_M, 1, function(x)  {  return(!any(x != 0 & v == 0))  }))
}


# Function Transforms a number of pathways in a reaction network into a latex table 
# that then is written to file. 
#
# net      - network that is transformed to table 
# pw       - matrix describing the pathways (1 pathway = 1 row)           
#          - 
# filename - latex table is written to this file
# style=1  - standard style (others not available yet)
# add_info - 1 (numeric) : means that counting is added as additional information
#          - c()         : nothing is added
#          -             : Data frame which has the same number of rows as pathway
#                        : matrix has rows

pa_pathways_to_ltable <- function(filename, pw, net, add_info=1, style=1, longtable=F) {
    con <- file(filename, "w")
    N <- jrnf_calculate_stoich_mat(net)

    # Standard argument for add_info (==1) is interpreted as having to replace it
    # by a data frame that adds a number for each pathway as single string. 
    if(is.numeric(add_info))
        add_info <- pa_h_add_count_column(nrow(pw))

    # fill up spaces (adds spaces to string to make the table code look smoother)
    fus <- function(x, s=25) {
        if(nchar(x) >= s)
            return(x)
        else 
            return(paste(c(x, rep(" ", s-nchar(x))), collapse=""))
    }

    # insert function writes string to file...
    i <- function(...)  {  writeLines(paste(list(...), collapse=""), con)  }

    # Returns string for the main layout of table 
    # <multiplicity> & <educts> & => & products ...
    get_layout_head <- function() {  
        b <- "l r c l"

        if(is.data.frame(add_info))
            b <- paste(c(b, rep(" c", ncol(add_info)), " "), collapse="")
            
        return(b)  
    }
 
    # draws (write to file) the header
    draw_header <- function()  {
        if(is.data.frame(add_info) && 
           (ncol(add_info) != 1 || names(add_info)[[1]] != "EMPTY")) {
            x <- " & & & "
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

    # builds (return as string) the specific reaction with specific multiplicity
    # (if multiplicity is negative it means direction is reversed)
    build_rea <- function(rea, mul) {
        educts_s <- jrnf_educts_to_string(net, rea, jrnf_species_name_to_latex_rep, "\\,+\\,")
        products_s <- jrnf_products_to_string(net, rea, jrnf_species_name_to_latex_rep, "\\,+\\,") 
        x <- ""

        if(mul > 1)
            x <- paste(as.character(mul), "X & $", fus(educts_s), "$ & $\\rightarrow $ &  $",
                       fus(products_s), " $ ", sep="")        
        else if (mul == 1) 
            x <- paste("  & $", fus(educts_s), "$ & $\\rightarrow $ &  $",
                       fus(products_s), " $ ", sep="")        
        else if(mul < -1)
            x <- paste(as.character(-mul), "X & $", fus(products_s), "$ & $\\rightarrow $ &  $",
                       fus(educts_s), " $ ", sep="")        
        else if (mul == -1) 
            x <- paste("  & $", fus(products_s), "$ & $\\rightarrow $ &  $",
                       fus(educts_s), " $ ", sep="")        
        return(x)
    }

    # builds (return as string) the additional columns / information for pathway
    # 'l'. If first is numberic, the data (multirow) is written, else empty columns. 
    build_info <- function(l, first=F) {
        x <- ""

        if(is.data.frame(add_info)) {
            for(k in 1:ncol(add_info)) {
                if(is.numeric(first))
                    x <- paste(x, "& ", fus(paste("\\multirow{", as.character(first) ,
                                                  "}{*}{", as.character(add_info[l,k]),
                                                  "}", sep="")), 
                               " ", sep="")
                else
                    x <- paste(x, "& ", fus(" "), " ", sep="")
            }
        }
 
        return(x)
    }

    # The function draws the entire pathway <j>. 
    draw_pw <- function(j) {
        # first list all the reactions in pathway
        x <- which(pw[j,] != 0)
        first <- length(x)+1
        accum_s <- N %*% pw[j,]

        for(k in x) {
            # For each reaction write reaction part followed by additional info
            # part. All additional info is only written the first time and afterwards
            # only entry entries are written (multicolumn). 
            i(paste(build_rea(k, pw[j,k]), build_info(j, first), "\\\\*", sep=""))
            first <- F
        }

        # now draw effective reaction
        # left hand side (generate string)
        lhs <- jrnf_side_to_string(which(accum_s < 0),
                                   abs(accum_s[accum_s < 0]),
                                   net, jrnf_species_name_to_latex_rep, "\\,+\\,")
        
        # right hand side (generate string)
        rhs <- jrnf_side_to_string(which(accum_s > 0),
                                   accum_s[accum_s > 0],
                                   net, jrnf_species_name_to_latex_rep, "\\,+\\,")
        
        i("\\cline{1-4}")
        i(paste("net: & $", lhs, " $ & $\\rightarrow $ &  $", rhs, " $ ", build_info(j, first), "\\\\*", sep="")) 
        i("\\hline")
    }

    # Code for actual table layout. Different parts when <longtable> is selected.
    con <- file(filename, "w")
    i("% created by pa_pathways_to_ltable!")
    if(!longtable) {
        i("% please don't forget to include \\usepackage{multirow}.")
        i("\\begin{table}[h]")
        i("%\\caption[table title]{long table caption}")
        i("%\\label{tab:mytabletable}")
        i("\\begin{tabular}{ ", get_layout_head(), " }")
        i("\\hline")
        draw_header()
        for(k in 1:nrow(pw)) 
            draw_pw(k)
        i("\\end{tabular}")
        i("\\end{table}")
    } else {
        i("% (longtable version - make sure longtable package is active)")
        i("% please don't forget to include \\usepackage{multirow}.")
        i("% also the following redefinition might be necessary:")
        i("% % from http://tex.stackexchange.com/questions/52100/longtable-multirow-problem-with-cline-and-nopagebreak")
        i("%\\makeatletter")
        i("%\\def\\@cline#1-#2\\@nil{%")
        i("%\\omit")
        i("%\\@multicnt#1%")
        i("%\\advance\\@multispan\\m@ne")
        i("%\\ifnum\\@multicnt=\\@ne\\@firstofone{&\\omit}\\fi")
        i("%\\@multicnt#2%")
        i("%\\advance\\@multicnt-#1%")
        i("%\\advance\\@multispan\\@ne")
        i("%\\leaders\\hrule\\@height\\arrayrulewidth\\hfill")
        i("%\\cr")
        i("%\\noalign{\\nobreak\\vskip-\\arrayrulewidth}}")
        i("%\\makeatother")
        i("\\begin{center}")
        i("\\begin{longtable}{ ", get_layout_head(), " }")
        i("\\caption[table title]{long table caption} \\label{tab:mylabel} \\\\")
        i("\\hline")
        draw_header()
        i("\\hline")
        i("\\endfirsthead")
        i("\\caption[]{(continued)}\\\\")
        i("\\hline")
        draw_header()
        i("\\hline")
        i("\\endhead")
        i("\\hline")
        i(paste("\\multicolumn{", as.character(3+ncol(add_info)) ,
                "}{c}{(continued on next page)}"))
        i("\\endfoot")
        i("\\hline")
        i("\\endlastfoot")
        for(k in 1:nrow(pw)) 
            draw_pw(k)
        i("\\end{longtable}")
        i("\\end{center}")
    }
    close(con)
}

