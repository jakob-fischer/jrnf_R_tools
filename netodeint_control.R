# author: jakob fischer (jakob@automorph.info)
# date: 6. April 2013
# description: 
# Scripts for automating the solving of reaction networks with the odeint_rnet
# tool. The functions in this file allow to generate scripts which then can 
# be submitted to the "biosys" batch system and solve a big ensemble of reaction
# equation at once.

source("jrnf_network.R")


# The file first analyses the netfile with jrnf_create_pnn_file. The results
# ("pfile.dat", "nfile.dat") are written to the same directory in which the
# <netfile> (jrnf-format) is in. After that for some (ordered) pairs of 
# boundary points inside of the biggest connected part and for every pair of
# concentrations in bvalues_l-list folders are created. In this folder ensemble_s
# subfolders are created that contain a concentration file and a network file that
# can be used for simulating the network. In the directory containing netfile 
# no_scripts of batch-submittable scripts are generated named "run_<n>.sh". Also a 
# script named "submit_all.sh" is generated.
#
#

netodeint_setup <- function(netfile, bvalues_l, no_scripts, ensemble_s,
                            sampling="spath", sampling_par=c(), sampling_sym=TRUE,
                            odeint_p="~/apps/odeint_rnet", Tmax=100000, deltaT=0.1, wint=10000) {
    # Save old and set new working directory
    scripts <- as.character()
    path_old <- getwd()
    path <- unlist(strsplit(netfile, "/", fixed=TRUE))
    setwd(paste(path[-length(path)], collapse="/"))

    cat("loading network\n")

    # load network / calculate topologic properties
    net <- jrnf_read(path[length(path)])

    cat("topological analysis \n")

    jrnf_create_pnn_file(net, "pfile.csv", "nfile.csv")


    # sample pairs of boundary conditions
    # TODO
    if(sampling_sym) {
       bvalues_l <- unique(bvalues_l)

       for(x in bvalues_l)
           bvalues_l[[length(bvalues_l)+1]] <- c(x[2], x[1])
    }


    pfile <- read.csv("pfile.csv")
    bids_l <- list()

    if(sampling == "all") {
        for(i in 1:nrow(pfile)) 
            if(pfile$from[i] < pfile$to[i])
                bids_l[[length(bids_l)+1]] <- c(pfile$from[i], pfile$to[i])  
    } else if (sampling == "spath") {
        for(i in sort(unique(pfile$shortest_path))) {
            if(is.finite(i) && i != 0) {
                l <- which(pfile$shortest_path == i & pfile$from < pfile$to)
                if(length(l) > sampling_par)
                    l <- sample(l, sampling_par)

                for(i in l)    
                    bids_l[[length(bids_l)+1]] <- c(pfile$from[i], pfile$to[i])  
            }
        }
    } else {
        cat("No valid sampling specified!\n")
        return() 
    }

 
    cat("creating directory structure and initial files\n")

    # create directory structure
    for(b in bids_l) {
        for(v in bvalues_l) {
            ff <- paste(c("b", as.character(b[1]), as.character(b[2]), "v", 
                          as.character(v[1]), as.character(v[2])), collapse="_")

            system(paste("mkdir", ff))  


            for(i in 1:ensemble_s) {
                fff <- paste(c(ff, as.character(i)), collapse="/")
                system(paste(c("mkdir", fff), collapse=" "))

               #cat("b=", b, "\n")
               #cat("v=", v, "\n")

                # create network-files and (initial) concentration files
                jrnf_create_initial(net, paste(fff, "run.con", sep="/"), 
                                    paste(fff, "net.jrnf", sep="/"), bc_id=b, bc_v=v)

                # add script command
                cmd <-odeint_p
                scripts <- c(scripts, paste(odeint_p, " simsim net=", fff, "/net.jrnf con=", fff, "/run.con deltaT=", as.character(deltaT), " Tmax=", as.character(Tmax), " wint=", as.character(wint) , sep=""))
            }
        }

    }
    

    cat("create batch-script-files\n")
    sel <- sample(no_scripts, length(scripts), replace=TRUE)

    for(i in 1:no_scripts) {
        con <- file(paste("bscript_", as.character(i), sep=""), "w")
        writeLines(scripts[sel == i], con)
        close(con)
    }
    
    cat("done\n")
    
    # Restore working directory
    setwd(path_old)
}
