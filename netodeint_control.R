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

netodeint_setup <- function(netfile, bvalues_l, no_scripts, sample_no,
                            sampling="spath", sampling_par=c(), sampling_sym=TRUE) {
    # Save old and set new working directory
    scripts <- as.character()
    path_old <- getwd()
    path <- unlist(strsplit(netfile, "/", fixed=TRUE))
    setwd(paste(path_old[-length(path_old)], collapse="/"))

    # load network / calculate topologic properties
    net <- jrnf_read(netfile, "pfile.dat", "nfile.dat")
    jrnf_create_pnn_file()

    # sample pairs of boundary conditions
    bids_l <- list(c(1,2), c(2,3), c(6,7))  # example

    # create directory structure
    for(b in bids_l) {
        for(v in bvalues_l) {
            ff <- paste(c("b", as.character(b[1]), as.character(b[2]), "v", 
                          as.character(v[1]), as.character(v[2])), collapse="/")

            system(paste("mkdir", ff))  


            for(i in 1:sample_no) {
                fff <- paste(c(ff, as.character(i)), collapse="/")
                system(paste(c("mkdir", fff), collapse=" "))

                # create network-files and (initial) concentration files
                jrnf_create_initial(net, paste(fff, "run.con", sep="/"), 
                                    paste(fff, "net.jrnf", sep="/"), bc_id=b, bc_v=v)

                # add script command
            }
        }

    }
    
    
    # Restore working directory
    setwd(path_old)
}
