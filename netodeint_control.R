# author: jakob fischer (jakob@automorph.info)
# description: 
# Scripts for automating the solving of reaction networks with the odeint_rnet
# tool. The functions in this file allow to generate scripts which then can 
# be submitted to the "biosys" batch system and solve a big ensemble of reaction
# equation at once.
#
# Parts of this this file / modules functionality are shared with "simulation_builder.R" 
# and "simulation_builder_....R", this could be unified in future. Nevertheless the
# functionality should be maintained as it was used for generate data for the publication
# Thermodynamics of Random Reaction Networks; Fischer, Kleidon, Dittrich; 2015
# TODO some documentation + cleanup

sourced_netodeint_control <- T

if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")

if(!exists("sourced_jrnf_network_io"))
    source("jrnf_network_io.R")


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


netodeint_setup <- function(netfile, bvalues_l, no_scripts, ensemble_s,
                            sampling="spath", sampling_par=c(), sampling_sym=TRUE,
                            odeint_p="~/apps/odeint_rnet", Tmax=100000, deltaT=0.1, v=1, 
                            wint=1000, b_seed=c(), script_lead="bscript_", zero_E=FALSE, bignet=FALSE) {
    # Save old and set new working directory
    scripts <- as.character()        # Vector of script entries / odeint_rnet calls
    scripts_level <- as.integer()    # difficulty of each scripts call
    path_old <- getwd()
    path <- unlist(strsplit(netfile, "/", fixed=TRUE))
    cat("path=", path, "\n")
    cat("newpath=", paste(path[-length(path)], collapse="/"), "\n")
    setwd(paste(path[-length(path)], collapse="/"))

    cat("loading network\n")

    # load network / calculate topologic properties
    net <- jrnf_read(path[length(path)])

        
    # sample pairs of boundary conditions
    # TODO
    if(sampling_sym) {
        bvalues_l <- unique(bvalues_l)

        for(x in bvalues_l)
            bvalues_l[[length(bvalues_l)+1]] <- c(x[2], x[1])
    }



    if(!file.exists("net_energies.jrnf")) {
        cat("sampling energies...\n")
        net <- jrnf_sample_energies(net, v=v, zero=zero_E)
        cat("writing energies in netfile.\n")
        jrnf_write("net_energies.jrnf", net)

        if(bignet) {
            cat("doing bignet sampling\n")

            jrnf_create_pnn_file(net, NA, "nfile.csv")

            bids_l <- list()
            N <- nrow(net[[1]])  
 
            for(i in 1:sampling_par) {
                a <- sample(1:N, 1)
                b <- sample(1:N, 1)

                while(a == b)
                    b <- sample(1:N, 1)

                bids_l[[length(bids_l)+1]] <- c(min(a,b), max(a,b))
            }


        } else {
            cat("topological analysis \n")

            jrnf_create_pnn_file(net, "pfile.csv", "nfile.csv")

            pfile <- read.csv("pfile.csv")
            bids_l <- list()

            bids_l <- calculate_bids_l(pfile, sampling, sampling_par, b_seed)
        }

        cat("creating directory structure and initial files\n")

        # create directory structure
        cat(paste("having", as.character(length(bids_l)), "boundary ids\n"))

        for(b in bids_l) {
            cat(".")
            ff_mid <- paste(c("b", as.character(b[1]), "_" ,as.character(b[2])), collapse="")
            system(paste("mkdir", ff_mid))

            for(v in 1:length(bvalues_l)) {
                ff <- paste(c("b", as.character(b[1]), "_", as.character(b[2]), "/v", 
                            as.character(bvalues_l[[v]][1]), "_" ,as.character(bvalues_l[[v]][2])), collapse="")

                system(paste("mkdir", ff))  


                for(i in 1:ensemble_s) {
                    fff <- paste(c(ff, as.character(i)), collapse="/")
                    system(paste(c("mkdir", fff), collapse=" "))

                    # create network-files and (initial) concentration files
                    if(i == 1 && v == 1) { # Create network-file only if first in ensemble and boundary
                        netref <- paste(ff_mid, "/net.jrnf", sep="/")
                        jrnf_create_initial(net, paste(fff, "run.con", sep="/"), 
                                            netref, bc_id=b, bc_v=bvalues_l[[v]])
                    } else { # 
                        jrnf_create_initial(net, paste(fff, "run.con", sep="/"), 
                                            NA, bc_id=b, bc_v=bvalues_l[[v]])
                    }

                    # add script command
                    cmd <-odeint_p
                    scripts <- c(scripts, paste(odeint_p, " simsim net=", ff_mid, "/net.jrnf con=", fff, "/run.con deltaT=", as.character(deltaT), " Tmax=", as.character(Tmax), " wint=", as.character(wint) , sep=""))
                }

            }
        }

    } else {
        cat("net_energies.jrnf already exists - extending found structure.\n")
	lf <- list.files()
        lf <- lf[file.info(lf)$isdir]
        bids_l <- list()

        for(x in lf) 
            if(substring(x,1,1) == "b"){
                x <- substring(x,2)
                bids_l[[length(bids_l)+1]] <- as.integer(unlist(strsplit(x, "_")))
            }



        for(b in bids_l) {
            cat(".")
            ff_mid <- paste(c("b", as.character(b[1]), "_" ,as.character(b[2])), collapse="")
            
            for(v in 1:length(bvalues_l)) {
                ff <- paste(c("b", as.character(b[1]), "_", as.character(b[2]), "/v", 
                            as.character(bvalues_l[[v]][1]), "_" ,as.character(bvalues_l[[v]][2])), collapse="")

                system(paste("mkdir", ff))  


                for(i in 1:ensemble_s) {
                    fff <- paste(c(ff, as.character(i)), collapse="/")
                    system(paste(c("mkdir", fff), collapse=" "))

                    jrnf_create_initial(net, paste(fff, "run.con", sep="/"), 
                                            NA, bc_id=b, bc_v=bvalues_l[[v]])
                    
                    # add script command
                    cmd <-odeint_p
                    scripts <- c(scripts, paste(odeint_p, " simsim net=", ff_mid, "/net.jrnf con=", fff, "/run.con deltaT=", as.character(deltaT), " Tmax=", as.character(Tmax), " wint=", as.character(wint) , sep=""))
                }

            }
        }
    }
    

    cat("create batch-script-files\n")
    sel <- sample(no_scripts, length(scripts), replace=TRUE)
    
    for(i in 1:no_scripts) {
        con <- file(paste(script_lead, as.character(i), ".sh", sep=""), "w")
        writeLines("#!/bin/sh",con)
        writeLines(paste("#$ -N ", script_lead, as.character(i) ,sep=""), con)
        writeLines("#$ -j y", con)
        writeLines("#$ -cwd", con)

        writeLines((scripts)[sel == i], con)
        close(con)
    }
    
    cat("done\n")
    
    # Restore working directory
    setwd(path_old)
}



# Executed in the directory ..
#
# TODO comment and implement

netodeint_collect_results <- function() {
    save_wd <- getwd()   # just to be save
    df <- data.frame(b1=numeric(), b2=numeric(), v1=numeric(), v2=numeric(),
                     flow1=numeric(), flow2=numeric(), ep_max=numeric(), sp_df=I(list()),
                     re_df=I(list()), last_time=numeric(), last_msd=numeric())

    bdir_v <- list.dirs(recursive=FALSE)

    for(bp in bdir_v) {
        l <- strsplit(bp, "b")[[1]][2]
        b1 <- as.numeric(strsplit(l, "_")[[1]][1])
        b2 <- as.numeric(strsplit(l, "_")[[1]][2])
        setwd(bp)
        net <- jrnf_read("net.jrnf")

        vdir_v <- list.dirs(recursive=FALSE)
        for(vp in vdir_v) {
            l <- strsplit(vp, "v")[[1]][2]
            v1 <- as.numeric(strsplit(l, "_")[[1]][1])
            v2 <- as.numeric(strsplit(l, "_")[[1]][2])
            setwd(vp)

            edir_v <- list.dirs(recursive=FALSE)
            for(ep in edir_v) {
                i <- as.numeric(strsplit(ep,"/")[[1]][2])
                setwd(ep)
                
                cat("doing b1=", b1, "b2=", b2, "v1=", v1, "v2=", v2, "i=", i, "\n")
                run <- read.csv("run.con")

                last_time_ <- run[nrow(run),1]
                last_msd_  <- run[nrow(run),2]
                last_con <- data.frame(con=as.numeric(run[nrow(run), 3:ncol(run)]))
                
                last_flow <- jrnf_calculate_flow(net, last_con$con)

                flowb1 <- jrnf_calculate_flow_dif(b1, net, last_flow$flow_effective)
                flowb2 <- jrnf_calculate_flow_dif(b2, net, last_flow$flow_effective)

                #cat("flowb1 = ", flowb1, "\n")
                #cat("flowb2 = ", flowb2, "\n")

                df <- rbind(df,
                              data.frame(b1=as.numeric(b1), b2=as.numeric(b2), v1=as.numeric(v1), v2=as.numeric(v2),
                              flow1=as.numeric(flowb1), flow2=as.numeric(flowb2), ep_tot=as.numeric(sum(last_flow$entropy_prod)), 
                              sp_df=I(list(last_con)), re_df=I(list(last_flow)), last_time=last_time_, last_msd=last_msd_))

                setwd("..")
            }
 
            setwd("..")
        }

        setwd("..")
    }

    setwd(save_wd)
    results <- df
    save(results, file="results.Rdata")
}


#
# Above command has to be already executed we are only reloading those samples that are not finished
#

netodeint_recollect_results <- function() {
    save_wd <- getwd()   # just to be save
    load("results.Rdata")

    df <-     data.frame(b1=numeric(), b2=numeric(), v1=numeric(), v2=numeric(),
                         flow1=numeric(), flow2=numeric(), ep_max=numeric(), sp_df=I(list()),
                         re_df=I(list()), last_time=numeric(), last_msd=numeric())


    nfinished <- which(results$last_time < 40000 & (results$last_msd > 1e-20 | results$last_msd < 1e-55)) 
    if(length(nfinished) == 0) {
        cat("already finished all the runs!\n")
        return()
    } else {
        cat(nfinished, "  - are not yet done!\n")
    }

    bdir_v <- list.dirs(recursive=FALSE)

    for(bp in bdir_v) {
        l <- strsplit(bp, "b")[[1]][2]
        b1 <- as.numeric(strsplit(l, "_")[[1]][1])
        b2 <- as.numeric(strsplit(l, "_")[[1]][2])
        setwd(bp)
        net <- jrnf_read("net.jrnf")

        vdir_v <- list.dirs(recursive=FALSE)
        for(vp in vdir_v) {
            l <- strsplit(vp, "v")[[1]][2]
            v1 <- as.numeric(strsplit(l, "_")[[1]][1])
            v2 <- as.numeric(strsplit(l, "_")[[1]][2])
            setwd(vp)

            edir_v <- list.dirs(recursive=FALSE)
            for(ep in edir_v) {
                if(!((nrow(df)+1) %in% nfinished)) {
                    df <- rbind(df, results[nrow(df)+1,])
                    cat(".")
                } else {
                    i <- as.numeric(strsplit(ep,"/")[[1]][2])
                    setwd(ep)
                
                    cat("doing b1=", b1, "b2=", b2, "v1=", v1, "v2=", v2, "i=", i, "\n")
                    run <- read.csv("run.con")

                    last_time_ <- run[nrow(run),1]
                    last_msd_  <- run[nrow(run),2]
                    last_con <- data.frame(con=as.numeric(run[nrow(run), 3:ncol(run)]))
                
                    last_flow <- calculate_flow(net, last_con$con)

                    flowb1 <- calculate_flow_dif(b1, net, last_flow$flow_effective)
                    flowb2 <- calculate_flow_dif(b2, net, last_flow$flow_effective)

                    #cat("flowb1 = ", flowb1, "\n")
                    #cat("flowb2 = ", flowb2, "\n")

                    df <- rbind(df,
                                data.frame(b1=as.numeric(b1), b2=as.numeric(b2), v1=as.numeric(v1), v2=as.numeric(v2),
                                    flow1=as.numeric(flowb1), flow2=as.numeric(flowb2), ep_tot=as.numeric(sum(last_flow$entropy_prod)), 
                                    sp_df=I(list(last_con)), re_df=I(list(last_flow)), last_time=last_time_, last_msd=last_msd_))

                    setwd("..")
                }
            }
 
            setwd("..")
        }

        setwd("..")
    }

    setwd(save_wd)
    results <- df
    save(results, file="results.Rdata")
}





calculate_bids_l <- function(pfile, sampling, sampling_par, b_seed) {
    bids_l <- list()

    tmp <- .Random.seed
    if(!is.null(b_seed)) {
        .Random.seed <<- b_seed
    }

    if(sampling == "all") {
        for(i in 1:nrow(pfile)) 
            if(pfile$from[i] < pfile$to[i])
                bids_l[[length(bids_l)+1]] <- c(pfile$from[i], pfile$to[i])  
    } else if (sampling == "spath") {
        for(i in sort(unique(pfile$shortest_path))) {   
            if(is.finite(i) && i != 0) {
                l <- which(pfile$shortest_path == i & pfile$from < pfile$to)

                if(length(l) > sampling_par) 
                    x <- sample(l, sampling_par)
                else 
                    x <- l

                for(y in x)    
                    bids_l[[length(bids_l)+1]] <- c(pfile$from[y], pfile$to[y])  
            }
        }
    } else {
        cat("No valid sampling specified!\n")
        return() 
    }


    .Random.seed <<- tmp

    return(bids_l)
}


create_pfile_bignet <- function(calc_sp_mul = TRUE) {
    net <- jrnf_read("net.jrnf")

    lf <- list.files()
    lf <- lf[file.info(lf)$isdir]
    bids_l <- list()

    for(x in lf) 
        if(substring(x,1,1) == "b"){
            x <- substring(x,2)
            bids_l[[length(bids_l)+1]] <- as.integer(unlist(strsplit(x, "_")))
        }

    jrnf_create_pfile_bignet(net, bids_l, "pfile.csv", calc_sp_mul) 
}
