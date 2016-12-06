# author: jakob fischer (jakob@automorph.info)
# description: 
# Scripts for automating the solving of reaction networks with the odeint_rnet
# tool. The functions in this file allow to generate scripts which then can 
# be submitted to a batch system and solve a big ensemble of reaction
# equation at once.
#
# Parts of this this file / modules functionality are shared with "simulation_builder.R" 
# and "simulation_builder_....R", this could be unified in future. Nevertheless the
# functionality should be maintained as it was used for generate data for the publication
# Thermodynamics of Random Reaction Networks; Fischer, Kleidon, Dittrich; 2015

sourced_netodeint_control <- T

if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")

if(!exists("sourced_jrnf_network_io"))
    source("jrnf_network_io.R")



# Create a initial file for a certain jrnf network which can 
# then be used as a starting point for solving an ode. In default every 
# concentration is initialized with random values drawn from a gaussian
# distribution, added to one. 
# Additional 'bc_id' can be given a vector of species id names or a vector
# of (1-indexed) ids. If this is done, then 'bc_v' should be set to a 
# numeric vector of the same size indicating the respective concentrations.
#
# parameters:
# <jrnf_network>  - Reaction network
# <init_fn>       - Filename of the file in that the initial concentration is written
# <network_fn>    - Name of the while to which the reaction network is written
# <bc_id>         - Vector of boundary species ids or names
# <bc_v>          - Boundary species concentrations
# <kb_T>          - Energy scale (Boltzmann Constant and temperature) 
# <setBCE0>       - Setting the boundary condition species energy to zero?

jrnf_create_initial <- function(jrnf_network, init_fn, network_fn=NA, bc_id=NA, bc_v=NA, kB_T=1, setBCE0=TRUE) {
    jrnf_species <- jrnf_network[[1]]
    jrnf_reactions <- jrnf_network[[2]]    

    # ensure bc_id (if given) is numeric
    if(!is.na(bc_id) && !is.numeric(bc_id)) 
        for(i in 1:length(bc_id)) 
            bc_id[i] <- which(jrnf_species$name == bc_id[i])

    # create data frame with one row 
    df <- data.frame(time=as.numeric(0),msd=as.numeric(0))
    df[as.vector(jrnf_species$name)] <- abs(rnorm(length(jrnf_species$name), mean=mean(bc_v), sd=sqrt(mean((bc_v- mean(bc_v))^2))))

    if(!is.na(bc_id) && !is.na(bc_v)) {
        if(length(bc_id) != length(bc_v)) {
            cat("Error in jrnf_create_initial; bc_id,bc_v length missmatch!")
            return() 
        }

        df[1,bc_id+2] <- bc_v 
    } 

    # and write 
    write.csv(df, init_fn, row.names=FALSE)

    # Don't write network if no filename is given.
    if(is.na(network_fn)) 
        return()

    # Write the network with the boundary species set constant.
    if(!is.na(network_file) && !is.na(bc_id))
        jrnf_network[[1]]$constant[bc_id] <- TRUE;

    # If <setBCE0> also set their energy to zero (and recalculate rates).
    if(!is.na(network_file) && !is.na(bc_id) && setBCE0) {
        jrnf_network[[1]]$energy[bc_id] <- 0
        jrnf_network <- jrnf_calculate_rconst(jrnf_network, kB_T)
    }

    jrnf_write(network_file, jrnf_network)
}


# This function calculates information on the network topology using igraph and 
# writes it to two files: 
# - One file (pfile) contains information specific to pairs of nodes
# shortest path, number of different shortest paths, ...
# - The second file (nfile) contains information specific to nodes, like degree,
# betweenness and whether the node belongs to the biggest cluster of the network.
# This two outputs can be enabled separately by giving a filename to write them 
# to for <pfile> or <nfile>. The topological analysis is done on the directed 
# graph one gets from 'jrnf_to_directed_network'.

jrnf_create_pnn_file <- function(jrnf_network, pfile=NA, nfile=NA) {
    N <- nrow(jrnf_network[[1]])
    g <- jrnf_to_undirected_network(jrnf_network)
    g_s <- simplify(g, remove.multiple=TRUE, remove.loops=FALSE)

    sps_g <- shortest.paths(g)

    if(!is.na(pfile)) {
        df_1 <- data.frame(from=as.numeric(rep(0, N*N)), to=as.numeric(rep(0, N*N)), 
                           shortest_path=as.numeric(rep(0, N*N)), sp_multiplicity=as.numeric(rep(0, N*N)), 
                           sp_multiplicity_s=as.numeric(rep(0, N*N)),
                           stringsAsFactors=FALSE)
        df_1$from <- floor((1:(N*N)-1)/N) + 1
        df_1$to <- (1:(N*N)-1) %% N + 1
            
        df_1$shortest_path <- as.vector(sps_g)
        #for(i in 1:(N*N)) {
        #    df_1$shortest_path[i] <- sps_g[floor((i-1)/N) + 1, (i-1) %% N + 1]
        #}
        df_1$sp_multiplicity <- 0
        df_1$sp_multiplicity_s <- 0


        for(k in 1:N) {
            if(k %% 10 == 1)
                cat(".")
            sp <- get.all.shortest.paths(g, from=k, mode="out")$res
            sp_s <- get.all.shortest.paths(g_s, from=k, mode="out")$res                     

            sp_x <- c()
            for(x in sp)
                sp_x <- c(sp_x, x[length(x)])

            adf <- as.data.frame(table(sp_x))
            df_1$sp_multiplicity[(k-1)*N+as.numeric(as.vector(adf$sp_x))] <- adf$Freq


            sp_s_x <- c()
            for(x in sp_s)
                sp_s_x <- c(sp_s_x, x[length(x)])

            adf <- as.data.frame(table(sp_s_x))
            df_1$sp_multiplicity_s[(k-1)*N+as.numeric(as.vector(adf$sp_s_x))] <- adf$Freq
        }
        cat("\n")

        write.csv(df_1, pfile, row.names=FALSE)
    }


    if(!is.na(nfile)) {
        df_2 <- data.frame(node=as.numeric(1:N), deg_in=as.numeric(rep(NA,N)), deg_out=as.numeric(rep(NA,N)), deg_all=as.numeric(rep(NA,N)),
                           deg_s_in=as.numeric(rep(NA,N)), deg_s_out=as.numeric(rep(NA,N)), deg_s_all=as.numeric(rep(NA,N)), betweenness=as.numeric(rep(NA,N)), betweenness_s=as.numeric(rep(NA,N)), main=as.logical(rep(NA,N)), main_w=as.logical(rep(NA,N)), stringsAsFactors=FALSE) 

        df_2$deg_in <- as.numeric(degree(g, mode="in"))
        df_2$deg_out <- as.numeric(degree(g, mode="out"))
        df_2$deg_all <- as.numeric(degree(g, mode="all"))
        df_2$deg_s_in <- as.numeric(degree(g_s, mode="in"))
        df_2$deg_s_out <- as.numeric(degree(g_s, mode="out"))
        df_2$deg_s_all <- as.numeric(degree(g_s, mode="all"))
        df_2$betweenness <- as.numeric(betweenness(g))
        df_2$betweenness_s <- as.numeric(betweenness(g_s))

        strong_clusters <- clusters(g, mode="strong")
        weak_clusters <- clusters(g, mode="weak")
        df_2$main <- (strong_clusters$membership == which.max(strong_clusters$csize))
        df_2$main_w <- (weak_clusters$membership == which.max(weak_clusters$csize))

        write.csv(df_2, nfile, row.names=FALSE)
    }

}



# a faster variant of pfile generation. Only calculates thes pairs that are given in 
# b_list structure...
# TODO

jrnf_create_pfile_bignet <- function(jrnf_network, b_list, pfile, calc_sp_mul=TRUE) {
    L <- length(b_list)
    g <- jrnf_to_undirected_network(jrnf_network)
    g_s <- simplify(g, remove.multiple=TRUE, remove.loops=FALSE)

    df_1 <- data.frame(from=as.numeric(rep(0, L)), to=as.numeric(rep(0, L)), 
                       shortest_path=as.numeric(rep(0, L)), sp_multiplicity=as.numeric(rep(0, L)), 
                       sp_multiplicity_s=as.numeric(rep(0, L)),
                       stringsAsFactors=FALSE)

    for(x in 1:length(b_list)) {
        from <- b_list[[x]][1]
        to <- b_list[[x]][2]
        df_1$from[x] <- from
        df_1$to[x] <- to
            
        cat(".")

        df_1$shortest_path[x] <- as.numeric(shortest.paths(g, from, to))
        df_1$sp_multiplicity[x] <- length(get.all.shortest.paths(g, from=from, to=to, mode="all")$res)        
        df_1$sp_multiplicity_s[x] <- length(get.all.shortest.paths(g_s, from=from, to=to, mode="all")$res)
    }        
    cat("\n")

    write.csv(df_1, pfile, row.names=FALSE)
}


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
                            odeint_p=sb_odeint_path, Tmax=1000, deltaT=0.1, v=1, 
                            wint=100, b_seed=c(), script_lead="bscript_", zero_E=FALSE, bignet=FALSE) {
    # Save old and set new working directory
    scripts <- as.character()        # Vector of script entries / odeint_rnet calls
    scripts_level <- as.integer()    # difficulty of each scripts call
    path_old <- getwd()
    path <- unlist(strsplit(netfile, "/", fixed=TRUE))
    cat("path=", path, "\n")
    cat("newpath=", paste(path[-length(path)], collapse="/"), "\n")

    if(length(path) != 1)
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

                flowb1 <- jrnf_calculate_concentration_change(net, last_flow$flow_effective)[b1]
                flowb2 <- jrnf_calculate_concentration_change(net, last_flow$flow_effective)[b2]

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

                    flowb1 <- jrnf_calculate_concentration_change(net, last_flow$flow_effective)[b1]
                    flowb2 <- jrnf_calculate_concentration_change(net, last_flow$flow_effective)[b2]

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
