# author: jakob fischer (jakob@automorph.info)
# description: 
# TODO

sourced_simulation_builder_ecol <- T


# Generator loads network from <netfile> and generates <N_energies> sets of energies 
# defining dynamics. For each of these networks (with energies) and each of the boundary
# value sets <bvalues_l> between the boundary species <bids> initial condition for 
# <N_runs> runs are generated and a line for calculation added in the script file <run.sh>
#
# parameters:
# <netfile>     - path to network file
# <bvalues_l>   - list of two element vectors containing the boundary values 
# <bids>        - two element vector with ids of the boundary species
# <N_energies>  - number of energy sets that are drawn for the network
# <N_runs>      - number of runs done for every energy set + boundary value set
#
# TODO: - Add parameter for generation of multiple script-files  
#
# TODO: Maybe this method can be written in terms of the next one (sb_generator_ecol_mul) 
#       with little code / effort? This would remove redundancies...


sb_generator_ecol <- function(netfile, bvalues, cvalues, N_energies, N_runs, 
                              odeint_p="~/apps/jrnf_int", Tmax=10000, wint=50, flat_energies=F, write_log=T, N_runs_eq=c()) {
    if(is.null(N_runs_eq)) N_runs_eq <- N_runs

    # Save old and set new working directory
    scripts <- as.character()        # Vector of script entries / odeint_rnet calls
    path_old <- getwd()              # Save to restore properly
    path <- unlist(strsplit(netfile, "/", fixed=TRUE))
    cat("path=", path, "\n")
    cat("newpath=", paste(path[-length(path)], collapse="/"), "\n")
    if(length(path) > 1)
        setwd(paste(path[-length(path)], collapse="/"))

    wlog_s <- "write_log"
    if(!write_log)
        wlog_s <- ""

    cat("loading network\n")

    # load network / calculate topologic properties
    net <- jrnf_read(path[length(path)])
    con <- jrnf_analyze_ecosystem_constituents(net)
  
    # make net reversible
    net[[2]]$reversible <- rep(T, nrow(net[[2]]))

    cat("topological analysis \n")
    jrnf_create_pnn_file(net, "pfile.csv", "nfile.csv")
    pfile <- read.csv("pfile.csv")
    
    # Calculate stoichiometric matrix            
    N <- jrnf_calculate_stoich_mat(net)

    # sample <N_energies> different sets of energies
    for(i in 1:N_energies) {
        # create new directory and enter
        system(paste("mkdir ", i, sep=""))
        setwd(as.character(i))
        if(!(!flat_energies && i == 1 && N_energies == 1)) 
            net <- jrnf_ae_draw_energies(net, i==1 && flat_energies)
        net <- jrnf_calc_reaction_r(net, 1)
     
        cat("writing energies in netfile.\n")
        jrnf_write("net_energies.jrnf", net)

        # Generate simulations without driving force (hv) by removing all photochemical reactions first
        for(c in cvalues) {
            net_red <- list(net[[1]][-nrow(N),], net[[2]][ N[nrow(N),]==0, ])
            N_red <- N[-nrow(N), N[nrow(N),]==0] 
            jrnf_write("net_reduced.jrnf", net_red)
            ff <- paste(c("v0_", "c", as.character(c)), collapse="")
            system(paste("mkdir", ff)) 
            setwd(ff) 
            if(N_runs_eq != 0)
            for(j in 1:N_runs_eq) {
                s_sample <- sample(ncol(N_red))
                s_sign <- sample(c(-1,1), size=ncol(N_red), replace=T)
                initial_con <- rep(c, length(net[[1]]$name)-1)

                for(k in 1:length(s_sign)) {
                    rea <- N_red[,s_sample[k]]*s_sign[k]

                    l <- min(-((initial_con - c/100)/rea)[rea < 0])
                    initial_con <- pmax(0, initial_con + rea*l)
                }

                df <- data.frame(time=as.numeric(0),msd=as.numeric(0))
                df[as.vector(net_red[[1]]$name)] <- initial_con
                write.csv(df, paste(j, ".con", sep=""), row.names=FALSE)

                scripts <- c(scripts, 
                             paste(odeint_p, " simulate solve_implicit ", wlog_s, " net=", i, "/net_reduced.jrnf con=", i, "/", ff, "/", j, 
                                   ".con wint=", as.character(wint), " Tmax=", as.character(Tmax), sep=""))
            }

            setwd("..")
        }

        # Inner loop (create simulations with different boundary values and differenc initial concentrations
        for(v in bvalues) for(c in cvalues) {
            ff <- paste(c("v", as.character(v), "_", "c", as.character(c)), collapse="")
            system(paste("mkdir", ff)) 
            setwd(ff) 

            for(j in 1:N_runs) {
                # It is not (trivially) possible to just preset the elementary components 
                # concentrations and derive species concentrations by that. 
                # To have random (but consistent) initial conditions all species' concentrations
                # are initialized with value c and afterwards randomly choosen reactions are 
                # applied as far as possible (no concentration should be below c/100).

                s_sample <- sample(ncol(N_red))
                s_sign <- sample(c(-1,1), size=ncol(N_red), replace=T)
                initial_con <- rep(c, length(net[[1]]$name)-1)

                for(k in 1:length(s_sign)) {
                    rea <- N_red[,s_sample[k]]*s_sign[k]

                    l <- min(-((initial_con - c/100)/rea)[rea < 0])
                    initial_con <- pmax(0, initial_con + rea*l)

                    # now apply reaction "rea" until one species has concentration of c/100
                }

                df <- data.frame(time=as.numeric(0),msd=as.numeric(0))
                df[as.vector(net[[1]]$name)] <- c(initial_con, v)
                write.csv(df, paste(j, ".con", sep=""), row.names=FALSE)

                scripts <- c(scripts, 
                             paste(odeint_p, " simulate solve_implicit ", wlog_s, " net=", i, "/net_energies.jrnf con=", i, "/", ff, "/", j, 
                                   ".con wint=", as.character(wint), " Tmax=", as.character(Tmax), sep=""))  
            }
            setwd("..")
        }

        setwd("..")
    }

    # write script file (file doing all the simulation when executed)
    cat("create batch-script-files\n")    
    con <- file("run.sh", "w")
    writeLines("#!/bin/sh",con)
    writeLines("#$ -N run.sh", con)
    writeLines("#$ -j y", con)
    writeLines("#$ -cwd", con)
    writeLines(scripts, con)
    close(con)
    cat("done\n")
    
    # Restore working directory
    setwd(path_old)
}


# Variant of sb_generator_ecol that does not load one network and assign <N_energies> 
# of energies, but is given a <generator>-function from which <N_nets> networks are 
# generated and each assigned <N_energies> of energies. 
#
# parameters:
# <netfile>     - path to network file
# <bvalues>   - vector containing the boundary values 
# <N_nets>      -
# <N_energies>  - number of energy sets that are drawn for the network
# <N_runs>      - number of runs done for every energy set + boundary value set


sb_generator_ecol_mul <- function(path_s, netgen, bvalues, cvalues, N_nets, N_energies, N_runs, 
                              odeint_p="~/apps/jrnf_int", Tmax=10000, wint=50, flat_energies=F, write_log=T, N_runs_eq=0, limit_AE=T, TlogMax=1e7) {
    if(is.null(N_runs_eq)) N_runs_eq <- N_runs
    Tmax_ <- Tmax

    if(Tmax > TlogMax && write_log)
        Tmax_ = TlogMax

    # Save old and set new working directory
    scripts <- as.character()        # Vector of script entries / odeint_rnet calls
    path_old <- getwd()              # Save to restore properly
    path <- unlist(strsplit(path_s, "/", fixed=TRUE))
    cat("path=", path, "\n")
    if(length(path) > 0)
        setwd(paste(path, collapse="/"))

    wlog_s <- "write_log"
    if(!write_log)
        wlog_s <- ""

    next_dir <- 1
    for(l in 1:N_nets) {
        cat("creating network structure!\n")
        net <- netgen()
        net_O <- net
        # make net reversible
        net[[2]]$reversible <- rep(T, nrow(net[[2]]))
        # Calculate stoichiometric matrix            
        N <- jrnf_calculate_stoich_mat(net)

        # sample <N_energies> different sets of energies
        for(i in 1:N_energies) {
            # create new directory and enter
            system(paste("mkdir ", next_dir, sep=""))
            setwd(as.character(next_dir))

            net <- jrnf_ae_draw_energies(net, i==1 && flat_energies, limit_AE)
            net <- jrnf_calc_reaction_r(net, 1)
     
            cat("writing energies in netfile.\n")
            jrnf_write("net_energies.jrnf", net)

            # Generate simulations without driving force (hv) by removing all photochemical reactions first
            for(c in cvalues) {
                net_red <- list(net[[1]][-nrow(N),], net[[2]][ N[nrow(N),]==0, ])
                N_red <- N[-nrow(N), N[nrow(N),]==0] 
                jrnf_write("net_reduced.jrnf", net_red)
                ff <- paste(c("v0_", "c", as.character(c)), collapse="")
                system(paste("mkdir", ff)) 
                setwd(ff) 
                if(N_runs_eq != 0)
                for(j in 1:N_runs_eq) {
                    s_sample <- sample(ncol(N_red))
                    s_sign <- sample(c(-1,1), size=ncol(N_red), replace=T)
                    initial_con <- rep(c, length(net[[1]]$name)-1)

                    for(k in 1:length(s_sign)) {
                        rea <- N_red[,s_sample[k]]*s_sign[k]

                        l <- min(-((initial_con - c/100)/rea)[rea < 0])
                        initial_con <- pmax(0, initial_con + rea*l)
                    }

                    df <- data.frame(time=as.numeric(0),msd=as.numeric(0))
                    df[as.vector(net_red[[1]]$name)] <- initial_con
                    write.csv(df, paste(j, ".con", sep=""), row.names=FALSE)


                    
                    scripts <- c(scripts, 
                                 paste(odeint_p, " simulate solve_implicit net=", next_dir, "/net_reduced.jrnf con=", next_dir, "/", ff, "/", j, 
                                           ".con Tmax=", as.character(1e-8), " wint=", as.character(3), sep=""),
                                 paste(odeint_p, " simulate solve_implicit ", wlog_s, " net=", next_dir, "/net_reduced.jrnf con=", next_dir, "/", ff, "/", j, 
                                           ".con Tmax=", as.character(Tmax_), " wint=", as.character(wint), sep=""))
                    if(Tmax > Tmax_) 
                        scripts <- c(scripts, 
                                 paste(odeint_p, " simulate solve_implicit net=", next_dir, "/net_reduced.jrnf con=", next_dir, "/", ff, "/", j, 
                                           ".con Tmax=", as.character(Tmax), " wint=", as.character(wint), sep=""))

                }

                setwd("..")
            }

            # Inner loop (create simulations with different boundary values and differenc initial concentrations
            for(v in bvalues) for(c in cvalues) {
                ff <- paste(c("v", as.character(v), "_", "c", as.character(c)), collapse="")
                system(paste("mkdir", ff)) 
                setwd(ff) 

                for(j in 1:N_runs) {
                    # It is not (trivially) possible to just preset the elementary components 
                    # concentrations and derive species concentrations by that. 
                    # To have random (but consistent) initial conditions all species' concentrations
                    # are initialized with value c and afterwards randomly choosen reactions are 
                    # applied as far as possible (no concentration should be below c/100).

                    s_sample <- sample(ncol(N_red))
                    s_sign <- sample(c(-1,1), size=ncol(N_red), replace=T)
                    initial_con <- rep(c, length(net[[1]]$name)-1)

                    for(k in 1:length(s_sign)) {
                        rea <- N_red[,s_sample[k]]*s_sign[k]

                        l <- min(-((initial_con - c/100)/rea)[rea < 0])
                        # TODO It occured that the argument of min was empty, check if that is really possible / how to handle it?
                        if(!is.finite(l)) {
                            cat("Error randomizing reaction initials: ", jrnf_reaction_to_string(net, k), "\n")
                            cat("rea=", rea, "\n")
                            print(N_red)
                            print(N)
                            return(list(net, net_O))
                        }
                        initial_con <- pmax(0, initial_con + rea*l)

                        # now apply reaction "rea" until one species has concentration of c/100
                    }

                    df <- data.frame(time=as.numeric(0),msd=as.numeric(0))
                    df[as.vector(net[[1]]$name)] <- c(initial_con, v)
                    write.csv(df, paste(j, ".con", sep=""), row.names=FALSE)

                    scripts <- c(scripts,  
                                 paste(odeint_p, " simulate solve_implicit net=", next_dir, "/net_energies.jrnf con=", next_dir, "/", ff, "/", j, 
                                           ".con Tmax=", as.character(1e-8), " wint=", as.character(3), sep=""),
                                 paste(odeint_p, " simulate solve_implicit ", wlog_s, " net=", next_dir, "/net_energies.jrnf con=", next_dir, "/", ff, "/", j, 
                                           ".con Tmax=", as.character(Tmax_), " wint=", as.character(wint), sep=""))

                    if(Tmax > Tmax_)
                        scripts <- c(scripts,  
                                 paste(odeint_p, " simulate solve_implicit net=", next_dir, "/net_energies.jrnf con=", next_dir, "/", ff, "/", j, 
                                           ".con Tmax=", as.character(Tmax), " wint=", as.character(wint), sep=""))
                    }

                setwd("..")
            }

            next_dir <- next_dir + 1
            setwd("..")
        }
    }

    # write script file (file doing all the simulation when executed)
    cat("create batch-script-files\n")  
    con <- file("run.sh", "w")
    writeLines("#!/bin/sh",con)
    writeLines("#$ -N run.sh", con)
    writeLines("#$ -j y", con)
    writeLines("#$ -cwd", con)
    writeLines(scripts, con)
    close(con)
    cat("done\n")
    
    # Restore working directory
    setwd(path_old)
}



# Function collects results of simulations defined by "sb_generator" after
# simulation was done using the generated script ("run.sh"). Function has 
# to be executed with working directory being the directory containing the
# initial network file! Results are saved in data frame results which is
# stored in the file "results.Rdata".
#
# TODO: add comments
# TODO: check adaption

sb_collect_results_ecol_legacy <- function() {
    save_wd <- getwd()   # just to be save
    df <- data.frame(Edraw=numeric(), Rdraw=numeric(), c=numeric(), v=numeric(), 
                     flow=numeric(), ep_max=numeric(), state_hash=numeric(), sp_df=I(list()), 
                     re_df=I(list()), last_time=numeric(), last_msd=numeric())
    results_nets <- list() # append all reaction networks in order of Edraw to this list

    Edir_v <- list.dirs(recursive=FALSE)

    for(bp in Edir_v) {
        setwd(bp)
        net <- jrnf_read("net_energies.jrnf")
        results_nets[[length(results_nets)+1]] <- net
        N <- jrnf_calculate_stoich_mat(net)
        if(file.exists("net_reduced.jrnf"))
            net_red <- jrnf_read("net_reduced.jrnf")
        else 
            net_red <- NA

        bp <- as.numeric(strsplit(bp,"/")[[1]][2])

        vdir_v <- list.dirs(recursive=FALSE)
        for(vp in vdir_v) {
            l <- strsplit(vp, "v")[[1]][2]
            v <- as.numeric(strsplit(l, "_")[[1]][1])
            c <- as.numeric(strsplit(strsplit(l, "_")[[1]][2], "c")[[1]][2])
            setwd(vp)

            edir_v <- list.files(recursive=FALSE)
            for(ep in edir_v) {
                i <- as.numeric(strsplit(ep,"\\.")[[1]][1])
                
                cat("v=", v, "c=", c, "Edraw=", as.numeric(bp), "Rdraw=", i, "\n")
                run <- read.csv(ep)

                last_time_ <- run[nrow(run),1]
                last_msd_  <- run[nrow(run),2]
                last_con <- data.frame(con=as.numeric(run[nrow(run), 3:ncol(run)]))

                # v == 0 means simulation is not driven and was done with net_red 
                # add concentration of hv==0 for calculation of rates
                if(v == 0 && nrow(net[[1]]) == nrow(net_red[[1]])+1)
                    last_con <- rbind(last_con, 0)

                last_flow <- jrnf_calculate_flow(net, last_con$con)

                # set flow for photochemical reactions to zero if v==0
                if(v == 0 && nrow(net[[1]]) == nrow(net_red[[1]])+1) 
                   last_flow[N[which(net[[1]]$name == "hv"),] != 0,] <- 0

                cc <- jrnf_calculate_concentration_change(net, last_flow$flow_effective)
                flowb <- -cc[nrow(net[[1]])]
                err_cc <- max(abs(cc[-nrow(net[[1]])]))
                err_cc_rel <- err_cc/flowb
                state_hash <- sb_v_to_hash(last_flow$flow_effective, 1e-20)

                df <- rbind(df,
                            data.frame(Edraw=as.numeric(bp),Rdraw=as.numeric(i), 
                                       v=as.numeric(v), c=as.numeric(c), flow=as.numeric(flowb), state_hash=as.numeric(state_hash),
                                       relaxing_sim=as.logical(v == 0), 
                                       err_cc=err_cc, err_cc_rel=err_cc_rel,
                                       ep_tot=as.numeric(sum(last_flow$entropy_prod)), 
                                       sp_df=I(list(last_con)), re_df=I(list(last_flow)),
                                       last_time=last_time_, last_msd=last_msd_))
            }
   
            setwd("..")
        }
        setwd("..")
    }

    setwd(save_wd)
    results <- df
    results_nets <- results_nets[order(unique(results$Edraw))]
    save(results, results_nets, file="results.Rdata")
}


# Function conducts a broad elementary mode analysis for results object generated
# from sb_collect_results function. All simulations were done with same network
# and same boundary species but directions of reactions may differ. Thus the set
# of elementary modes differs. But some elementary modes may be present in different
# simulations. To track the change in occurence + rates of these modes / pathways
# a global list / matrix is created. In this list coefficients of reactions in pathways are
# negative if the dynamic of the network in that the pathway was found favours the 
# backward direction.
#
# TODO: check adaption

sb_em_analysis_ecol_legacy <- function(res, res_nets, c_max=4, do_precise=T) {
    N_sp <- length(res$sp_df[[1]]$con)
    N_re <- length(results$re_df[[1]]$flow_effective)
   
    # matrix where combined elementary modes are saved
    em_matrix <- matrix(0, ncol=N_re+1, nrow=0)
    
    # list of data frames for species specific, reaction specific and elementary mode expansion data
    sp <- res$sp_df
    re <- res$re_df
    em_ex <- list()

    # 
    C_matrix <- matrix(0, ncol=c_max, nrow=nrow(res))
    C_sum <- rep(0, nrow(res))
    sh_p <- rep(0, nrow(res))
    sh_p_d <- rep(0, nrow(res))
    ep_l <- rep(0,nrow(res))
    ep_nl <- rep(0,nrow(res))
    err_rel_max <- rep(0,nrow(res))
    err <- rep(0,nrow(res))
    em_no <- rep(0, nrow(res))



    # which reactions of network are linear
    is_linear <- (unlist(lapply(res_nets[[1]][[2]]$educts_mul, sum)) == 1) 
    # species degree
    deg <- as.vector(degree(jrnf_to_undirected_network(res_nets[[1]])))    

    # helper function
    # calculating elementary modes
    # TODO Extract all general functionality of this method to sb_h_calc_em_info
    #      (that is also used by sb_stability_analysis_eco)
    calc_em_info <- function(rev, cutoff) {
        rates_rev <- abs(rev) 
        net_rev <- jrnf_reverse_reactions(net, rev)
      
        cdif_r <- jrnf_calculate_concentration_change(net_rev, rates_rev)

        # add pseudoreaction generating hv
        n_hv <- nrow(net[[1]])     # species id of hv
        net_rev[[2]] <- rbind(net_rev[[2]], data.frame(reversible=factor(c(FALSE)), 
                                       c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                                       activation=as.numeric(c(0)),educts=I(list(c())), educts_mul=I(list(c())),
                                       products=I(list(n_hv)), products_mul=I(list(1))))
        rates_rev <- c(rates_rev, -cdif_r[n_hv])

        if(do_precise)
            x <- pa_decompose(jrnf_calculate_stoich_mat(net_rev), rates_rev, branch_all=T, cutoff=cutoff)
        else
            x <- pa_analysis(net_rev, rates_rev, 0.01, 0.01, T)
        
        
        # rename results and sort decreasing with fraction of v explained...
        x_em <- x[[1]]
        x_rates <- x[[2]]
        x_sum <- apply(x_em, 1, sum)
        o <- order(x_rates*x_sum, decreasing=T)

        x_em <- matrix(x_em[o,], nrow=nrow(x_em))
        x_rates <- x_rates[o] 
        x_sum <- x_sum[o]

        em_no[i] <<- nrow(x_em)
        # TODO pa_decompose (in contrary to pa_analysis) does not calculate errors
        err_rel_max[i] <<- NA
        err[i] <<- NA
        #err_rel_max[i] <- x[[3]]
        #err[i] <- x[[4]]

        # construct vector containing fraction of rate explained by pathway (with rate)
        exp_r <- rep(0, nrow(x_em))
        exp_d <- rep(0, nrow(x_em))
        

        em_id <- rep(0, nrow(x_em))

        for(l in 1:nrow(x_em)) {
            # 
            exp_r[l] <- sum(x_em[l,]*x_rates[l])/sum(rates_rev)
            exp_d[l] <- x_em[l,length(rates_rev)]*x_rates[l] / rates_rev[length(rates_rev)]

            # make copy from elementary mode
            em <- x_em[l,]   
                
            # now negate those elements (reactions) that are reverted
            em <- em*(c(sign(rev), 1))

            m <- which(apply(em_matrix, 1, identical, em))
            if(length(m) == 0) {    # add em to em_matrix
                em_matrix <<- rbind(em_matrix, em)
                em_id[l] <- nrow(em_matrix)
            } else {
                em_id[l] <- m[1]
            }
        }

        return(data.frame(id=em_id, rate=x_rates, exp_r=exp_r, exp_r_cum=cumsum(exp_r), 
                          exp_d=exp_d, exp_d_acc=cumsum(exp_d)))
    }


    # MAIN LOOP
    #

    for(i in 1:nrow(res)) {
        cat("============================================================\n")
        cat(" i=", i, "  v=", res$v[i], "  c=", res$c[i], "\n")  
        net <- res_nets[[res$Edraw[i]]]
        # calculating directed network
        gc()
        g <- jrnf_to_directed_network_d(net, res$re_df[[i]]$flow_effective)
        cat(".")         
        
        # build data frame refering to expansion of simulation's rate by elementary modes
        if(!res$relaxing_sim[i] & (res$err_cc_rel[i] < 0.1 | is.na(res$err_cc_rel[i])))
            em_ex[[i]] <- calc_em_info(res$re_df[[i]]$flow_effective, res$err_cc[i]/1e2)
        else
            em_ex[[i]] <- NA

        # calculation of cycles 
        for(j in 1:c_max) 
            C_matrix[i,j] <- get_n_cycles_directed(g,j)[[1]]

        C_sum[i] <- sum(C_matrix[i,])

        
        # Entropy production of linear part and of nonlinear part
        ep_l[i] <- sum(res$re_df[[i]]$entropy_prod[is_linear])
        ep_nl[i] <- sum(res$re_df[[i]]$entropy_prod[!is_linear])

        # reaction specific data
        re[[i]]$is_linear <- is_linear

        # species specific data
        sp[[i]]$degree <- deg
        sp[[i]]$mu <- net[[1]]$energy + log(res$sp_df[[i]]$con)
        sp[[i]]$mean_ep <- associate_reaction_to_species(net, res$sp_df[[i]]$entropy_prod)
    }

    
    # Copying back from local objects to <res> object
    res$ep_linear <- ep_l
    res$ep_nonlinear <- ep_nl

    res$re_df <- re
    res$sp_df <- sp
    res$em_ex <- em_ex

    res$C1 <- C_matrix[,1]
    res$C2 <- C_matrix[,2]
    res$C3 <- C_matrix[,3]
    res$C4 <- C_matrix[,4]
    res$C_sum <- C_sum

    res$em_no <- em_no
    res$em_err_max <- err_rel_max
    res$em_err_sum <- err


    results_em <- res
    em_m <- em_matrix
    save(results_em, em_m, file="results_em.Rdata")

    return(list(res, em_matrix))
}




# Function collects results of simulations defined by "sb_generator" after
# simulation was done using the generated script ("run.sh"). Function has 
# to be executed with working directory being the directory containing the
# initial network file! Results are saved in data frame results which is
# stored in the file "results.Rdata". This is the redone function that 
# was done for including all intermediate state in the results object.
# It's results objects are not compatible with the old format!
#
# TODO: add comments
# TODO: check adaption

sb_collect_results_ecol <- function(col_dynamics=T) {
    save_wd <- getwd()   # just to be save
    df <- data.frame(Edraw=numeric(), Rdraw=numeric(), c=numeric(), v=numeric(), 
                     flow=numeric(), ep_max=numeric(), state_hash=numeric(), sp_df=I(list()), 
                     re_df=I(list()), time=numeric(), msd=numeric(), is_last=logical())
    results_nets <- list() # append all reaction networks in order of Edraw to this list

    Edir_v <- list.dirs(recursive=FALSE)

    for(bp in Edir_v) {
        setwd(bp)
        net <- jrnf_read("net_energies.jrnf")
        results_nets[[length(results_nets)+1]] <- net
        N <- jrnf_calculate_stoich_mat(net)
        if(file.exists("net_reduced.jrnf"))
            net_red <- jrnf_read("net_reduced.jrnf")
        else 
            net_red <- NA

        bp <- as.numeric(strsplit(bp,"/")[[1]][2])

        vdir_v <- list.dirs(recursive=FALSE)
        for(vp in vdir_v) {
            l <- strsplit(vp, "v")[[1]][2]
            v <- as.numeric(strsplit(l, "_")[[1]][1])
            c <- as.numeric(strsplit(strsplit(l, "_")[[1]][2], "c")[[1]][2])
            setwd(vp)

            edir_v <- list.files(recursive=FALSE)
            for(ep in edir_v) {
                i <- as.numeric(strsplit(ep,"\\.")[[1]][1])
                
                cat("v=", v, "c=", c, "Edraw=", as.numeric(bp), "Rdraw=", i, "\n")
                run <- read.csv(ep)

                sel_r <- 1:nrow(run)
                if(!col_dynamics) 
                    sel_r <- nrow(run)

                for(t in sel_r) {
                    time_ <- run[t,1]
                    msd_  <- run[t,2]
                    is_last_ <- t == nrow(run)
                    con <- data.frame(con=as.numeric(run[t, 3:ncol(run)]))

                    # v == 0 means simulation is not driven and was done with net_red 
                    # add concentration of hv==0 for calculation of rates
                    if(v == 0 && nrow(net[[1]]) == nrow(net_red[[1]])+1)
                        con <- rbind(con, 0)

                    flow <- jrnf_calculate_flow(net, con$con)

                    # set flow for photochemical reactions to zero if v==0
                    if(v == 0 && nrow(net[[1]]) == nrow(net_red[[1]])+1) 
                        flow[N[which(net[[1]]$name == "hv"),] != 0,] <- 0

                    cc <- jrnf_calculate_concentration_change(net, flow$flow_effective)
                    flowb <- -cc[nrow(net[[1]])]
                    err_cc <- max(abs(cc[-nrow(net[[1]])]))
                    err_cc_rel <- err_cc/flowb
                    state_hash <- sb_v_to_hash(flow$flow_effective, 1e-20)

                    df <- rbind(df,
                                data.frame(Edraw=as.numeric(bp),Rdraw=as.numeric(i), 
                                           v=as.numeric(v), c=as.numeric(c), flow=as.numeric(flowb), state_hash=as.numeric(state_hash),
                                           relaxing_sim=as.logical(v == 0), 
                                           err_cc=err_cc, err_cc_rel=err_cc_rel,
                                           ep_tot=as.numeric(sum(flow$entropy_prod)), 
                                           sp_df=I(list(con)), re_df=I(list(flow)),
                                           time=time_, msd=msd_, is_last=is_last_))

                    if(t == nrow(run) && col_dynamics && nrow(run) > 1) {
                        s <- (nrow(df)-nrow(run)+1):(nrow(df)-1)
                        df$err_cc[s] <- df$err_cc[t]
                        df$err_cc_rel[s] <- NA
                    }
                }
            }
   
            setwd("..")
        }
        setwd("..")
    }

    setwd(save_wd)
    results <- df
    results_nets <- results_nets[order(unique(results$Edraw))]
    save(results, results_nets, file="results.Rdata")
}




sb_h_calc_em_info <- function(net, v, em_matrix, ex_sp, em_cutoff=0) {
    # first add pseudoreactions (with rates to balance change of exchanged species - ex_sp)
    cdif_r <- jrnf_calculate_concentration_change(net, v)

    i <- length(ex_sp)

    net[[2]] <- rbind(net[[2]], data.frame(reversible=factor(rep(FALSE, i)), 
                          c=as.numeric(rep(1, i)), k=as.numeric(rep(1, i)),
                          k_b=as.numeric(rep(0, i)), activation=as.numeric(rep(0, i)),
                          educts=I(rep(I(list(c())), i)), educts_mul=I(rep(I(list(c())), i)),
                          products=I(as.list(ex_sp)), products_mul=I(rep(I(list(1)), i))))
    v_ <- c(v, -cdif_r[ex_sp])


    # now reverse all reactions with negative rates 
    rates_rev <- abs(v_) 
    net_rev <- jrnf_reverse_reactions(net, v_)

    x <- pa_decompose(jrnf_calculate_stoich_mat(net_rev), rates_rev, branch_all=T, cutoff=1e-4*max(abs(cdif_r)))
    #x <- pa_analysis_A(net_rev, rates_rev, 0, 0, T) 

    # rename results and sort decreasing with fraction of v explained...
    x_em <- x[[1]]
    x_rates <- x[[2]]
    x_sum <- apply(x_em, 1, sum)
    o <- order(x_rates*x_sum, decreasing=T)

    x_em <- matrix(x_em[o,], nrow=nrow(x_em))
    x_rates <- x_rates[o] 
    x_sum <- x_sum[o]

    
    # construct vector containing fraction of rate explained by pathway (with rate)
    df <- data.frame(id=rep(0, nrow(x_em)), rate=x_rates, 
                     exp_r=rep(0, nrow(x_em)), exp_d=rep(0, nrow(x_em)))

    for(l in 1:nrow(x_em)) {
        # 
        df$exp_r[l] <- sum(x_em[l,]*x_rates[l])/sum(rates_rev)
        # Next line only works for hv driven system in steady state... TODO fix
        df$exp_d[l] <- x_em[l,length(rates_rev)]*x_rates[l] / rates_rev[length(rates_rev)]

        em <- x_em[l,1:ncol(em_matrix)]  # cut current elementary mode (to width of em_matrix) 
        em <- em*(sign(v))   # transform to original network (net)

        m <- which(apply(em_matrix, 1, identical, em))
        if(length(m) == 0) {    # add em to em_matrix
            em_matrix <- rbind(em_matrix, em)
            df$id[l] <- nrow(em_matrix)
        } else 
            df$id[l] <- m[1]
    }

    df$exp_r_cum <- cumsum(df$exp_r)
    df$exp_d_acc <- cumsum(df$exp_d)

    # TODO pa_decompose (in contrary to pa_analysis) does not calculate errors
    # (thus err_rel_max and err are NA    
    return(list(em_ex=df, em_matrix=em_matrix, em_no=nrow(x_em), 
                err_rel_max=NA, err=NA))
}



# Function conducts a broad elementary mode analysis for results object generated
# from sb_collect_results function. All simulations were done with same network
# and same boundary species but directions of reactions may differ. Thus the set
# of elementary modes differs. But some elementary modes may be present in different
# simulations. To track the change in occurence + rates of these modes / pathways
# a global list / matrix is created. In this list coefficients of reactions in pathways are
# negative if the dynamic of the network in that the pathway was found favours the 
# backward direction.
#
# TODO: check adaption

sb_em_analysis_ecol <- function(res, res_nets, c_max=4, do_precise=T) {
    # Elementary mode structure depends on network (res_nets) which can have different
    # sizes if multiple networks were drawn randomly. Thus there is a extended
    # (pseudoreactions added) network and a list of elementary modes for each network
    # in res_nets.

    net_ext <- list()
    em_matrix <- list()
    ex_state_hashes <- list()
    N_sp <- c()
    N_re <- c()

    for(net in res_nets) {
        net_x <- pa_extend_net(net, rep(0,nrow(net[[2]])), 0, T)[[1]]
        net_ext[[length(net_ext)+1]] <- net_x 
        em_matrix[[length(em_matrix)+1]] <- matrix(0, ncol=nrow(net_x[[2]]), nrow=0)
        ex_state_hashes[[length(ex_state_hashes)+1]] <- character()
        N_sp <- c(N_sp, nrow(net[[1]]))
        N_re <- c(N_sp, nrow(net[[1]]))
    }
    
    # list of data frames for species specific, reaction specific and elementary mode expansion data
    sp <- res$sp_df
    re <- res$re_df
    em_ex <- list()

    # 
    C_matrix <- matrix(0, ncol=c_max, nrow=nrow(res))
    C_sum <- rep(0, nrow(res))
    sh_p <- rep(0, nrow(res))
    sh_p_d <- rep(0, nrow(res))
    ep_l <- rep(0,nrow(res))
    ep_nl <- rep(0,nrow(res))
    err_rel_max <- rep(0,nrow(res))
    err <- rep(0,nrow(res))
    em_no <- rep(0, nrow(res))

 
    # helper function one, 
    add_em <- function(rev, cutoff, i_net) {
        if(is.na(cutoff))
            cutoff <- 0
        #print(res_nets[[i_net]])
        cat("length(ex_state_hashes)= ", length(ex_state_hashes[[i_net]]), "\n")

        r_ext <- c(rev, jrnf_calculate_concentration_change(res_nets[[i_net]], rev))    
        x_hash <- sb_v_to_hash_s(r_ext, cutoff) 
        if(!(x_hash %in% ex_state_hashes[[i_net]])) {
            cat("added hash: ", x_hash, "\n")
            ex_state_hashes[[i_net]] <<- c(ex_state_hashes[[i_net]], x_hash)

            net_rev <- jrnf_reverse_reactions(net_ext[[i_net]], r_ext)
            rates_rev <- abs(r_ext)        

            if(do_precise)
                x_em <- pa_decompose(jrnf_calculate_stoich_mat(net_rev), rates_rev, branch_all=T)
            else
                x_em <- pa_analysis(net_rev, rates_rev, 0.05, 0.5, T, F)[[1]]
            
            for(l in 1:nrow(x_em))
                x_em[l,] <- x_em[l,]*(sign(r_ext))
                
            em_matrix[[i_net]] <<- rbind(em_matrix[[i_net]], x_em)

            x_keep <- !duplicated(matrix(em_matrix[[i_net]],nrow=nrow(em_matrix[[i_net]])))
            em_matrix[[i_net]] <<- matrix(em_matrix[[i_net]][x_keep,], ncol=ncol(em_matrix[[i_net]]))   
        }
    }
    



    # helper function
    # calculating elementary modes
    # TODO Extract all general functionality of this method to sb_h_calc_em_info
    #      (that is also used by sb_stability_analysis_eco)
    calc_em_info <- function(rev, cutoff, i_net) {
        if(is.na(cutoff))
            cutoff <- 0
        r_ext <- c(rev, jrnf_calculate_concentration_change(res_nets[[i_net]], rev))    

        net_rev <- jrnf_reverse_reactions(res_nets[[i_net]], r_ext)
        rates_rev <- abs(r_ext)  
        
        M <- em_matrix[[i_net]]

        comp_em <- !as.logical((M > 0) %*% (r_ext < 0) | (M < 0) %*% (r_ext > 0))
          
        abs_M_comp <- abs(matrix(M[comp_em,], ncol=ncol(M)))
        x <- pa_calc_coefficients(abs_M_comp, rates_rev, F)   
        coef <- x$coef
        err_rel_max[i] <<- x$score_b   # TODO maybe calculate score only from non-pseudo reactions?

        abs_sum <- apply(matrix(abs_M_comp[,1:length(rev)], nrow=nrow(abs_M_comp)), 1, sum)
        o <- order(coef*abs_sum, decreasing=T)

        em_id <- which(comp_em)[o]
        exp_r <- (coef*abs_sum)[o]/sum(rates_rev[1:length(rev)])
        exp_d <- (coef*abs_M_comp[,ncol(abs_M_comp)])[o]/abs(rates_rev[length(rates_rev)])
        
        em_no[i] <<- sum(coef != 0)
        err[i] <<- 1-sum(exp_r)

        return(data.frame(id=em_id, rate=coef[o], exp_r=exp_r, exp_r_cum=cumsum(exp_r), 
                          exp_d=exp_d, exp_d_acc=cumsum(exp_d)))
    }


    # PRE LOOP
    # (calculate extended rate vector and it's state hash - using err_cc as cutoff)
    # (if state hash hasn't been decomposed before, decompose current one and add
    #  found pathways to the pathway list)

    for(i in 1:nrow(res)) {
        cat("============================================================\n")
        cat(" i=", i, "  v=", res$v[i], "  c=", res$c[i], "\n")  
        add_em(res$re_df[[i]]$flow_effective, res$err_cc[i], res$Edraw[i])
    }

    # MAIN LOOP

    for(i in 1:nrow(res)) {
        cat("============================================================\n")
        cat(" i=", i, "  v=", res$v[i], "  c=", res$c[i], "\n")  
        net <- res_nets[[res$Edraw[i]]]
        # calculating directed network
        gc()
        g <- jrnf_to_directed_network_d(net, res$re_df[[i]]$flow_effective)
        cat(".")         

        # which reactions of network are linear
        is_linear <- (unlist(lapply(net[[2]]$educts_mul, sum)) == 1) 
        # species degree
        deg <- as.vector(degree(jrnf_to_undirected_network(net)))    
        
        # build data frame refering to expansion of simulation's rate by elementary modes
        if(!res$relaxing_sim[i] & (res$err_cc_rel[i] < 0.1 | is.na(res$err_cc_rel[i])))
            em_ex[[i]] <- calc_em_info(res$re_df[[i]]$flow_effective, res$err_cc[i], res$Edraw[i])
        else
            em_ex[[i]] <- NA

        # calculation of cycles 
        for(j in 1:c_max) 
            C_matrix[i,j] <- get_n_cycles_directed(g,j)[[1]]

        # Entropy production of linear part and of nonlinear part
        ep_l[i] <- sum(res$re_df[[i]]$entropy_prod[is_linear])
        ep_nl[i] <- sum(res$re_df[[i]]$entropy_prod[!is_linear])

        # reaction specific data
        re[[i]]$is_linear <- is_linear

        # species specific data
        sp[[i]]$degree <- deg
        sp[[i]]$mu <- net[[1]]$energy + log(res$sp_df[[i]]$con)
        sp[[i]]$mean_ep <- associate_reaction_to_species(net, res$sp_df[[i]]$entropy_prod)
    }

    
    # Copying back from local objects to <res> object
    res$ep_linear <- ep_l
    res$ep_nonlinear <- ep_nl

    res$re_df <- re
    res$sp_df <- sp
    res$em_ex <- em_ex

    for(k in 1:c_max) 
        res[paste("C", as.character(k), sep="")] <- C_matrix[,k]

    res$C_sum <- apply(C_matrix, 1, sum)

    res$em_no <- em_no
    res$em_err_max <- err_rel_max
    res$em_err_sum <- err

    results_em <- res
    em_m <- em_matrix
    save(results_em, em_m, file="results_em.Rdata")

    return(list(res, em_matrix))
}


#
#
# TODO the function has to be changed to the new convention. We expects a individual network
#      for each Edraw in a list called <res_nets>
#
# deleted em_exp_r_cross because it is not compatible with having different networks data inside 
# one results object (res_nets). One has to subset res_nets for a certain network / Edraw first.

sb_em_cross_analysis_ecol <- function(res_nets, em_m, res_em) {
    em_der <- list()           # 

    for(i in 1:length(res_nets)) {
        cat("i=", i, "\n")
        N <- jrnf_calculate_stoich_mat(res_nets[[i]])
        if(nrow(em_m[[i]]) != 0)
            em_der[[i]] <- pa_em_derive(em_m[[i]][,1:ncol(N)], res_nets[[i]], 3, ignore_hv=T)
        else
            em_der[[i]] <- c()
    }

    # Over all explained fraction of rates (quality of pathway expansion)
    res_em$exp_r_sum <- rep(0, nrow(res_em))         
    # Explained fraction of most important pathway
    res_em$exp_r_max <- rep(0, nrow(res_em))        
    # Estimation on number of pathway necessary to explain 90% of rates
    res_em$need_em_90 <- rep(0, nrow(res_em))       

    # Cycle number (weighted number of cycing - through pathways)
    res_em$cycles_r <- rep(0, nrow(res_em))
    # Weighted number of reactions (in pathways)
    res_em$reactions_r <- rep(0, nrow(res_em))
    # Weighted number of species (participating in Pathways)
    res_em$species_r <- rep(0, nrow(res_em))
    # Information entropy in pathway distribution (explained fractions)			
    res_em$informationE <- rep(0, nrow(res_em))
    # Fraction of pathways by explained fraction that are in steady state			
    res_em$em_steadystate_r <- rep(0, nrow(res_em))
    # Fraction of pathways that is changing concentration (not steady state)		 
    res_em$em_con_ch_r <- rep(0, nrow(res_em))

    for(i in 1:nrow(res_em)) 
        if(is.data.frame(res_em$em_ex[[i]])) {
            i_net <- res_em$Edraw[i]
            res_em$exp_r_sum[i] <- sum(res_em$em_ex[[i]]$exp_r)
            res_em$exp_r_max[i] <- max(res_em$em_ex[[i]]$exp_r)

            xp_r <- sort(res_em$em_ex[[i]]$exp_r, decreasing=T)
            xp <- which(cumsum(xp_r) > 0.9)
            if(length(xp) == 0)
                res_em$need_em_90[i] <- NA
            else
                res_em$need_em_90[i] <- xp[1]

            # normalize
            id <- res_em$em_ex[[i]]$id
            r <- res_em$em_ex[[i]]$exp_r
            res_em$cycles_r[i] <- sum(em_der[[i_net]]$C_s_sum[id]*r)
            res_em$reactions_r[i] <- sum(em_der[[i_net]]$Re_s[id]*r)
            res_em$species_r[i] <-  sum(em_der[[i_net]]$Sp_no[id]*r)
            res_em$informationE[i] <- -sum(r*log(r), na.rm=T)

            res_em$em_steadystate_r[i] <- sum(r[(em_der[[i_net]]$Ex_s[id] == 0)])
            res_em$em_con_ch_r[i] <- sum(r[(em_der[[i_net]]$Ex_s[id] != 0)])
        }

    results=results_em_cross <- res_em
    save(results_em_cross, em_der, file="results_em_cross.Rdata")
    return(list(results_em_cross=results_em_cross, em_derive=em_der))
}


# 
#
#

sb_cross_analysis_ecol <- function(res_nets, res) {
    N <- jrnf_calculate_stoich_mat(res_nets[[1]])
    # species degree
    deg <- as.vector(degree(jrnf_to_undirected_network(res_nets[[1]])))   

    res$flow_energy <- rep(0, nrow(res))      # flow (of energy to hv)
    res$energy <- rep(0, nrow(res))           # complete energy in system (except of hv species)
    res$diseq <- rep(0, nrow(res))            # difference between energy of this SS and eq. SS
    res$core_sp_no <- rep(NA, nrow(res))      # number of species in core (species with mu higher 
                                              # than the reference state - eq. SS)
    res$core_hash <- rep(NA, nrow(res))       # hash of vector indicating core species
    res$core_diseq_f <- rep(NA, nrow(res))    # fraction of disequilibrium contained in network core
    res$core_mass_f <- rep(NA, nrow(res))     # fraction of mass contained in network core
    res$equilibrium_ref <- rep(NA, nrow(res)) # id of equivalent equilibrium simulation for out of equilibrium sim.
    
    # Matrix containing chemical potential difference for every reaction
    dmu_cross <- matrix(0, ncol=ncol(N), nrow=nrow(res))

    # First calculate 'equilibrium_ref' for non-equilibrium rows and all the data for
    # the equilibrium rows
    for(i in 1:nrow(res)) {
        net <- res_nets[[res$Edraw[i]]]
        cat("-")
        if(i %% 1000 == 0) cat("\n", i/nrow(res))
        res$sp_df[[i]]$degree <- deg
        res$sp_df[[i]]$mu <- net[[1]]$energy + log(res$sp_df[[i]]$con)
        res$sp_df[[i]]$mu[!is.finite(res$sp_df[[i]]$mu)] <- 0
        res$sp_df[[i]]$mean_ep <- associate_reaction_to_species(net, res$re_df[[i]]$entropy_prod)

        con <- res$sp_df[[i]]$con
        mu <- res$sp_df[[i]]$mu
        x <- con*mu
           
        res$energy[i] <- sum(x[-length(x)])

        if(res$relaxing_sim[i]) {     # If it is a relaxing simulation calculate energy
            con <- res$sp_df[[i]]$con
            mu <- res$sp_df[[i]]$mu
            x <- con*mu 
            res$energy[i] <- sum(x[-length(x)])
            res$core_diseq_f[i] <- 1
            res$core_mass_f[i] <- 1
        } else {
            s <- which(res$relaxing_sim & res$Edraw == res$Edraw[i] &
                       res$Rdraw == res$Rdraw[i] & res$c == res$c[i])
            if(length(s) == 0)
                s <- which(res$relaxing_sim & res$Edraw == res$Edraw[i] &
                           res$c == res$c[i])[1]

            res$equilibrium_ref[i] <- s
        }
    }
     
    for(i in which(!res$relaxing_sim)) { 
        cat(".")
        if(i %% 1000 == 0) cat("\n", i/nrow(res))
        dmu_cross[i,] <- res$sp_df[[i]]$mu %*% N
 
        # load data on mu (chemical potential) and con (concentration)
        eqref <- res$equilibrium_ref[i]
        con <- res$sp_df[[i]]$con
        mu <- res$sp_df[[i]]$mu 
        x <- con*mu
        x2 <- res$sp_df[[eqref]]$con*res$sp_df[[eqref]]$mu
        
        core_sp <- mu > res$sp_df[[ eqref ]]$mu
        core_sp[length(core_sp)] <- F    # Ensure 'hv' is not in core
        res$sp_df[[i]]$core_sp <- core_sp

        res$flow_energy[i] <- res$flow[i]*mu[length(mu)]
        res$diseq[i] <- res$energy[i] - res$energy[eqref]
        res$core_sp_no[i] <- sum(core_sp) 
                                             # than the reference state - eq. SS)
        res$core_hash[i] <-  sb_v_to_hash(core_sp, 0.1)
        res$core_diseq_f[i] <- sum(x[core_sp]-x2[core_sp])/res$diseq[i]
        res$core_mass_f[i] <- sum(con[core_sp])/sum(con[-length(con)])
    }
     
    results_cross <- res

    save(results_cross, dmu_cross, file="results_cross.Rdata")
    return(list(results_cross=results_cross, dmu_cross=dmu_cross))
}
