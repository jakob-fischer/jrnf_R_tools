# author: jakob fischer (jakob@automorph.info)
# description: 
# Legacy code for functions in simulation_builder_ecol.R. If those functions
# are thoroughly checked all this (file) can be removed.

sourced_simulation_builder_ecol_legacy <- T



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
