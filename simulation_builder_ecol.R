# author: jakob fischer (mail@jakobfischer.eu)
# description: 
# Module for creating simulations of artificial ecosystems with a wide range of
# parameters (different Energies, different driving force, different 
# concentration). "Artificial ecosystems" reffers to artificial reaction network
# that have no matter exchange with their environment but are drivent to a
# steady state in thermodynamic disequilibrium by reactions similar to photo-
# chemical reactions. This file contains code for generating network files
# and initial concentrations (in the file system), preparing the call to the 
# ode solver and methods for collecting the results. Code for evaluation allows
# to automatically conclude a pathway analysis and can also be used for 
# results of evoveled networks generated with simulation_builder_evol.R.

sourced_simulation_builder_ecol <- T


if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")


# Function loads the network from <netfile> and returns a generator functor that 
# can be used at "netgen" parameter of "sb_generator_ecol_mul". 

sb_build_generator <- function(netfile) {
    net <- jrnf_read(netfile)
    return(function()  {  return(net)  })
} 


# Variant of sb_generator_ecol that does not load one network and assign <N_energies> 
# of energies, but is given a generator-function (<netgen>) from which <N_nets> 
# networks are generated and each assigned <N_energies> of energies. 
#
# parameters:
# <path_s>        - Path to where simulation files will be put
# <netgen>        - Function that generates a net (randomly) and returns it 
# <bvalues>       - Vector containing the values for the driving force (hv-species) 
# <N_nets>        - Number of nets that are generated using <netgen> function
# <N_energies>    - Number of energy sets that are drawn for each the network
# <N_runs>        - Number of simulations done for every network+energy set + 
#                   boundary value set
# <odeint_p>      - Path to ode integrator (jrnf_int)
# <Tmax>          - Time units up to that the networks are simulated
# <wint>          - Number of times between 0 and <Tmax> output is written to file
# <flat_energies> - If set, one time (of <N_energies> times) the network is taken 
#                   with trivial energies (mu = 0, EA = 1)
# <write_log>     - Results are written logarithmically equidistant
# <N_runs_eq>     - Number of times for each network + energies + concentration the
#                   system is relaxed without driving force for reference.
# <limit_AE>      - Maximum activation energy (is forwarded as parameter to 
#                   "jrnf_ae_draw_energies").

sb_generator_ecol_mul <- function(path_s, netgen, bvalues, cvalues, N_nets, N_energies, N_runs, 
                              odeint_p=sb_odeint_path, Tmax=10000, wint=50, 
                              flat_energies=F, write_log=T, N_runs_eq=c(), limit_AE=T) {
    # Standard value for equilibrium runs is the same as the normal runs
    if(is.null(N_runs_eq)) 
        N_runs_eq <- N_runs

    Tmax_ <- Tmax
    if(Tmax > sb_max_step_size  && write_log)
        Tmax_ = sb_max_step_size 

    # Save old and set new working directory
    scripts <- as.character()        # Vector of script entries / odeint_rnet calls
    path_old <- getwd()              # Save to restore properly
    path <- unlist(strsplit(path_s, "/", fixed=TRUE))
    cat("path=", path, "\n")
    if(length(path) > 0)
        setwd(paste(path, collapse="/"))

    # 
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

            # if no flat energies use first network with given energies!
            if(flat_energies || i > 1)  
                net <- jrnf_ae_draw_energies(net, i==1 && flat_energies, limit_AE)
            net <- jrnf_calculate_rconst(net, 1)  

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
# initial network file! Results are saved in data frame "results" and the list
# "results_nets" which are stored in the file "results.Rdata". This is the 
# redone function that was done for including all intermediate state in the 
# results object. It's results objects are not compatible with the old format!
#
# The parameter <col_dynamics> indicates if the dynamics of the simulation is
# going to be analyzed. Then a entry / row in the results is created for EACH
# time step that was written for the simulations. Else only the last entry is
# collected.

sb_collect_results_ecol <- function(col_dynamics=T) {
    save_wd <- getwd()   # just to be save
    df <- data.frame(Edraw=numeric(), Rdraw=numeric(), c=numeric(), v=numeric(), 
                     flow=numeric(), ep_max=numeric(), state_hash=character(), sp_df=I(list()), 
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
                    state_hash <- sb_v_to_hash_s(flow$flow_effective, 1e-20)

                    df <- rbind(df,
                                data.frame(Edraw=as.numeric(bp),Rdraw=as.numeric(i), 
                                           v=as.numeric(v), c=as.numeric(c), flow=as.numeric(flowb), state_hash=as.character(state_hash),
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


# Function conducts a broad elementary mode analysis for results object generated
# from sb_collect_results function. All simulations with same <Edraw> / <res_nets> 
# were done with same network and same boundary species but directions of 
# reactions may differ. Thus the set of elementary modes differs. But some 
# elementary modes may be present in different simulations. To track the change 
# in occurence + rates of these modes / pathways a list (matrix) is created for
# each network. In this matrices coefficients of reactions in pathways are
# negative if the dynamic of the network in that the pathway was found favours the 
# backward direction.
#
# parameters:
# <res_nets>   - results-nets object (list of networks)
# <res>        - results object 
# <c_max>      - Maximum cycle size for determining cycles
# <do_precies> - Calculate ALL pathways (only recommended for VERY small networks)
# <param>      - Parameters for pathway calculation (if do_precise == F)
# <param$fext> - Cutoff for explained fraction of rates / turnover
# <param$pmin> - Cutoff for combination / branching of pathways
#                (see pa_analysis in "pathway_analysis.R" for details)

sb_em_analysis_ecol <- function(res_nets, res, c_max=4, do_precise=F, param=list(fext=0.05, pmin=0.5)) {
    # Elementary mode structure depends on network (res_nets) which can have different
    # sizes if multiple networks were drawn randomly. Thus there is a extended
    # (pseudoreactions added) network and a list of elementary modes for each network
    # in res_nets.

    net_ext <- list()          # list of networks extended by in- / outflow reactions
    em_matrix <- list()        # list of pathway matrices for each extended network
    ex_state_hashes <- list()  # list of vectors of already handled states' hashes

    for(net in res_nets) {
        net_x <- pa_extend_net(net, rep(0,nrow(net[[2]])), 0, T)[[1]]
        net_ext[[length(net_ext)+1]] <- net_x 
        em_matrix[[length(em_matrix)+1]] <- matrix(0, ncol=nrow(net_x[[2]]), nrow=0)
        ex_state_hashes[[length(ex_state_hashes)+1]] <- character()
    }
    
    # List of data frames for species specific, reaction specific and elementary 
    # mode expansion data
    sp <- res$sp_df
    re <- res$re_df
    em_ex <- list()

    # Number of cycles of different length
    C_matrix <- matrix(0, ncol=c_max, nrow=nrow(res))
    # Total number of cycles
    C_sum <- rep(0, nrow(res))
    # Entropy production in linear reactions
    ep_l <- rep(0,nrow(res))
    # Entropy production in nonlinear reactions
    ep_nl <- rep(0,nrow(res))
    # Error of pathway decomposition
    err_rel_max <- rep(0,nrow(res))
    err <- rep(0,nrow(res))
    em_no <- rep(0, nrow(res))

 
    # Helper function which takes steady state vector <rev>, a <cutoff> and an id 
    # of a network <i_net> to determine through a hash if this state was already 
    # seen. If not, pathways are calculated and collected for this steady state. 
    add_em <- function(rev, cutoff, i_net) {
        if(is.na(cutoff))
            cutoff <- 0

        # create and hash steady state for extended network net_ext[[x]]
        r_ext <- c(rev, jrnf_calculate_concentration_change(res_nets[[i_net]], rev))    
        x_hash <- sb_v_to_hash_s(r_ext, cutoff) 
        if(!(x_hash %in% ex_state_hashes[[i_net]])) {
            cat("added hash: ", x_hash, "\n")
            ex_state_hashes[[i_net]] <<- c(ex_state_hashes[[i_net]], x_hash)
            # Pathway decomposition is done with all rates positive
            net_rev <- jrnf_reverse_reactions(net_ext[[i_net]], r_ext)
            rates_rev <- abs(r_ext)        

            # Calculate pathways
            if(do_precise)
                x_em <- pa_decompose(jrnf_calculate_stoich_mat(net_rev), rates_rev, branch_all=T)
            else
                x_em <- pa_analysis(net_rev, rates_rev, param$fext, param$pmin, T, F)[[1]]
             
            for(l in 1:nrow(x_em))
                x_em[l,] <- x_em[l,]*(sign(r_ext))

            # have to remove degenerate cases (if found because of numerical problems)
            not_zero = apply(x_em, 1, sum) != 0
            x_em = matrix(x_em[not_zero,], ncol=ncol(x_em))
              
            # Add found pathways to existing ones for this network... 
            em_matrix[[i_net]] <<- rbind(em_matrix[[i_net]], x_em)

            # ... and remove duplicates
            x_keep <- !duplicated(matrix(em_matrix[[i_net]],nrow=nrow(em_matrix[[i_net]])))
            em_matrix[[i_net]] <<- matrix(em_matrix[[i_net]][x_keep,], ncol=ncol(em_matrix[[i_net]]))   
        }
    }

    # Helper function that calculates the pathway decomposition specific data 
    # frames using the previously calculated pathways.
    # 
    # parameters:
    # <rev>   - Steady state vector (without exchange reactions)
    # <mu>    - Chemical potential vector 
    # <i_net> - Id of the respective network
    #
    # TODO maybe calculate score only from non-pseudo reactions?

    calc_em_info <- function(rev, mu, i_net) {
        # Calculate state vector that is compatible with extended network net_ext[[i]]
        r_ext <- c(rev, jrnf_calculate_concentration_change(res_nets[[i_net]], rev))    
        rates_rev <- abs(r_ext)  

        M <- em_matrix[[i_net]]    # shortcut

        # Find compatible pathways
        comp_em <- !as.logical((M > 0) %*% (r_ext < 0) | (M < 0) %*% (r_ext > 0))
          
        # Build subset of pathways and map to positive coefficients
        abs_M_comp <- abs(matrix(M[comp_em,], ncol=ncol(M)))
        # calculate coefficients
        x <- pa_calculate_coef(abs_M_comp, rates_rev, F)   
        coef <- x$coef
        err_rel_max[i] <<- x$score_b   # quality of pathway decomposition  

        # Order pathways with decreasing explained fractions
        abs_sum <- apply(matrix(abs_M_comp[,1:length(rev)], nrow=nrow(abs_M_comp)), 1, sum)
        o <- order(coef*abs_sum, decreasing=T)

        # Calculate explained fraction (only considering reactions of original net)
        exp_r <- (coef*abs_sum)[o]/sum(rates_rev[1:length(rev)])
 
        # calculate dissipation through coefficient of hv inflow reaction
        id_hv_ex_rea <- nrow(res_nets[[i_net]][[2]])+which(res_nets[[i_net]][[1]]$name == "hv")
        exp_d <- (coef*abs_M_comp[,id_hv_ex_rea])[o]/sum(coef*abs_M_comp[,id_hv_ex_rea])
        
        em_no[i] <<- sum(coef != 0)
        err[i] <<- 1-sum(exp_r)

        return(data.frame(id=which(comp_em)[o], rate=coef[o], exp_r=exp_r, exp_r_cum=cumsum(exp_r), 
                          exp_d=exp_d, exp_d_acc=cumsum(exp_d)))
    }


    # PRE LOOP
    # (Check if state already occured. Do pw decomposition if not.)

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
        sp[[i]]$mean_ep <- jrnf_associate_reaction_to_species(net, res$sp_df[[i]]$entropy_prod)

        # build data frame refering to expansion of simulation's rate by elementary modes
        if(!res$relaxing_sim[i] & (res$err_cc_rel[i] < 0.1 | is.na(res$err_cc_rel[i])))
            em_ex[[i]] <- calc_em_info(res$re_df[[i]]$flow_effective, sp[[i]]$mu, res$Edraw[i])
        else
            em_ex[[i]] <- NA

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

    # Save to file / return
    results_em <- res
    em_m <- em_matrix
    save(results_em, em_m, file="results_em.Rdata")
    return(list(results_em=res, em_m=em_matrix))
}


# This function can be called after pathways have been calculated (sb_em_analysis_ecol).
# It first derives further properties of the pathways (pa_em_derive) and then uses
# this properties and the coefficients of the pathway decomposition to calculate further
# properties of the respective steady states.
#
# deleted em_exp_r_cross because it is not compatible with having different networks data inside 
# one results object (res_nets). One has to subset res_nets for a certain network / Edraw first.
#
# parameters:
# <res_nets>  - results-nets object (list of networks)
# <res_em>    - results object (results_em) 
# <em_m>      - List of elementary modes / pathway matrices (em_m)
# <c_max>     - Maximum cycle size for determining cycles


sb_em_cross_analysis_ecol <- function(res_nets, res_em, em_m, c_max=4) {
    em_der <- list()           # 

    for(i in 1:length(res_nets)) {
        cat("i=", i, "\n")
        N <- jrnf_calculate_stoich_mat(res_nets[[i]])
        if(nrow(em_m[[i]]) != 0)
            em_der[[i]] <- pa_em_derive(em_m[[i]][,1:ncol(N)], res_nets[[i]], c_max, ignore_hv=T)
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
        # Makes only sense if pathway information is available for entry
        if(is.data.frame(res_em$em_ex[[i]])) {
            i_net <- res_em$Edraw[i] # shortcut to network id
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
 
    # Write back / to file
    results_em_cross <- res_em
    em_derive <- em_der
    save(results_em_cross, em_derive, file="results_em_cross.Rdata")
    return(list(results_em_cross=results_em_cross, em_derive=em_derive))
}


# This function executes a further analysis for simulation results (not related
# to found pathways). Espespecially it calculates quantities that can only be
# determined if one can compare an out of equilibrium state with a reference 
# state in thermodynamic equilibrium.
#
# parameters:
# <res_nets>  - results-nets object (list of networks)
# <res>       - results object 

sb_cross_analysis_ecol <- function(res_nets, res) {
    # species degree
    deg <- list()
    for(net in res_nets)
        deg[[length(deg)+1]] <- as.vector(degree(jrnf_to_undirected_network(res_nets[[1]])))   

    res$flow_energy_hv <- rep(0, nrow(res))   # flow (of energy to hv)
    res$flow_energy_tot <- rep(0, nrow(res))  # flow of energy
    res$energy <- rep(0, nrow(res))           # complete energy in system (except of hv species)
    res$diseq <- rep(0, nrow(res))            # difference between energy of this SS and eq. SS
    res$core_sp_no <- rep(NA, nrow(res))      # number of species in core (species with mu higher 
                                              # than the reference state - eq. SS)
    res$core_hash <- character(nrow(res))     # hash of vector indicating core species
    res$core_diseq_f <- rep(NA, nrow(res))    # fraction of disequilibrium contained in network core
    res$core_mass_f <- rep(NA, nrow(res))     # fraction of mass contained in network core
    res$equilibrium_ref <- rep(NA, nrow(res)) # id of equivalent equilibrium simulation for out of equilibrium sim.

    # First calculate some additional data for all rows (energy flow). Then 
    # calculate 'equilibrium_ref' for non-equilibrium rows and all the data for
    # the equilibrium rows.
    for(i in 1:nrow(res)) {
        net <- res_nets[[res$Edraw[i]]]
        # Function should not take too long to run compared to others, but want
        # to have some user feedback.
        cat("-")
        if(i %% 1000 == 0) cat("\n", i/nrow(res))

        # calculate additional species specific data
        res$sp_df[[i]]$degree <- deg[[ res$Edraw[i] ]]
        res$sp_df[[i]]$mu <- net[[1]]$energy + log(res$sp_df[[i]]$con)
        res$sp_df[[i]]$mu[!is.finite(res$sp_df[[i]]$mu)] <- 0
        res$sp_df[[i]]$mean_ep <- jrnf_associate_reaction_to_species(net, res$re_df[[i]]$entropy_prod)

        # shortcuts
        con <- res$sp_df[[i]]$con
        mu <- res$sp_df[[i]]$mu
        x <- con*mu
        hv_sel <- "hv" == res_nets[[res$Edraw[i]]][[1]]$name   
        res$energy[i] <- sum(x[!hv_sel])
 
        # calculate concentration change from flow and derive energy flow
        cc <- jrnf_calculate_concentration_change(net, res$re_df[[i]]$flow_effective)
        res$flow_energy_hv[i] <- sum(mu[hv_sel]*res$flow[i])
        res$flow_energy_tot[i] <- -sum(mu*cc)
        
        if(res$relaxing_sim[i]) {     
            # If it is a relaxing simulation disequilibrium and mass fraction are 
            # 1 by definition.
            res$core_diseq_f[i] <- 1
            res$core_mass_f[i] <- 1
        } else {
            # Try finding best fit for equilibrium reference by loosening 
            # constraints step by step.
            s <- which(res$relaxing_sim & res$Edraw == res$Edraw[i] &
                       res$Rdraw == res$Rdraw[i] & res$c == res$c[i] &
                       res$is_last)
            if(length(s) == 0)
                s <- which(res$relaxing_sim & res$Edraw == res$Edraw[i] &
                           res$c == res$c[i] & res$is_last)[1]

            if(length(s) == 0)
                s <- 0

            res$equilibrium_ref[i] <- s
        }
    }
     
    # Now calculate all the data for out of equilibrium states which have an
    # equilibrium reference state (simulation).
    for(i in which(!res$relaxing_sim & res$equilibrium_ref != 0)) { 
        cat(".")
        if(i %% 1000 == 0) cat("\n", i/nrow(res))
 
        # load data on mu (chemical potential) and con (concentration)
        eqref <- res$equilibrium_ref[i]
        con <- res$sp_df[[i]]$con
        mu <- res$sp_df[[i]]$mu 
        x <- con*mu
        x2 <- res$sp_df[[eqref]]$con*res$sp_df[[eqref]]$mu
        hv_sel <- "hv" == res_nets[[res$Edraw[i]]][[1]]$name
        
        # "Core" species are those that have higher energy than in their
        # respective refference state.     
        core_sp <- (mu > res$sp_df[[ eqref ]]$mu) & !hv_sel
        res$sp_df[[i]]$core_sp <- core_sp

        # Disequilibrium is distance of core species from reference state in
        # terms of energy
        res$diseq[i] <- res$energy[i] - res$energy[eqref]
        res$core_sp_no[i] <- sum(core_sp) 
                          
        # Hash core species and calculate fraction of mass + dissequilibrium in core.                  
        res$core_hash[i] <-  sb_v_to_hash_s(core_sp, 0.1)
        res$core_diseq_f[i] <- sum(x[core_sp]-x2[core_sp])/res$diseq[i]
        res$core_mass_f[i] <- sum(con[core_sp])/sum(con[!hv_sel])
    }
     
    results_cross <- res

    save(results_cross, file="results_cross.Rdata")
    return(list(results_cross=results_cross))
}
