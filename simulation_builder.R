# author: jakob fischer (jakob@automorph.info)
# description: 
# Script for generating files and directory structure of artificial chemistry
# simulations with thermodynamic constraints. (simulation builder = "sb_")
# This file contains the core functionality to build and evaluate wide ranging
# simulations - entire directories with multiple changing parameters. The 
# functions to plot (and partially evaluate) the results is in "simulation_builder_plot.R".
#
# Parts of this this file / modules functionality are shared with "netodeint_control.R", this
# could be unified in future.

sourced_simulation_builder <- T

if(!exists("sourced_simulation_builder_plot"))
    source("simulation_builder_plot.R")

if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")

if(!exists("sourced_art_ecosystem_gen"))
    source("art_ecosystem_gen.R")

if(!exists("sourced_pathway_analysis"))
    source("pathway_analysis.R")


# Given a steady state vector <v> and an error range <zero_range> the function
# calculates a hash encoding the reaction directions of the steady state vector.
# Sign of each reaction can be negative, positive or 
# zero (-<zero_range>,+<zero_range>). The hash should be unique as long as the 
# number of reactions (size of steady state vector) is below 60

sb_v_to_hash <- function(v, zero_range) {
    k <- 0

    for(i in 1:length(v)) {
        if(v[i] < -zero_range)
            j <- 0
        else if(v[i] > zero_range)
            j <- 1
        else 
            j <- 2

        k <- k + j * (3 ** (i-1))
    }

    return(k)
}



assemble_em_from_flow <- function(results_ss, em, net) {
    df <- data.frame(delta_b=numeric(), flow=numeric(), C_sum=numeric(), em_id=numeric(), em_exp_r=numeric(), 
                     em_exp_d=numeric(), C_f_r=numeric(), C_f_d=numeric(), ep_lin=numeric(), ep_nonl=numeric())

    der <- pa_em_derive(em[,1:nrow(net[[2]])], net)

    for(i in 1:nrow(results_ss)) {
        em_ex <- results_ss$em_ex[[i]]
        N <- nrow(em_ex)
        C_f_r <- 0
        C_f_d <- 0

        if(nrow(em_ex) != 0) {
            C_f_r <- sum(em_ex$exp_r * der$C_sum[em_ex$id])
            C_f_d <- sum(em_ex$exp_d * der$C_sum[em_ex$id])
        }

        df <- rbind(df,
                    data.frame(delta_v= as.numeric(rep(max(results_ss$v1[i], results_ss$v2[i]), N)),
                               flow=as.numeric(rep(results_ss$flow[i], N)),
                               C_sum=as.numeric(rep(results_ss$C_sum[i], N)),
                               em_id=as.numeric(em_ex$id),
                               em_exp_r=as.numeric(em_ex$exp_r), 
                               em_exp_d=as.numeric(em_ex$exp_d),
                               C_f_r=as.numeric(rep(C_f_r, N)),
                               C_f_d=as.numeric(rep(C_f_d, N)),
                               ep_lin=as.numeric(rep(results_ss$ep_linear[i]/(results_ss$ep_linear[i]+results_ss$ep_nonlinear[i]), N)),                           ep_nonl=as.numeric(rep(results_ss$ep_nonlinear[i]/(results_ss$ep_linear[i]+results_ss$ep_nonlinear[i]), N))))

    }
  
    return(df)
}




sb_extend_results <- function(results) {
    state_hash <- rep(0, nrow(results))
    no_states <- rep(0, nrow(results))

    for(i in 1:nrow(results))
        state_hash[i] <- sb_v_to_hash(results$re_df[[i]]$flow_effective, results$flow[i]/1e5)

    for(ed in unique(results$Edraw)) {
        sel <- which(results$Edraw == ed)
        count <- length(unique(state_hash[sel]))
        no_states[sel] <- count
    }

    results$state_hash <- state_hash
    results$no_states <- no_states
    return(results)
}



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


sb_generator <- function(netfile, bvalues_l, bids, N_energies, N_runs, odeint_p="~/apps/jrnf_int", Tmax=10000, deltaT=1) {
    # Save old and set new working directory
    scripts <- as.character()        # Vector of script entries / odeint_rnet calls
    path_old <- getwd()
    path <- unlist(strsplit(netfile, "/", fixed=TRUE))
    cat("path=", path, "\n")
    cat("newpath=", paste(path[-length(path)], collapse="/"), "\n")
    if(length(path) > 1)
        setwd(paste(path[-length(path)], collapse="/"))

    cat("loading network\n")

    # load network / calculate topologic properties
    net <- jrnf_read(path[length(path)])

    # make net reversible
    net[[2]]$reversible <- rep(T, nrow(net[[2]]))

    cat("topological analysis \n")
    jrnf_create_pnn_file(net, "pfile.csv", "nfile.csv")
    pfile <- read.csv("pfile.csv")
  
    # sample <N_energies> different sets of energies
    for(i in 1:N_energies) {
        # create new directory and enter
        system(paste("mkdir ", i, sep=""))
        setwd(as.character(i))
   

        # sample energies and save to net_energies.jrnf
        net <- jrnf_sample_energies(net)

        # as we are having only one set of boundary species we are setting their 
        # boundary energies to zero first to save then

        net[[1]]$energy[bids[[1]]] <- 0
        net[[1]]$energy[bids[[2]]] <- 0
        net <- jrnf_calc_reaction_r(net, 1)
     
        cat("writing energies in netfile.\n")
        jrnf_write("net_energies.jrnf", net)


        # Inner loop (create simulations with different boundary values

        for(v in 1:length(bvalues_l)) {
            ff <- paste(c("v", as.character(bvalues_l[[v]][1]), "_" ,as.character(bvalues_l[[v]][2])), collapse="")
            system(paste("mkdir", ff)) 
            setwd(ff) 
       
            for(j in 1:N_runs) {
                jrnf_create_initial(net, paste(j, ".con", sep=""), "../net.jrnf", bc_id=bids, bc_v=bvalues_l[[v]])

                scripts <- c(scripts, 
                             paste(odeint_p, " simulate solve_implicit net=", i, "/net.jrnf con=", i, "/", ff, "/", j, ".con deltaT=", as.character(deltaT), " Tmax=", as.character(Tmax), sep=""))
            }

            setwd("..")
        }


        # back to start directory
        setwd(path_old)
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

sb_collect_results <- function(b1, b2) {
    save_wd <- getwd()   # just to be save
    df <- data.frame(Edraw=numeric(), Rdraw=numeric(), b1=numeric(), b2=numeric(), 
                     v1=numeric(), v2=numeric(), flow1=numeric(), flow2=numeric(), flow=numeric(),
                     err_cc=numeric(), err_cc_rel=numeric(),
                     ep_max=numeric(), sp_df=I(list()), re_df=I(list()), 
                     last_time=numeric(), last_msd=numeric())

    Edir_v <- list.dirs(recursive=FALSE)


    for(bp in Edir_v) {
        setwd(bp)
        net <- jrnf_read("net_energies.jrnf")
        bp <- as.numeric(strsplit(bp,"/")[[1]][2])

        vdir_v <- list.dirs(recursive=FALSE)
        for(vp in vdir_v) {
            l <- strsplit(vp, "v")[[1]][2]
            v1 <- as.numeric(strsplit(l, "_")[[1]][1])
            v2 <- as.numeric(strsplit(l, "_")[[1]][2])
            setwd(vp)

            edir_v <- list.files(recursive=FALSE)
            for(ep in edir_v) {
                i <- as.numeric(strsplit(ep,"\\.")[[1]][1])
                
                cat("doing b1=", b1, "b2=", b2, "v1=", v1, "v2=", v2, "Edraw=", as.numeric(bp), "Rdraw=", i, "\n")
                run <- read.csv(ep)

                last_time_ <- run[nrow(run),1]
                last_msd_  <- run[nrow(run),2]
                last_con <- data.frame(con=as.numeric(run[nrow(run), 3:ncol(run)]))

               
                last_flow <- jrnf_calculate_flow(net, last_con$con)

                cc <- jrnf_calculate_concentration_change(net, last_flow$flow_effective)
                flowb1 <- cc[b1]
                flowb2 <- cc[b2]
                err_cc <- max(abs(cc[c(-b1,-b2)]))
                flow <- 0.5*abs(flowb1-flowb2)
                err_cc_rel <- err_cc/flow

                #cat("flowb1 = ", flowb1, "\n")
                #cat("flowb2 = ", flowb2, "\n")

                df <- rbind(df,
                              data.frame(Edraw=as.numeric(bp), Rdraw=as.numeric(i), b1=as.numeric(b1), 
                              b2=as.numeric(b2), v1=as.numeric(v1), v2=as.numeric(v2),
                              flow1=as.numeric(flowb1), flow2=as.numeric(flowb2), flow=flow, err_cc=err_cc,
                              err_cc_rel=err_cc_rel,
                              ep_tot=as.numeric(sum(last_flow$entropy_prod)), 
                              sp_df=I(list(last_con)), re_df=I(list(last_flow)), last_time=last_time_, last_msd=last_msd_))
            }
 
            setwd("..")
        }

        setwd("..")
    }

    setwd(save_wd)
    results <- df
    save(results, file="results.Rdata")
}



# Function conducts a broad elementary mode analysis for results object generated
# from sb_collect_results function. All simulations were done with same network
# and same boundary species but directions of reactions may differ. Thus the set
# of elementary modes differs. But some elementary modes may be present in different
# simulations. To track the change in occurence + rates of these modes / pathways
# a global list / matrix is created. In this list coefficients of reactions in pathways are
# negative if the dynamic of the network in that the pathway was found favours the 
# backward direction.

sb_em_analysis <- function(res, net, c_max=5) {
    N_sp <- length(res$sp_df[[1]]$con)
    N_re <- length(results$re_df[[1]]$flow_effective)
   
    # matrix where combined elementary modes are saved / two reactions added -> inflow / outflow pseudoreactions
    em_matrix <- matrix(0, ncol=N_re+2, nrow=0)
    
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
    is_linear <- (unlist(lapply(net[[2]]$educts_mul, sum)) == 1) 
    # species degree
    deg <- as.vector(degree(jrnf_to_undirected_network(net)))    

    # flow through system (has to be calculated before iteration because it's a parameter for pa_analysis)
    res$flow <- 0.5*abs(res$flow1-res$flow2)
    res$err_cc_rel <- res$err_cc/res$flow

    for(i in 1:nrow(res)) {
        cat("============================================================\n")
        cat(" i=", i, "  b1=", res$b1[i], "  b2=", res$b2[i], "  v1=", res$v1[i], "  v2=", res$v2[i], "\n")  

        # calculating directed network
        gc()
        g <- jrnf_to_directed_network_d(net, res$re_df[[i]]$flow_effective)
        cat(".")         

        # calculating elementary modes
        # TODO
        rev <- res$re_df[[i]]$flow_effective
        rates_rev <- abs(res$re_df[[i]]$flow_effective) 
        net_rev <- jrnf_reverse_reactions(net, rev)
      
        cdif_r <- jrnf_calculate_concentration_change(net_rev, rates_rev)

        # add reactions to balance growth / decrease
        for(k in c(res$b1[i], res$b2[i])) 
            # if species' concentration increases one pseudoreaction has to be included to remove it ("X -> ")
            if(cdif_r[k] > 0) {   
	        net_rev[[2]] <- rbind(net_rev[[2]], data.frame(reversible=factor(c(FALSE)), 
                                      c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                                      activation=as.numeric(c(0)),educts=I(list(k)), educts_mul=I(list(1)),
                                      products=I(list(c())), products_mul=I(list(c()))))
                rates_rev <- c(rates_rev, cdif_r[k])
            # if species' concentration decreases one pseudoreaction is included to add it ("-> X ")
            } else { 
	         net_rev[[2]] <- rbind(net_rev[[2]], data.frame(reversible=factor(c(FALSE)), 
                                       c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                                       activation=as.numeric(c(0)),educts=I(list(c())), educts_mul=I(list(c())),
                                       products=I(list(k)), products_mul=I(list(1))))
                rates_rev <- c(rates_rev, -cdif_r[k])
            }

        x <- pa_analysis(net_rev, rates_rev, 0, 0)

        # rename results and sort decreasing with fraction of v explained...
        x_em <- x[[1]]
        x_rates <- x[[2]]
        x_sum <- apply(x_em, 1, sum)
        o <- order(x_rates*x_sum, decreasing=T)


        x_em <- matrix(x_em[o,], ncol=ncol(x_em))
        x_rates <- x_rates[o] 
        x_sum <- x_sum[o]

        em_no[i] <- nrow(x_em)
        err_rel_max[i] <- x[[3]]
        err[i] <- x[[4]]

        # construct vector containing fraction of rate explained by pathway (with rate)
        exp_r <- rep(0, nrow(x_em))
        exp_d <- rep(0, nrow(x_em))
        

        em_id <- rep(0, nrow(x_em))

        if(nrow(x_em) > 0)
        for(l in 1:nrow(x_em)) {
            # 
            exp_r[l] <- sum(x_em[l,]*x_rates[l])/sum(rates_rev)
            exp_d[l] <- x_em[l,length(rates_rev)]*x_rates[l] / rates_rev[length(rates_rev)]

            # make copy from elementary mode
            em <- x_em[l,]   
                
            # now negate those elements (reactions) that are reverted
            em <- em*(c(sign(rev), 1, 1))

            m <- which(apply(em_matrix, 1, identical, em))
            if(length(m) == 0) {    # add em to em_matrix
                em_matrix <- rbind(em_matrix, em)
                em_id[l] <- nrow(em_matrix)
            } else {
                em_id[l] <- m[1]
            }
        }

        # build data frame refering to expansion of simulation's rate by elementary modes
        em_ex[[i]] <- data.frame(id=em_id, rate=x_rates, exp_r=exp_r, exp_r_cum=cumsum(exp_r), 
                                 exp_d=exp_d, exp_d_acc=cumsum(exp_d))

        # calculation of cycles 
        for(j in 1:c_max) 
            C_matrix[i,j] <- get_n_cycles_directed(g,j)[[1]]

        C_sum[i] <- sum(C_matrix[i,])

        
        # calculating shortest paths
        sh_p[i] <- matrix(shortest.paths(g), ncol=N_sp, nrow=N_sp)[res$b1[i], res$b2[i]]
        
        if(res$v1[i] > res$v2[i])
            sh_p_d[i] <- matrix(shortest.paths(g, mode="out"), ncol=N_sp, nrow=N_sp)[res$b1[i], res$b2[i]]
        else
            sh_p_d[i] <- matrix(shortest.paths(g, mode="out"), ncol=N_sp, nrow=N_sp)[res$b2[i], res$b1[i]]

        # Entropy production of linear part and of nonlinear part

        ep_l[i] <- sum(res$re_df[[i]]$entropy_prod[is_linear])
        ep_nl[i] <- sum(res$re_df[[i]]$entropy_prod[!is_linear])

        # reaction specific data
        re[[i]]$is_linear <- is_linear

        # species specific data
        sp[[i]]$degree <- deg
        sp[[i]]$mu <- net[[1]]$energy + log(res$sp_df[[i]]$con)
        sp[[i]]$mean_ep <- associate_reaction_to_species(net, res$sp_re[[i]]$entropy_prod)
    }

    
    # Copying back from local objects to <res> object
    res$ep_linear <- ep_l
    res$ep_nonlinear <- ep_nl

    res$re_df <- re
    res$sp_df <- sp
    res$em_ex <- em_ex
    #res <- cbind(res, em_ex)

    res$reliability <- as.numeric(abs(abs(res$flow*log(res$v1/res$v2)) - res$ep_tot)/res$ep_tot > 0.1) + 1  

    res$shortest_path <- sh_p
    res$shortest_path_d <- sh_p_d

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
# TODO: - copy and pasted from sb_general_ecol   - adapt and document
# TODO: - instead of having different energies one might have different initial conditions


sb_generator_ecol <- function(netfile, bvalues, cvalues, N_energies, N_runs, 
                              odeint_p="~/apps/jrnf_int", Tmax=10000, deltaT=100) {
    # Save old and set new working directory
    scripts <- as.character()        # Vector of script entries / odeint_rnet calls
    path_old <- getwd()              # Save to restore properly
    path <- unlist(strsplit(netfile, "/", fixed=TRUE))
    cat("path=", path, "\n")
    cat("newpath=", paste(path[-length(path)], collapse="/"), "\n")
    if(length(path) > 1)
        setwd(paste(path[-length(path)], collapse="/"))

    cat("loading network\n")

    # load network / calculate topologic properties
    net <- jrnf_read(path[length(path)])
    con <- jrnf_analyze_ecosystem_constituents(net)
  
    # make net reversible
    net[[2]]$reversible <- rep(T, nrow(net[[2]]))

    cat("topological analysis \n")
    jrnf_create_pnn_file(net, "pfile.csv", "nfile.csv")
    pfile <- read.csv("pfile.csv")
    
    # Calculate reduced (no hv) stoichiometric matrix            
    N_red <- jrnf_calculate_stoich_mat(net)
    N_red <- N_red[-nrow(N_red),] 

  
    # sample <N_energies> different sets of energies
    for(i in 1:N_energies) {
        # create new directory and enter
        system(paste("mkdir ", i, sep=""))
        setwd(as.character(i))
        
         net <- jrnf_calc_reaction_r(net, 1)
     
        cat("writing energies in netfile.\n")
        jrnf_write("net_energies.jrnf", net)


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
                             paste(odeint_p, " simulate solve_implicit net=", i, "/net_energies.jrnf con=", i, "/", ff, "/", j, ".con deltaT=", as.character(deltaT), " Tmax=", as.character(Tmax), " wint=500", sep=""))
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

sb_collect_results_ecol <- function() {
    save_wd <- getwd()   # just to be save
    df <- data.frame(Edraw=numeric(), Rdraw=numeric(), c=numeric(), v=numeric(), 
                     flow=numeric(), ep_max=numeric(), sp_df=I(list()), 
                     re_df=I(list()), last_time=numeric(), last_msd=numeric())


    Edir_v <- list.dirs(recursive=FALSE)

    for(bp in Edir_v) {
        setwd(bp)
        net <- jrnf_read("net_energies.jrnf")
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
                
                last_flow <- jrnf_calculate_flow(net, last_con$con)
                #flowb <- calculate_flow_dif(nrow(net[[1]]), net, last_flow$flow_effective)
            
                cc <- jrnf_calculate_concentration_change(net, last_flow$flow_effective)
                flowb <- -cc[nrow(net[[1]])]
                err_cc <- max(abs(cc[-nrow(net[[1]])]))
                err_cc_rel <- err_cc/flowb

                df <- rbind(df,
                            data.frame(Edraw=as.numeric(bp),Rdraw=as.numeric(i), 
                                       v=as.numeric(v), c=as.numeric(c), flow=as.numeric(flowb),  
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
    save(results, file="results.Rdata")
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

sb_em_analysis_ecol <- function(res, net, c_max=5) {
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
    is_linear <- (unlist(lapply(net[[2]]$educts_mul, sum)) == 1) 
    # species degree
    deg <- as.vector(degree(jrnf_to_undirected_network(net)))    


    for(i in 1:nrow(res)) {
        cat("============================================================\n")
        cat(" i=", i, "  v=", res$v[i], "  c=", res$c[i], "\n")  

        # calculating directed network
        gc()
        g <- jrnf_to_directed_network_d(net, res$re_df[[i]]$flow_effective)
        cat(".")         

        # calculating elementary modes
        # TODO
        rev <- res$re_df[[i]]$flow_effective
        rates_rev <- abs(res$re_df[[i]]$flow_effective) 
        net_rev <- jrnf_reverse_reactions(net, rev)
      
        cdif_r <- jrnf_calculate_concentration_change(net_rev, rates_rev)

        # add pseudoreaction generating hv
        n_hv <- nrow(net[[1]])     # species id of hv
        net_rev[[2]] <- rbind(net_rev[[2]], data.frame(reversible=factor(c(FALSE)), 
                                       c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                                       activation=as.numeric(c(0)),educts=I(list(c())), educts_mul=I(list(c())),
                                       products=I(list(n_hv)), products_mul=I(list(1))))
        rates_rev <- c(rates_rev, -cdif_r[n_hv])

        x <- pa_analysis(net_rev, rates_rev, 0, 0)

        # rename results and sort decreasing with fraction of v explained...
        x_em <- x[[1]]
        x_rates <- x[[2]]
        x_sum <- apply(x_em, 1, sum)
        o <- order(x_rates*x_sum, decreasing=T)

        x_em <- matrix(x_em[o,], nrow=nrow(x_em))
        x_rates <- x_rates[o] 
        x_sum <- x_sum[o]

        em_no[i] <- nrow(x_em)
        err_rel_max[i] <- x[[3]]
        err[i] <- x[[4]]

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
                em_matrix <- rbind(em_matrix, em)
                em_id[l] <- nrow(em_matrix)
            } else {
                em_id[l] <- m[1]
            }
        }

        # build data frame refering to expansion of simulation's rate by elementary modes
        em_ex[[i]] <- data.frame(id=em_id, rate=x_rates, exp_r=exp_r, exp_r_cum=cumsum(exp_r), 
                                 exp_d=exp_d, exp_d_acc=cumsum(exp_d))

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
        sp[[i]]$mean_ep <- associate_reaction_to_species(net, res$sp_re[[i]]$entropy_prod)
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



sb_em_cross_analysis <- function(net, em, res_em) {
    N <- jrnf_calculate_stoich_mat(net)
    em_der <- pa_em_derive(em[,1:ncol(N)], net, 4)

    em_exp_r_cross <- matrix(0, ncol=nrow(em), nrow=nrow(res_em)) 
    exp_r_sum <- rep(0, nrow(res_em))
    exp_r_max <- rep(0, nrow(res_em))
    need_em_90 <- rep(0, nrow(res_em))

    cycles_r <- rep(0, nrow(res_em))
    reactions_r <- rep(0, nrow(res_em))
    species_r <- rep(0, nrow(res_em))
    informationE <- rep(0, nrow(res_em))
    
 
    for(i in 1:nrow(res_em)) {
        em_exp_r_cross[i, res_em$em_ex[[i]]$id] <- res_em$em_ex[[i]]$exp_r
        exp_r_sum[i] <- sum(res_em$em_ex[[i]]$exp_r)
        exp_r_max[i] <- max(res_em$em_ex[[i]]$exp_r)

        xp_r <- sort(res_em$em_ex[[i]]$exp_r, decreasing=T)
        xp <- which(cumsum(xp_r) > 0.9)
        if(length(xp) == 0)
            need_em_90[i] <- NA
        else
            need_em_90[i] <- xp[i]

        # normalize
        r <- em_exp_r_cross[i,] / sum(em_exp_r_cross[i,])
        cycles_r[i] <- sum(em_der$C_s_sum*r)
        reactions_r[i] <- sum(em_der$Re_s*r)
        species_r[i] <-  sum(em_der$Sp_no*r)

        informationE[i] <- -sum(r*log(r), na.rm=T)
           
    }   

    results_cross <- res_em
    results_cross$exp_r_sum <- exp_r_sum
    results_cross$exp_r_max <- exp_r_max
    results_cross$need_em_90 <- need_em_90

    results_cross$cycles_r <- cycles_r
    results_cross$reactions_r <- reactions_r
    results_cross$species_r <- species_r
    results_cross$informationE <- informationE

    em_derive=em_der

    save(results_cross, em_exp_r_cross, em_derive, file="results_cross.Rdata")
 
    return(list(results_cross=results_cross, em_exp_r_cross=em_exp_r_cross, em_derive=em_der))
}
