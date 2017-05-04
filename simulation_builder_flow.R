# author: jakob fischer (mail@jakobfischer.eu)
# description: 
# Similar to "simulation_builder_ecol.R" this module generates a gib number of
# simulations in the file system and evaluates them after they have been 
# executed. But the systems generated from this module here are driven by flow
# of matter instead of photochemical reactions.

sourced_simulation_builder_flow <- T


# Generator loads network from <netfile> and generates <N_energies> sets of energies 
# defining dynamics. For each of these networks (with energies) and each of the boundary
# value sets <bvalues_l> between the boundary species <bids> initial condition for 
# <N_runs> runs are generated and a line for calculation added in the script file <run.sh>
#
# parameters:
# <netfile>     - Path to network file
# <bvalues_l>   - List of two element vectors containing the boundary values 
# <bids>        - Two element vector with ids of the boundary species
# <N_energies>  - Number of energy sets that are drawn for the network
# <N_runs>      - Number of runs done for every energy set + boundary value set
# <odeint_p>    - Path to ode integrator (jrnf_int)
# <Tmax>        - Time up to which to ode is solved.
# <wint>        - Number of writes up to the <Tmax>.
# <deltaT>      - Initial time step for integrator.

sb_generator <- function(netfile, bvalues_l, bids, N_energies, N_runs,
                         odeint_p=sb_odeint_path, Tmax=10000, wint=50, deltaT=0.1) {
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
        net <- jrnf_calculate_rconst(net, 1)
     
        cat("writing energies in netfile.\n")
        jrnf_write("net_energies.jrnf", net)

        # Inner loop (create simulations with different boundary values
        for(v in 1:length(bvalues_l)) {
            # create directory
            ff <- paste(c("v", as.character(bvalues_l[[v]][1]), "_" ,
                        as.character(bvalues_l[[v]][2])), collapse="")
            system(paste("mkdir", ff)) 
            setwd(ff) 
       
            # create initial files for multiple runs
            for(j in 1:N_runs) {
                jrnf_create_initial(net, paste(j, ".con", sep=""), "../net.jrnf", 
                                    bc_id=bids, bc_v=bvalues_l[[v]])

                scripts <- c(scripts, 
                             paste(odeint_p, " simulate solve_implicit net=", i, 
                                   "/net.jrnf con=", i, "/", ff, "/", j, 
                                   ".con deltaT=", as.character(deltaT), " Tmax=", 
                                   as.character(Tmax), " wint=", as.character(wint), 
                                   sep=""))
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

sb_collect_results <- function(b1, b2) {
    save_wd <- getwd()   # just to be save
    df <- data.frame(Edraw=numeric(), Rdraw=numeric(), b1=numeric(), b2=numeric(), 
                     v1=numeric(), v2=numeric(), flow1=numeric(), flow2=numeric(), flow=numeric(),
                     err_cc=numeric(), err_cc_rel=numeric(),
                     ep_max=numeric(), sp_df=I(list()), re_df=I(list()), 
                     last_time=numeric(), last_msd=numeric())

    # Collect folders' names - Each corresponds to one simulation of the same
    # network with a different set of energies.
    Edir_v <- list.dirs(recursive=FALSE)


    for(bp in Edir_v) {
        setwd(bp)
        # Load networks (because they contain energy set).
        net <- jrnf_read("net_energies.jrnf")
        bp <- as.numeric(strsplit(bp,"/")[[1]][2])

        # Collect subfolders' names - they contain the boundary species' concentration.
        vdir_v <- list.dirs(recursive=FALSE)
        for(vp in vdir_v) {
            # Extract concentrations.
            l <- strsplit(vp, "v")[[1]][2]
            v1 <- as.numeric(strsplit(l, "_")[[1]][1])
            v2 <- as.numeric(strsplit(l, "_")[[1]][2])
            setwd(vp)

            # Every file corresponds to one simulation run.
            edir_v <- list.files(recursive=FALSE)
            for(ep in edir_v) {
                i <- as.numeric(strsplit(ep,"\\.")[[1]][1])
                
                # Read file
                cat("doing b1=", b1, "b2=", b2, "v1=", v1, "v2=", v2, "Edraw=", as.numeric(bp), "Rdraw=", i, "\n")
                run <- read.csv(ep)

                last_time_ <- run[nrow(run),1]
                last_msd_  <- run[nrow(run),2]
                last_con <- data.frame(con=as.numeric(run[nrow(run), 3:ncol(run)]))

                # Calculate flow (steady state rates)
                last_flow <- jrnf_calculate_flow(net, last_con$con)

                # Use concentration change to calculate exchange flow through 
                # both boundary species.
                cc <- jrnf_calculate_concentration_change(net, last_flow$flow_effective)
                flowb1 <- cc[b1]
                flowb2 <- cc[b2]
                # Estimate absolute error from "concentration change" of non 
                # exchange species and relative error by comparing this to
                # the flow through the system.
                err_cc <- max(abs(cc[c(-b1,-b2)]))
                flow <- 0.5*abs(flowb1-flowb2)
                err_cc_rel <- err_cc/flow

                # Add run to data frame
                df <- rbind(df,
                              data.frame(Edraw=as.numeric(bp), Rdraw=as.numeric(i),
                              b1=as.numeric(b1), b2=as.numeric(b2), 
                              v1=as.numeric(v1), v2=as.numeric(v2),
                              flow1=as.numeric(flowb1), flow2=as.numeric(flowb2), 
                              flow=flow, err_cc=err_cc, err_cc_rel=err_cc_rel,
                              ep_tot=as.numeric(sum(last_flow$entropy_prod)), 
                              sp_df=I(list(last_con)), re_df=I(list(last_flow)), 
                              last_time=last_time_, last_msd=last_msd_))
            }
 
            setwd("..")
        }

        setwd("..")
    }

    # Save results in original directory
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

sb_em_analysis <- function(res, net, c_max=4) {
    N_sp <- length(res$sp_df[[1]]$con)
    N_re <- length(results$re_df[[1]]$flow_effective)
   
    # matrix where combined elementary modes are saved / two reactions added -> inflow / outflow pseudoreactions
    em_matrix <- matrix(0, ncol=N_re+2, nrow=0)
    
    # list of data frames for species specific, reaction specific and elementary mode expansion data
    sp <- res$sp_df
    re <- res$re_df
    em_ex <- list()

    # Cycles (matrix and sum)
    C_matrix <- matrix(0, ncol=c_max, nrow=nrow(res))  
    C_sum <- rep(0, nrow(res))
    # Shortest pathways
    sh_p <- rep(0, nrow(res))
    sh_p_d <- rep(0, nrow(res)) # directed
    # Entropy production (linear / nonlinear reactions)
    ep_l <- rep(0,nrow(res))
    ep_nl <- rep(0,nrow(res))
    # Error
    err_rel_max <- rep(0,nrow(res))
    err <- rep(0,nrow(res))
    # Number of pathways / elementary modes
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

        # Calculating directed network
        gc()
        g <- jrnf_to_directed_network_d(net, res$re_df[[i]]$flow_effective)
        cat(".")         

        # Calculate reverse network with all rates positive
        rev <- res$re_df[[i]]$flow_effective
        rates_rev <- abs(res$re_df[[i]]$flow_effective) 
        net_rev <- jrnf_reverse_reactions(net, rev)
      
        cdif_r <- jrnf_calculate_concentration_change(net_rev, rates_rev)

        # Add reactions to balance growth / decrease
        for(k in c(res$b1[i], res$b2[i])) 
            # If species' concentration increases one pseudoreaction has to be included to remove it ("X -> ")
            if(cdif_r[k] > 0) {   
	        net_rev[[2]] <- rbind(net_rev[[2]], data.frame(reversible=factor(c(FALSE)), 
                                      c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                                      activation=as.numeric(c(0)),educts=I(list(k)), educts_mul=I(list(1)),
                                      products=I(list(c())), products_mul=I(list(c()))))
                rates_rev <- c(rates_rev, cdif_r[k])
            # If species' concentration decreases one pseudoreaction is included to add it ("-> X ")
            } else { 
	         net_rev[[2]] <- rbind(net_rev[[2]], data.frame(reversible=factor(c(FALSE)), 
                                       c=as.numeric(c(1)), k=as.numeric(c(1)),k_b=as.numeric(c(0)), 
                                       activation=as.numeric(c(0)),educts=I(list(c())), educts_mul=I(list(c())),
                                       products=I(list(k)), products_mul=I(list(1))))
                rates_rev <- c(rates_rev, -cdif_r[k])
            }

        # Calculate pathways and coefficients
        x_em <- pa_analysis(net_rev, rates_rev, param$fext, param$pmin, T, F)[[1]]
        x_rates <- pa_calculate_coef(x_em, rates_rev)$coef

        # Sort decreasing with fraction of v explained...
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
            # Calculate explained rate and explained dissipation
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
        sp[[i]]$mean_ep <- jrnf_associate_reaction_to_species(net, res$sp_df[[i]]$entropy_prod)
    }

    
    # Copying back from local objects to <res> object
    res$ep_linear <- ep_l
    res$ep_nonlinear <- ep_nl

    res$re_df <- re
    res$sp_df <- sp
    res$em_ex <- em_ex

    # Additional score to determine error
    res$reliability <- as.numeric(abs(abs(res$flow*log(res$v1/res$v2)) - res$ep_tot)/res$ep_tot > 0.1) + 1  

    res$shortest_path <- sh_p
    res$shortest_path_d <- sh_p_d

    for(k in 1:c_max) 
        res[paste("C", as.character(k), sep="")] <- C_matrix[,k]

    res$C_sum <- apply(C_matrix, 1, sum)

    res$em_no <- em_no
    res$em_err_max <- err_rel_max
    res$em_err_sum <- err

    # Save and return data
    results_em <- res
    em_m <- em_matrix
    save(results_em, em_m, file="results_em.Rdata")

    return(list(res, em_matrix))
}



