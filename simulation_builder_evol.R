# author: jakob fischer (jakob@automorph.info)
# description: 
# TODO

sourced_simulation_builder_evol <- T


#  Function generates a set of initial conditions, prepares everything and calls the 
#  solver for which the binary is located in "~/apps/jrnf_int". Temporary files are
#  created and deleted in the current directory. The solver can be found on
#  https://github.com/jakob-fischer/jrnf_int
#
# parameters:
# <net>    - network containing information on association of species to organisms - inorganic part
# <con_a>  - initial concentration (before shuffling) of anorganic species
# <con_o>  - initial concentration (before shuffling) of species related to organisms
# <con_hv> - concentration of driving / hv species
# <Tmax>   - time to which is simulated
# <wint>   - number of intervalls that the results is written to file

sb_solve_ae_org <- function(net, con_a, con_o, con_hv, Tmax, wint) {
    
    init <- rep(con_o, nrow(net[[1]]))
    init[net$assoc$sp == 0] <- con_a       # set concentration for anorganic part
    init[net[[1]]$name == "hv"] <- con_hv  # for hv / energy
    N <- jrnf_calculate_stoich_mat(net)   # calculate stoich matrix 
    N[net[[1]]$name == "hv",] <- 0        # and set row for "hv" to zero

    # shuffle concentrations by randomly applying different reactions to (almost) exhaustion
    s_sample <- sample(ncol(N), ncol(N)*3, replace=T)
    s_sign <- sample(c(-1,1), size=ncol(N)*3, replace=T)
    
    for(k in 1:length(s_sign)) {
        rea <- N[,s_sample[k]]*s_sign[k]
        l <- min(-((init - con_o/100)/rea)[rea < 0])
        init <- pmax(0, init + rea*l)
    }

    # now write network file and initial concentration file to filesystem
    jrnf_write("X34tmp_net.jrnf", net)
    df <- data.frame(time=as.numeric(0),msd=as.numeric(0))
    df[as.vector(net[[1]]$name)] <- init         
    write.csv(df, "X34tmp_ini.con", row.names=FALSE)

    # call the c++-program that solves the ode
    call_cmd <- function(cmd) {
        cat("CMD: ", cmd, "\n")
        system(cmd)

    }

    odeint_p <- "~/apps/jrnf_int"
    # first make 5 small linear steps and then <wint> logarithmic spaced steps
    # TODO check if linear steps are necessary or if the tool odeint_p has to be
    # adapted.
    call_cmd(paste(odeint_p, " simulate solve_implicit  net=X34tmp_net.jrnf con=X34tmp_ini.con wint=", as.character(5), " Tmax=", as.character(1e-8), sep=""))
    call_cmd(paste(odeint_p, " simulate solve_implicit write_log net=X34tmp_net.jrnf con=X34tmp_ini.con wint=", as.character(wint), " Tmax=", as.character(Tmax), sep=""))

    # load the dynamics / final concentrations / rates
    run <- read.csv("X34tmp_ini.con")

    # cleanup
    system("rm X34tmp_net.jrnf", ignore.stderr=T, ignore.stdout=T) 
    system("rm X34tmp_ini.con", ignore.stderr=T, ignore.stdout=T)

    # return final concentrations / rates / evaluation of error!
    con <- as.numeric(as.matrix(run)[nrow(run), -(1:2)])
    flow <- jrnf_calculate_flow(jrnf_calc_reaction_r(net), con)
    cc <- jrnf_calculate_concentration_change(net, flow$flow_effective)
    # To estimate the error / closenes to steady state the change of 
    # concentration is calculated with the last rate vector. As hv is kept 
    # constant its concentration change is through "missing" inflow reactions
    # and thuss discarded.
    # Relative error through concentration change scaled by concentration
    err_cc_con <- abs(cc) / con
    err_cc_con[net[[1]]$name == "hv"] <- 0   
    # relative error concentration change scaled by maximal effective reaction rate
    err_cc_flow <- abs(cc) / max(flow$flow_effective)
    err_cc_flow[net[[1]]$name == "hv"] <- 0

    return(list(concentration=con, flow=flow, err_cc_con=max(err_cc_con), err_cc_flow=max(err_cc_flow)))
}



# 
#
#
#

sb_ae_evo_score <- function(net, result, eval_worst=(function (r) which.min(r[[1]]$flux_org)), do_pw_ana=T) {
    subsum <- function(a, b) {
        max_b <- max(b)
        x <- rep(0, max_b)
        for(i in 1:max_b)
            x[i] <- sum(a[b == i])
        return(x)
    }

    # Calculates the flux that one part (specific organism, anorganic) of the 
    # network has with the other parts. For this the network part is separated
    # and all effective in- and out-fluxes are transformed to elementary 
    # components and all of them are summed up.
    #
    #
    # 
    subflux <- function(N_rev, rates_rev, re_assoc, i_i) {
        comp <- net$composition
        comp[net[[1]]$name == "hv",1] <- 1
        x <- c()
 
        for(i in i_i) {
            dx <- abs(matrix(N_rev[,re_assoc == i], nrow=nrow(N_rev)) %*% rates_rev[re_assoc == i])
            x <- c(x, sum(t(comp) %*% dx))
        }
        return(x)
    }

    # conclude pathway analysis
    hv_id <- which(net[[1]]$name == "hv")

    x <- pa_extend_net(net, result$flow$flow_effective, invwhich(hv_id, nrow(net[[1]])), T)
    net_ex <- x[[1]]
    rates_ex <- x[[2]]

    net_rev <- jrnf_reverse_reactions(net_ex, rates_ex)
    N_rev <- jrnf_calculate_stoich_mat(net_rev)
    rates_rev <- abs(rates_ex)
    # TODO
    # res_ana <- pa_analysis(net_rev, rates_rev, 0.001, 0.7)

    if(do_pw_ana) {
        res_ana <- pa_analysis(net_rev, rates_rev, 0.01, 0.7)
        a <- pa_calc_coefficients(res_ana[[1]], rates_rev, T)

    } else {
        res_ana <- list(1) # pa_analysis(net_rev, rates_rev, 0.01, 0.7)
        a <- list(coef=0)  #pa_calc_coefficients(res_ana[[1]], rates_rev, T)
    }

    pathways <- list()
    pathways$net_ex <- net_ex
    pathways$rates_ex <- rates_ex
    pathways$em <- res_ana[[1]]
    pathways$em_coef <- a$coef

    if(do_pw_ana) {
        for(l in 1:nrow(pathways$em))
            pathways$em[l,] <- pathways$em[l,]*(sign(rates_ex))
    }

    # for each organism calculate exchange fluxes
    # and calculate explained fraction of steady state 
    # and especially fraction of exchange and cummulative steady state
    rates <- list()     

    # sum of rates (total, organisms, anorganic)
    # The pseudoreaction adding hv is accounted as anorganic. 
    re_assoc_ex <- c(net$assoc$re, 0)
    rates$rates_tot <- sum(rates_rev)                            # total accumulated rate
    rates$rates_org <- subsum(rates_rev, re_assoc_ex)            # vector containing accumulated rate of all
                                                                # reactions associated with different organisms
    rates$rate_anorg <- rates$rates_tot - sum(rates$rates_org)  # accumulated rate of anorganic part
 
    # sum of fluxes (total, organisms, anorganic)
    rates$flux_anorg <- subflux(N_rev, rates_rev, re_assoc_ex, 0)
    rates$flux_org <- subflux(N_rev, rates_rev, re_assoc_ex, 1:max(re_assoc_ex))
    rates$flux_tot <- sum(rates$flux_org) + rates$flux_anorg
    rates$flux_anorg_f <- rates$flux_anorg/rates$flux_tot

    # fraction of rates (organisms) - in terms of total rates
    rates$rates_org_f <- rates$rates_org / rates$rates_tot
    rates$flux_org_f <- rates$flux_org / rates$flux_tot   # fraction a organism has on total flux
    rates$org_cycling_f <- rates$rates_org / rates$flux_org   # cycling ratio defined by ratio of rates by flux
                     # Because thermodynamically rates has to be driven by fluxes a high fraction implies there
                     # has to be a lot of cycling.


    # for each organism calculate mass (fraction)
    # especially mass fraction in elementary components
    # also especially that of the private part
    weight <- list()
    weight$con <- apply(result$concentration * net$composition, 1, sum)
    weight$total <- sum(weight$con)
    weight$org <- subsum(weight$con, net$assoc$sp)
    weight$anorg <- weight$total - sum(weight$org) 
    weight$org_f <- weight$org/weight$total
    weight$anorg_f <- weight$anorg / weight$total

    # 
    result$eval_pathways <- pathways
    result$eval_weight <- weight
    result$eval_rates <- rates
    result$eval_worst_id <- eval_worst(list(result))
    
    cat(" weight (fraction) of org.:", weight$org_f, "\n")
    cat(" flux (fraction) of org.:", rates$flux_org_f, "\n")
    cat(" rates(fraction) of org.:", rates$rates_org_f, "\n")
    cat(" cycling fraction of org.:", rates$org_cycling_f, "\n")
    cat(" worst is organism with id ", result$eval_worst_id, "\n")
    cat("-----------------------------------------------------------\n")


    return(result)
}


# helper function tha merges a list of results and calculates all interesting
# parameters for directly observing the parameters.
#

sb_ae_evo_merge <- function(res_eval) {
    weight_anorg_f <- c()
    flux_anorg_f <- c()
    flux_f <- c()
    rates_f <- c()
    weight_f <- c()
    cycling_f <- c()
    worst_id <- c()

    for(x in res_eval) {
        weight_anorg_f <- c(weight_anorg_f, x$eval_weight$anorg_f)
        flux_anorg_f <- c(flux_anorg_f, x$eval_rates$flux_anorg_f)
        flux_f <- c(flux_f, x$eval_rates$flux_org_f)
        rates_f <- c(rates_f, x$eval_rates$rates_org_f)
        weight_f <- c(weight_f, x$eval_weight$org_f)
        cycling_f <- c(cycling_f, x$eval_rates$org_cycling_f)
        worst_id <- c(worst_id, x$eval_worst_id)
    }

    return( data.frame(weight_anorg_f=mean(weight_anorg_f), flux_anorg_f=mean(flux_anorg_f),
                       flux_f_max=max(flux_f), rates_f_max=max(rates_f), 
                       weight_f_max=max(weight_f), cycling_f_max=max(cycling_f), 
                       flux_f_mean=mean(flux_f), rates_f_mean=mean(rates_f),
                       weight_f_mean=mean(weight_f), rates_f_mean=mean(rates_f), 
                       unique_worst=(length(unique(worst_id)) == 1) ) )
}


# Function evolves a artificial ecosstem consisting of an anorganic network <net_ac>
# and a number of organisms. Evolution in this case refers to a very primitive 
# optimization / sampling of the organisms and not to a sophisticated method 
# containing mutation of genetic information. In every step the combined network
# of 
#
# parameters:
# <net_ac> - anorganic network 
# <no_o>   - number of organisms
# <no_gen> - number of generations
# <eval_worst>  
#          - score function that i used by sb_ae_evo_score to evaluate which
#            of the organisms can be dropped in the next simulation step. 
#            (see above)

sb_ae_org_evolve <- function(net_ac, no_o, no_gen, eval_worst, no_eva=1, do_pw_ana=T) {
    # Submethod calling <sb_solve_ae_org> to calculate initial conditions and
    # simulate the differential equation. The parameters are hardcoded now 
    # but that might be changed later. 
    # Most important parameter (initial concentration before shuffeling for 
    #                           anorganic part 10 and for organisms 1e-3)
    solve <- function(net) {
        res <- sb_solve_ae_org(net, 10, 1e-3, 10, 1e7, 50)
        return(res)
    }

    # create extended network <ext> and add all but one organisms 
    ext <- net_ac
    if(no_o != 1)
        for(i in 1:(no_o-1))
            ext <- jrnf_ae_add_organism(ext)

    # List of extended networks and corresponding simulations
    net_list <- list()
    res_list <- list()
    dyn <- list()
    dyn$worst_id <- c()
    dyn$weight_anorg_f <- c()
    dyn$flux_anorg_f <- c()
    dyn$flux_f_max <- c()
    dyn$rates_f_max <- c()
    dyn$weight_f_max <- c()
    dyn$cycling_f_max <-c()
    dyn$flux_f_mean <- c()
    dyn$rates_f_mean <- c()
    dyn$weight_f_mean <- c()
    dyn$cycling_f_mean <-c()

    # Do simulation loom <no_gen> times
    for(i in 1:no_gen) {
        cat("-----------------------------------------------------\n")
        cat("---  DOING GENERATION ", i, "\n")
        cat("-----------------------------------------------------\n")
 
        ext <- jrnf_ae_add_organism(ext)   # add one additional organism

        if(nrow(ext[[2]]) != length(ext$assoc$re)) {
            cat("BAD LENGTH MISMATCH in ASSOC RE (AA):\n")
            return(ext)
        }

        res <- list()
        res_eval <- list()

        for(k in 1:no_eva) {
            res[[i]] <- solve(ext)                  # solve
            res_eval[[i]] <- sb_ae_evo_score(ext, res[[i]], eval_worst, do_pw_ana)        # evaluate result
            cat("EVALUATED!\n")
        }

        x <- sb_ae_evo_merge(res_eval)  # TODO implement this function
        dyn$worst_id <- c(dyn$worst_id, eval_worst(res_eval))

        dyn$weight_anorg_f <- c(dyn$weight_anorg_f, x$weight_anorg_f)
        dyn$flux_anorg_f <- c(dyn$flux_anorg_f, x$flux_anorg_f)
        dyn$flux_f_max <- c(dyn$flux_f_max, x$flux_f_max)
        dyn$rates_f_max <- c(dyn$rates_f_max, x$rates_f_max)
        dyn$weight_f_max <- c(dyn$weight_f_max, x$weight_f_max)
        dyn$cycling_f_max <-c(dyn$cycling_f_max, x$cycling_f_max)
        dyn$flux_f_mean <- c(dyn$flux_f_mean, x$flux_f_mean)
        dyn$rates_f_mean <- c(dyn$rates_f_mean, x$rates_f_mean)
        dyn$weight_f_mean <- c(dyn$weight_f_mean, x$weight_f_mean)
        dyn$cycling_f_mean <-c(dyn$cycling_f_mean, x$cycling_f_mean)

        net_list[[length(net_list)+1]] <- ext         # list collecting data
        res_list[[length(res_list)+1]] <- res_eval    # list collecting data

        # kill off worst organism
        ext <- jrnf_ae_remove_organism(ext, dyn$worst_id[length(dyn$worst_id)])

        if(nrow(ext[[2]]) != length(ext$assoc$re)) {
            cat("BAD LENGTH MISMATCH in ASSOC RE (BB):\n")
            return(ext)
        }
    } 

    return(list(net_list=net_list, res_list=res_list, dyn_data=dyn))
}



# Function generate a results object for a network evolution (sb_ae_org_evolve)
# object. The results object is compatible with the results generated from
# sb_collect_results_ecol - this means als consecutive methods in 
# simulation_builder_ecol.R can then be applied th the resulting object.
#
# Time reffers to the time of the last timestep of the respective simulation in
# this case. The evolution generation i implicitly contained in Edraw.
#
# TODO complete (copy and paste of sb_collect_results_ecol at the moment)s

sb_evol_build_results <- function(res) {
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


# Similar as above, function generate a results object for a network evolution 
# (sb_ae_org_evolve) object. The results object is compatible with the results 
# generated from sb_collect_results_ecol - this means als consecutive methods in 
# simulation_builder_ecol.R can then be applied th the resulting object. But in 
# comparison to above function, this one looks only at the anorganic core of
# the network. As this is unchanged all the states are using the same network
# (Edraw = constant). In this case the generation is implicitly contained in
# Rdraw.
#
# TODO complete (copy and paste of sb_collect_results_ecol at the moment)

sb_evol_build_results_core <- function(res) {
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
