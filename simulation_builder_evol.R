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

sb_ae_evo_score <- function(net, result, eval_worst=(function (r) which.min(r$flux_org)), do_pw_ana=T) {
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
        cat("re_assoc=", re_assoc, "  - len(..)=", length(re_assoc), "\n")
        cat("i_i=", i_i, "\n")
        cat("dim(N_rev)=", dim(N_rev), "\n")
 
        for(i in i_i) {
            dx <- abs(matrix(N_rev[,re_assoc == i], nrow=nrow(N_rev)) %*% rates_rev[re_assoc == i])
            x <- c(x, sum(t(comp) %*% dx))
        }
        return(x)
    }

    # conclude pathway analysis
    hv_id <- which(net[[1]]$name == "hv")

    cat("net - size: ", nrow(net[[1]]), " - ", nrow(net[[2]]), "\n")

    x <- pa_extend_net(net, result$flow$flow_effective, invwhich(hv_id, nrow(net[[1]])), T)
    net_ex <- x[[1]]
    rates_ex <- x[[2]]


    cat("net_ex - size: ", nrow(net_ex[[1]]), " - ", nrow(net_ex[[2]]), "\n")

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
    weight$org_f = weight$org/weight$total

    # 
    result$eval_pathways <- pathways
    result$eval_weight <- weight
    result$eval_rates <- rates
    result$eval_worst_id <- eval_worst(result)
    
    cat(" weight (fraction) of org.:", weight$org_f, "\n")
    cat(" flux (fraction) of org.:", rates$flux_org_f, "\n")
    cat(" rates(fraction) of org.:", rates$rates_org_f, "\n")
    cat(" cycling fraction of org.:", rates$org_cycling_f, "\n")
    cat(" worst is organism with id ", result$eval_worst_id, "\n")
    cat("-----------------------------------------------------------\n")


    return(result)
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

sb_ae_org_evolve <- function(net_ac, no_o, no_gen, eval_worst, do_pw_ana=T) {
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
    dyn$weight_f_anorg <- c()
    dyn$flux_anorg <- c()
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
        res <- solve(ext)                  # solve
        res_eval <- sb_ae_evo_score(ext, res, eval_worst, do_pw_ana)        # evaluate result
        cat("EVALUATED!\n")

        dyn$worst_id <- c(dyn$worst_id, res_eval$eval_worst_id)
        dyn$weight_f_anorg <- c(dyn$weight_f_anorg, res_eval$eval_weight$anorg/res_eval$eval_weight$total)
        dyn$flux_anorg <- c(dyn$flux_anorg, res_eval$eval_rates$flux_anorg)
        dyn$flux_f_max <- c(dyn$flux_f_max, max(res_eval$eval_rates$flux_org_f))
        dyn$rates_f_max <- c(dyn$rates_f_max, max(res_eval$eval_rates$rates_org_f))
        dyn$weight_f_max <- c(dyn$weight_f_max, max(res_eval$eval_weight$org_f))
        dyn$cycling_f_max <-c(dyn$cycling_f_max, max(res_eval$eval_rates$org_cycling_f))
        dyn$flux_f_mean <- c(dyn$flux_f_mean, mean(res_eval$eval_rates$flux_org_f))
        dyn$rates_f_mean <- c(dyn$rates_f_mean, mean(res_eval$eval_rates$rates_org_f))
        dyn$weight_f_mean <- c(dyn$weight_f_mean, mean(res_eval$eval_weight$org_f))
        dyn$cycling_f_mean <-c(dyn$cycling_f_mean, mean(res_eval$eval_rates$org_cycling_f))

        net_list[[length(net_list)+1]] <- ext         # list collecting data
        res_list[[length(res_list)+1]] <- res_eval    # list collecting data

        #return(list(ext, res_eval))

        # kill off worst organism
        ext <- jrnf_ae_remove_organism(ext, res_eval$eval_worst_id)

        if(nrow(ext[[2]]) != length(ext$assoc$re)) {
            cat("BAD LENGTH MISMATCH in ASSOC RE (BB):\n")
            return(ext)
        }
    } 

    return(list(net_list=net_list, res_list=res_list, dyn_data=dyn))
}

