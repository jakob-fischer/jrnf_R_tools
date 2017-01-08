# author: jakob fischer (jakob@automorph.info)
# description: 
# Methods for evolving a reaction network consisting of an anorganic core and 
# multiple organisms. Each time the network is simulated for a fixed time that
# should allow it to reach steady state in most cases. Afterwards the 
# different organisms are scored with a choosen score and the worst one is 
# replaced with a new, randomly drawn organism (subnetwork).

sourced_simulation_builder_evol <- T

#   EVAL / FITNESS FUNCTIONS
# A fitness function is necessary to determine what organisms are successful.
# Its parameter is a list of results object (from "sb_solve_ae_org") for which
# "sb_ae_evo_score" already has been called. The eval function is called for
# the individual results (one element list) as well as all simulations of one
# generation (if multiple runs are done). The return value is either the worst
# organism (that will be deleted) or "NA" if none can be determined. Best 
# practice is to check if results are converged ("x$converged_proper" / see below)
# before using any results. One can use one of the fitness functions below or
# write their own (and pass it to sb_ae_org_evolve).

# build eval X - constructor for various eval functions
# As most eval functions differ just in one term this functions gets a function
# f as parameter that calculates the fitness of the organisms of one simulation.
# <f> takes a results object for which sb_ae_evo_score has been called and that
# is guaranteed to have converged.

sb_build_eval <- function(f) {
    sb_eval <- function (r) {
        l = length(r[[1]]$eval_weight$org_f)
        a = rep(NA, l) 
     
        # for all organisms calculate the maximum fitness over all simulations
        for(x in r)
            if(x$converged_proper)
                a <- pmax(a, f(x), na.rm=T)

        # if there was no simulation (value) for at least one organism return NA
        if(any(is.na(a)))
            return(NA)

        # else return the least fit one
        return(which.min(a))
    }

    return(sb_eval)
}


# EVAL FUNCTION: fitness by organism's mean fraction on the systems mass.

sb_eval_weight <- sb_build_eval(function(x) x$eval_weight$org_f)


# EVAL FUNCTION: fitness by organism's mean fraction on the systems steady state rate.

sb_eval_rates <- sb_build_eval(function(x) x$eval_rates$rates_org_f)


# EVAL FUNCTION: fitness by organism's mean fraction on the systems mass times
# its mean fraction on thes system's internal fluxes.

sb_eval_weight_X_flux <- sb_build_eval(function(x) x$eval_weight$org_f*x$eval_rates$flux_org_f)


# EVAL FUNCTION: fitness by organism's mean fraction on the systems mass times
# its mean fraction on the system's reaction rates.

sb_eval_weight_X_rates <- sb_build_eval(function(x) x$eval_weight$org_f*x$eval_rates$rates_org_f)



# Function generates a set of initial conditions, prepares everything and calls the 
# solver for which the binary is located in <sb_odeint_path> (see 
# "simulation_builder.R". Temporary files are created and deleted in the current 
# directory. The solver can be found on https://github.com/jakob-fischer/jrnf_int
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
    N <- jrnf_calculate_stoich_mat(net)    # calculate stoich matrix 
    N[net[[1]]$name == "hv",] <- 0         # and set row for "hv" to zero

    # shuffle concentrations by randomly applying different reactions to (almost) exhaustion
    s_sample <- sample(ncol(N), ncol(N)*3, replace=T)
    s_sign <- sample(c(-1,1), size=ncol(N)*3, replace=T)
    
    # As <con_o> is normaly set around one thousandth of <con_a> it is taken as
    # lower bound. All reactions (s_sample*s_sign) are applied until so far that
    # no species has a concentration lower than <con_a>
    for(k in 1:length(s_sign)) {
        rea <- N[,s_sample[k]]*s_sign[k]
        l <- min(-((init - con_o)/rea)[rea < 0])
        init <- pmax(0, init + rea*l)
    }

    # now write network file and initial concentration file to filesystem
    jrnf_write("X34tmp_net.jrnf", net)
    df <- data.frame(time=as.numeric(0),msd=as.numeric(0))
    df[as.vector(net[[1]]$name)] <- init         
    write.csv(df, "X34tmp_ini.con", row.names=FALSE)

    # call the c++-program (through shell) that solves the ode
    call_cmd <- function(cmd) {
        cat("CMD: ", cmd, "\n")
        system(cmd)
    }

    # strings with parameters for ode solver
    sim_log_s <- " simulate solve_implicit write_log net=X34tmp_net.jrnf con=X34tmp_ini.con wint="
    sim_lin_s <- " simulate solve_implicit net=X34tmp_net.jrnf con=X34tmp_ini.con wint="

    # first make 5 small linear steps and then <wint> logarithmic spaced steps
    call_cmd(paste(sb_odeint_path, sim_lin_s, as.character(5), " Tmax=", as.character(1e-8), sep=""))
    call_cmd(paste(sb_odeint_path, sim_log_s, as.character(wint), " Tmax=", as.character(min(Tmax, 1e7)), sep="")) 
    # after 1e7 do linear steps not bigger than 1e7 if necessary
    if(Tmax > 1e7)
        call_cmd(paste(sb_odeint_path, sim_lin_s, as.character(as.integer(Tmax/1e7)), " Tmax=", as.character(Tmax), sep=""))

    # load the dynamics / final concentrations / rates
    run <- read.csv("X34tmp_ini.con")

    # cleanup
    system("rm X34tmp_net.jrnf", ignore.stderr=T, ignore.stdout=T) 
    system("rm X34tmp_ini.con", ignore.stderr=T, ignore.stdout=T)

    # return final concentrations / rates / evaluation of error!
    time <- as.numeric(as.matrix(run)[nrow(run), 1])
    msd <- as.numeric(as.matrix(run)[nrow(run), 2])
    con <- as.numeric(as.matrix(run)[nrow(run), -(1:2)])
    flow <- jrnf_calculate_flow(jrnf_calculate_rconst(net), con)
    cc <- jrnf_calculate_concentration_change(net, flow$flow_effective)

    # To estimate the error / closenes to steady state the change of 
    # concentration is calculated with the last rate vector. As hv is kept 
    # constant its concentration change is through "missing" inflow reactions
    # and thuss discarded.
    # Relative error through concentration change scaled by concentration
    err_cc <- abs(cc)
    flow_b <- err_cc[net[[1]]$name == "hv"]
    err_cc[net[[1]]$name == "hv"] <- 0
    err_cc_con <- abs(cc) / con
    err_cc_con[net[[1]]$name == "hv"] <- 0   

    # Relative error concentration change scaled by maximal effective reaction rate
    err_cc_flow <- abs(cc) / max(flow$flow_effective)
    err_cc_flow[net[[1]]$name == "hv"] <- 0

    # Create and return results object / for being considered as converged the 
    # relative error estimated through concentration change has to be below 0.1.
    return(list(concentration=con, flow=flow, err_cc=max(err_cc), err_cc_con=max(err_cc_con), 
                err_cc_flow=max(err_cc_flow), flow_b=flow_b, 
                converged_proper=as.logical(max(err_cc_con) < 0.1), msd=msd, time=time))
}



# Function calculates score for an individual results object of network evolution.
#
# <net>        - Network that was simulated
# <result>     - Results object (returned from "sb_solve_ae_org")
# <eval_worst> - Eval function (see above, "sb_eval_...")
# <do_pw_ana>  - Do pathway analysis (not recommendet + not thorroghly tested)

sb_ae_evo_score <- function(net, result, eval_worst=(function (r) which.min(r[[1]]$flux_org)), do_pw_ana=F) {
    # Function takes some numeric vector <a> and an integer class vector <id> of
    # same length. For all integer values 1 ... max(id) the sum of all associated
    # entries in <a> is returned.
    subsum <- function(a, id) {
        max_id <- max(id)
        x <- rep(0, max_id)
        for(i in 1:max_id)
            x[i] <- sum(a[id == i])
        return(x)
    }

    # Calculates the flux that one part (specific organism, anorganic) of the 
    # network has with the other parts. For this the network part is separated
    # and all effective in- and out-fluxes are transformed to elementary 
    # components and all of them are summed up.
    #
    # parameters:
    # <N_rev>     - stoichiometric matrix
    # <rates_rev> - reaction rates
    # <re_assoc>  - association of reactions to different types / organisms
    # <i_i>       - Ids of organisms for which the exchange flux will be 
    #               calculated. Output is of same size (+ same order) as this
    #               vector.
    subflux <- function(N_rev, rates_rev, re_assoc, i_i) {
        comp <- net$composition
        comp[net[[1]]$name == "hv",1] <- 1  # "hv" corresponds to one mass unit.
                # If pseudo inflow reaction for "hv" is included in <N_rev> 
                # photoreactions in the anorganic part doesn't lead to exchange,
                # but photoreactions in organisms does.
        x <- c()
 
        for(i in i_i) {
            # Exchange flux in terms of network species.
            dx <- abs(matrix(N_rev[,re_assoc == i], nrow=nrow(N_rev)) %*% rates_rev[re_assoc == i])
            # Transform to elementary components and add up.
            x <- c(x, sum(t(comp) %*% dx))
        }
        return(x)
    }

    # conclude pathway analysis
    hv_id <- which(net[[1]]$name == "hv")

    # Extend network with hv-inflow-reaction and reverse reactions (if necessary).
    x <- pa_extend_net(net, result$flow$flow_effective, invwhich(hv_id, nrow(net[[1]])), T)
    net_ex <- x[[1]]
    rates_ex <- x[[2]]

    net_rev <- jrnf_reverse_reactions(net_ex, rates_ex)
    N_rev <- jrnf_calculate_stoich_mat(net_rev)
    rates_rev <- abs(rates_ex)

    # If active, calculate pathways and coefficient. Else, create empty objects.
    if(do_pw_ana) {
        res_ana <- pa_analysis(net_rev, rates_rev, 0.01, 0.7)
        a <- pa_calculate_coef(res_ana[[1]], rates_rev, T)

    } else {
        res_ana <- list(1) 
        a <- list(coef=0)  
    }

    # Collect pathwar relevant data in one list / object.
    pathways <- list()
    pathways$net_ex <- net_ex
    pathways$rates_ex <- rates_ex
    pathways$em <- res_ana[[1]]
    pathways$em_coef <- a$coef

    # Revert reaction directions that were reversed before calculating pathways.
    if(do_pw_ana) 
        for(l in 1:nrow(pathways$em))
            pathways$em[l,] <- pathways$em[l,]*(sign(rates_ex))

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
    rates$rate_anorg <- rates$rates_tot - sum(rates$rates_org)   # accumulated rate of anorganic part
 
    # sum of fluxes (total, organisms, anorganic)
    rates$flux_anorg <- subflux(N_rev, rates_rev, re_assoc_ex, 0)
    rates$flux_org <- subflux(N_rev, rates_rev, re_assoc_ex, 1:max(re_assoc_ex))
    rates$flux_tot <- sum(rates$flux_org) + rates$flux_anorg
    rates$flux_anorg_f <- rates$flux_anorg/rates$flux_tot

    # fraction of rates / fluxes (organisms) - in terms of total rates
    rates$rates_org_f <- rates$rates_org / rates$rates_tot
    rates$flux_org_f <- rates$flux_org / rates$flux_tot     # fraction a organism has on total flux
    rates$org_cycling_f <- rates$rates_org / rates$flux_org # cycling ratio defined by ratio of rates by flux
             # Because thermodynamically rates has to be driven by fluxes a high fraction implies there
             # has to be a lot of cycling.

    # For each organism calculate mass (fraction). In this case this implicitly 
    # means the mass fraction of the private part. 
    weight <- list()
    weight$con <- apply(result$concentration * net$composition, 1, sum)
    weight$total <- sum(weight$con)
    weight$org <- subsum(weight$con, net$assoc$sp)
    weight$anorg <- weight$total - sum(weight$org) 
    weight$org_f <- weight$org/weight$total
    weight$anorg_f <- weight$anorg / weight$total

    # Add calculated values / objects to <result$eval_...>
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


# Helper function tha merges a list of results and calculates all interesting
# parameters for directly observing the parameters of this generation.

sb_ae_evo_merge <- function(res_eval) {
    converged_proper <- c()
    weight_anorg_f <- c()
    flux_anorg_f <- c()
    flux_f <- c()
    rates_f <- c()
    weight_f <- c()
    cycling_f <- c()
    worst_id <- c()

    for(x in res_eval) {
        converged_proper <- c(converged_proper, x$converged_proper)
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
                       cycling_f_mean=mean(cycling_f),
                       flux_f_mean=mean(flux_f), rates_f_mean=mean(rates_f),
                       weight_f_mean=mean(weight_f), rates_f_mean=mean(rates_f), 
                       converged_proper=any(converged_proper),
                       unique_worst=(length(unique(worst_id)) == 1) && any(converged_proper)) )
}


# Function evolves a artificial ecosstem consisting of an anorganic network <net_ac>
# and a number of organisms. Evolution in this case refers to a very primitive 
# optimization / sampling of the organisms and not to a sophisticated method 
# containing mutation of genetic information. In every step the combined network
# of 
#
# parameters:
# <net_ac>     - anorganic network 
# <no_o>       - number of organisms
# <no_gen>     - number of generations
# <eval_worst> - score function that i used by sb_ae_evo_score to evaluate which
#                of the organisms can be dropped in the next simulation step. 
# <no_eva>     - number of evaluations (simulations per generation)
# <do_pw_ana>  - conclude pathway analysis while simulating (not recommended)
# <param>      - Parameters for simulation (solving the ODE). If one does not
#                want to use the standard parameter all of the following have
#                to be set:
# <param$con_a>   - Initial (before shuffling) concentration of anorganic species
# <param$con_o>   - Initial concentration of organic species
# <param$con_hv>  - Concentration of enegy / hv pseudospecies.
# <param$Tmax>    - Time up to which the ODE / network is solved.
# <param$wint>    - Number of times the ODE is writen to file while solving. 

sb_ae_org_evolve <- function(net_ac, no_o, no_gen, eval_worst, no_eva=1, do_pw_ana=F,
                             param=list(con_a=10, con_o=1e-2, con_hv=10, Tmax=1e7, wint=50)) {
    # Submethod calling <sb_solve_ae_org> to calculate initial conditions and
    # simulate the differential equation. The parameters are hardcoded now 
    # but that might be changed later. 
    # Most important parameter (initial concentration before shuffeling for 
    #                           anorganic part 10 and for organisms 1e-3)
    solve <- function(net) {
        res <- sb_solve_ae_org(net, param$con_a, param$con_o, param$con_hv, param$Tmax, param$wint)
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
    dyn$innovated <- c()
    dyn$converged_proper <- c()
    dyn$worst_id <- c()
    dyn$unique_worst <- c()
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

        res <- list()
        res_eval <- list()

        for(k in 1:no_eva) {
            res[[k]] <- solve(ext)                  # solve
            res_eval[[k]] <- sb_ae_evo_score(ext, res[[k]], eval_worst, do_pw_ana)        # evaluate result
            cat("EVALUATED!\n")
        }

        # Merge / collect information for dynamics
        x <- sb_ae_evo_merge(res_eval) 

        # Save old worst organism and determine new worst
        old_worst <- dyn$worst_id[length(dyn$worst_id)]
        if(is.null(old_worst))
            old_worst <- 1      
        new_worst <- eval_worst(res_eval)

        if(is.na(new_worst))
            new_worst <- old_worst

        # Extend the dynamical information 
        dyn$worst_id <- c(dyn$worst_id, new_worst)
        dyn$innovated <- c(dyn$innovated, new_worst != old_worst)
        dyn$converged_proper <- c(dyn$converged_proper, x$converged_proper)
        dyn$unique_worst <- c(dyn$unique_worst, x$unique_worst)
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
    } 

    # First generation is marked as innovated because it is interesting in itself.
    dyn$innovated[1] = T

    return(list(net_list=net_list, res_list=res_list, dyn_data=dyn, param=param))
}



# Function generate a results object for a network evolution (sb_ae_org_evolve)
# object. The results object is compatible with the results generated from
# sb_collect_results_ecol - this means als consecutive methods in 
# simulation_builder_ecol.R can then be applied th the resulting object.
#
# Time reffers to the time of the last timestep of the respective simulation in
# this case. The evolution generation i implicitly contained in Edraw.

sb_evol_build_results <- function(res) {
    if(is.null(res$param)) # Old version without param / add standard values that were used
        res$param <- list(con_a=10, con_o=1e-2, con_hv=10, Tmax=1e7, wint=50)

    df <- data.frame(Edraw=numeric(), Rdraw=numeric(), c=numeric(), v=numeric(), 
                     flow=numeric(), ep_max=numeric(), state_hash=character(), sp_df=I(list()), 
                     re_df=I(list()), time=numeric(), msd=numeric(), is_last=logical(), 
                     od_run_worst_id=numeric(), od_innovated=logical(), 
                     od_converged_proper=logical(), od_worst_id=numeric(),
                     od_unique_worst=logical(), od_weight_anorg_f=numeric(),
                     od_flux_anorg_f=numeric(), od_flux_f_max=numeric(),
                     od_rates_f_max=numeric(), od_weight_f_max=numeric(),
                     od_cycling_f_max=numeric(), od_flux_f_mean=numeric(),
                     od_rates_f_mean=numeric(), od_weight_f_mean=numeric(),
                     od_cycling_f_mean=numeric(), od_weight_rates_f_mean=numeric(),
                     od_weight_rates_f_max=numeric())
    # shortcuts
    net_l = res$net_list;
    res_l = res$res_list
    dyn = res$dyn_data
    param = res$param
    v <- param$con_hv
    c <- param$con_a

    for(bp in 1:length(net_l)) {
        # Generate entry in results object for each simulation of each generation.
        # Because different simulations of same generations are done with same
        # network they also have the same <Edraw> number.
        net <- net_l[[bp]]
        N <- jrnf_calculate_stoich_mat(net)

        for(i in 1:length(res_l[[bp]])) {
            cat("v=", v, "c=", c, "Edraw=", as.numeric(bp), "Rdraw=", i, "\n")
           
            xres <- res_l[[bp]][[i]]
            state_hash <- sb_v_to_hash_s(xres$flow$flow_effective, 1e-20)
            con <- data.frame(con=xres$concentration)

            xres$err_cc <- 0  # set zero as it is used to estimate error / calculate direction hash in em-analysis

            df <- rbind(df,
                        data.frame(Edraw=as.numeric(bp),Rdraw=as.numeric(i), 
                                   v=as.numeric(v), c=as.numeric(c), 
                                   flow=as.numeric(xres$flow_b), 
                                   state_hash=as.character(state_hash),
                                   relaxing_sim=as.logical(F), 
                                   err_cc=xres$err_cc, err_cc_rel=xres$err_cc_con,
                                   ep_tot=as.numeric(sum(xres$flow$entropy_prod)), 
                                   sp_df=I(list(con)), re_df=I(list(xres$flow)),
                                   time=xres$time, msd=xres$msd, is_last=as.logical(T),
                                   od_run_worst_id=as.numeric(xres$eval_worst_id),
                                   od_innovated=as.logical(dyn$innovated[bp]), 
                                   od_converged_proper=as.logical(dyn$converged_proper[bp]), 
                                   od_worst_id=as.numeric(dyn$worst_id[bp]),
                                   od_unique_worst=as.logical(dyn$unique_worst[bp]), 
                                   od_weight_anorg_f=as.numeric(dyn$weight_anorg_f[bp]),
                                   od_flux_anorg_f=as.numeric(dyn$flux_anorg_f[bp]), 
                                   od_flux_f_max=as.numeric(dyn$flux_f_max[bp]),
                                   od_rates_f_max=as.numeric(dyn$rates_f_max[bp]), 
                                   od_weight_f_max=as.numeric(dyn$weight_f_max[bp]),
                                   od_cycling_f_max=as.numeric(dyn$cycling_f_max[bp]), 
                                   od_flux_f_mean=as.numeric(dyn$flux_f_mean[bp]),
                                   od_rates_f_mean=as.numeric(dyn$rates_f_mean[bp]),  
                                   od_weight_f_mean=as.numeric(dyn$weight_f_mean[bp]),
                                   od_cycling_f_mean=as.numeric(dyn$cycling_f_mean[bp]),
                                   od_weight_rates_f_mean=as.numeric(mean(xres$eval_weight$org_f * xres$eval_rates$rates_org_f)),
                                   od_weight_rates_f_max=as.numeric(max(xres$eval_weight$org_f * xres$eval_rates$rates_org_f))))
        }
    }

    results <<- df
    results <<- results[!duplicated(results$Edraw),]
    results_nets <<- net_l
    system("mkdir -p joint")
    save(results, results_nets, file="joint/results.Rdata")
}


# Function extracts the anorganic subnet identified by assoc$re and assoc$sp.
# This is nontrivial because the function maintains all the metainformation
# and updates it correctly. 

sb_get_anorganic_subnet <- function(net) {
    net[[1]] <- net[[1]][net$assoc$sp == 0,]
    net[[2]] <- net[[2]][net$assoc$re == 0,]
    net$para$org$next_id <- 1
    net$composition <- matrix(net$composition[net$assoc$sp == 0,], 
                              ncol=ncol(net$composition))
    net$assoc$sp <- net$assoc$sp[net$assoc$sp == 0]
    net$assoc$re <- net$assoc$re[net$assoc$re == 0]

    return(net)
}


# Similar as above, function generates a results object for a network evolution 
# (sb_ae_org_evolve) object. The results object is compatible with the results 
# generated from sb_collect_results_ecol - this means als consecutive methods in 
# simulation_builder_ecol.R can then be applied th the resulting object. In 
# contrary to above function, this one looks only at the anorganic core of
# the network. As this is unchanged, all the states are using the same network
# (Edraw = constant). In this case the generation is implicitly contained in
# Rdraw. (Not quite sure how to handle the case if one makes multiple simulations
# / evaluations for the same generation, here. For starters just analyze the first
# simulation.) 
# The function saves the resulting object under "core/results.Rdata" relative
# to the current path.

sb_evol_build_results_core <- function(res) {
    if(is.null(res$param)) # Old version without param / add standard values that 
                           # were used back then.
        res$param <- list(con_a=10, con_o=1e-3, con_hv=10, Tmax=1e7, wint=50)

    df <- data.frame(Edraw=numeric(), Rdraw=numeric(), c=numeric(), v=numeric(), 
                     flow=numeric(), ep_max=numeric(), state_hash=character(), sp_df=I(list()), 
                     re_df=I(list()), time=numeric(), msd=numeric(), is_last=logical(), 
                     od_run_worst_id=numeric(), od_innovated=logical(), 
                     od_converged_proper=logical(), od_worst_id=numeric(),
                     od_unique_worst=logical(), od_weight_anorg_f=numeric(),
                     od_flux_anorg_f=numeric(), od_flux_f_max=numeric(),
                     od_rates_f_max=numeric(), od_weight_f_max=numeric(),
                     od_cycling_f_max=numeric(), od_flux_f_mean=numeric(),
                     od_rates_f_mean=numeric(), od_weight_f_mean=numeric(),
                     od_cycling_f_mean=numeric(), od_weight_rates_f_mean=numeric(),
                     od_weight_rates_f_max=numeric())
    # shortcuts
    net_l = res$net_list;
    res_l = res$res_list
    dyn_l = res$dyn_data
    param = res$param
    v <- param$con_hv
    c <- param$con_a
    
    net_a <- sb_get_anorganic_subnet(res$net_list[[1]])
    N <- jrnf_calculate_stoich_mat(net_a)
    sel_sp <- 1:nrow(N)
    sel_re <- 1:ncol(N)    

    # Add row to results data frame for each simulation in each generation
    for(i in 1:length(res_l)) 
        for(xres in res_l[[i]]) {

            if(!is.na(xres$eval_worst_id) && xres$eval_worst_id == dyn_l$worst_id[i]) {
                cat("v=", v, "c=", c, "Edraw=", 1, "Rdraw=", i, "\n")
                state_hash <- sb_v_to_hash_s(xres$flow$flow_effective[sel_sp], 1e-20)
                con <- data.frame(con=xres$concentration[sel_sp])

                xres$err_cc <- 0  # set zero as it is used to estimate error / calculate 
                                  # direction hash in em-analysis

                df <- rbind(df,
                             data.frame(Edraw=as.numeric(1),Rdraw=as.numeric(i), 
                                        v=as.numeric(v), c=as.numeric(c), 
                                        flow=as.numeric(xres$flow_b), 
                                        state_hash=as.character(state_hash),
                                        relaxing_sim=as.logical(F), 
                                        err_cc=xres$err_cc, err_cc_rel=xres$err_cc_con,
                                        ep_tot=as.numeric(sum(xres$flow$entropy_prod[sel_re])), 
                                        sp_df=I(list(con)), re_df=I(list(xres$flow[sel_re,])),
                                        time=xres$time, msd=xres$msd, is_last=as.logical(T),
                                        od_run_worst_id=as.numeric(xres$eval_worst_id),
                                        od_innovated=as.logical(dyn_l$innovated[i]), 
                                        od_converged_proper=as.logical(dyn_l$converged_proper[i]), 
                                        od_worst_id=as.numeric(dyn_l$worst_id[i]),
                                        od_unique_worst=as.logical(dyn_l$unique_worst[i]), 
                                        od_weight_anorg_f=as.numeric(dyn_l$weight_anorg_f[i]),
                                        od_flux_anorg_f=as.numeric(dyn_l$flux_anorg_f[i]), 
                                        od_flux_f_max=as.numeric(dyn_l$flux_f_max[i]),
                                        od_rates_f_max=as.numeric(dyn_l$rates_f_max[i]), 
                                        od_weight_f_max=as.numeric(dyn_l$weight_f_max[i]),
                                        od_cycling_f_max=as.numeric(dyn_l$cycling_f_max[i]), 
                                        od_flux_f_mean=as.numeric(dyn_l$flux_f_mean[i]),
                                        od_rates_f_mean=as.numeric(dyn_l$rates_f_mean[i]),  
                                        od_weight_f_mean=as.numeric(dyn_l$weight_f_mean[i]),
                                        od_cycling_f_mean=as.numeric(dyn_l$cycling_f_mean[i]),
                                        od_weight_rates_f_mean=as.numeric(mean(xres$eval_weight$org_f * xres$eval_rates$rates_org_f)),
                                        od_weight_rates_f_max=as.numeric(max(xres$eval_weight$org_f * xres$eval_rates$rates_org_f))))
             }
         }

    # Write to common location in global environment and into subdirectory "core".
    results <<- df
    results <<- results[!duplicated(results$Rdraw),]
    results_nets <<- list(net_a)
    system("mkdir -p core")
    save(results, results_nets=results_nets, file="core/results.Rdata")
}
