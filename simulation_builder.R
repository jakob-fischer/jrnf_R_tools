# author: jakob fischer (jakob@automorph.info)
# description: 
# Script for generating files and directory structure of artificial chemistry
# simulations with thermodynamic constraints. (simulation builder = "sb_")
# The functionality has been distributed in multiple files:
# simulation_builder_flow.R - Build simulation driven by in- and outflow
# simulation_builder_ecol.R - Ecosystems with constant mass driven by energy
# simulation_builder_evol.R - Evolving Ecosystem (containing organisms)
# simulaiton_builder_plot.R - Plot methods for different types of simulations
#
# Parts of this this file / modules functionality are shared with "netodeint_control.R", this
# could be unified in future.

sourced_simulation_builder <- T


if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")

if(!exists("sourced_art_ecosystem_gen"))
    source("art_ecosystem_gen.R")

if(!exists("sourced_pathway_analysis"))
    source("pathway_analysis.R")

if(!exists("sourced_simulation_builder_flow"))
    source("simulation_builder_flow.R")

if(!exists("sourced_simulation_builder_ecol"))
    source("simulation_builder_ecol.R")

if(!exists("sourced_simulation_builder_evol"))
    source("simulation_builder_evol.R")


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


# This function can be used instead of <sb_v_to_hash> if the number of reactions
# exceeds 60. It uses a string as hash instead of a long integer. This implies 
# that it needs one byte for each reaction of memory. 
# TODO: Check usage pattern (especially in this file. Maybe it is better to 
#       use sb_v_to_hash_s as standard case to avoid problems if somone uses
#       the toolkit with bigger networks.

sb_v_to_hash_s <- function(v, zero_range) {
    f <- function(x) {  if(x > zero_range) return("+") else if (x < -zero_range) return("-") else return("0")  }
    return(paste(sapply(v, f), collapse=""))
}


# This function probably only has legacy function.
# TODO: Check an remove if save

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


# Same as above (<assemble_em_from_flow>)
# TODO check usage and remove if possible

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



# Function to estimate the stability of the steady state assuming the system
# got to it. For this the jacobi matrix is calculated for each concentration
# using the <sa> object generated from <jrnf_build_linear_stability_analyzer>.
# Then the spectrum of the matrix is calculated. Then a object is returned that
# contains lists with all jacobi-matrices, all spectrums, all maximal eigenvalues
# and all minimal eigenvalues (for all rows in results <res>)

sb_lin_stab_analysis_ecol <- function(res_nets, res) {
    mat <- list()
    ew <- list()
    ew_max <- c()
    ew_min <- c()    

    for(i in 1:nrow(res)) {
        cat(".")
        sa <- jrnf_build_linear_stability_analyzer(res_nets[[res$Edraw[i]]])
        x <- sa$calculate(res$sp_df[[i]]$con)
        mat[[length(mat)+1]] <- x
        y <- eigen(x, symmetric=F)$values
        ew[[length(ew)+1]] <- y
        rel <- abs(Re(y)) > 1e-10
        ew_max[length(ew_max)+1] <- max(Re(y[rel]))
        ew_min[length(ew_min)+1] <- min(Re(y[rel]))
    }
    cat("\n")

    return(list(jacobi=mat, jacobi_spec=ew, re_max=ew_max, re_min=ew_min))
}



#
#

sb_reduce_last_save <- function() {
    results <<- results[results$is_last,]
    save(results, results_nets, file="results_red.Rdata")
}


