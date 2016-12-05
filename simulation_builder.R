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
# Parts of this this files / modules functionality are shared with "netodeint_control.R", this
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


# Constant that defines the maximum step size for calls of integrator (jrnf_int).
# Logarithmic stepping can be used up to this time point and after that linear
# steps are used.

sb_max_step_size <- 1e7

# Constant that defines the local path to the integrator for odes "jrnf_int".
# (program is available on github: https://github.com/jakob-fischer/jrnf_int)

sb_odeint_path <- "~/apps/jrnf_int"


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
# that it needs one byte of memory for each reaction. 

sb_v_to_hash_s <- function(v, zero_range) {
    f <- function(x) {  if(x > zero_range) return("+") else if (x < -zero_range) return("-") else return("0")  }
    return(paste(sapply(v, f), collapse=""))
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


# Function reduces the results object (in the global environment) to those that 
# are the last of their simulation (have "is_last" flag) and saves it in the 
# current directory as "results_red.Rdata".

sb_reduce_last_save <- function() {
    results <<- results[results$is_last,]
    save(results, results_nets, file="results_red.Rdata")
}


# Function calculates how contribution of different pathways to a steady state
# changes for different rows in a results_em / results_em_cross object. 
# Because an entire results object can contain data reffering to different 
# network objects one has to subset the results first. This function takes
# the subsetted results object, and selected network and elementary modes
# objects. One can consider a subset of all elementary modes by using the
# parameter <sub_em>. In all cases the matrix for explained fraction is calculated
# normalized by row as well as non-normalized.
#
# parameters:
# <results> - subetted results object (Edraw has to be identical for all entries)
# <net>     - associated network object
# <em>      - associated matrix of elementary modes
# <sub_em>  - optional, a boolean vector indicating a subset of elementary modes

sb_calc_pw_cross <- function(results, net, em, sub_em=c()) {
    if(is.null(sub_em))
        sub_em <- rep(T, nrow(em))
    x <- n <- r <- matrix(0, nrow(results), sum(sub_em))

    
    for(i in 1:nrow(results)) 
        if(is.list(results$em_ex[[i]])) {
            # calculate explained rate and coefficient for full set of pathways
            exp_r <- rate <- rep(0, nrow(em))
            exp_r[results$em_ex[[i]]$id] <- results$em_ex[[i]]$exp_r
            rate[results$em_ex[[i]]$id] <- results$em_ex[[i]]$rate
            # now take the relevant subset
            x[i,] <- exp_r[sub_em]
            n[i,] <- x[i,] / sum(x[i,])
            r[i,] <- rate[sub_em]
        } else {
            x[i,] <- NA
            n[i,] <- NA
            r[i,] <- NA      
        }    
    return(list(x=x, x_norm=n, rates=r))
}


# This function calls sb_calc_pw_cross to calculate how contribution of different 
# pathways to a steady state changes for different rows in a results_em / 
# results_em_cross object. The difference of this function is that it subsets 
# the results object and selects net and em...
#
# parameters:
# <results> - results object 
# <net>     - associated list of networks
# <em>      - associated list of matrices of elementary modes
# <i>       - Which network (Edraw) to analyze
# <sub_em>  - optional, a boolean vector indicating a subset of elementary modes

sb_calc_pw_cross_i <- function(results, results_net, em, i, sub_em=c()) {
    sel <- results$Edraw == i
    return(sb_calc_pw_cross(results[sel,] , results_net[[i]], em[[i]], sub_em))
}


# Given a matrix <m> that indicates contribution / explained fraction of different
# elementary modes to different simulation runs, the function determines the <N>
# most important elementary modes.

sb_pw_cross_sub_N <- function(m, N) {
    a <- apply(m, 2, function(x) max(x, na.rm=T))
    return(head(order(a, decreasing=T), N))
}


# Given a matrix <m> that indicates contribution / explained fraction of different
# elementary modes to different simulation runs, the function determines the most
# important elementary nodes that have explained fraction higher than <f>.

sb_pw_cross_sub_f <- function(m, f) {
    a <- apply(m, 2, function(x) max(x, na.rm=T))
    r <- which(a > f)
    return(r)
}
