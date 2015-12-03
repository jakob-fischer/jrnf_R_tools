# author: jakob fischer (jakob@automorph.info)
# date: 3rd December 2015
# description: 

sourced_potential_analysis <- T

if(!exists("sourced_jrnf_network"))
    source("jrnf_network.R")

if(!exists("sourced_pathway_analysis"))
    source("pathway_analysis.R")

# Thermodynamic constants
k_B <- 1.381E-23  # Boltzmann constant (J / K)
R <- 8.314        # Ideal gas constant (J/(K mol))
p0 <- 100000      # Standard pressure (Pa | N/m²)


# Load coefficients of nasa polynomials
nasa_polys <- read.csv("nasa_polynomials_200_1000.csv", header=F)


# Sum function. All NA values are removed. If all values of the vector are NA
# the function returns NA instead of 0 (standard behaviour)
msum <- function(x) {
    if(all(is.na(x)))
        return(NA)

    return(sum(x, na.rm=T))
}


# Multiplicates the matrix 'm' with the vector 'v' 
# If a zero entry in the matrix is multiplied with a NA entry in the vector a 
# value of 0 is taken in the corresponding part of the sum 
mmul <- function(m, v) {
    if(ncol(m) != length(v)) {
        cat("mmul: dimension mismatch!\n")
        return(0)
    }

    res <- rep(0, nrow(m))

    for(i in 1:nrow(m))
        res[i] <- sum((m[i,]*v)[m[i,] != 0])
    
    return(res)
}


# Calculates the chemical potential "mu" for the species named "sp_name" at
# temperature "T" and number density "N" (the last one given in particles/cm^3) 

get_mu <- function(sp_name, T, N) {
    # first find the species in the list of nasa polynomials coefficients 
    sp_id <- which(nasa_polys$V1 == sp_name)

    # check it was identified and read the constants
    if(length(sp_id) != 1)
        return(NA)

    a1 <- nasa_polys$V2[sp_id] 
    a2 <- nasa_polys$V3[sp_id]
    a3 <- nasa_polys$V4[sp_id]
    a4 <- nasa_polys$V5[sp_id]
    a5 <- nasa_polys$V6[sp_id]
    a6 <- nasa_polys$V7[sp_id]
    a7 <- nasa_polys$V8[sp_id]


    # now calculate h0 and s0 as in paper by Venot
    h0 <- a1*R*T + a2/2*R*T**2 + a3/3*R*T**3 + a4/4*R*T**4 + a5/5*R*T**5 + a6*R
    s0 <- a1*R*log(T) + a2*R*T + a3/2*R*T**2 + a4/3*R*T**3 + a5/4*R*T**4 + a7*R

    # Calculate (partial) pressure P and from that the chemical potential mu
    p <- N*k_B*T*100**3              # conversion because of different units 'cm' -> 'm' 
    
    mu <- h0 - T*s0 + R*T*log(p/p0)

    # If mu is -Inf because p==0 set it to NA
    if(!is.finite(mu))
        return(NA)

    return(mu)
} 


# Use the concentration profile 'profile' and other atmospheric parameter atparam 
# (especially temperature profile) 'atparam' to calculate the profile of chemical
# potential for all species.

calculate_mu_profile <- function(profile, atparam) {
    # Concentration profile and atmospheric parameters have to match
    if(nrow(profile) != nrow(atparam)) {
        cat("calc_mu_matrix: Number of rows does not match!\n")
        return()
    }

    if(any(profile$z != atparam$z)) {
        cat("calc_mu_matrix: Profile does not match!\n")
        return()   
    }

    # clone and set all entries to zero
    mu_m <- profile   
    mu_m[,-1] <- 0

    # Calculate chemical potential for each species and each height.
    for(i in 1:nrow(mu_m)) {
        T <- atparam$T[i]
        for(j in 2:ncol(mu_m))
            mu_m[i,j] <- get_mu(names(profile)[j], T, profile[i,j])
    }

    return(mu_m)
}


# Calculates the density profile of the atmosphere from the profile dataframe 
# containing the concentration profiles of all species 
calculate_density_profile <- function(profile) {
    return(as.numeric(apply(data.matrix(profile[,-1]),1, sum)))
}

# Integrates the matrix of chemical potentials and returns them as a vector.
# User can specify whether the specific density of the species is used for weighting ('sp') 
# or if the density of all species is used.
integrate_mu <- function(mu_m, profile, atparam, sp=FALSE, Tr=FALSE) {
    acc <- c()

    if(Tr) {
        if(sp) {
            # Use species concentration profile for weighting
            for(i in 2:ncol(mu_m))
                acc <- c(acc, msum(mu_m[,i]*profile[,i]/atparam$T)/sum(profile[,i]))
        } else {
            # Use density / all species concentration profile for weighting
            for(i in 2:ncol(mu_m))
                acc <- c(acc, msum(mu_m[,i]*atparam$density/atparam$T)/sum(atparam$density))
        }  

    } else {
        if(sp) {
            # Use species concentration profile for weighting
            for(i in 2:ncol(mu_m))
                acc <- c(acc, msum(mu_m[,i]*profile[,i])/sum(profile[,i]))
        } else {
            # Use density / all species concentration profile for weighting
            for(i in 2:ncol(mu_m))
                acc <- c(acc, msum(mu_m[,i]*atparam$density)/sum(atparam$density))
        }    
    }

    return(acc)
}


# Assigns the chemical potentials (ordered as indicated in the 'mu_names' parameter)
# to the species in the order of the jrnf reaction network 'net'. If 'hv' and 'M' are
# not included in 'mu_names' the values of 'v_hv' and 'v_M' are used.
jrnf_assign_mu <- function(net, mu, mu_names, v_hv=1e20, v_M=0) {
    # All chemical potentials are initialized with NA
    mu_new <- rep(NA, nrow(net[[1]]))

    # Assign all those potentials that are included in 'mu_names'
    for(i in 1:length(mu_new)) {
        x <- which(mu_names == net[[1]]$name[i])

        if(length(x) == 1)
            mu_new[i] <- mu[x]
    }

    # If 'hv's chemical potential is still NA use 'v_hv'
    x <- which("hv" == net[[1]]$name)
    if(length(x) == 1)
        if(is.na(mu_new[x]))
            mu_new[x] <- v_hv;

    # If 'M's chemical potential is still NA use 'v_M'
    x <- which("M" == net[[1]]$name)
    if(length(x) == 1)
        if(is.na(mu_new[x]))
            mu_new[x] <- v_M;

    return(mu_new)
}


# Function checks if potentials are compatible with the direction of the 
# reactions.
# (The function is a legacy function and in principle all functionality containted
# shoud be also containedin the function "calculate_reaction_energetics" that returns
# on big dataframe with all results.

check_potentials <- function(net, rates, E) {
    N <- jrnf_calculate_stoich_mat(net)
    se <- c()
 
    d_mu <- mmul(t(N), E)
    se <- d_mu*rates < 0

    cat("NA = ", length(which(is.na(se))), "\n")
    cat("T  = ", length(which(TRUE == se )), "\n")
    cat("F  = ", length(which(FALSE == se)), "\n")

    return(which(se == FALSE))
}



# 
# Methods test potentials for plausibility with reaction directions
#

calculate_reactions_energetics <- function(net, mu, re_rates=c()) {
    no_reas <- nrow(net[[2]])

    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)
    N <- N_out - N_in
    hv_id <- which(net[[1]]$name == "hv")
    if(length(hv_id) != 1) {
        cat("ERROR: did not find unique hv species!\n")
        return(0)
    }

    # set chemical potential of hv to 0 for convenient calculation of mu_in and mu_out
    mu[hv_id] <- 0

    mu_in <- mmul(t(N_in), mu)
    mu_out <- mmul(t(N_out), mu)
    mu_eff <- mu_out - mu_in

    has_hv <- N[hv_id,] != 0

    dir_match <- (mu_out < mu_in) | (N[hv_id,] < 0)

    undecided <- which(is.na(dir_match))
    mismatch <- which(dir_match == FALSE)

    cat("===================================================================================\n")
    cat("From totally", length(dir_match), "reaction directions", sum(dir_match == TRUE, na.rm=T), "are correct ")
    cat(length(mismatch), "are incorrect and", length(undecided), "are undecided.\n")
    if(length(re_rates) != 0)
        cat("Weighted with reaction's rates mismatch is", sum(re_rates[mismatch])/sum(re_rates), 
            "and undecided", sum(re_rates[undecided])/sum(re_rates), ".\n")
    cat("From ", length(which(has_hv)), "reactions with hv there are", length(which(has_hv & mu_eff < 0)),
        "that are directly dissipating and", length(which(has_hv & is.na(mu_eff))), "with NA\n")
    cat("There are", length(which(has_hv & N[hv_id,] > 0)), "reactions that produce hv!",
        length(which(has_hv & N[hv_id,] > 0 & mu_eff < 0)), "of them are dissipating and",
        length(which(has_hv & N[hv_id,] > 0 & (dir_match == F | is.na(dir_match)))), "of them are mismatching or NA.\n")
    cat("===================================================================================\n")

    return(data.frame(mu_in=mu_in, mu_out=mu_out, mu_eff=mu_eff, dir_match=dir_match, has_hv=has_hv))
}


#
# Methods test potentials for plausibility with pathways directions
# (in its final version this function should also check the potentials
#  - give feedback on the potentials quality) 
#  
#

calculate_pathways_energetics <- function(net, mu, ems, em_rates, re_rates=c()) {
    re_en <- calculate_reactions_energetics(net, mu, re_rates)
    N <- jrnf_calculate_stoich_mat(net)

    if(ncol(ems) != nrow(re_en)) {
        cat("calculate_pathway_energetics: ems / re_ens dimension mismatch!\n")
        return(0)
    }

    hv_id <- which(net[[1]]$name == "hv")
    if(length(hv_id) != 1) {
        cat("ERROR: did not find unique hv species!\n")
        return(0)
    }

    all_species_present <- rep(F, nrow(ems))           #
    interacting_species_present <- rep(F, nrow(ems))   #
    
    # Thermodynamic data on pathways derived from per reaction data
    hv_in <- rep(0, nrow(ems))
    dis <- rep(0, nrow(ems))  
    match_all_reactions <- rep(F, nrow(ems))

    # Thermodynamic data on pathways derived directly 
    # (intermediate species with missing / inconsistent potentials cancel out)
    hv_in_D <- rep(0, nrow(ems))
    delta_mu_D <- rep(0, nrow(ems))   
    match_interactions <- rep(F, nrow(ems))

    # Also topological properties of the reaction pathways are calculated
    # (similar to the function pa_em_derive)    
    
    # calculate species bilance in matrix form (every row contains change of
    # species concentration through pathway)
    bilance_sp <- t(N %*% t(ems))
    
    ems_hv <- ems
    ems_hv[,!re_en$has_hv] <- 0
    bilance_sp_hv <- t(N %*% t(ems_hv))
         
    # number of species taking part in elementary mode (not counting hv!)
    Sp_no <- apply(bilance_sp[,-hv_id] != 0, 1, sum)
    # number of reactions taking part in elementary mode
    Re <- apply(abs(ems),1, sum)
    # number of reactions (counting each only once)
    Re_s <- apply(ems != 0, 1, sum)
    # Number of input from environment (not counting hv!)
    In <- apply(bilance_sp[,-hv_id], 1, function(x)  {  x[x>0] <- 0;  return(-sum(x))  })
    # Number of input (counting each species only once)
    In_s <- apply(bilance_sp[,-hv_id], 1, function(x)  {  x[x>0] <- 0; x[x<0] <- 1; return(sum(x))  })
    # Number of output to environment (not counting hv!)
    Out <- apply(bilance_sp[,-hv_id], 1, function(x)  {  x[x<0] <- 0;  return(sum(x))  })
    # Number of output to environment (each species only once) 
    Out_s <- apply(bilance_sp[,-hv_id], 1, function(x)  {  x[x<0] <- 0; x[x>0] <- 1; return(sum(x))  })

    # Number of hv-pseudospecies that are consumed by pathway
    Hv <- -bilance_sp[,hv_id]

    # set hv in bilance matrices to zero
    bilance_sp[,hv_id] <- 0
    bilance_sp_hv[,hv_id] <- 0

    for(i in 1:nrow(ems)) {
        #cat(".")

        # calculate thermodynamic properties from reaction data
        hv_in[i] <- sum((ems[i,]*re_en$mu_eff)[re_en$mu_eff > 0 & ems[i,] != 0])
        dis[i] <- -sum((ems[i,]*re_en$mu_eff)[re_en$mu_eff < 0  & ems[i,] != 0])
        match_all_reactions[i] <- all(re_en$dir_match[which(ems[i,]!=0)])

        # calculate the thermodynamic properties of the pathways directly    
        hv_in_D[i] <- max(0, sum((bilance_sp_hv[i,]*mu)[bilance_sp_hv[i,]!=0]))
        delta_mu_D[i] <- sum((bilance_sp[i,]*mu)[bilance_sp[i,]!=0])
        all_species_present[i] <- all(is.finite( re_en$mu_eff[ems[i,] != 0] ))         

        # If one of the interacting species (bilance_sp / bilance_sp_hv) are NA
        # matching direction can not be decided (-> NA)
        b <- bilance_sp[i,]-bilance_sp_hv[i,]
        match_interactions[i] <- sum((b*mu)[b != 0]) < 0
        if(length(which(b != 0)) == 0)   # special case (only dissipating hv)
            match_interactions[i] <- T
        interacting_species_present[i] <- all(is.finite(c(mu[bilance_sp[i,] != 0],
                                                           mu[bilance_sp_hv[i,] != 0]))) 
    }

    # quantities that can be calculated after hv_in and dis are known
    delta_mu <- hv_in - dis
    eff <- pmin(pmax(delta_mu/hv_in,0),1)  

    dis_D <- hv_in_D - delta_mu_D
    eff_D <- pmin(pmax(delta_mu_D/hv_in_D,0),1)  

    cat("===================================================================================\n")
    cat("From", nrow(ems), "pathways on a reaction level", length(which(match_all_reactions)),
        "match on reaction level with", length(which(!match_all_reactions)), "mismatches and", 
        length(which(is.na(match_all_reactions))), "are NA!\n")
    cat("On pathway level", length(which(match_interactions)),
        "match", length(which(!match_interactions)), "mismatches and", 
        length(which(is.na(match_interactions))), "are NA!\n")

    cat("Weighted with pathway's rates mismatch is", sum(em_rates[!match_interactions], na.rm=T)/sum(em_rates), 
        "and undecided", sum(em_rates[is.na(match_interactions)])/sum(em_rates), ".\n")

    cat("===================================================================================\n")

    return(data.frame(hv_in, delta_mu, dis, eff, match_all_reactions, all_species_present, hv_in_D,
                       delta_mu_D, dis_D, eff_D, match_interactions, interacting_species_present,
                       Sp_no, Re, Re_s, In, In_s, Out, Out_s, Hv))
}



# This function creates 
#
#  * The data layout of genomes / fitness cases is like this 
#  * (N is the number of chemical species and M the number of reactions):
#  * - The first N values (floats) contains the chemical potential of the chemical species.
#  * - The next M values are the activation energies of the species.
#  * - The last M values are the effective reaction rates of the reactions.
#  
# In this case the reaction directions + the priors of the chemical potentials are taken
# to calculate the fitness. This means the actual rates, the activation energies and the 
# temperature of the system are not significant.

create_energy_evolver_initial <- function(net, rates, mu_in, net_file, mutate_f_file, fitness_f_file, sample_file, initial_file, mutate_all=T) {

    hv_id <- which(net[[1]]$name == "hv")
    if(length(hv_id) != 1) {
        cat("ERROR: did not find unique hv species!\n")
        return(0)
    }
    
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)
    N <- N_out - N_in

    re_sel <- N_in[hv_id,] == 0  &  N_out[hv_id,] == 0  &  N[hv_id,] == 0
    net <- list(net[[1]], net[[2]][re_sel,])
    rates <- rates[re_sel]

    no_sp <- nrow(net[[1]])
    no_re <- nrow(net[[2]])

    # write reaction network
    jrnf_write(net_file, net)

    # prepare initial state
    mu_in_ <- mu_in
    mu_in_[is.na(mu_in)] <- runif(sum(is.na(mu_in)), mean(mu_in, na.rm=T), sd(mu_in, na.rm=T)) 
    initial <- c(mu_in_, rep(1, no_re), sign(rates))
    write.table(initial, file=initial_file, row.names=F, col.names=F)

    # prepare sample state
    mu_in_ <- mu_in
    mu_in_[is.na(mu_in)] <- 0 
    sample <- c(mu_in_, rep(10, no_re), rates)
    write.table(sample, file=sample_file, row.names=F, col.names=F)

    # prepare the mutate flag ("1" means position is mutated, "0" means it stays unchanged)
    mutate_f <- c(is.na(mu_in), rep(0,2*no_re))

    if(mutate_all)
        mutate_f <- c(rep(1,no_sp), rep(0,2*no_re))
    write.table(mutate_f, file=mutate_f_file, row.names=F, col.names=F)


    # the fitness flag indicate which entries are used for evaluating the fitness 
    # by comparing the sign of the fluxes (="1") and which are used precisely (="2").
    # Entries not relevant for the fitness are indicated by the value "0".

    fitness_f <- c(rep(2, no_sp), rep(0, no_re), rep(1, no_re))
    fitness_f[which(is.na(mu_in))] <- 0
    write.table(fitness_f, file=fitness_f_file, row.names=F, col.names=F)
}

