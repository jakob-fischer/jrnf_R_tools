source("jrnf_network.R")
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
integrate_mu <- function(mu_m, profile, atparam, sp=FALSE) {
    acc <- c()

    if(sp) {
        # Use species concentration profile for weighting
        for(i in 2:ncol(mu_m))
            acc <- c(acc, msum(mu_m[,i]*profile[,i])/sum(profile[,i]))
    } else {
        # Use density / all species concentration profile for weighting
        for(i in 2:ncol(mu_m))
            acc <- c(acc, msum(mu_m[,i]*atparam$density)/sum(atparam$density))
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



check_ems <- function(em, em_expand, mm_rea) {
    x <- rep(F, nrow(em))

    for(i in mm_rea) {
        x <- x | em[,i]
    }

    cat(length(which(x)), " elementary modes of ", length(x), " have a direction mismatch (f=", length(which(x))/length(x), "\n")
    cat("The explained fraction is ", sum(em_expand$exp_f[x]), " of totally ",  sum(em_expand$exp_f),"\n")

    return(which(x))
}


# 
# Methods test potentials for plausibility with reaction directions
#

calculate_reactions_energetics <- function(net, mu) {
    N_in <- jrnf_calculate_stoich_mat_in(net)
    N_out <- jrnf_calculate_stoich_mat_out(net)
    hv_id <- which(net[[1]]$name == "hv")
    if(length(hv_id) != 1) {
        cat("ERROR: did not find unique hv species!\n")
        return(0)
    }

    mu[hv_id] <- 0

    mu_in <- rep(0, nrow(net[[2]]))
    mu_out <- rep(0, nrow(net[[2]]))
    mu_hv <- rep(0, nrow(net[[2]])) 
    dis <- rep(0, nrow(net[[2]]))

    for(i in 1:nrow(net[[2]])) {
        ed <- which(N_in[,i] != 0)
        pro <- which(N_out[,i] != 0)

        # simple outflow or inflow reaction
        if(length(ed) == 0 && length(pro) == 1) {
            mu_in[i] <- mu[pro]*N_out[pro,i]
        } else if(length(ed) == 1 && length(pro) == 0) {
            mu_out[i] <- mu[ed]*N_in[ed,i]

        # photochemical reaction
        } else if(hv_id %in% ed) {
            #ed <- ed[-which(ed == hv_id)]
            mu_hv[i] <- sum((N_out-N_in)[,i]*mu)   
        } else if(hv_id %in% pro) {
            mu_hv[i] <- sum((N_out-N_in)[,i]*mu) 
            cat("WARNING: reaction ", i, " is producing hv!\n")
            #return(0)

        # normal reaction (not driven externally)
        } else {  
            dis[i] <- sum((N_in-N_out)[,i]*mu)           
        }

        if(mu_hv[i] < 0) {
            dis[i] <- -mu_hv[i]
            mu_hv[i] <- 0
        }     

    }

    return(data.frame(mu_in=mu_in, mu_out=mu_out, mu_hv=mu_hv, dis=dis))
}


#
# Methods test potentials for plausibility with pathways directions
#

calculate_pathways_energetics <- function(re_en, ems, em_rates) {
    mu_in <- rep(0, nrow(ems))
    mu_out <- rep(0, nrow(ems))
    mu_hv <- rep(0, nrow(ems))
    dis <- rep(0, nrow(ems))
    hv_eff <- rep(0,nrow(ems))

    for(i in 1:nrow(ems)) {
        mu_in[i] <- sum(re_en$mu_in*ems[i,])
        mu_out[i] <- sum(re_en$mu_out*ems[i,])
        mu_hv[i] <- sum(re_en$mu_hv*ems[i,])
        dis[i] <- sum(re_en$dis*ems[i,])

        if(dis[i] < mu_hv[i])
            hv_eff[i] <- 1-dis[i]/mu_hv[i] 
    }

    mu_exchange <- mu_out-mu_in
    f_exchange <- abs(mu_exchange*em_rates)/sum(abs(mu_exchange*em_rates))
    f_hv <- pmax(0,mu_hv*em_rates)/sum(pmax(0,mu_hv*em_rates))
    f_dis <- (dis*em_rates)/sum(dis*em_rates)
     

    return(data.frame(mu_in=mu_in, mu_out=mu_out, mu_hv=mu_hv, dis=dis, mu_exchange=mu_exchange, hv_eff=hv_eff, 
                       f_exchange=f_exchange, f_hv=f_hv, f_dis=f_dis))
}

