# author: jakob fischer (jakob@automorph.info)
# description: 
# Script for generating files and directory structure of artificial chemistry
# simulations with thermodynamic constraints. (simulation builder = "sb_")
# This file contains the functilonality to plot (and partially evaluate) the results.
# The core functionality is in "simulation_builder.R".

sourced_simulation_builder_plot <- T

library(ggplot2)
library(gridExtra)

if(!exists("sourced_simulation_builder"))
    source("simulation_builder.R")


#
#
#

sb_plot_em_spectre <- function(em_cross, x_data, x_name="xname", em_name=c(), filename="em_spectre.eps") { 
    a <- data.frame(x_data=numeric(), 
                    exp_r=numeric(),
                    em=factor())
 
    for(i in 1:ncol(em_cross)) {
        a <- rbind(a, 
                   data.frame(x_data=as.numeric(x_data), 
                              exp_r=as.numeric(em_cross[,i]),
                              em=as.factor(rep(em_name[i], length(x_data)))))
    }

    postscript(filename, width=10, height=6)
    
    gg<-ggplot(a, aes(x=x_data, y=exp_r, colour=em, group=em, shape=em))  + 
           scale_y_continuous(limits=c(1e-3,1)) + 
           ylab("exp. fraction") + 
           xlab(x_name) +
           geom_point(size=2.2) + 
           theme_bw(22) +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position="bottom")

    print(gg)
    
    dev.off()

    return(gg)
}



# function calculates how contribution of different pathways to a steady state
# changes for different rows in a results_em / results_em_cross object. 
# Because an entire results object can contain data reffering to different 
# network objects one has to subset the results first. This function takes
# the subsetted results object, and selected network and elementary modes
# objects. One can only consider a subset of all elementary modes by using the
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

    
    for(i in 1:nrow(results)) {
        # calculate explained rate and coefficient for full set of pathways
        exp_r <- rate <- rep(0, nrow(em))
        exp_r[results$em_ex[[i]]$id] <- results$em_ex[[i]]$exp_r
        rate[results$em_ex[[i]]$id] <- results$em_ex[[i]]$rate
        # now take the relevant subset
        x[i,] <- exp_r[sub_em]
        n[i,] <- x[i,] / sum(x[i,])
        r[i,] <- rate[sub_em]
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
    a <- apply(m, 2, max)
    return(tail(order(a), N))
}


# Given a matrix <m> that indicates contribution / explained fraction of different
# elementary modes to different simulation runs, the function determines the most
# important elementary nodes that have explained fraction higher than f


sb_pw_cross_sub_f <- function(m, f) {
    a <- apply(m, 2, max)
    r <- which(a > f)
    return(r)
}


#
#

sb_pw_cross_plot <- function(net, em, cross, x_data, x_name="xname", sub_s=c(), filename="pw_cross_%.eps") {
    if(length(grep("%", filename)) == 0) {
        cat("WARNING: sb_pw_cross_plot:\n")
        cat("filename should contain placeholder '%' - output may be faulty!\n") 
    }

    if(is.null(sub_s)) 
        sub_s <- 1:nrow(em)

    if(is.logical(sub_s))
        sub_s <- which(sub_s)


    # First print spectrum of selected / subseted pathways
    em_name <- 1:ncol(cross)
    y <- sb_plot_em_spectre(cross[,sub_s], x_data, x_name, em_name[sub_s], filename=sub("%", "spec", filename))

    # Now plot all (selected) pathways in separate files
    for(i in sub_s) {
        postscript(sub("%", as.character(i), filename), width=5, height=5, family="serif")
        x <- jrnf_plot_pathway(net, em[i,1:nrow(net[[2]])], layout_f=layout.lgl, lim_plot=T)
        dev.off()
    }

    return(y)
}




sb_write_network_pathways <- function(sel=0) {
    jrnf_network_to_ltable("network.tex", net)

    sel_em <- which(apply(em_exp_r_cross, 2, max) >= sel)
    sel_em_info <- data.frame(cycles=em_derive$C_s_sum[sel_em] ,EMPTY=paste("(", as.character(sel_em), ")", sep=""))
    pa_pathways_to_ltable("pathways.tex", em_m[sel_em,1:nrow(net[[2]])], net, sel_em_info)
}


#
#

evaluate_Edraw_plot <- function(Edraw, my_results=c(), my_em=c(), my_net=c(), filename="test.eps") {
    if(is.null(my_em))
        my_em <- em

    if(is.null(my_net))
        my_net <- net

    if(is.null(my_results))
        my_results <- results_em

    my_results_X <- sb_extend_results(my_results)
    my_results_X <- my_results_X[my_results_X$Edraw == Edraw,]

    x <- assemble_em_from_flow(my_results_X, my_em, my_net)
    

    gg1 <- ggplot(x, aes(x=flow, y=em_exp_r, colour=as.factor(em_id), group=as.factor(em_id), lty=as.factor(em_id), size=2)) + 
           scale_x_log10() + 
           scale_y_log10() + 
           xlab("flow") +
           ylab("fraction") +
           geom_line(size=1.2) + 
           theme_bw(28) + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position="bottom", legend.title=element_blank())

    gg2 <- ggplot(x, aes(x=flow, y=C_sum, size=2)) + 
           scale_x_log10() + 
           xlab("flow") +
           ylab("cycles") +
           geom_line(size=1.2) + 
           theme_bw(28) + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position=c(1,.7))

    gg3 <- ggplot(x, aes(x=flow, y=C_f_d, size=2)) + 
           scale_x_log10() + 
           xlab("flow") +
           ylab("cycle factor") +
           geom_line(size=1.2) + 
           theme_bw(28) + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position=c(1,.7))

    y <- data.frame(flow=c(x$flow, x$flow), dis_f=c(x$ep_lin, x$ep_nonl), type=c(rep("linear", nrow(x)), rep("nonlinear", nrow(x))))


    gg4 <- ggplot(y, aes(x=flow, y=dis_f, colour=as.factor(type), group=as.factor(type), lty=as.factor(type), size=2)) + 
           scale_x_log10() + scale_y_log10() +
           xlab("flow") +
           ylab("fraction") +
           geom_line(size=1.2) + 
           theme_bw(28) + 
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position="bottom", legend.title=element_blank())


    setEPS()
    postscript(filename, width=14, height=9)
    grid.arrange(gg1, gg2, gg3, gg4, ncol=2)
    dev.off()


#    return(list(gg1, gg2, gg3, gg4)) 

}



plot_core_sp <- function(res_cross, filename="core_sp.eps") { 
    setEPS()
    postscript(filename, width=10, height=6)
    
    gg<-ggplot(res_cross, aes(x=v, y=core_sp_no, colour=as.factor(c), group=as.factor(c), shape=as.factor(c)))  + 
          scale_x_log10() + scale_y_continuous() + 
           ylab("#core species") + 
           xlab("v (hv concentration)") +
           geom_point(size=2.2) + 
           theme_bw(22) +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position="bottom")

    print(gg)
    
    dev.off()

    return(gg)
}


plot_flow_v <- function(res_cross, filename="flow_v.eps") { 
    setEPS()
    postscript(filename, width=10, height=6)
    
    gg<-ggplot(res_cross, aes(x=v, y=flow, colour=as.factor(c), group=as.factor(c), shape=as.factor(c)))  + 
          scale_x_log10() + scale_y_log10() + 
           ylab("flow") + 
           xlab("hv concentration") +
           geom_point(size=2.2) + 
           theme_bw(22) +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position="bottom")

    print(gg)
    
    dev.off()

    return(gg)
}



plot_flow_v2 <- function(res_cross, filename="flow_v2.eps") { 
    setEPS()
    postscript(filename, width=10, height=6)
    
    gg<-ggplot(res_cross, aes(x=v2, y=flow, colour=as.factor(Edraw), group=as.factor(Edraw), shape=as.factor(Edraw)))  + 
          scale_x_log10() + scale_y_log10() + 
           ylab("flow") + 
           xlab("boundary concentration") +
           geom_point(size=2.2) + 
           theme_bw(22) +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position="bottom")

    print(gg)
    
    dev.off()

    return(gg)
}



plot_cycles_v <- function(res_cross, filename="cycles_v.eps") { 
    setEPS()
    postscript(filename, width=10, height=6)
    
    gg<-ggplot(res_cross, aes(x=v, y=cycles_r, colour=as.factor(c), group=as.factor(c), shape=as.factor(c)))  + 
          scale_x_log10() + scale_y_continuous() + 
           ylab("cycles") + 
           xlab("hv concentration") +
           geom_point(size=2.2) + 
           theme_bw(22) +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position="bottom")

    print(gg)
    
    dev.off()

    return(gg)
}


plot_cycles_v2 <- function(res_cross, filename="cycles_v2.eps") { 
    setEPS()
    postscript(filename, width=10, height=6)
    
    gg<-ggplot(res_cross, aes(x=v2, y=cycles_r, colour=as.factor(Edraw), group=as.factor(Edraw), shape=as.factor(Edraw)))  + 
          scale_x_log10() + scale_y_continuous() + 
           ylab("cycles") + 
           xlab("boundary concentration") +
           geom_point(size=2.2) + 
           theme_bw(22) +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position="bottom")

    print(gg)
    
    dev.off()

    return(gg)
}



plot_shannon_v <- function(res_cross, filename="shannon_v.eps") { 
    setEPS()
    postscript(filename, width=10, height=6)
    
    gg<-ggplot(res_cross, aes(x=v, y=informationE, colour=as.factor(c), group=as.factor(c), shape=as.factor(c)))  + 
          scale_x_log10() + scale_y_continuous() + 
           ylab("inf. entropy") + 
           xlab("hv concentration") +
           geom_point(size=2.2) + 
           theme_bw(22) +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position="bottom")

    print(gg)
    
    dev.off()

    return(gg)
}


plot_shannon_v2 <- function(res_cross, filename="shannon_v2.eps") { 
    setEPS()
    postscript(filename, width=10, height=6)
    
    gg<-ggplot(res_cross, aes(x=v2, y=informationE, colour=as.factor(Edraw), group=as.factor(Edraw), shape=as.factor(Edraw)))  + 
          scale_x_log10() + scale_y_continuous() + 
           ylab("inf. entropy") + 
           xlab("boundary concentration") +
           geom_point(size=2.2) + 
           theme_bw(22) +
           theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position="bottom")

    print(gg)
    
    dev.off()

    return(gg)
}

