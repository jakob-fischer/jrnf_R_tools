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


plot_em_spectre <- function(em_cross, x_data, x_name="xname", sel=c(), filename="em_spectre.eps") { 
    if(is.null(sel))
        sel <- rep(T, length(x_data))

    sel <- which(sel)

    a <- data.frame(x_data=numeric(), 
                    exp_r=numeric(),
                    em=factor())
 
    for(i in 1:ncol(em_cross)) {
        a <- rbind(a, 
                   data.frame(x_data=as.numeric(x_data[sel]), 
                              exp_r=as.numeric(em_cross[sel,i]),
                              em=as.factor(rep(i, length(sel)))))
    }

    setEPS()
    postscript(filename, width=10, height=6)
    
    gg<-ggplot(a, aes(x=x_data, y=exp_r, colour=em, group=em, shape=em))  + 
          scale_x_log10() + scale_y_continuous(limits=c(1e-2,1)) + 
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

