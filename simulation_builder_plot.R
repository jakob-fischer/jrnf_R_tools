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


# Function plots evolution of different pathway related properties in the 
# current directory. They are either plotted as a function of Edraw (<Edyn> == T)
# or as a function of Rdraw (<Edyn> == F).
#
# The <Edyn> parameter determines if dynamics is defined by increasing <Edraw>
# in the results object as it is for networks evolving their structure or 
# alternatively through <Rdraw> as for an anorganic core network where the 
# network itself stays constant through all generations.
#
# The function does not take more arguments but assumes that the results of a
# simulation are loaded in the global enviroment (as it is put there by 
# "sb_em_cross_analysis_ecol").

sb_plot_evol_pw <- function(Edyn=T) {
    # Subfunction to plot into file <fn>. Y-axis is named <yname> and the data 
    # plotted is <ydata>. Length of <ydata> has to match the number of rows in 
    # <results_em_cross>. Y-axis is logscale if <logY> is set.
    plt <- function(fn, yname, ydata, width=7, height=5, logY=F) {
        if(Edyn)
            a <- data.frame(b=results_em_cross$Edraw, c=ydata)
        else
            a <- data.frame(b=results_em_cross$Rdraw, c=ydata)
        postscript(fn, width=width, height=height)

        scale <- scale_y_continuous()
        if(logY) 
            scale <- scale_y_log10()


        gg <- ggplot(a, aes(x=b, y=c))  + 
            scale_x_continuous() +
            scale +
            ylab(yname) + 
            xlab("generation") +
            geom_point(size=2.2, color="grey") + 
            geom_point(data=a[results_em_cross$od_innovated,], aes(x=b, y=c), size=2.3, color="blue") +
            theme_bw(22) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position="bottom")
        
         print(gg)
         dev.off()
    }

    rec <- results_em_cross
    # Don't call plt if ALL values of need_em_90 are NA (it's surely a bad 
    # simulation, but we still want to plot the other values and not abort
    # because of an error).
    if(!all(is.na(results_em_cross$need_em_90)))
        plt("evol_need_em_90.eps", "pw for 90%", rec$need_em_90, logY=T)
    plt("evol_exp_r_max.eps", "max exp. r",  rec$exp_r_max)
    plt("evol_cycles_r.eps", "cycle number", rec$cycles_r)
    plt("evol_species_r.eps", "species number", rec$species_r)
    plt("evol_reactions_r.eps", "reaction number", rec$reactions_r)
    plt("evol_em_steadystate_r.eps", "steady state fraction", rec$em_steadystate_r)
    plt("evol_em_con_ch_r.eps", "concentration change fraction", rec$em_con_ch_r)
    plt("evol_informationE.eps", "inf. E. of pw dist.", rec$informationE)
}


# Generic plot function that can be applied in the "core" directory generated
# from "sb_evol_build_results_core". The function plots all those parameters that
# where saved while evolving (and ARE NOT specific to the core but to the entire
# simulation / network). In the subdirectory "pw_evol" the importance of the 5 
# most relevant pathways as function of generation / "Rdraw" is plotted using
# "pw_evol". Additionally the next five most relevant pathways are plotted.
#
# The function does not take more arguments but assumes that the results of a
# simulation are loaded in the global enviroment (as it is put there by 
# "sb_em_cross_analysis_ecol").

sb_plot_evol_core <- function() {
    # Subfunction to plot into file <fn>. Y-axis is named <yname> and the data 
    # plotted is <ydata>. Length of <ydata> has to match the number of rows in 
    # <results_em_cross>. Y-axis is logscale if <logY> is set.
    plt <- function(fn, yname, ydata, width=7, height=5, logY=F) {
        a <- data.frame(b=results_em_cross$Rdraw, c=ydata)
        postscript(fn, width=width, height=height)

        scale <- scale_y_continuous()
        if(logY) 
            scale <- scale_y_log10()


        gg <- ggplot(a, aes(x=b, y=c))  + 
            scale_x_continuous() +
            scale +
            ylab(yname) + 
            xlab("generation") +
            geom_point(size=2.2, color="grey") + 
            geom_point(data=a[results_em_cross$od_innovated,], aes(x=b, y=c), size=2.3, color="blue") +
            theme_bw(22) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 legend.key = element_rect(color = "white", fill="white"), legend.position="bottom")
        
         print(gg)
         dev.off()
    }

    rec <- results_em_cross


    plt("evol_weight_f_mean.eps", "mean weight fraction", rec$od_weight_f_mean)
    plt("evol_weight_f_max.eps", "max weight fraction", rec$od_weight_f_max)
    plt("evol_weight_f_anorg.eps", "anorg. weight fraction", rec$od_weight_anorg_f, logY=T)
    plt("evol_rates_f_mean.eps", "mean rates fraction", rec$od_rates_f_mean)
    plt("evol_rates_f_max.eps", "max rates fraction", rec$od_rates_f_max)
    plt("evol_flux_f_mean.eps", "mean flux fraction", rec$od_flux_f_mean)
    plt("evol_flux_f_max.eps", "max flux fraction", rec$od_flux_f_max)

    plt("evol_flux_anorg_f.eps", "anorg. flux fraction", rec$od_flux_anorg_f)
    plt("evol_cycling_f_mean.eps", "mean cycling", rec$od_cycling_f_mean)
    plt("evol_cycling_f_max.eps", "max cycling", rec$od_cycling_f_max)

    plt("evol_weight_rates_f_mean_DIR.eps", "mean weight*rates", rec$od_weight_rates_f_mean)
    plt("evol_weight_rates_f_max_DIR.eps", "mean weight*rates", rec$od_weight_rates_f_max)
    plt("evol_weight_rates_f_mean.eps", "mean weight*rates", rec$od_weight_f_mean*rec$od_rates_f_mean)
    plt("evol_weight_flux_f_mean.eps", "mean weight*flux", rec$od_weight_f_mean*rec$od_flux_f_mean)

    plt("evol_worst_id.eps", "worst id", rec$od_worst_id)
    plt("evol_unique_worst.eps", "unique worst", as.numeric(rec$od_unique_worst))
    plt("evol_flow.eps", "flow (hv)", rec$flow, logY=T)

    # plot change of important parameters in main directory:
    sb_plot_evol_pw(F)

    system("mkdir -p pw_evol")
    setwd("pw_evol")
    for(Edraw in unique(rec$Edraw)) {
        # calculate pathway evolution of core in subdirectory
        # plot diagramms of 10 most important pathways
        fn <- paste("E", as.character(Edraw), "_pw_cross_%.eps", sep="")

        # Select network and pathway matrix...
        net <- results_nets[[Edraw]]
        em <- em_m[[Edraw]]
        # ...and calculate contribution cross matrix.
        a <- sb_calc_pw_cross(rec, net, em)
       
        # Select 5 most relevant pathways...
        sub_x <- sb_pw_cross_sub_N(a$x_norm, 5)
        # ... and the following five.
        sub_y <- sb_pw_cross_sub_N(a$x_norm, 10)[-(1:5)]

        # Plot pathway spectre
        sb_pw_cross_plot(net, em, a$x_norm, rec$Rdraw, x_name="generation", sub_s=sub_x, filename=fn)
  
        # plot pw 6-10 (that are not plotted with the spectre but might be interesting)
        for(i in sub_y) {
            fn_ <- sub("%", as.character(i), fn)
            postscript(fn_, width=5, height=5, family="serif")
            x <- jrnf_plot_pathway(net, em[i,1:nrow(net[[2]])], layout_f=layout.lgl, lim_plot=T)
            dev.off()
        }
    }
    setwd("..")
}


# Function plots the full evolved networks. For each state (value of Edraw) the
# ten most significant pathways are plotted into a separate folder. Explained
# fraction and explained dissipation are written in the filename. Additionally
# information derived from the pathway decomposition is plotted with
# "sb_plot_evol_pw".
#
# The function does not take more arguments but assumes that the results of a
# simulation are loaded in the global enviroment (as it is put there by 
# "sb_em_cross_analysis_ecol").

sb_plot_evol <- function() {
    rec <- results_em_cross

    # create a directory for each network / state 
    for(i in 1:nrow(rec)) 
        if(is.list(rec$em_ex[[i]])) { 
            # For each generation / network generate an extra subdirectory
            dir = paste("pw_details_E", as.character(results_em_cross$Edraw[i]), "_R", results_em_cross$Rdraw[i], sep="")
            system(paste("mkdir -p", dir))
            setwd(dir)

            Edraw <- rec$Edraw[i]
            net <- results_nets[[Edraw]]

            # plot the ten most important pathways
            # put information on explained fraction and explained dissipation in filename
            for(k in 1:min(10,nrow(rec$em_ex[[i]]))) {
                id <- rec$em_ex[[i]]$id[k]
                fn = paste("pw", as.character(k), "_id", as.character(id), "_er", rec$em_ex[[i]]$exp_r[k],
                           "_ed", rec$em_ex[[i]]$exp_d[k], ".eps", sep="")

                postscript(fn, width=5, height=5, family="serif")
                x <- jrnf_plot_pathway(net, em_m[[Edraw]][id,1:nrow(net[[2]])], layout_f=layout.lgl, lim_plot=T)
                dev.off()
            }

            # Plot explained fraction of rates and explained fraction of 
            # dissipation (decreasing for all wathways). Both are plotted in the
            # same plot as a function of the ordered pathway id. (Only first 50 
            # pathways are plotted!)

            K <- rec$em_ex[[i]]
            J <- data.frame(x=rep(1:nrow(K), 2),
                            f=c(K$exp_r, K$exp_d),
                            type=as.factor(c(rep("rates", nrow(K)), rep("dissipation", nrow(K))))) 

            J <- J[J$x<51,]

            # Plot with linear x and logarithmic y axis
            postscript("exp_f_linlog.eps", width=7, height=5)
            gg <- ggplot(J, aes(x=x, y=f, colour=type, group=type, shape=type))  + 
                scale_x_continuous() +
                scale_y_log10() +
                ylab("fraction") + 
                xlab("pathway no") +
                geom_point(size=3) + 
                theme_bw(22) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.key = element_rect(color = "white", fill="white"), legend.position="bottom")
        
             print(gg)
             dev.off()

            # Plot with both axis logarithmic
            postscript("exp_f_loglog.eps", width=7, height=5)
            gg <- ggplot(J, aes(x=x, y=f, colour=type, group=type, shape=type))  + 
                scale_x_log10() +
                scale_y_log10() +
                ylab("fraction") + 
                xlab("pathway no") +
                geom_point(size=3) + 
                theme_bw(22) +
                theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     legend.key = element_rect(color = "white", fill="white"), legend.position="bottom")
        
            print(gg)
            dev.off()
        
            setwd("..")
        }

    # plot change of important parameters in main directory:
    sb_plot_evol_pw(T)
}



# Function plots how pathway decomposition changes with some variable <x_data>. 
# Number of pathways (rows in <em_cross>) should be limited to 5.
#
# parameters:
# <em_cross>  - Explained rate fraction for same pathways in different simulations. 
#               Each row represents one simulation, each column one pathway.
# <x_data>    - Data for x-axis (pathway decomposition is plotted as a function of this)
# <x_name>    - X-axis label text
# <em_name>   - Names used for identifying different pathways (rows in em_cross)
# <filename>  - Name of output file
# <logY>      - Logarithmic y-scale?
# <logX>      - Logarithmic x-scale?

sb_plot_em_spectre <- function(em_cross, x_data, x_name="xname", em_name=c(), 
                               filename="em_spectre.eps", logY=F, logX=F) { 
    minY <- 1e-2
    a <- data.frame(x_data=numeric(), # x-axis
                    exp_r=numeric(),  # y-axis
                    em=factor())      # colour / group / shape
 
    for(i in 1:ncol(em_cross)) {
        a <- rbind(a, 
                   data.frame(x_data=as.numeric(x_data), 
                              exp_r=as.numeric(em_cross[,i]),
                              em=as.factor(rep(em_name[i], length(x_data)))))
    }

    a <- a[a$exp_r >= minY,]

    postscript(filename, width=7, height=5)

    # Decide on linear or logarithmic scale
    scaleY <- scale_y_continuous()
            
    if(logY)
        scaleY <- scale_y_log10(limits=c(minY,1)) 

    scaleX <- scale_x_continuous()
           
    if(logX)
        scaleX <- scale_x_log10() 


    gg<-ggplot(a, aes(x=x_data, y=exp_r, colour=em, group=em, shape=em))  + 
           scaleY + scaleX +
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


# This function uses "sb_plot_em_spectre" for plotting the dependency of the
# pathway decomposition of an external parameter. Thus, the number of pathways 
# (entries in <sub_s>) should be limited to 5. The function also makes separate
# plots of the topology of these 5 pathways.
#
# parameters:
# <net>      - The reaction network (non extended)
# <em>       - Set of elementary modes / pathways
# <cross>    - Explained pathways for same pathways in different simulations. 
#              Each row represents one simulation, each column one pathway.
# <x_data>   - Data for x-axis (pathway decomposition is plotted as a function of this)
# <x_name>   - X-axis label text
# <sub_s>    - Subset of interesting pathways (to plot) 
# <filename> - Name of the postscript output files. Needs to contain the 
#              character "%" which will be replaced with different values
#              for the different plots.
# <logX>     - Logarithmic x-scale?

sb_pw_cross_plot <- function(net, em, cross, x_data, x_name="xname", sub_s=c(), filename="pw_cross_%.eps", logX=F) {
    if(length(grep("%", filename)) == 0) {
        cat("WARNING: sb_pw_cross_plot:\n")
        cat("filename should contain placeholder '%' - output may be faulty!\n") 
    }

    # Extend net.
    y <- pa_extend_net(net, em[1,1:nrow(net[[2]])], unique_dir=T)
    net_x <- y[[1]] 

    # Ensure that <sub_s> contains indices of pathways
    if(is.null(sub_s)) 
        sub_s <- 1:nrow(em)

    if(is.logical(sub_s))
        sub_s <- which(sub_s)

    # First print spectrum of selected / subseted pathways
    em_name <- 1:ncol(cross)
    fn <- sub("%", "spec", filename)
    y <- sb_plot_em_spectre(cross[,sub_s], x_data, x_name, em_name[sub_s], filename=sub(".eps", "_log.eps", fn), logY=T, logX=logX) 
    y <- sb_plot_em_spectre(cross[,sub_s], x_data, x_name, em_name[sub_s], filename=sub(".eps", "_lin.eps", fn), logY=F, logX=logX)

    # Now plot all (selected) pathways in separate files
    for(i in sub_s) {
        postscript(sub("%", as.character(i), filename), width=5, height=5, family="serif")
        x <- jrnf_plot_pathway(net_x, em[i,], layout_f=layout.lgl, lim_plot=T)
        dev.off()
    }

    return(y)
}


# Simple plot function that plots the number of core species in dependency of
# driving force <v>.
#
# parameters:
# <res_cross> - Results object (from "sb_cross_analysis_ecol"!)
# <filename>  - Filename of output file

plot_core_sp <- function(res_cross, filename="core_sp.eps") { 
    setEPS()
    postscript(filename, width=7, height=5)
    
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


# Simple plot function that plots flow as a function as function of driving 
# force strength <v>.
#
# parameters:
# <res_cross> - Results object
# <filename>  - Filename of output file

plot_flow_v <- function(res_cross, filename="flow_v.eps") { 
    setEPS()
    postscript(filename, width=7, height=5)
    
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


#          LEGACY FUNCTION
#
# TODO Functions below is not functional any more because of some reengineering 
#      that was done. Nevertheless it is kept here as it gives valuable insight
#      in how results of pathway analysis may be layouted and how a 4 panel plot
#      can be created with ggplot.

evaluate_Edraw_plot <- function(Edraw, my_results=c(), my_em=c(), my_net=c(), filename="test.eps") {
    if(is.null(my_em))
        my_em <- em

    if(is.null(my_net))
        my_net <- net

    if(is.null(my_results))
        my_results <- results_em_cross

    my_results_X <- sb_extend_results(my_results)
    my_results_X <- my_results_X[my_results_X$Edraw == Edraw,]

    #x <- sb_assemble_em_from_flow(my_results_X, my_em, my_net)
    

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

    return(list(gg1, gg2, gg3, gg4)) 
}

