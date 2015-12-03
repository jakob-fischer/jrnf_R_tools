# author: jakob fischer (jakob@automorph.info)
# description: 
# Script for generating files and directory structure of artificial chemistry
# simulations with thermodynamic constraints. (simulation builder = "sb_")
# This file contains the core functionality to build and evaluate wide ranging
# simulations - entire directories with multiple changing parameters. The 
# functions to plot (and partially evaluate) the results is in "simulation_builder_plot.R".

sourced_simulation_builder <- T

if(!exists("sourced_simulation_builder_plot"))
    source("simulation_builder_plot.R")
