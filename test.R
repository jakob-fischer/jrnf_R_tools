library(igraph)
source("jrnf_network.R")

enthalpy <- read.csv("enthalpy.csv")[,2:4]


kasting <- load_jrnf("networks/kasting.jrnf")
kasting_dn <- jrnf_to_directed_network(kasting)
kasting_dn_dd <- degree.distribution(kasting_dn, mode="in", cumulative=TRUE)
kasting_ud <- as.undirected(kasting_dn)
kasting_ud_dd <- degree.distribution(kasting_ud, mode="in", cumulative=TRUE)
kasting_dn_s <- simplify(kasting_dn, remove.multiple=TRUE, remove.loops=TRUE)
kasting_dn_s_dd <- degree.distribution(kasting_dn_s, mode="in", cumulative=TRUE)
kasting_ud_s <- simplify(kasting_ud, remove.multiple=TRUE, remove.loops=TRUE)
kasting_ud_s_dd <- degree.distribution(kasting_ud_s, mode="in", cumulative=TRUE)

yung_demore <- load_jrnf("networks/yung_demore.jrnf")
yung_demore_dn <- jrnf_to_directed_network(yung_demore)
yung_demore_dn_dd <- degree.distribution(yung_demore_dn, mode="in", cumulative=TRUE)
yung_demore_ud <- as.undirected(yung_demore_dn)
yung_demore_ud_dd <- degree.distribution(yung_demore_ud, mode="in", cumulative=TRUE)
yung_demore_dn_s <- simplify(yung_demore_dn, remove.multiple=TRUE, remove.loops=TRUE)
yung_demore_dn_s_dd <- degree.distribution(yung_demore_dn_s, mode="in", cumulative=TRUE)
yung_demore_ud_s <- simplify(yung_demore_ud, remove.multiple=TRUE, remove.loops=TRUE)
yung_demore_ud_s_dd <- degree.distribution(yung_demore_ud_s, mode="in", cumulative=TRUE)


ref_er_kasting_ud_s <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/erN31M164.jrnf")))
ref_er_kasting_ud_s_dd <- degree.distribution(ref_er_kasting_ud_s, mode="in", cumulative=TRUE)
ref_ba_kasting_ud_s <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/baN31M164.jrnf")))
ref_ba_kasting_ud_s_dd <- degree.distribution(ref_ba_kasting_ud_s, mode="in", cumulative=TRUE)
ref_ws_kasting_ud_s <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/wsN31M164.jrnf")))
ref_ws_kasting_ud_s_dd <- degree.distribution(ref_ws_kasting_ud_s, mode="in", cumulative=TRUE)
ref_ps_kasting_ud_s <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/psN31M164.jrnf")))
ref_ps_kasting_ud_s_dd <- degree.distribution(ref_ps_kasting_ud_s, mode="in", cumulative=TRUE)

ref_er_kasting_ud <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/erN31M176m.jrnf")))
ref_er_kasting_ud_dd <- degree.distribution(ref_er_kasting_ud, mode="in", cumulative=TRUE)
ref_ba_kasting_ud <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/baN31M176m.jrnf")))
ref_ba_kasting_ud_dd <- degree.distribution(ref_ba_kasting_ud, mode="in", cumulative=TRUE)
ref_ws_kasting_ud <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/wsN31M176m.jrnf")))
ref_ws_kasting_ud_dd <- degree.distribution(ref_ws_kasting_ud, mode="in", cumulative=TRUE)
ref_ps_kasting_ud <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/psN31M176m.jrnf")))
ref_ps_kasting_ud_dd <- degree.distribution(ref_ps_kasting_ud, mode="in", cumulative=TRUE)

ref_er_yung_demore_ud_s <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/erN280M893.jrnf")))
ref_er_yung_demore_ud_s_dd <- degree.distribution(ref_er_yung_demore_ud_s, mode="in", cumulative=TRUE)
ref_ba_yung_demore_ud_s <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/baN280M893.jrnf")))
ref_ba_yung_demore_ud_s_dd <- degree.distribution(ref_ba_yung_demore_ud_s, mode="in", cumulative=TRUE)
ref_ws_yung_demore_ud_s <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/wsN280M893.jrnf")))
ref_ws_yung_demore_ud_s_dd <- degree.distribution(ref_ws_yung_demore_ud_s, mode="in", cumulative=TRUE)
ref_ps_yung_demore_ud_s <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/psN280M893.jrnf")))
ref_ps_yung_demore_ud_s_dd <- degree.distribution(ref_ps_yung_demore_ud_s, mode="in", cumulative=TRUE)

ref_er_yung_demore_ud <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/erN280M1185m.jrnf")))
ref_er_yung_demore_ud_dd <- degree.distribution(ref_er_yung_demore_ud, mode="in", cumulative=TRUE)
ref_ba_yung_demore_ud <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/baN280M1185m.jrnf")))
ref_ba_yung_demore_ud_dd <- degree.distribution(ref_ba_yung_demore_ud, mode="in", cumulative=TRUE)
ref_ws_yung_demore_ud <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/wsN280M1185m.jrnf")))
ref_ws_yung_demore_ud_dd <- degree.distribution(ref_ws_yung_demore_ud, mode="in", cumulative=TRUE)
ref_ps_yung_demore_ud <- as.undirected(jrnf_to_directed_network(load_jrnf("ref_nws/psN280M1185m.jrnf")))
ref_ps_yung_demore_ud_dd <- degree.distribution(ref_ps_yung_demore_ud, mode="in", cumulative=TRUE)










std_analysis1_jrnf <- function(rn, name, enthalpy=NULL) {
    rn_dn <- jrnf_to_directed_network(rn)
    rn_dn_s <- simplify(rn_dn, remove.multiple=TRUE, remove.loops=TRUE)
    rn_ud <- as.undirected(rn_dn, mode="each")
    rn_ud_s <- simplify(rn_ud, remove.multiple=TRUE, remove.loops=TRUE)



}


std_analysis2_jrnf <- function(rn, name, enthalpy=NULL) {
    rn_red_sp <- jrnf_get_s_con_subnet(rn, jrnf_get_s_con_subnet(rn), "s")
    rn_red_re <- jrnf_get_s_con_subnet(rn, jrnf_get_s_con_subnet(rn), "r") 
}

# After one plottet a 
#





#ebc <- edge.betweenness.community(rn_kasting_dn_scred)
#enthalpy_kastings <- get_fenthalpy_species(rn_kasting_dn_scred, enthalpy)

#postscript("kasting.ps", fonts=c("serif", "Palatino"))
#plot(rn_kasting_dn_scred)
#par(mar=c(0,0,0,0))
#plot(ebc, rn_kasting_dn_s_scred)
#dev.off()

#V(rn_kasting_dn_s_scred)$
#plot(ebc, rn_kasting_dn_s_scred)
#V(rn_kasting_dn_s_scred)$color <- get_color(h, enthalpy_kastings)


#postscript("kasting_nodes_by_enthalpy.ps", fonts=c("serif", "Palatino"))
#plot(rn_kasting_dn_s_scred)
#dev.off()





postscript("scaling_ud.ps", fonts=c("serif", "Palatino"))
par(mfrow = c(2,2))

plot(kasting_ud_dd, log="xy", xlab="degree", ylab="cum. frequency", col=1, main="Kasting-M.")
lines(ref_er_kasting_ud_dd, col="red")
lines(ref_ws_kasting_ud_dd, col="blue")
lines(ref_ba_kasting_ud_dd, col="green")
lines(ref_ps_kasting_ud_dd, col="orange")
legend(x=1, y=0.5, legend="red - Erdos Renyi\n blue - Watts-Strogatz\n green - Barabasi-Albert\n orange - Pan Sinha")

plot(kasting_ud_s_dd, log="xy", xlab="degree", ylab="cum. frequency", col=1, main="Kasting-M. (simple)")
lines(ref_er_kasting_ud_s_dd, col="red")
lines(ref_ws_kasting_ud_s_dd, col="blue")
lines(ref_ba_kasting_ud_s_dd, col="green")
lines(ref_ps_kasting_ud_dd, col="orange")

plot(yung_demore_ud_dd, log="xy", xlab="degree", ylab="cum. frequency", col=1, main="Yung DeMore")
lines(ref_er_yung_demore_ud_dd, col="red")
lines(ref_ws_yung_demore_ud_dd, col="blue")
lines(ref_ba_yung_demore_ud_dd, col="green")
lines(ref_ps_kasting_ud_dd, col="orange")

plot(yung_demore_ud_s_dd, log="xy", xlab="degree", ylab="cum. frequency", col=1, main="Yung DeMore (simple)")
lines(ref_er_yung_demore_ud_dd, col="red")
lines(ref_ws_yung_demore_ud_dd, col="blue")
lines(ref_ba_yung_demore_ud_dd, col="green")
lines(ref_ps_kasting_ud_dd, col="orange")
dev.off()
