source("jrnf_network.R")

enthalpy <- read.csv("enthalpy.csv")[,2:4]



postscript("formation_enthalpy_hist.ps", fonts=c("serif", "Palatino"))
par(ps=22)
hist_dHF <- hist(enthalpy$dHf, breaks=20, col=rainbow(20), xlab="Formation Enthalpy dH0  (kcal/mol)", main="Formation Enthalpies")
dev.off()



postscript("standard_entropy_hist.ps", fonts=c("serif", "Palatino"))
par(ps=22)
hist_S0 <- hist(enthalpy$S0, breaks=20, col=rainbow(20), xlab="Standard Enthalpy S0 (cal/molK)", main="Standard Entropies")
dev.off()



kasting_s <- load_jrnf("networks/kasting_s.jrnf")
kasting_s_dn <- jrnf_to_directed_network(kasting_s)
kasting_s_dn_dd <- degree.distribution(kasting_s_dn, mode="in", cumulative=TRUE)

E(kasting_s_dn)$weight <- 1
kasting_s_dn <- simplify(kasting_s_dn, remove.loops=TRUE)
kasting_s_dn_dH0 <- get_fenthalpy_species(V(kasting_s_dn)$name, enthalpy)
V(kasting_s_dn)$color <- get_color_by_hist(hist_dHF, kasting_s_dn_dH0)

postscript("kasting_net_ent.ps", fonts=c("serif", "Palatino"))
par(ps=18)
plot(kasting_s_dn)
dev.off()


yung_demore <- load_jrnf("networks/yung_demore.jrnf")
yung_demore_dn <- jrnf_to_directed_network(yung_demore)
yung_demore_dd <- degree.distribution(yung_demore_dn, cumulative=TRUE)

postscript("kasting_net_ent_dd.ps", fonts=c("serif", "Palatino"))
par(ps=22)
plot(yung_demore_dd, main="Data on AC of Yung + DeMore", xlab="degree", ylab="acc. degree distribution", log="xy")
dev.off()

