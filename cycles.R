# author: jakob fischer (jakob@automorph.info)
# date: 7. May 2014
# description: 
# Functions to count the occurences of cycles in directed graphs.  

library(igraph)


# Example-graph for testing
x <- data.frame( V1 = c(1,1,2,3,3,2,4,4,5,6), V2 = c(2,2,3,2,1,4,3,5,6,6) )
g_x <- graph(t(as.matrix(x)))


# The function calculates the number of cycles of length <n> in the directed 
# graph <g> and the number of occurences of each edge in this cycles by using
# the graph subisomorphism method of igraph. Because the subisomorphism is 
# defined in terms of nodes multiple edges between two nodes don't imply the
# cycle being counted multiple times.

get_n_cycles_directed_A <- function(g, n) {
    gc()                              # Try using garbace gollector because function tends to segfault (igraph)
    v_count <- rep(0, length(V(g)))   # vector for counting occurence of edges in cycles (initialized 0)

    # if n==1 just count loops
    if(n == 1) {
        loops <- which(is.loop(g))
        for(i in loops) { 
           v <- get.edge(g, E(g)[i])[1]
           v_count[v] <- v_count[v] + 1 
        }

        return(list(length(loops), v_count))   
    }

    # find subisomorphisms and count occurence of nodes
    si <- graph.get.subisomorphisms.vf2(g, graph.ring(n, directed=T))
    for(i in si) 
        v_count[i+1] <- v_count[i+1] + 1
    
    # first element of list is number of subisomorphisms, second number of vertices occurence
    # has to be divided by n because there are n equivalent subisomorphisms
    return(list(length(si)/n, v_count/n))
}


# get_n_cycles_directed_B fixes the problem of not counting cycles with multiple
# edges multiple times by replacing each edge by two edges with one "virtual" 
# node and using get_n_cycles_directed_A with n=2*n

get_n_cycles_directed_B <- function(g, n) {
    gc()                         # Try using garbace gollector because function tends to segfault (igraph)
    vertex_no <- length(V(g))    # number of vertices
    edge_no <- length(E(g))      # and edges
    el_new <- matrix(0, edge_no*2, 2)   # matrix containing edges of transformed network

    nxt_new <- vertex_no+1              # index of next vertex 

    # split edges into two with "virtual" node inbetween
    for(i in 1:edge_no) {
        el_new[i*2-1,1] <- get.edge(g,i)[1]
        el_new[i*2-1,2] <- nxt_new

        el_new[i*2,1] <- nxt_new
        el_new[i*2,2] <- get.edge(g,i)[2]
        
        nxt_new <- nxt_new + 1
    }

    # Go back to graph representation and use (original) get_n_cycles_directed_A function
    g_new <- graph(el_new)
    x <- get_n_cycles_directed_A(g_new, n*2)
    return(list(x[[1]], x[[2]][1:vertex_no]))
}


# Standard method for calculating cycle occurences in directed graphs (see above)

get_n_cycles_directed <- get_n_cycles_directed_B
