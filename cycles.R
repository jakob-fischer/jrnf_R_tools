# author: jakob fischer (jakob@automorph.info)
# date: 7. May 2014
# description: 
# Functions to count the occurences of cycles in directed graphs.

sourced_cycles <- T

library(igraph)


# Example-graph for testing
x <- data.frame( V1 = c(1,1,2,3,3,2,4,4,5,6), V2 = c(2,2,3,2,1,4,3,5,6,6) )
g_x <- graph(t(as.matrix(x)))


# The function calculates the number of cycles of length <n> in the directed 
# graph <g> and the number of occurences of each edge in this cycles by using
# the graph subisomorphism method of igraph. Because the subisomorphism is 
# defined in terms of nodes, multiple edges between two nodes don't imply the
# cycle being counted multiple times.

get_n_cycles_directed_A <- function(g, n, list_cycles=F) {
    gc()                              # Try using garbace gollector because function tends to segfault (igraph)
    v_count <- rep(0, length(V(g)))   # vector for counting occurence of nodes in cycles (initialized 0)
    res <- list()

    # if n==1 just count loops
    if(n == 1) {
        loops <- which(is.loop(g))        # which edges are loops?

        for(i in loops) {                 
           v <- get.edge(g, E(g)[i])[1]   # For each, the first node is the only one
           v_count[v] <- 1                # Count only once for each node (to be consistent with n>1) 
        }

        res <- list(sum(v_count != 0), v_count)
        
        # Listing cycles is simple here, just find the vertices that have nonzero
        # loops (v_count != 0) and add one row for each of them in cycles matrix.
        if(list_cycles) {
            u_loops_s <- which(v_count != 0)
            cycles <- matrix(0, ncol=length(V(g)), nrow=length(u_loops_s))  
            
            if(length(u_loops_s) != 0)
                for(i in 1:length(u_loops_s))
                    cycles[i,u_loops_s[i]] <- 1

            res[[3]] <- cycles
        }
    # handles the general case of n > 1
    } else {
        # find subisomorphisms of directed graph with ring of length n
        si <- graph.get.subisomorphisms.vf2(g, graph.ring(n, directed=T))
        si <- subisomorphism_rm_permutation(si)  # remove permutations

        for(i in si)
            v_count[i] <- v_count[i] + 1

        # first element of list is number of subisomorphisms, second is the number of vertices occurence,
        # both have to be divided by n because there are n equivalent subisomorphisms
        res <- list(length(si), v_count)

        if(list_cycles) {
            cycles <- matrix(0, ncol=length(V(g)), nrow=length(si))  
            
            if(length(si) != 0)
                for(i in 1:length(si))
                    cycles[ i,si[[i]] ] <- 1

            res[[3]] <- cycles
        }
    }

    return(res)
}


# get_n_cycles_directed_B fixes the problem of not counting cycles with multiple
# edges multiple times by replacing each edge by two edges with one "virtual" 
# node and using get_n_cycles_directed_A with n=2*n
#
# TODO There is a faster implementation below. If it is sufficiently tested this 
# one can be removed!

get_n_cycles_directed_B <- function(g, n, list_cycles=F) {
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
    res <- get_n_cycles_directed_A(graph(el_new), n*2, list_cycles)
    res[[2]] <- res[[2]][1:vertex_no]
    if(list_cycles) res[[3]] <- res[[3]][,1:vertex_no]

    return(res)
}



# This is a version of "get_n_cycles_directed_B" that is optimized for speed. For this,
# instead of introducing additional "in between" nodes to avoid multiple edges all the
# multiple edges are removed and their effect on the number of cycles is just calculated
# (multiplied) after the subisomorphism algorithm has been invoked.

get_n_cycles_directed_B_fast <- function(g, n, list_cycles=F) {
    gc()                              # Try using garbace gollector because function tends to segfault (igraph)
    v_count <- rep(0, length(V(g)))   # vector for counting occurence of nodes in cycles (initialized 0)
    res <- list()

    if(n == 1) {
        loops <- which(is.loop(g))        # which edges are loops?

        for(i in loops) {                 
           v <- get.edge(g, E(g)[i])[1]   # For each, the first node is the only one
           v_count[v] <- v_count[v] + 1   # 
        }

        res <- list(sum(v_count), v_count)
        
        if(list_cycles) {
            cycles <- matrix(0, ncol=length(V(g)), nrow=sum(v_count != 0))  
            cycles_m <- rep(0, nrow(cycles))
            
            u_loops_s <- which(v_count != 0)
            if(length(u_loops_s) != 0)
                for(i in 1:length(u_loops_s)) {
                    cycles[i,u_loops_s[i]] <- 1
                    cycles_m[i] <- v_count[u_loops_s[i]]
                }

            res[[3]] <- cycles
            res[[4]] <- cycles_m
        }
    } else {
        m <- jrnf_graph_to_amatrix(g)                          # calculate graph adjacency matrix
        g_s <- simplify(g, remove.multiple=T, remove.loops=T)  # remove multiple edges and loops

        # find subisomorphisms and count occurence of nodes
        si <- graph.get.subisomorphisms.vf2(g, graph.ring(n, directed=T))
        si <- subisomorphism_rm_permutation(si)  # remove permutations
  
        cycles_m <- c()    

        for(i in si) { 
            k <- 1  # number of variants to traverse the multiple edges for this subisomorphism / cycle

            for(j in 1:length(i))    # just multiply all the values of the adjacency matrix in the cycle 
                if(j != length(i))            
                    k <- k * m[i[j], i[j+1]]
                else
                    k <- k * m[i[length(i)], i[1]]    # edge between first and last node

            v_count[i] <- v_count[i] + k    # accumulate number of cycles individual nodes are part of
            cycles_m <- c(cycles_m, k)
        }

        res <- list(sum(cycles_m), v_count)

        if(list_cycles) {
            cycles <- matrix(0, ncol=length(V(g)), nrow=length(cycles_m)) 

            if(length(si) != 0)
                for(i in 1:length(si))
                    cycles[ i,si[[i]] ] <- 1

            res[[3]] <- cycles  
            res[[4]] <- cycles_m
        }  
    }

    return(res)
}



# Standard method for calculating cycle occurences in directed graphs (see above)
get_n_cycles_directed <- get_n_cycles_directed_B_fast
