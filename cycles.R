# author: jakob fischer (jakob@automorph.info)
# description: 
# Functions to count the occurences of cycles in directed graphs.

sourced_cycles <- T

library(igraph)


# Example-graph for testing
cycles.R_x <- data.frame( V1 = c(1,1,2,3,3,2,4,4,5,6), V2 = c(2,2,3,2,1,4,3,5,6,6) )
cycles.R_g_x <- graph(t(as.matrix(cycles.R_x)))


# The function calculates the number of cycles of length <n> in the directed 
# graph <g> and the number of occurences of each edge in this cycles by using
# the graph subisomorphism method of igraph. Because the subisomorphism is 
# defined in terms of nodes, multiple edges between two nodes don't imply the
# cycle being counted multiple times. The suffix "_V" indicates thus that the
# cycles are identified in respect to the vertices. If <list_cycles> is set
# then the third part of the returned list is a matrix containing for which 
# each row indicates the vertices taking part in one cycle.

get_n_cycles_directed_V <- function(g, n, list_cycles=F) { 
    gc()                              # Try using garbace gollector because function tends to segfault (igraph)
    v_count <- rep(0, length(V(g)))   # vector for counting occurence of nodes in cycles (initialized 0)

    # if n==1 just count loops
    if(n == 1) {
        loops <- which(is.loop(g))        # which edges are loops?

        for(i in loops) {                 
           v <- ends(g, i, names=F)[1,1]           # For each, the first node is the only one
           v_count[v] <- 1                # Count only once for each node (to be consistent with n>1) 
        }

        # Listing cycles is simple here, just find the vertices that have nonzero
        # loops (v_count != 0) and add one row for each of them in cycles matrix.
        if(list_cycles) {
            u_loops_s <- which(v_count != 0)
            cycles <- matrix(u_loops_s, ncol=1)  
            cycles_m <- rep(1, length(u_loops_s))
        }
    # handles the general case of n > 1
    } else {
        # find subisomorphisms of directed graph with ring of length n
        si <- graph.get.subisomorphisms.vf2(g, graph.ring(n, directed=T))
        si <- subisomorphism_rm_permutation(si)  # remove permutations

        for(i in si)
            v_count[i] <- v_count[i] + 1

        if(list_cycles) {
            if(length(si) != 0)
                cycles <- do.call(rbind, si)
            else
                cycles <- matrix(0, ncol=n, nrow=0)    
        
            cycles_m <- rep(1, nrow(cycles))
        }
    }

    if(list_cycles)
        return(list(count=sum(v_count)/n, v_incidence=v_count, cycles=cycles, c_mul=cycles_m))
    else
        return(list(count=sum(v_count)/n, v_incidence=v_count))
}


# get_n_cycles_directed_B fixes the problem of not counting cycles with multiple
# edges multiple times by replacing each edge by two edges with one "virtual" 
# node and using get_n_cycles_directed_A with n=2*n
#
# This is a version of "get_n_cycles_directed_B" that is optimized for speed. For this,
# instead of introducing additional "in between" nodes to avoid multiple edges all the
# multiple edges are removed and their effect on the number of cycles is just calculated
# (multiplied) after the subisomorphism algorithm has been invoked.

get_n_cycles_directed_E <- function(g, n, list_cycles=F) {
    gc()                              # Try using garbace gollector because function tends to segfault (igraph)
    v_count <- rep(0, length(V(g)))   # vector for counting occurence of nodes in cycles (initialized 0)

    if(n == 1) {
        loops <- which(is.loop(g))        # which edges are loops?

        for(i in loops) {                 
           v <- ends(g, i, names=F)[1,1]           # For each, the first node is the only one
           v_count[v] <- v_count[v] + 1   # 
        }
        
        if(list_cycles) {
            cycles <- matrix(0, ncol=length(V(g)), nrow=sum(v_count != 0))  
            cycles_m <- rep(0, nrow(cycles))
            
            u_loops_s <- which(v_count != 0)
            cycles <- matrix(u_loops_s, ncol=1)  
            cycles_m <- rep(1, length(u_loops_s))
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

        if(list_cycles) {
            if(length(si) != 0)
                cycles <- do.call(rbind, si)
            else
                cycles <- matrix(0, ncol=n, nrow=0)   
        }  
    }

    if(list_cycles)
        return(list(count=sum(v_count)/n, v_incidence=v_count, cycles=cycles, c_mul=cycles_m))
    else
        return(list(count=sum(v_count)/n, v_incidence=v_count))
}



# Standard method for calculating cycle occurences in directed graphs (see above)
get_n_cycles_directed <- get_n_cycles_directed_E
