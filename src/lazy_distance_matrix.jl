# Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
# Microsoft Software License Terms for Microsoft Quantum Development Kit Libraries
# and Samples. See LICENSE in the project root for license information.

"""
    LazyDistance

An array-like structure that represents unused edges with weight 1, and used edges with weight nv(g)
"""
struct LazyDistance <: AbstractMatrix{Int}
    graph :: AbstractGraph
    # Store used edges in another graph
    # Assuming that LightGraphs has an efficient datastructure to store this kind of information,
    # so I reuse it.
    used_graph :: SimpleGraph 

    LazyDistance(g :: AbstractGraph) = new(g, SimpleGraph(nv(g)))
end

Base.size(ld :: LazyDistance) = size(ld.used_graph)

setused!(ld :: LazyDistance, e) = add_edge!(ld.used_graph, e)
setused!(ld :: LazyDistance, s :: Integer, d :: Integer) = add_edge!(ld.used_graph, s, d)

function Base.getindex(ld :: LazyDistance, s :: Integer, d :: Integer) :: Int
    if has_edge(ld.used_graph, s, d)
        return nv(ld.graph)
    else
        return has_edge(ld.graph, s, d)
    end
end