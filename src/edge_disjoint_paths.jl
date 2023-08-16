# Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
# Microsoft Software License Terms for Microsoft Quantum Development Kit Libraries
# and Samples. See LICENSE in the project root for license information.

"A pair of terminal vertices, the first vertex being the origin (from)."
const TerminalPair = Tuple{Int,Int}
const Path = Vector{Int}

"""Get the from field of a TerminalPair"""
function from(terminal::TerminalPair)
    terminal[1]
end
"""Get the to field of a TerminalPair"""
function to(terminal::TerminalPair)
    terminal[2]
end

"""
Remove the edges in a path from the graph.
"""
function rem_path!(graph, path::Path)::Nothing
    # Delete edges used along the path
    for i = 1:length(path)-1
        from = path[i]
        to = path[i+1]
        # In directed graph remove edges in both directions.
        rem_edge!(graph, from, to)
        rem_edge!(graph, to, from)
    end
end

"""
Find the shortest path between terminals using A*

The heuristic_dist function takes the source and destination vertices.
"""
function a_star(
    graph,
    terminal::TerminalPair;
    distmx::AbstractMatrix{T}=weights(graph),
    heuristic_dist::Function=(n, t) -> 0
)::Array{Edge{Int}} where {T}
    # Partially apply the destination vertex.
    heur_dist_dest = from -> heuristic_dist(from, to(terminal))
    LightGraphs.a_star(graph, from(terminal), to(terminal), distmx, heur_dist_dest)
end

"""
Greedily find shortest paths and construct a set of edge-disjoint paths.

O(√|E|)-approximation algorithm from Kolliopolis & Stein (2004) "Approximating disjoint-path problems
using packing integer programs"
O(√|E|)-approximation algorithms from Kleinberg & Tardos "Algorithm Design" §11.5
"""
function greedy_edp(
    graph,
    terminals;
    heuristic_dist::Function=(n, t) -> 0
)::Vector{Path}
    graph = copy(graph) # TODO: Get rid of copy
    terminals = Set(terminals)
    ed_paths = Path[]

    while true
        # Find shortest paths between all remaining terminal vertices
        paths = Path[]
        sizehint!(paths, length(terminals))
        for (i, terminal) in enumerate(terminals)
            shortest_path =
                a_star(graph, terminal, heuristic_dist=heuristic_dist)
            # Terminals cannot be connected so we delete them
            if isempty(shortest_path)
                delete!(terminals, terminal)
            else
                # Simplify path to vertex ids.
                vertex_path = [shortest_path[1].src; dst.(shortest_path)]
                push!(paths, vertex_path)
            end
        end

        # If we cannot connect terminals any more, stop.
        if isempty(paths)
            break
        end

        # Find which terminal pair was closest and route that.
        shortest_idx = getindex.(paths, 2) .|> length |> argmin
        shortest_path = paths[shortest_idx]
        push!(ed_paths, shortest_path)

        # Remove routed terminal from set of available terminals
        shortest_terminal = (shortest_path[1], shortest_path[end])
        @assert shortest_terminal ∈ terminals
        delete!(terminals, shortest_terminal)

        rem_path!(graph, shortest_path)
    end

    return ed_paths
end

"""
Faster greedy algorithm for finding edge-disjoint paths.

The bounded version is a O(√|E|)-approximation algorithm on mesh
according to Kleinberg's dissertation "Approximation Algorithms for Disjoint
Paths Problems" Theorem 5.1.1.
Here we cannot use the bound because all paths must be routed at some point.
"""
function greedier_edp(
    graph,
    terminals;
    heuristic_dist::Function=(n, t) -> 0,
    rng=Random.default_rng()
)::Vector{Path}
    distmx = LazyDistance(graph)
    ed_paths = Path[]
    terminals = collect(terminals)

    # Proceed through terminals in random order
    Random.shuffle!(rng, terminals)
    for terminal in terminals
        shortest_path = a_star(graph, terminal; distmx=distmx, heuristic_dist=heuristic_dist)
        # Check if path is longer than graph size, if so we have reused an edge
        distance = 0
        for edge in shortest_path
            distance += distmx[edge.src, edge.dst]
        end
        # Terminals cannot be connected so we do not use them
        if isempty(shortest_path) || distance ≥ nv(graph)
            continue
        else
            # Simplify path to vertex ids.
            vertex_path = [shortest_path[1].src; dst.(shortest_path)]
            push!(ed_paths, vertex_path)

            for edge in shortest_path
                setused!(distmx, edge)
            end
        end
    end

    return ed_paths
end
