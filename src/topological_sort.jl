# Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
# Microsoft Software License Terms for Microsoft Quantum Development Kit Libraries
# and Samples. See LICENSE in the project root for license information.

struct DAGTopologicalVisitor{T <: Integer}
    dag :: LightGraphs.SimpleDiGraph{T}
    visited :: Vector{Bool}
    frontier :: Set{T}
end

function DAGTopologicalVisitor(dag :: LightGraphs.SimpleDiGraph{T}) where T <: Integer
    visited = fill(false, LightGraphs.nv(dag))
    frontier = Set{T}()

    for v in LightGraphs.vertices(dag)
        if isfrontier(dag, visited, v)
            push!(frontier, v)
        end
    end
    DAGTopologicalVisitor(dag, visited, frontier)
end

function isfrontier(dag :: LightGraphs.SimpleDiGraph{T}, visited :: Vector{Bool}, v :: T) where T <: Integer
    # Check if all inneighbors have been visited
    return all(getindex(visited, LightGraphs.inneighbors(dag, v)))
end

isfrontier(dagvisitor :: DAGTopologicalVisitor{T}, v :: T) where T = isfrontier(dagvisitor.dag, dagvisitor.visited, v)
function visit!(dagvisitor :: DAGTopologicalVisitor{T}, v :: T) where T
    pop!(dagvisitor.frontier, v)
    dagvisitor.visited[v] = true

    for neighbor in LightGraphs.outneighbors(dagvisitor.dag, v)
        if isfrontier(dagvisitor, neighbor)
            push!(dagvisitor.frontier, neighbor)
        end
    end
end

"Topologically sort the the vertices in the given DAG"
function topological_sort(dag :: LightGraphs.SimpleDiGraph{T}) :: Vector{T} where T <: Integer
    sorted = Vector{T}()
    sizehint!(sorted, LightGraphs.nv(dag))

    visitor = DAGTopologicalVisitor(dag)

    while !isempty(visitor.frontier)
        current = first(visitor.frontier)
        push!(sorted, current)
        visit!(visitor, current)
    end
    if length(sorted) != nv(dag)
        throw(ArgumentError("The graph has cycles."))
    end
    return sorted
end

"Find all vertices in the DAG that do not have unvisited predecessors"
function layers(
    dag :: LightGraphs.SimpleDiGraph{T};
) :: Vector{Vector{T}} where T <: Integer
    layers = Vector{Vector{T}}()
    visitor = DAGTopologicalVisitor(dag)

    while !isempty(visitor.frontier)
        # Visit all vertices in the current frontier
        layer = collect(visitor.frontier)
        push!(layers, layer)
        for v in layer
            visit!(visitor, v)
        end
    end
    if sum(length, layers) != nv(dag)
        throw(ArgumentError("The graph has cycles."))
    end
    return layers
end