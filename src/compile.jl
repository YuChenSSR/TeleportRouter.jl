# Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
# Microsoft Software License Terms for Microsoft Quantum Development Kit Libraries
# and Samples. See LICENSE in the project root for license information.

"""
Sample unique random pairs from the given collection

Pairs is rounded down in the case of odd length.
"""
function random_pairs(collection)::Array{Tuple{Int,Int}}
    shuffled = Random.shuffle(collection)
    if isodd(length(shuffled))
        pop!(shuffled)
    end
    return Tuple.(Iterators.partition(shuffled, 2))
end

"""Convert linear index into Cartesian index for the grid"""
function grid_coord(id::Int, width::Int)::Tuple{Int,Int}
    (((id - 1) % width) + 1, ((id - 1) ÷ width) + 1)
end

"""Convert an ID to an ID on the scaled grid with ancilla."""
function ancilla_grid_id(id::Int, width::Int)::Int
    scaled_width = 2 * (width + 2) - 1
    x = 2 * ((id - 1) % width) + 3 # Add two columns left of every id
    y = 2 * ((id - 1) ÷ width) + 3 # Add two rows above every id
    return x + (y - 1) * scaled_width
end

function snake_id(id::Int, width::Int)
    if iseven((id - 1) ÷ width)
        return id
    else
        y = (id - 1) ÷ width
        x = width - ((id - 1) % width)
        return x + y * width
    end
end

"""
Construct a graph with outgoing directed edges from sources onto grid
and ingoing directed edges only from the grid and boundary qubits.
We insert ancilla between data qubits.
The boundary will add one row and column of logical qubits (with ancilla separating
it from the bulk of the grid).
"""
function ancilla_grid(
    sources::Vector{Int},
    destinations::Vector{Int},
    width::Int,
    height::Int,
)::DiGraph{Int}
    function id(id::Int)::Int
        ancilla_grid_id(id, width)
    end
    # Add two columns around bulk to add boundary, and intersperse with ancilla
    scaled_width = 2 * (width + 2) - 1
    scaled_height = 2 * (height + 2) - 1
    graph = DiGraph(LightGraphs.grid([scaled_width, scaled_height]))

    scaled_sources = id.(sources) # Scale to new graph
    scaled_destinations = id.(destinations)

    # Verify input
    badsource = findfirst(scaled_sources .> nv(graph))
    if !isnothing(badsource)
        throw(ArgumentError("Source $(sources[badsource]) is out of bounds ($(scaled_sources[badsource]) ∉ $scaled_width×$scaled_height)"))
    end
    baddest = findfirst(destinations .> nv(graph))
    if !isnothing(baddest)
        throw(ArgumentError("Source $(destinations[baddest]) is out of bounds ($(scaled_destinations[baddest]) ∉ $scaled_width×$scaled_height)"))
    end

    # Add boundary vertices as destinations
    # Rows
    append!(scaled_destinations, 1:2:scaled_width)
    append!(scaled_destinations, nv(graph):-2:(nv(graph)-scaled_width))
    # Columns
    append!(scaled_destinations, 1:scaled_width:nv(graph))
    append!(scaled_destinations, scaled_width:scaled_width:nv(graph))

    for vertex in vertices(graph)
        # Either x or y is even
        if any(iseven.(grid_coord(vertex, scaled_width)))
            continue
        end
        if vertex in scaled_sources
            # Collect neighborhood because we need to mutate it.
            for neighbor in collect(outneighbors(graph, vertex))
                # Remove incoming edge
                rem_edge!(graph, neighbor, vertex)
                # Remove outgoing edges left and right (can only control from top & bottom)
                if abs(neighbor - vertex) == 1
                    rem_edge!(graph, vertex, neighbor)
                end
            end
        elseif vertex in scaled_destinations
            for neighbor in collect(inneighbors(graph, vertex))
                # Remove outgoing edge
                rem_edge!(graph, vertex, neighbor)
                # Remove ingoing edges top and bottom (can only target from left & right)
                if abs(neighbor - vertex) != 1
                    rem_edge!(graph, neighbor, vertex)
                end
            end
        else
            for neighbor in collect(all_neighbors(graph, vertex))
                # Remove incoming and outgoing edge (not participating)
                rem_edge!(graph, vertex, neighbor)
                rem_edge!(graph, neighbor, vertex)
            end
        end
    end
    return graph
end

"""Give the list of logical boundary qubits in a given graph

Note that the width should be include ancilla and boundary"""
function boundary_qubits(width::Int, height::Int)::Vector{Int}
    nrnodes = width * height
    qs = collect(1:2:width)
    append!(qs, height*width:-2:(height-1)*width)
    append!(qs, 1:2*width:nrnodes)
    append!(qs, width:2*width:nrnodes)
    return unique(qs)
end

"""
Given disjoint CNOTs, find an edge-disjoint paths schedule to perform them on the grid.

One time step of a schedule is a set of edge-disjoint paths can be performed
in the same time step using constant time.
The schedule consists of multiple such time steps
"""
function apply_disjoint_cnots(
    cnots::Vector{Tuple{Int,Int}},
    graph;
    greedier=false,
    heuristic_dist=(n, y) -> 0
)::Vector{Path}
    @assert length(Set(Iterators.flatten(cnots))) == 2 * length(cnots) "cnots aren't disjoint"
    if greedier
        greedier_edp(
            graph,
            cnots;
            heuristic_dist=heuristic_dist
        )
    else
        greedy_edp(
            graph,
            cnots;
            heuristic_dist=heuristic_dist
        )
    end
end

"""Apply a cnot to a boundary qubit for a boundary operation
"""
function apply_boundary_ops(
    qubits::Vector{Int},
    boundary::Vector{Int},
    graph,
)::Vector{Path}
    if isempty(qubits)
        return Path[]
    end

    # Modify graph in-place
    add_vertices!(graph, 2)
    t = nv(graph)
    s = nv(graph) - 1

    for qubit in qubits
        add_edge!(graph, s, qubit)
    end
    for boundary_node in boundary
        add_edge!(graph, boundary_node, t)
    end

    maxflow, F = maximum_flow(graph, s, t)
    @assert maxflow > 0 "At least one boundary operation must be applied"
    # F cannot contain length 2 cycles because then the flow capacity will cancel out

    #undo graph modifications
    rem_vertex!(graph, t)
    rem_vertex!(graph, s)


    # Reconstruct paths from flow matrix F
    paths = Vector{Path}(undef, maxflow)
    Succ = F .== 1 # precompute boolean successor matrix
    for i in 1:maxflow
        path = [s]
        while path[end] != t
            cur = path[end]
            succ = findfirst(Succ[cur, :])
            Succ[cur, succ] = false
            append!(path, succ)
        end

        paths[i] = path[2:end-1]
    end

    return paths
end

"""Compute a schedule to apply the operations given in ops"

    Every element of the returned vector is a set of layers of operations. Such a set
    needs to be surrounded by appropriate Hadamard gates that add a bit of time cost.
    One set of operations can be executed in constant time. Therefore, this is sufficient
    information to upper bound the total time cost of running the operations in ops.

    `width` and `height` specify the number of logical data qubits in each dimension
"""
function apply_ops(
    ops::Vector{Op},
    dag::SimpleDiGraph,
    width::Int,
    height::Int;
    greedier=false,
    snake=false
)::Vector{Vector{Vector{MappedOp}}}
    @assert width >= 2
    @assert height >= 2
    schedule = Vector{Vector{Vector{MappedOp}}}()
    scaled_width = ancilla_grid_id(width, width)

    for op in ops
        for qubit in op.qubits
            @assert qubit ≤ width * height "Grid too small for op $op"
        end
    end

    # Define a distance heuristic for A* to use
    # that always underestimates the true distance.
    # Use Manhatten distance (L1-distance) as heuristic
    function heuristic_dist(from::Int, to::Int)::Int
        x0, y0 = grid_coord(from, scaled_width)
        x1, y1 = grid_coord(to, scaled_width)
        abs(x0 - x1) + abs(y0 - y1)
    end

    # Since the EDP schedule is unidirectional and we have Hadamards anyway,
    # we do not care about distinguishing CX, cxX, and cz
    binop_names = ("CX", "cxX", "cz")
    boundary_names = ("tz", "tx", "sz", "sx")

    # Translate qubit ids to snake allocation
    if snake
        ops = map(ops) do op
            qubits = snake_id.(op.qubits, width)
            Op(op.id, op.op, qubits)
        end
    end

    dag_layers = layers(dag)
    total_layers = length(dag_layers)
    layer_counter = 1
    for layer in dag_layers
        # @debug "Processing layer $layer_counter/$total_layers"
        # Each layer is schedule in a number of time steps
        # Surrounding one such layer are appropriate Hadamard gates
        layer_schedule = Vector{Vector{MappedOp}}()

        layer_ops = ops[layer]
        binops = filter(x -> x.op ∈ binop_names, layer_ops)
        bounops = filter(x -> x.op ∈ boundary_names, layer_ops)
        unops = filter(x -> x.op ∉ binop_names && x.op ∉ boundary_names, layer_ops)

        # Find an EDP schedule for the boundary and binary operations
        while !isempty(binops) || !isempty(bounops)
            sources = [binop.qubits[1] for binop in binops]
            append!(sources, [bounop.qubits[1] for bounop in bounops])
            destinations = [binop.qubits[2] for binop in binops]
            graph = ancilla_grid(sources, destinations, width, height)
            boundary = boundary_qubits(2 * (width + 2) - 1, 2 * (height + 2) - 1)

            # Apply boundary ops, which are all unary operations
            # Change qubit ids to grid with ancilla
            scaled_qubits = [ancilla_grid_id(bounop.qubits[1], width) for bounop in bounops]
            bounop_paths = apply_boundary_ops(scaled_qubits, boundary, graph)
            if !isempty(scaled_qubits)
                @assert !isempty(bounop_paths) "At least 1 boundary path must be found"
            end

            # Remove used edges from graph
            for path in bounop_paths
                for e in zip(path, path[2:end])
                    rem_edge!(graph, e[1], e[2])
                    rem_edge!(graph, e[2], e[1])
                end
            end

            # Find associated bounops and build schedule
            timestep = Vector{MappedOp}()
            bounop_idxs = indexin(first.(bounop_paths), scaled_qubits)
            selected_bounops = bounops[bounop_idxs]
            # Remove scheduled bounops
            deleteat!(bounops, sort(bounop_idxs))
            append!(timestep, MappedOp.(selected_bounops, bounop_paths))

            # Apply as many binops as possibly (greedily). Could be 0 depending on bounops
            scaled_terminals = TerminalPair[
                (
                    ancilla_grid_id(binop.qubits[1], width),
                    ancilla_grid_id(binop.qubits[2], width)
                ) for binop in binops]
            edp = apply_disjoint_cnots(
                scaled_terminals, graph;
                greedier=greedier,
                heuristic_dist=heuristic_dist)

            # Find associated binops and append to schedule
            terminal_pairs = [(path[1], path[end]) for path in edp]
            binop_idxs = indexin(terminal_pairs, scaled_terminals)
            selected_binops = binops[binop_idxs]
            # Remove scheduled binops
            deleteat!(binops, sort(binop_idxs))
            append!(timestep, MappedOp.(selected_binops, edp))

            push!(layer_schedule, timestep)
        end

        # Schedule all unary operations immediately
        unop_schedule = [MappedOp(unop, [ancilla_grid_id(unop.qubits[1], width)]) for unop in unops]
        if length(layer_schedule) == 0
            push!(layer_schedule, unop_schedule)
        else
            append!(layer_schedule[1], unop_schedule)
        end

        push!(schedule, layer_schedule)
        layer_counter += 1
    end

    return schedule
end

"""Compute a CNOT schedule
"""
function apply_cnot_circuit(
    ops::Vector{Op},
    dag::SimpleDiGraph,
    width::Int,
    height::Int;
    greedier=false,
    snake=false
)::Vector{Vector{MappedOp}}
    @assert width >= 2
    @assert height >= 2
    schedule = Vector{Vector{Vector{MappedOp}}}()
    scaled_width = ancilla_grid_id(width, width)

    for op in ops
        for qubit in op.qubits
            @assert qubit ≤ width * height "Grid too small for op $op"
        end
    end

    # Define a distance heuristic for A* to use
    # that always underestimates the true distance.
    # Use Manhatten distance (L1-distance) as heuristic
    function heuristic_dist(from::Int, to::Int)::Int
        x0, y0 = grid_coord(from, scaled_width)
        x1, y1 = grid_coord(to, scaled_width)
        abs(x0 - x1) + abs(y0 - y1)
    end

    # Translate qubit ids to snake allocation
    if snake
        ops = map(ops) do op
            qubits = snake_id.(op.qubits, width)
            Op(op.id, op.op, qubits)
        end
    end

    visitor = DAGTopologicalVisitor(dag)

    schedule = Vector{Vector{MappedOp}}()
    while !isempty(visitor.frontier)
        frontier = collect(visitor.frontier)
        current_ops = ops[frontier] # Get the ops associated with vertices in DAG

        sources = [binop.qubits[1] for binop in current_ops]
        destinations = [binop.qubits[2] for binop in current_ops]
        graph = ancilla_grid(sources, destinations, width, height)

        scaled_terminals = TerminalPair[
            (
                ancilla_grid_id(op.qubits[1], width),
                ancilla_grid_id(op.qubits[2], width)
            ) for op in current_ops]
        edp = apply_disjoint_cnots(
            scaled_terminals, graph;
            greedier=greedier,
            heuristic_dist=heuristic_dist)

        # Find associated binops and append to schedule
        terminal_pairs = [(path[1], path[end]) for path in edp]
        binop_idxs = indexin(terminal_pairs, scaled_terminals)
        selected_binops = current_ops[binop_idxs]
        selected_frontier = frontier[binop_idxs]
        # visit selected binops
        for v in selected_frontier
            visit!(visitor, v)
        end

        push!(schedule, MappedOp.(selected_binops, edp))
    end

    return schedule
end

function schedule_counts(
    schedule::Vector{Vector{Vector{MappedOp}}}
)
    hadamards = 0
    vdps = 0
    edps = 0
    nr_pqubits = maximum(qubit for layer in schedule for step in layer for op in step for qubit in op.qubits)

    prev_hadamard = zeros(Bool, nr_pqubits)

    for layer in schedule
        # Check if any gate in layer needs Hadamard
        needed_hadamard = false
        hadamard_qubits = Set(q for step in layer
                              for op in step if op.op in HADAMARD_OPS
                              for q in op.qubits)
        for q in eachindex(prev_hadamard)
            if !xor(prev_hadamard[q], q in hadamard_qubits)
                # Hadamards cancel or no Hadamard
                continue
            end

            needed_hadamard = true
            prev_hadamard[q] = !prev_hadamard[q]
        end
        if needed_hadamard
            hadamards += 1
        end
        # Hadamard after final layer can be corrected by rotating pauli measurements

        # Iterate through EDP steps
        for step in layer
            # # Check if any EDP/VDP was necessary (maybe only unops)
            # if all(length(op[2]) == 1 for op in step)
            #     # Unop layer can be executed at no time cost
            #     # because it can be merged with subsequent layers
            #     continue
            # end

            # Check if EDP step was VDP
            is_edp = false
            visited = Set{Int}()
            for op in step, node in op.path
                if node ∈ visited
                    is_edp = true
                    break
                end
                push!(visited, node)
            end
            if is_edp
                edps += 1
            else
                vdps += 1
            end
        end
    end
    return (hadamards=hadamards, vdps=vdps, edps=edps)
end