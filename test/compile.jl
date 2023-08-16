# Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
# Microsoft Software License Terms for Microsoft Quantum Development Kit Libraries
# and Samples. See LICENSE in the project root for license information.

@testset "ancilla graph id" begin
    @test TeleportRouter.ancilla_grid_id(1, 2) == 17
    @test TeleportRouter.ancilla_grid_id(2, 2) == 19
    @test TeleportRouter.ancilla_grid_id(3, 2) == 31
    @test TeleportRouter.ancilla_grid_id(4, 2) == 33
end

@testset "ancilla graph" begin
    small_grid = TeleportRouter.ancilla_grid([1], [4], 2, 2)
    # Should be 7x7 (we add boundary rows/columns around and ancilla).
    @test length(vertices(small_grid)) == 49
    terminal_vertex = TeleportRouter.ancilla_grid_id(1, 2)
    width = 7
    @test !has_edge(small_grid, terminal_vertex, terminal_vertex + 1)
    @test has_edge(small_grid, terminal_vertex, terminal_vertex + width) # only control vertically
    @test isempty(inneighbors(small_grid, terminal_vertex))
    @test length(outneighbors(small_grid, terminal_vertex)) == 2
    terminal_vertex = TeleportRouter.ancilla_grid_id(4, 2)
    @test !has_edge(small_grid, terminal_vertex - width, terminal_vertex)
    @test has_edge(small_grid, terminal_vertex - 1, terminal_vertex) # only target horizontally
    @test isempty(outneighbors(small_grid, terminal_vertex))
    @test length(inneighbors(small_grid, terminal_vertex)) == 2
    # Central ancilla qubit
    @test Set(inneighbors(small_grid, 25)) ∩ Set(outneighbors(small_grid, 25)) == Set([24, 26, 25 + width, 25 - width])
    #Nonparticipating data qubits
    terminal_vertex = TeleportRouter.ancilla_grid_id(2, 2)
    @test isempty(all_neighbors(small_grid, terminal_vertex))
    terminal_vertex = TeleportRouter.ancilla_grid_id(3, 2)
    @test isempty(all_neighbors(small_grid, terminal_vertex))

    # Terminal out of bounds
    @test_throws ArgumentError TeleportRouter.ancilla_grid([1], [100], 2, 2)

    # Boundary qubits
    @test Set(TeleportRouter.boundary_qubits(7, 7)) ==
          Set([1, 3, 5, 7, 15, 21, 29, 35, 43, 45, 47, 49])
end


"""Remove edges to the boundary"""
function rem_boundary!(graph, width, height)
    logical_topleft = TeleportRouter.ancilla_grid_id(1, width)
    scaled_width = (TeleportRouter.ancilla_grid_id(1 + width, width) - logical_topleft) ÷ 2
    boundary_topleft = logical_topleft - scaled_width - 1
    logical_botright = TeleportRouter.ancilla_grid_id(width * height, width)
    boundary_botright = logical_botright + scaled_width + 1

    function rem_edges!(id)
        for neighbor in collect(all_neighbors(graph, id))
            rem_edge!(graph, id, neighbor)
            rem_edge!(graph, neighbor, id)
        end
    end
    rem_edges!.(boundary_topleft:boundary_topleft+scaled_width-3)
    rem_edges!.(boundary_botright:-1:boundary_botright-scaled_width+3)
    rem_edges!.(boundary_topleft:scaled_width:boundary_botright)
    rem_edges!.(boundary_botright:-scaled_width:boundary_topleft)
end

function path_exists(graph, path)::Bool
    all(has_edge(graph, e) for e in zip(path, path[2:end]))
end

@testset "disjoint cnots" begin
    width = 3
    height = 2
    sources = [1, 2, 4]
    destinations = [3, 5, 6]
    graph = TeleportRouter.ancilla_grid(sources, destinations, width, height)
    scaled_terminals = [(TeleportRouter.ancilla_grid_id(terminal[1], width),
        TeleportRouter.ancilla_grid_id(terminal[2], width))
                        for terminal in zip(sources, destinations)]
    rem_boundary!(graph, width, height)


    edp = apply_disjoint_cnots(scaled_terminals, graph)
    @test length(edp) == 1
    @test (edp[1][1], edp[1][end]) ∈ scaled_terminals

    width = 2
    height = 2
    sources = [1, 4]
    destinations = [2, 3]
    graph = TeleportRouter.ancilla_grid(sources, destinations, width, height)
    scaled_terminals = [(TeleportRouter.ancilla_grid_id(terminal[1], width),
        TeleportRouter.ancilla_grid_id(terminal[2], width))
                        for terminal in zip(sources, destinations)]
    rem_boundary!(graph, 2, 2)
    edp = apply_disjoint_cnots(scaled_terminals, graph)
    @test length(edp) == 2
    @test (edp[1][1], edp[1][end]) ∈ scaled_terminals
    @test (edp[2][1], edp[2][end]) ∈ scaled_terminals
end

@testset "apply boundary ops" begin
    graph = DiGraph(8)
    edges = [
        (1, 2), (1, 3), (1, 4), (2, 3), (2, 5),
        (2, 6), (3, 4), (3, 6), (4, 7), (5, 6),
        (5, 8), (6, 7), (6, 8), (7, 3), (7, 8)
    ]
    for e in edges
        add_edge!(graph, e[1], e[2])
    end

    graph_copy = copy(graph)
    paths = TeleportRouter.apply_boundary_ops([1], [8], graph)
    @test graph_copy == graph # not mutated
    @test length(paths) == 1
    @test paths[1][1] == 1
    @test paths[1][end] == 8
    @test path_exists(graph, paths[1])

    paths = TeleportRouter.apply_boundary_ops([1, 2], [8], graph)
    @test length(paths) == 1
    @test path_exists(graph, paths[1])

    paths = TeleportRouter.apply_boundary_ops([1, 2], [8, 7], graph)
    @test length(paths) == 2
    @test Set([paths[1][1], paths[2][1]]) == Set([1, 2])
    @test Set([paths[1][end], paths[2][end]]) == Set([8, 7])
    for path in paths
        @test path_exists(graph, path)
    end
end

@testset "compile operations" begin
    @testset "2×2 grid" begin
        op = Op(1, "CX", [1, 2])
        ops = [op]

        dependent_on = fill(Vector{Int}(), length(ops))
        dag = TeleportRouter.dag_circuit(ops, dependent_on)

        width = 2
        height = 2

        schedule = apply_ops(ops, dag, width, height)
        @test length(schedule) == 1
        @test Op(schedule[1][1][1]) == op

        op2 = Op(2, "cxX", [3, 4])
        push!(ops, op2)
        dag = TeleportRouter.dag_circuit(ops, dependent_on)
        schedule = apply_ops(ops, dag, width, height)
        @test length(schedule) == 1 # All ops belong to the same layer in the dag
        @test length(schedule[1]) == 1 # The ops can be performed in one EDP using boundary
        @test Set(ops) == Set(Op.([schedule[1][1][1], schedule[1][1][2]]))

        push!(dependent_on, [1]) # op2 depends on op1
        dag = TeleportRouter.dag_circuit(ops, dependent_on)
        schedule = apply_ops(ops, dag, width, height)
        @test length(schedule) == 2 # The two ops belong to different layers in the dag
        @test Op(schedule[1][1][1]) == op
        @test Op(schedule[2][1][1]) == op2

        op3 = Op(3, "tx", [1])
        push!(ops, op3)
        dag = TeleportRouter.dag_circuit(ops, dependent_on)
        schedule = apply_ops(ops, dag, width, height)
        @test length(schedule) == 2 # The two ops belong to different layers in the dag
        # Check if T gate was routed to boundary
        boundary = Set(TeleportRouter.boundary_qubits(width, height))
        # Find T op
        timestep2 = schedule[2][1]
        idx = findfirst(x -> x.op == "tx", timestep2)
        txOp = timestep2[idx]
        @test last(txOp.path) ∈ boundary
        # Make a dummy copy of the graph and check if the paths are valid
        graph = TeleportRouter.ancilla_grid([op3.qubits[1], op2.qubits[1]], [op2.qubits[2]], width, height)
        for scheduled_op in timestep2
            @test path_exists(graph, scheduled_op.path)
        end
    end
    @testset "3×3 grid" begin
        width = 3
        height = 3
        ops = [Op(1, "CX", [1, 7]), Op(2, "cxX", [2, 9])]
        dependent_on = fill(Vector{Int}(), length(ops))
        dag = TeleportRouter.dag_circuit(ops, dependent_on)
        schedule = apply_ops(ops, dag, width, height)
        @test length(schedule) == 1 # The ops belong to the same layer
        @test length(schedule[1]) == 1 # The ops can be executed in the same EDP
        scheduled_ops = Set(Op.(schedule[1][1]))
        @test scheduled_ops == Set(ops) # All ops were scheduled

        push!(ops, Op(3, "X", [4]))
        dag = TeleportRouter.dag_circuit(ops, dependent_on)
        schedule = apply_ops(ops, dag, width, height)
        @test length(schedule) == 1 # The ops belong to the same layer
        @test length(schedule[1]) == 1 # The ops can be executed in the same EDP
        scheduled_ops = Set(Op.(schedule[1][1]))
        @test scheduled_ops == Set(ops) # All ops were scheduled

        ops_l1 = copy(ops)
        op3 = Op(4, "CX", [4, 7])
        push!(ops, op3)
        dag = TeleportRouter.dag_circuit(ops, dependent_on)
        schedule = apply_ops(ops, dag, width, height)
        @test length(schedule) == 2 # The ops belong to two layers
        @test length(schedule[1]) == 1 # The ops in layer 1 can be executed in the same EDP
        scheduled_ops = Set(Op.(schedule[1][1]))
        @test scheduled_ops == Set(ops_l1) # All layer 1 ops were scheduled
        @test Op(schedule[2][1][1]) == op3 # And the layer 2 op was scheduled
    end
end

@testset "apply cnot circuit" begin
    @testset "2×2 grid" begin
        op = Op(1, "CX", [1, 2])
        ops = [op]

        dependent_on = fill(Vector{Int}(), length(ops))
        dag = TeleportRouter.dag_circuit(ops, dependent_on)

        width = 2
        height = 2

        schedule = apply_cnot_circuit(ops, dag, width, height)
        @test length(schedule) == 1
        @test Op(schedule[1][1]) == op

        op2 = Op(2, "cxX", [3, 4])
        push!(ops, op2)
        dag = TeleportRouter.dag_circuit(ops, dependent_on)
        schedule = apply_cnot_circuit(ops, dag, width, height)
        @test length(schedule) == 1 # All ops belong to the same layer in the dag
        @test length(schedule[1]) == 2 # The ops can be performed in one EDP using boundary
        @test Set(ops) == Set(Op.([schedule[1][1], schedule[1][2]]))

        push!(dependent_on, [1]) # op2 depends on op1
        dag = TeleportRouter.dag_circuit(ops, dependent_on)
        schedule = apply_cnot_circuit(ops, dag, width, height)
        @test length(schedule) == 2 # The two ops belong to different layers in the dag
        @test Op(schedule[1][1]) == op
        @test Op(schedule[2][1]) == op2
    end

    @testset "3×3 grid" begin
        width = 3
        height = 3
        ops = [Op(1, "CX", [1, 7]), Op(2, "cxX", [2, 9])]
        dependent_on = fill(Vector{Int}(), length(ops))
        dag = TeleportRouter.dag_circuit(ops, dependent_on)
        schedule = apply_cnot_circuit(ops, dag, width, height)
        @test length(schedule) == 1 # The ops belong to the same layer
        @test length(schedule[1]) == 2 # The ops can be executed in the same EDP
        scheduled_ops = Set(Op.(schedule[1]))
        @test scheduled_ops == Set(ops) # All ops were scheduled

        ops_l1 = copy(ops)
        op3 = Op(4, "CX", [4, 7])
        push!(ops, op3)
        dag = TeleportRouter.dag_circuit(ops, dependent_on)
        schedule = apply_cnot_circuit(ops, dag, width, height)
        @test length(schedule) == 2 # The ops belong to two layers
        @test length(schedule[1]) == 2 # The ops in layer 1 can be executed in the same EDP
        scheduled_ops = Set(Op.(schedule[1]))
        @test scheduled_ops == Set(ops_l1) # All layer 1 ops were scheduled
        @test Op(schedule[2][1]) == op3 # And the layer 2 op was scheduled
    end
end
