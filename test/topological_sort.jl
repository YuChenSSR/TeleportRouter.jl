# Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
# Microsoft Software License Terms for Microsoft Quantum Development Kit Libraries
# and Samples. See LICENSE in the project root for license information.

@testset "Topological Sort" begin
    @testset "topological_sort" begin
        dag = SimpleDiGraph(5)
        LightGraphs.add_edge!(dag, 1,4)
        LightGraphs.add_edge!(dag, 2,4)
        LightGraphs.add_edge!(dag, 2,3)
        LightGraphs.add_edge!(dag, 3,5)
        top_sorted = TeleportRouter.topological_sort(dag)
        @test Set(top_sorted) == Set(1:5)
        @test findfirst(x -> x == 1, top_sorted) < findfirst(x -> x == 4, top_sorted)
        @test findfirst(x -> x == 2, top_sorted) < findfirst(x -> x == 4, top_sorted)
        @test findfirst(x -> x == 2, top_sorted) < findfirst(x -> x == 3, top_sorted)
        @test findfirst(x -> x == 3, top_sorted) < findfirst(x -> x == 5, top_sorted)

        # Add cycle
        LightGraphs.add_edge!(dag, 4, 1)
        @test_throws ArgumentError TeleportRouter.topological_sort(dag)
    end
    @testset "layers" begin
        dag = SimpleDiGraph(5)
        LightGraphs.add_edge!(dag, 1,4)
        LightGraphs.add_edge!(dag, 2,4)
        LightGraphs.add_edge!(dag, 2,3)
        LightGraphs.add_edge!(dag, 3,5)
        layers = TeleportRouter.layers(dag)
        @test Set(layers[1]) == Set([1,2])
        @test Set(layers[2]) == Set([3,4])
        @test layers[3] == [5]

        LightGraphs.add_edge!(dag, 1,5)
        layers = TeleportRouter.layers(dag)
        @test Set(layers[1]) == Set([1,2])
        @test Set(layers[2]) == Set([3,4])
        @test layers[3] == [5]

        LightGraphs.add_edge!(dag, 2,5)
        layers = TeleportRouter.layers(dag)
        @test Set(layers[1]) == Set([1,2])
        @test Set(layers[2]) == Set([3,4])
        @test layers[3] == [5]

        LightGraphs.rem_edge!(dag, 3,5)
        layers = TeleportRouter.layers(dag)
        @test Set(layers[1]) == Set([1,2])
        @test Set(layers[2]) == Set([3,4,5])

        # Add cycle
        LightGraphs.add_edge!(dag, 4, 1)
        @test_throws ArgumentError TeleportRouter.layers(dag)

        # Get layers of a random graph
        nr_nodes = 10000
        random_dag = random_orientation_dag(SimpleGraph(nr_nodes, div(nr_nodes*(nr_nodes-1), 20)))
        @test sum(length, TeleportRouter.layers(random_dag)) == nr_nodes

        # Add cycle
        cycle_nodes = Random.randsubseq(vertices(random_dag), 0.01)
        rotated_cycle = [cycle_nodes[2:end]; cycle_nodes[1]]
        for (v1, v2) in zip(cycle_nodes, rotated_cycle)
            add_edge!(random_dag, v1, v2)
        end
        @test_throws ArgumentError TeleportRouter.layers(random_dag)
    end
end