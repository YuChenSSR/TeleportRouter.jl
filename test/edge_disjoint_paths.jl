# Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
# Microsoft Software License Terms for Microsoft Quantum Development Kit Libraries
# and Samples. See LICENSE in the project root for license information.

using Test

using TeleportRouter
using LightGraphs
import Random

@testset "Edge disjoint path tests" begin
    @testset "greedy" begin
        graph = grid([5, 1])
        disjoint_paths = TeleportRouter.greedy_edp(
            grid([5, 1]),
            Set([(1, 2), (3, 5)]);
        )
        @test length(disjoint_paths) == 2
        # don't care about order
        @test Set([[1, 2], [3, 4, 5]]) == Set(disjoint_paths)

        overlapping_paths = TeleportRouter.greedy_edp(
            grid([5, 1]),
            Set([(2, 3), (1, 5)]),
        )
        @test overlapping_paths == [[2, 3]]

        alternate_paths = TeleportRouter.greedy_edp(
            grid([5, 2]),
            Set([(2, 4), (1, 5)]),
        )
        @test length(alternate_paths) == 2
        @test alternate_paths[1] == [2, 3, 4]
        @test alternate_paths[2][1] == 1
        @test alternate_paths[2][end] == 5
        @test 3 ∉ alternate_paths[2]

        # Check if an edge-disjoint but not vertex-disjoint path is allowed.
        vertex_intersect = TeleportRouter.greedy_edp(
            grid([5, 3]),
            Set([(2, 12), (6, 10)]),
        )
        @test length(vertex_intersect) == 2
        # Cuts the graph in two with edges.
        @test Set(vertex_intersect) == Set([
            [2, 7, 12],
            6:10]) # EDP still allows for this path
    end

    @testset "greedier" begin
        @show disjoint_paths = TeleportRouter.greedier_edp(
            grid([5, 1]),
            Set([(1, 2), (3, 5)]),
        )
        @test length(disjoint_paths) == 2
        # don't care about order
        @test Set([[1, 2], [3, 4, 5]]) == Set(disjoint_paths)

        overlapping_paths = TeleportRouter.greedier_edp(
            grid([5, 1]),
            Set([(2, 3), (1, 5)]),
        )
        @test length(overlapping_paths) == 1 # Order is random so either is correct

        # Even if (1,5) is picked first, there is still a path from 2 to 4.
        # Fix RNG for test so that (1,5) is picked
        alternate_paths = TeleportRouter.greedier_edp(
            grid([5, 2]),
            Set([(2, 4), (1, 5)]);
            rng=Random.MersenneTwister(11)
        )
        @test length(alternate_paths) == 2
        @test alternate_paths[1] == [1, 2, 3, 4, 5]
        @test alternate_paths[2] == [2, 7, 8, 9, 4]
    end
end
