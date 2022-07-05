# Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
# Microsoft Software License Terms for Microsoft Quantum Development Kit Libraries
# and Samples. See LICENSE in the project root for license information.

using Test

using TeleportRouter
using LightGraphs
import Random

@testset "Does it work" begin
    @test 1 + 1 == 2
end

include("edge_disjoint_paths.jl")
include("compile.jl")
include("parse_circuit.jl")
include("topological_sort.jl")
