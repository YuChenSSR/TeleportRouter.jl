# Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
# Microsoft Software License Terms for Microsoft Quantum Development Kit Libraries
# and Samples. See LICENSE in the project root for license information.

@testset "Parse circuit" begin
    (ops, dependent_on) = open(joinpath(".","circuits","example.circuit.json"), "r") do io
        parse_circuit(io)
    end
    @test length(ops) == 5
end