# Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
# Microsoft Software License Terms for Microsoft Quantum Development Kit Libraries
# and Samples. See LICENSE in the project root for license information.

const GATE_COSTS = Dict(
    "swap" => 2
    , "cx" => 2
    , "CX" => 2
    , "cxX" => 2*3 + 2 # 2*2 Hadamard + cnot
    , "cz" => 2*3 + 2
    , "tz" => 0 # Assume performed on boundary
    , "tx" => 0
    , "sz" => 0
    , "sx" => 0
    , "mx" => 0
    , "mz" => 0
)

const HADAMARD_OPS = ["ccx", "cxX", "cz", "mx", "tx"]

function dag_cost(dag, ops)
    cost = zeros(Int, length(ops))
    for node in LightGraphs.topological_sort_by_dfs(dag)
        # Get maximum depth of preceding nodes
        node_depth = maximum(cost[pre] for pre in inneighbors(dag, node); init=0)
        # Add the cost of the current gate
        node_depth += GATE_COSTS[ops[node].op]
        cost[node] = node_depth
    end
    maximum(cost)
end