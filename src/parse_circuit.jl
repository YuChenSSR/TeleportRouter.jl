# Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
# Microsoft Software License Terms for Microsoft Quantum Development Kit Libraries
# and Samples. See LICENSE in the project root for license information.

abstract type AbstractOp end
struct Op <: AbstractOp
    id :: Int
    op :: String
    qubits :: Vector{Int}
end

struct MappedOp <: AbstractOp
    id :: Int
    op :: String
    qubits :: Vector{Int}
    path :: Vector{Int}
end
MappedOp(op :: Op, path :: Vector{Int}) = MappedOp(op.id, op.op, op.qubits, path)
Op(mappedOp :: MappedOp) = Op(mappedOp.id, mappedOp.op, mappedOp.qubits)

function parse_circuit(io:: IO) :: Tuple{Vector{Op}, Vector{Vector{Int}}}
    parsed = JSON.parse(io)

    # Assume ids are sequential
    ops = Vector{Op}()
    dependent_on = Vector{Vector{Int}}()
    for parsed_op in parsed
        op = Op(parsed_op["id"], parsed_op["op"], parsed_op["qubits"])
        if op.id < 1 || any(op.qubits .< 1)
            throw(ArgumentError("Index or qubit is not strictly positive. I expect 1-indexed values.\n$op"))
        end
        push!(ops, op)
        push!(dependent_on, Vector{Int}(get(parsed_op, "depends-on", [])))
    end

    return (ops, dependent_on)
end

"""
    Parse a JSON file from io to vector of MappedOps
"""
function parse_mapped_circuit(io:: IO) :: Vector{Vector{Vector{MappedOp}}}
    parsed = JSON.parse(io)

    # Assume ids are sequential
    ops = Vector{MappedOp}()
    dependent_on = Vector{Vector{Int}}()
    ops = map(parsed) do layer
        map(layer) do step
            map(step) do parsed_op
                op = MappedOp(parsed_op["id"], parsed_op["op"], parsed_op["qubits"], parsed_op["path"])
                if op.id < 1 || any(op.qubits .< 1) || any(op.path .< 1)
                    throw(ArgumentError("Index, qubits, or path are not strictly positive. I expect 1-indexed values.\n$op"))
                end
                op
            end
        end
    end

    return ops 
end

function dag_circuit(ops :: Vector{<:AbstractOp}, dependent_on :: Vector{Vector{Int}}) :: LightGraphs.SimpleDiGraph
    dag = LightGraphs.SimpleDiGraph(length(ops))
    predecessor = Dict{Int, Int}()

    for (i, op) in enumerate(ops)
        for qubit in op.qubits
            if haskey(predecessor, qubit)
                LightGraphs.add_edge!(dag, predecessor[qubit], i)
            end
            predecessor[qubit] = i
        end
    end

    for (i, dependents) in enumerate(dependent_on)
        for dependent in dependents
            LightGraphs.add_edge!(dag, dependent, i)
        end
    end

    return dag
end
