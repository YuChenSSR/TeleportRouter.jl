# Copyright (c) Microsoft Corporation. All rights reserved. Licensed under the
# Microsoft Software License Terms for Microsoft Quantum Development Kit Libraries
# and Samples. See LICENSE in the project root for license information.

using TeleportRouter
using JSON

function parse_circuit(str :: String)
    line_re = r"^([A-Za-z]+) q\[([0-9]+)\](?:,q\[([0-9]+)\])?;$"
    lines = split(str, "\n")[2:end-1] # Remove blank end line and qreg first line
    id = 1
    matched_lines = match.(line_re, lines)
    ops = map(matched_lines) do matched_line
        qubits = [parse(Int, matched_line[2])]
        if matched_line[3] !== nothing
            push!(qubits, parse(Int, matched_line[3]))
        end
        qubits .+= 1 #Use 1-indexing
        op = Op(id, matched_line[1], qubits)
        id += 1
        op
    end
    return ops
end
