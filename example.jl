import TeleportRouter


# Load in a circuit from JSON and parse it to Ops and a DAG
ops, idx = TeleportRouter.parse_circuit(open("circuits/example.circuit.json", "r"))
dag = TeleportRouter.dag_circuit(ops, idx)

println(dag)
println(ops)

# Run the compiler on the given ops on a 3 by 2 grid.
schedule = TeleportRouter.apply_ops(ops, dag, 3, 2)

println(schedule)
