""" Model as in exemplar runs:
    - Purification 
    - TEN and TEE 
    Used in Sanity Checks
"""

function run_measurement_only_dynamics_PURE!(state::MixedDestabilizer,
    get_random_pauli::Function,
    system,
    simulation)
"""this runs the mesaurement only dynamics. 
state            : state to be evolved. It is taken by reference (hence the ! in the name by convention).
get_random_pauli : a function that is called each time step and should generate the paulis to be measured.
system               : Defines the geometry of the system
simulation       : Defines the paramters of the simulation
"""

nt = simulation.nt
measure_ts = simulation.measure_ts

# get system size n (number of qubits)
n = system.nbits
# this will store the total entropy (con Neumann) as a function of time
entropies = Vector{Float64}([entanglement_entropy(state, 1:n, Val(:clip))])

# now perform nt measurement sweeps
for t in 1:nt
for tt in 1:n
# projectrand randomizes also the measurement outcomes
projectrand!(state,             # state to measure
get_random_pauli() # the model is fully defined by this function
)
end
if t > measure_ts[length(entropies)]
# calculate the full van-neumann entropy of the system
s = entanglement_entropy(state,     # state to compute entropy for
  1:n,       # subsystem as list of indices. Specifying full cell_sizetem gives the van-neumann entropy
  Val(:clip) # algorithm to use (see documentation)
  )
push!(entropies, s)
end
end

# print(stabilizerview(state))
return entropies 
end


function run_measurement_only_dynamics_TEN!(state::MixedDestabilizer,
    get_random_pauli::Function,
    system,
    simulation)
"""this runs the mesaurement only dynamics. 
state            : state to be evolved. It is taken by reference (hence the ! in the name by convention).
get_random_pauli : a function that is called each time step and should generate the paulis to be measured.
system               : Defines the geometry of the system
simulation       : Defines the paramters of the simulation
"""

nt = simulation.nt
measure_ts = simulation.measure_ts

# get system size n (number of qubits)
n = system.nbits
# this will store the total entropy (con Neumann) as a function of time
TENs = Vector{Float64}([entanglement_negativity!(state, system)])

# now perform nt measurement sweeps
for t in 1:nt
for tt in 1:n
# projectrand randomizes also the measurement outcomes
projectrand!(state,             # state to measure
get_random_pauli() # the model is fully defined by this function
)
end
if t > measure_ts[length(TENs)]
# calculate the full van-neumann entropy of the system
ten = entanglement_negativity!(state,     # state to compute TEN for
  system,       # Details of the geometry
  )
push!(TENs, ten)
end
end

# print(stabilizerview(state))
return TENs 
end


function run_measurement_only_dynamics_TEE!(state::MixedDestabilizer,
    get_random_pauli::Function,
    system,
    simulation)
"""this runs the mesaurement only dynamics. 
state            : state to be evolved. It is taken by reference (hence the ! in the name by convention).
get_random_pauli : a function that is called each time step and should generate the paulis to be measured.
system               : Defines the geometry of the system
simulation       : Defines the paramters of the simulation
"""

nt = simulation.nt
measure_ts = simulation.measure_ts

# get system size n (number of qubits)
n = system.nbits
# this will store the total entropy (con Neumann) as a function of time
TEEs = Vector{Float64}([entanglement_entropy_topo!(state, system)])

# now perform nt measurement sweeps
for t in 1:nt
for tt in 1:n
# projectrand randomizes also the measurement outcomes
projectrand!(state,             # state to measure
get_random_pauli() # the model is fully defined by this function
)
end
if t > measure_ts[length(TEEs)]
# calculate the full van-neumann entropy of the system
tee = entanglement_entropy_topo!(state,     # state to compute TEE for
  system,       # Details of the geometry of the system.
  )
push!(TEEs, tee)
end
end

# print(stabilizerview(state))
return TEEs 
end
