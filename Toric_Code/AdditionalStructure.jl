"""
Some auxiliary types used to clean up the code and remove the reliance on Global Constants!

Additional Types:
"""

#
#
#############################################################################################################################

# Some Useful Packages

using Base.Filesystem  # for creaating directories etc
using Random, Distributions  # random numbers
using HDF5  # hdf5 files
using QuantumClifford  # this is the stabilizer simulation package
using Plots # for plotting
using Formatting # string formatting
using LinearAlgebra # some useful matrices etc.

#
#
#############################################################################################################################

# New Types

struct VertexSquareLattice
    # Dimensions of lattice on the torus
    L::Int
    nbits::Int
    cell_num::Int
    # Projectors and Ranges used in TEN and TEE measurments, assuming ususal 
    subsys::Int
    cell_size::Int
    P_A
    P_B
    P_C
    S
    r_A
    r_B
    r_C
    r_AB
    r_BC
    r_ABC
end

struct EdgeSquareLattice
    # Dimensions of lattice on the torus
    L::Int
    nbits::Int
    cell_num::Int
    # Projectors and Ranges used in TEN and TEE measurments, assuming ususal 
    subsys::Int
    cell_size::Int
    P_A
    P_B
    P_C
    S
    r_A
    r_B
    r_C
    r_AB
    r_BC
    r_ABC
end

struct SimulationTime
    nt
    measure_ts
end

function Init_VertexSquareLattice(subsys::Int)
    
    # Dimensions
    L = 4*subsys
    k = 1 # Number of Qubits per unit cell
    cell_num = L*L
    nbits = k*cell_num
    
    # For TEN: define the projection onto the ususal multipartitions of the torus
    cell_size = k*L*subsys
    P_A = Diagonal(repeat(vcat(repeat([true], cell_size), repeat([false], 3*cell_size)), 2))
    P_B = Diagonal(repeat(vcat(repeat([false], cell_size), repeat([true], cell_size), repeat([false], 2*cell_size)), 2))
    P_C = Diagonal(repeat(vcat(repeat([false], 2*cell_size), repeat([true], cell_size), repeat([false], cell_size)), 2))
    S = kron([false true; true false], Diagonal(repeat([true], k*L*L)));

    # For TEE: refined qubit indices for diferent continous partitions
    r_A = 0*cell_size+1:1*cell_size
    r_B = 1*cell_size+1:2*cell_size
    r_C = 2*cell_size+1:3*cell_size
    r_AB = 0*cell_size+1:2*cell_size
    r_BC = 1*cell_size+1:3*cell_size
    r_ABC = 0*cell_size+1:3*cell_size

    return VertexSquareLattice(L, nbits, cell_num, subsys, cell_size, P_A, P_B, P_C, S, r_A, r_B, r_C, r_AB, r_BC, r_ABC)
end

function Init_EdgeSquareLattice(subsys::Int)
    
    # Dimensions
    L = 4*subsys
    k = 2 # Number of Qubits per unit cell
    cell_num = L*L
    nbits = k*cell_num
    
    # For TEN: define the projection onto the ususal multipartitions of the torus
    cell_size = k*L*subsys # Misleading NAME, be careful!
    P_A = Diagonal(repeat(vcat(repeat([true], cell_size), repeat([false], 3*cell_size)), 2))
    P_B = Diagonal(repeat(vcat(repeat([false], cell_size), repeat([true], cell_size), repeat([false], 2*cell_size)), 2))
    P_C = Diagonal(repeat(vcat(repeat([false], 2*cell_size), repeat([true], cell_size), repeat([false], cell_size)), 2))
    S = kron([false true; true false], Diagonal(repeat([true], k*L*L)));

    # For TEE: refined qubit indices for diferent continous partitions
    r_A = 0*cell_size+1:1*cell_size
    r_B = 1*cell_size+1:2*cell_size
    r_C = 2*cell_size+1:3*cell_size
    r_AB = 0*cell_size+1:2*cell_size
    r_BC = 1*cell_size+1:3*cell_size
    r_ABC = 0*cell_size+1:3*cell_size

    return EdgeSquareLattice(L, nbits, cell_num, subsys, cell_size, P_A, P_B, P_C, S, r_A, r_B, r_C, r_AB, r_BC, r_ABC)
end

"""
Some auxiliary functions used at various points. Augments the QunatumClifford.jl package and does not go beyond qubits!

Additional Functions:
"""

#
#
#############################################################################################################################

# define enum just for just for better readibility of code (Do not need PauliY in this case!)

@enum PauliType begin
    PauliX = 1
    PauliY = 2   
    PauliZ = 3
end

@enum StabType begin
    Star = 1
    Plaquette = 2
    e_hor = 3
    e_ver = 4
    m_hor = 5
    m_ver = 6
    f_hor = 7
    f_ver = 8
end

#
#
#############################################################################################################################

# Bit string opreations

function bit_string_i(i::Integer, n::Integer)
    """returns a bitstring (array of bool) of size n, that is false anywhere but at position i.
       This is convenient since paulis in the QuantumClifford package are defined by two bitstrings"""
    str = zeros(Bool, n);
    str[i] = true;
    return str
end

function bit_string_ij(i::Integer, j::Integer, n::Integer)
    """returns a bitstring (array of bool) of size n, that is false anywhere but at position i and j.
       This is convenient since paulis in the QuantumClifford package are defined by two bitstrings"""
    @assert(i != j);
    str = zeros(Bool, n);
    str[i] = true;
    str[j] = true;
    return str
end

function bit_string_ijkl(i::Integer, j::Integer, k::Integer, l::Integer, n::Integer)
    """returns a bitstring (array of bool) of size n, that is false anywhere but at position i, j, k and l.
       This is convenient since paulis in the QuantumClifford package are defined by two bitstrings"""
    @assert(i != j);
    @assert(i != k);
    @assert(i != l);
    @assert(j != k);
    @assert(j != l);
    @assert(k != l);
    str = zeros(Bool, n);
    str[i] = true;
    str[j] = true;
    str[k] = true;
    str[l] = true;
    return str
end

#
#
####################################################################################################################################################

# Some Measuremnts: Mostly wokring on a LxL toric code and predefined partitioning

function entanglement_negativity!(state::MixedDestabilizer,
                                system)
    """this calculates the TEN of the state on a torus, using the usual quadpartition. 
    state            : state to be evolved. It is taken by reference (hence the ! in the name by convention).
    system           : system object that has knows everyhting!
    """

    S = system.S
    P_A = system.P_A
    P_B = system.P_B
    P_C = system.P_C

    T = stab_to_gf2(stabilizerview(state))'

    e_A = rank(T'*(P_A)*S*(P_A)*T)/2
    e_B = rank(T'*(P_B)*S*(P_B)*T)/2
    e_C = rank(T'*(P_C)*S*(P_C)*T)/2
    e_AB = rank(T'*(P_A+P_B)*S*(P_A+P_B)*T)/2
    e_AC = rank(T'*(P_A+P_C)*S*(P_A+P_C)*T)/2
    e_BC = rank(T'*(P_B+P_C)*S*(P_B+P_C)*T)/2
    e_ABC = rank(T'*(P_A+P_B+P_C)*S*(P_A+P_B+P_C)*T)/2

    e_topo = -e_A-e_B-e_C+e_AB+e_AC+e_BC-e_ABC
    return -e_topo
end

function entanglement_entropy_topo!(state::MixedDestabilizer,
                            system)
    """this calculates the TEE of the state on a torus, using the usual quadpartition. 
    state            : state to be evolved. It is taken by reference (hence the ! in the name by convention).
    n                : total number of qubits.
    """

    r_A = system.r_A
    r_B = system.r_B
    r_C = system.r_C
    r_AB = system.r_AB
    r_BC = system.r_BC
    r_ABC = system.r_ABC

    e_A = entanglement_entropy(state,     # state to compute entropy for
                r_A,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )
    e_B = entanglement_entropy(state,     # state to compute entropy for
                r_B,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )
    e_C = entanglement_entropy(state,     # state to compute entropy for
                r_C,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )
    e_AB = entanglement_entropy(state,     # state to compute entropy for
                r_AB,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )
    e_AC = entanglement_entropy(state,     # state to compute entropy for
                [r_A; r_C],       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )
    e_BC = entanglement_entropy(state,     # state to compute entropy for
                r_BC,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )
    e_ABC = entanglement_entropy(state,     # state to compute entropy for
                r_ABC,       # subsystem as list of indices. Specifying full subsystem gives the van-neumann entropy
                Val(:rref) # algorithm to use (see documentation)
                )

    e_topo = -e_A-e_B-e_C+e_AB+e_AC+e_BC-e_ABC
    # print(e_topo)
    return -e_topo
end

#
#
####################################################################################################################################################

# Some Runs: Purification, Dev of TEN and Dev of TEE

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

#
#
####################################################################################################################################################

# Some Models: Be Careful that most of the above functions care that we are on a LxL torus with qubits at the edges, TC Setup!

function toric_code(system, stab_type_dist::DiscreteNonParametric)
    """this defines the kitaev model. For any other (qubit) model changing this function should in principle suffice."""

    ncells = system.cell_num

    # get a random cell to measure
    # note that cell going from 0 to ncells-1 is against typical julia convention, which indexes starting from 1
    random_cell = rand(0:(ncells-1)) 


    # get a random pauli type (not the conversion since rand(pauli_type_dist) generates a number between 1 and 3)
    stab_type::StabType = StabType(rand(stab_type_dist));
    
    # now return the respective pauli operator
    if stab_type == Star || stab_type == Plaquette 
        return tc_stab(stab_type, random_cell, system)
    else
        return ribs_stab(stab_type, random_cell, system)
    end
end

function classical_gauge_code(system)
    """this defines the Z_2 gauge theory on a 45deg rotated square lattice. For any other (qubit) model changing this function should in principle suffice.
    Amounts to making both star and plaquette into Z operators!
    No ribbons! Make sure the stab_type_dist takes that into account!
    """

    ncells = system.cell_num

    # get a random cell to measure
    # note that cell going from 0 to ncells-1 is against typical julia convention, which indexes starting from 1
    random_cell = rand(0:(ncells-1)) 
    
    # now return the respective pauli operator
    return classical_gauge_stab(random_cell, system)
end

#
#
##################################################################################################################################

# Helper functions for the Model definitions!

function tc_stab(t::StabType, i::Integer, system)
    """returns a four-body-pauli stabilier of the TC, acting on a site i
        Note the indexing descepancy, i is a site index (counting from 0)!"""
    
    n = system.nbits
    L = system.L
    if t == Star::StabType
        if i == 0 # The Periodic Boundary Condition
            bits = bit_string_ijkl(1, 2, 2*L*L, 2*L*L-1, 2*L*L)
        elseif i < L
            bits = bit_string_ijkl(2*i+1, 2*i+2, 2*(i+L*(L-1))+1, 2*(i-1)+2, 2*L*L)
        elseif mod(i, L) == 0
            bits = bit_string_ijkl(2*i+1, 2*i+2, 2*(i-L)+1, 2*(i+L-1)+2, 2*L*L)
        else
            bits = bit_string_ijkl(2*i+1, 2*i+2, 2*(i-L)+1, 2*(i-1)+2, 2*L*L)
        end
        return PauliOperator(0x0, bits, zeros(Bool, n));
    elseif t == Plaquette::StabType
        if i == L*L-1
            bits = bit_string_ijkl(2*L*L-1, 2*(L-1)+2, 2*L*(L-1)+1, 2*L*L, 2*L*L)
        elseif i >= L*(L-1)
            bits = bit_string_ijkl(2*i+1, 2*mod(i, L)+2, 2*(i+1)+1, 2*i+2, 2*L*L)
        elseif mod(i, L) == L-1
            bits = bit_string_ijkl(2*i+1, 2*(i+L)+2, 2*(i-L+1)+1, 2*i+2, 2*L*L)
        else
            bits = bit_string_ijkl(2*i+1, 2*(i+L)+2, 2*(i+1)+1, 2*i+2, 2*L*L)
        end
        return PauliOperator(0x0, zeros(Bool, n), bits);
    end
end

function ribs_stab(t::StabType, i::Integer, system)
    """returns a one-to-two-body-pauli short ribbon operator, acting on a site i
        Note the indexing descepancy, i is a site index (counting from 0)!"""
    
    n = system.nbits
    L = system.L
    if t == e_hor::StabType
        if i < L
            bits = bit_string_i(2*(i+L*(L-1))+1, 2*L*L)
        else
            bits = bit_string_i(2*(i-L)+1, 2*L*L)
        end
        return PauliOperator(0x0, zeros(Bool, n), bits);

    elseif t == e_ver::StabType
        if mod(i, L) == 0
            bits = bit_string_i(2*(i+L-1)+2, 2*L*L)
        else
            bits = bit_string_i(2*(i-1)+2, 2*L*L)
        end
        return PauliOperator(0x0, zeros(Bool, n), bits);

    elseif t == m_hor::StabType
        bits = bit_string_i(2*i+2, 2*L*L)
        return PauliOperator(0x0, bits, zeros(Bool, n));

    elseif t == m_ver::StabType
        bits = bit_string_i(2*i+1, 2*L*L)
        return PauliOperator(0x0, bits, zeros(Bool, n));

    elseif t == f_hor::StabType
        if i < L
            bits_z = bit_string_i(2*(i+L*(L-1))+1, 2*L*L)
        else
            bits_z = bit_string_i(2*(i-L)+1, 2*L*L)
        end
        bits_x = bit_string_i(2*i+2, 2*L*L)
        return PauliOperator(0x0, bits_x, bits_z);

    elseif t == f_ver::StabType
        if mod(i, L) == 0
            bits_z = bit_string_i(2*(i+L-1)+2, 2*L*L)
        else
            bits_z = bit_string_i(2*(i-1)+2, 2*L*L)
        end
        bits_x = bit_string_i(2*i+1, 2*L*L)
        return PauliOperator(0x0, bits_x, bits_z);
    end
end

function classical_gauge_stab(i::Integer, system)
    """ Uses the VertexLattice, so cell index is the same as the qubit index.
    There is only one stabiliser type, the plaquette! 
    Indexing of object from 0 and arrays from 1! """
    
    ncell = system.cell_num
    nbits = system.nbits
    L = system.L

    if i == ncell - 1
        bits_z = bit_string_ijkl(ncell, L, L*(L-1)+1, 1, nbits)
    elseif mod(i, L) == L - 1
        bits_z = bit_string_ijkl(i+1, i-L+2, i+L+1, i+2, nbits)
    elseif i >= ncell-L
        bits_z = bit_string_ijkl(i+1, i+2, mod(i, L)+1, mod(i, L)+2, nbits)
    else
        bits_z = bit_string_ijkl(i+1, i+2, i+L+1, i+L+2, nbits)
    end

    return PauliOperator(0x0, zeros(Bool, nbits), bits_z)
end

#
#
##################################################################################################################################

# Iterators and other

function iterate_measurements_only!(state::MixedDestabilizer, system, get_random_pauli::Function, full_steps::Integer)
    """ Evolves the state by a set number of full steps (intensive variable).
    state - the refernece to a state involved, it is mutable and will change
    system - details of the geometry
    get_random_pauli - defines the model
    full_steps - the number of steps in units of system size
    """

    nbits = system.nbits

    for i_n = 1:nbits
        for i_s = 1:full_steps
            projectrand!(state, get_random_pauli())
        end
    end
    
    return 0
end
