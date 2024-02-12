"""
Some auxiliary functions used at various points. Augments the QunatumClifford.jl package and does not go beyond qubits!

Additional Functions:
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

function topological_partitioning(subsys::Integer)
    """given a size of a quarter of a LxL torus this generates useful quantities needed for topological measurements
    """
    L = 4*subsys  # linear system size

    # For TEN: define the projection onto the ususal multipartitions of the torus
    P_A = Diagonal(repeat(vcat(repeat([true], 8*subsys*subsys), repeat([false], 24*subsys*subsys)), 2))
    P_B = Diagonal(repeat(vcat(repeat([false], 8*subsys*subsys), repeat([true], 8*subsys*subsys), repeat([false], 16*subsys*subsys)), 2))
    P_C = Diagonal(repeat(vcat(repeat([false], 16*subsys*subsys), repeat([true], 8*subsys*subsys), repeat([false], 8*subsys*subsys)), 2))
    P_D = Diagonal(repeat(vcat(repeat([false], 24*subsys*subsys), repeat([true], 8*subsys*subsys)), 2))
    S = kron([false true; true false], Diagonal(repeat([true], 2*L*L)));

    # For TEE: refined qubit indices for diferent continous partitions
    cell_size = 8*subsys*subsys
    r_A = 0*cell_size+1:1*cell_size
    r_B = 1*cell_size+1:2*cell_size
    r_C = 2*cell_size+1:3*cell_size
    r_AB = 0*cell_size+1:2*cell_size
    r_BC = 1*cell_size+1:3*cell_size
    r_ABC = 0*cell_size+1:3*cell_size
    
    return L, P_A, P_B, P_C, P_D, S, r_A, r_B, r_C, r_AB, r_BC, r_ABC
end 

function entanglement_negativity!(state::MixedDestabilizer,
                                n::Integer)
    """this calculates the TEN of the state on a torus, using the usual quadpartition. 
    state            : state to be evolved. It is taken by reference (hence the ! in the name by convention).
    n                : total number of qubits.

    The projectors onto the support of the regions are predefined: P_A, P_B, ...
    The "symplectic" 2n-by-2n matrix is also predefined: S.
    """

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
                            n::Integer)
    """this calculates the TEE of the state on a torus, using the usual quadpartition. 
    state            : state to be evolved. It is taken by reference (hence the ! in the name by convention).
    n                : total number of qubits.
    """

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
                                       nt::Integer,
                                       measure_ts::AbstractArray)
    """this runs the mesaurement only dynamics. 
       state            : state to be evolved. It is taken by reference (hence the ! in the name by convention).
       get_random_pauli : a function that is called each time step and should generate the paulis to be measured.
       nt               : the total number of time steps (one time step is defined as number-of-sites measurements
       measure_ts       : times at which entropy is measured)
       """
    
    # get system size n (number of qubits)
    n = size(stabilizerview(state))[2]
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
                                       nt::Integer,
                                       measure_ts::AbstractArray)
    """this runs the mesaurement only dynamics. 
       state            : state to be evolved. It is taken by reference (hence the ! in the name by convention).
       get_random_pauli : a function that is called each time step and should generate the paulis to be measured.
       nt               : the total number of time steps (one time step is defined as number-of-sites measurements
       measure_ts       : times at which entropy is measured)
       """
    
    # get system size n (number of qubits)
    n = size(stabilizerview(state))[2]
    # this will store the total entropy (con Neumann) as a function of time
    TENs = Vector{Float64}([entanglement_negativity!(state, n)])
    
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
                                     n,       # number of qubits on the square lattice torus.
                                     )
            push!(TENs, ten)
        end
    end
    
    # print(stabilizerview(state))
    return TENs 
end


function run_measurement_only_dynamics_TEE!(state::MixedDestabilizer,
                                       get_random_pauli::Function,
                                       nt::Integer,
                                       measure_ts::AbstractArray)
    """this runs the mesaurement only dynamics. 
       state            : state to be evolved. It is taken by reference (hence the ! in the name by convention).
       get_random_pauli : a function that is called each time step and should generate the paulis to be measured.
       nt               : the total number of time steps (one time step is defined as number-of-sites measurements
       measure_ts       : times at which entropy is measured)
       """
    
    # get system size n (number of qubits)
    n = size(stabilizerview(state))[2]
    # this will store the total entropy (con Neumann) as a function of time
    TEEs = Vector{Float64}([entanglement_entropy_topo!(state, n)])
    
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
            tee = entanglement_entropy_topo!(state,     # state to compute TEN for
                                     n,       # number of qubits on the square lattice torus.
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

function toric_code(L, stab_type_dist::DiscreteNonParametric)
    """this defines the kitaev model. For any other (qubit) model changing this function should in principle suffice."""

    ncells = L*L

    # get a random cell to measure
    # note that cell going from 0 to ncells-1 is against typical julia convention, which indexes starting from 1
    random_cell = rand(0:(ncells-1)) 


    # get a random pauli type (not the conversion since rand(pauli_type_dist) generates a number between 1 and 3)
    stab_type::StabType = StabType(rand(stab_type_dist));
    
    # now return the respective pauli operator
    if stab_type == Star || stab_type == Plaquette 
        return tc_stab(stab_type, random_cell, L)
    else
        return ribs_stab(stab_type, random_cell, L)
    end
end

function classical_gauge_code(L, stab_type_dist::DiscreteNonParametric)
    """this defines the Z_2 gauge theory on a 45deg rotated square lattice. For any other (qubit) model changing this function should in principle suffice.
    Amounts to making both star and plaquette into Z operators!
    No ribbons! Make sure the stab_type_dist takes that into account!
    """

    ncells = L*L

    # get a random cell to measure
    # note that cell going from 0 to ncells-1 is against typical julia convention, which indexes starting from 1
    random_cell = rand(0:(ncells-1)) 


    # get a random pauli type (not the conversion since rand(pauli_type_dist) generates a number between 1 and 3)
    stab_type::StabType = StabType(rand(stab_type_dist));
    
    # now return the respective pauli operator
    return classical_gauge_stab(stab_type, random_cell, L)
end

#
#
##################################################################################################################################

# Stabiliser Setup for Toric Code on a square lattice on a LxL torus

function tc_stab(t::StabType, i::Integer, L::Integer)
    """returns a four-body-pauli stabilier of the TC, acting on a site i
        Note the indexing descepancy, i is a site index (counting from 0)!"""
    
    n = 2*L*L
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

function ribs_stab(t::StabType, i::Integer, L::Integer)
    """returns a one-to-two-body-pauli short ribbon operator, acting on a site i
        Note the indexing descepancy, i is a site index (counting from 0)!"""
    
    n = 2*L*L
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

function classical_gauge_stab(t::StabType, i::Integer, L::Integer)
    """returns a four-body-pauli stabilier of the TC, acting on a site i
        Note the indexing descepancy, i is a site index (counting from 0)!"""
    
    n = 2*L*L
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
        return PauliOperator(0x0, zeros(Bool, n), bits);
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