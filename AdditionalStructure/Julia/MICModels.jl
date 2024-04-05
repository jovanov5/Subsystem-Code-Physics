""" Models: Measurement sets used to define measurement only dynamics.
"""

# define enum just for just for better readibility of code (Do not need PauliY in this case!)

@enum PauliType begin
    PauliX = 1
    PauliY = 2   
    PauliZ = 3
end

@enum StabTypeTC begin
    Star = 1
    Plaquette = 2
    Electric = 3
    Magnetic = 4
    Fermionic = 5
end

@enum StabTypeHC begin
    XX = 1
    YY = 2
    ZZ = 3
end

@enum Orientation begin
    Vertical = 1
    Horizontal = 2
end

#
#
####################################################################################################################################################

# Some Models: Each functions most likely work on each a specific type of system.

function toric_code(system::EdgeSquareLattice, stab_type_dist::DiscreteNonParametric)
    """this defines the TC model. For any other (qubit) model changing this function should in principle suffice."""

    ncells = system.cell_num

    # get a random cell to measure
    # note that cell going from 0 to ncells-1 is against typical julia convention, which indexes starting from 1
    random_cell = rand(0:(ncells-1)) 


    # get a random pauli type (not the conversion since rand(pauli_type_dist) generates a number between 1 and 3)
    stab_type::StabTypeTC = StabTypeTC(rand(stab_type_dist));
    
    # now return the respective pauli operator
    if stab_type == Star || stab_type == Plaquette 
        return tc_stab(stab_type, random_cell, system)
    else
        return ribs_stab(stab_type, random_cell, system)
    end
end

function toric_code_GS(system::EdgeSquareLattice)
    """this return one of the four coputational logic states.
    At the moment it is only working for (0,0)"""

    nbits = system.nbits
    state = z_polarised_state(system)
    L = system.L

    for cell_index = 0:L*L-1
        state, anticom_index, result = project!(state, tc_stab(Star, cell_index, system))
        # See how to force the projection result!
    end

    for cell_index = 0:L*L-1
        state, anticom_index, result = project!(state, tc_stab(Plaquette, cell_index, system)) 
        # See how to force the projection result!
    end


    return state
end

function toric_code_GS_MIXED(system::EdgeSquareLattice)
    """this return fully mixed state in the GS manifold."""

    nbits = system.nbits
    state = maximally_mixed_state(system)
    L = system.L

    for cell_index = 0:L*L-1
        state, anticom_index, result = project!(state, tc_stab(Star, cell_index, system)) 
        # See how to force the projection result!
    end

    for cell_index = 0:L*L-1
        state, anticom_index, result = project!(state, tc_stab(Plaquette, cell_index, system)) 
        # See how to force the projection result!
    end

    return state
end

function kitaev_code(system::VertexHoneyLattice, stab_type_dist::DiscreteNonParametric)
    """this defines the kitaev model. For any other (qubit) model changing this function should in principle suffice."""

    ncells = system.cell_num
    nbits = system.nbits
    L = system.L

    # get a random cell to measure
    # note that cell going from 0 to ncells-1 is against typical julia convention, which indexes starting from 1
    random_cell = rand(0:(ncells-1)) 


    # get a random pauli type (not the conversion since rand(pauli_type_dist) generates a number between 1 and 3)
    stab_type::StabTypeHC = StabTypeHC(rand(stab_type_dist));
    
    # now return the respective pauli operator    
    if stab_type == XX::StabTypeHC
        adjacent_cell = mod(random_cell+1, L) + div(random_cell, L) * L
        bits = bit_string_ij(2*random_cell+1, 2*adjacent_cell+2, nbits)
        return PauliOperator(0x0, bits, zeros(Bool, nbits))
    elseif stab_type == YY::StabTypeHC
        adjacent_cell = mod(random_cell, L) + mod(div(random_cell, L)+1, L) * L
        bits = bit_string_ij(2*random_cell+2, 2*adjacent_cell+1, nbits)
        return PauliOperator(0x0, bits, bits)
    else
        adjacent_cell = random_cell
        bits = bit_string_ij(2*random_cell+1, 2*adjacent_cell+2, nbits)
        return PauliOperator(0x0, zeros(Bool, nbits), bits)
    end
end

function classical_gauge_code(system::VertexSquareLattice)
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

function edge_picker(cell_index::Integer, edge_index::Integer)
    # returns the array index (from 1) of the edge DoF! Cell and Edge Indices are from 0!
    return 2*cell_index + edge_index + 1
end

function tc_stab(t::StabTypeTC, cell_index::Integer, system::EdgeSquareLattice)
    """returns a four-body-pauli stabilier of the TC, acting on a site i
        Note the indexing descepancy, i is a site index (counting from 0)!"""
    
    nbits = system.nbits
    L = system.L
    if t == Star::StabTypeTC
        adj_cell_left = mod(cell_index, L) + mod(div(cell_index, L)-1, L) * L
        adj_cell_down = mod(cell_index-1, L) + div(cell_index, L) * L
        bits = bit_string_ijkl(edge_picker(cell_index, 0), edge_picker(cell_index, 1), 
                        edge_picker(adj_cell_down, 1), edge_picker(adj_cell_left, 0), nbits)
        return PauliOperator(0x0, bits, zeros(Bool, nbits));
    elseif t == Plaquette::StabTypeTC
        adj_cell_right = mod(cell_index, L) + mod(div(cell_index, L)+1, L) * L
        adj_cell_up = mod(cell_index+1, L) + div(cell_index, L) * L
        bits = bit_string_ijkl(edge_picker(cell_index, 0), edge_picker(cell_index, 1), 
                        edge_picker(adj_cell_up, 0), edge_picker(adj_cell_right, 1), nbits)
        return PauliOperator(0x0, zeros(Bool, nbits), bits);
    end
end

function ribs_stab(t::StabTypeTC, cell_index::Integer, system::EdgeSquareLattice)
    """returns a one-to-two-body-pauli short ribbon operator, acting on a site i
        Note the indexing descepancy, i is a site index (counting from 0)!"""
    
    nbits = system.nbits
    L = system.L

    direction = Orientation(rand(1:2))

    if t == Electric::StabTypeTC
        if direction == Vertical::Orientation
            adj_cell_down = mod(cell_index-1, L) + div(cell_index, L) * L
            bits = bit_string_i(edge_picker(adj_cell_down, 1), nbits)
        else
            adj_cell_left = mod(cell_index, L) + mod(div(cell_index, L)-1, L) * L
            bits = bit_string_i(edge_picker(adj_cell_left, 0), nbits)
        end
        return PauliOperator(0x0, zeros(Bool, nbits), bits);
    elseif t == Magnetic::StabTypeTC
        if direction == Vertical::Orientation
            bits = bit_string_i(edge_picker(cell_index, 0), nbits)
        else
            bits = bit_string_i(edge_picker(cell_index, 1), nbits)
        end
        return PauliOperator(0x0, bits, zeros(Bool, nbits));
    else
        if direction == Vertical::Orientation
            adj_cell_down = mod(cell_index-1, L) + div(cell_index, L) * L
            bits_z = bit_string_i(edge_picker(adj_cell_down, 1), nbits)
            bits_x = bit_string_i(edge_picker(cell_index, 0), nbits)
        else
            adj_cell_left = mod(cell_index, L) + mod(div(cell_index, L)-1, L) * L
            bits_z = bit_string_i(edge_picker(adj_cell_left, 0), nbits)
            bits_x = bit_string_i(edge_picker(cell_index, 1), nbits)
        end
        return PauliOperator(0x0, bits_x, bits_z);
    end
end

function classical_gauge_stab(i::Integer, system::VertexSquareLattice)
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
