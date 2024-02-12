""" Models: Measurement sets used to define measurement only dynamics.
"""

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
####################################################################################################################################################

# Some Models: Each functions most likely work on each a specific type of system.

function toric_code(system::EdgeSquareLattice, stab_type_dist::DiscreteNonParametric)
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

function tc_stab(t::StabType, i::Integer, system::EdgeSquareLattice)
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

function ribs_stab(t::StabType, i::Integer, system::EdgeSquareLattice)
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
