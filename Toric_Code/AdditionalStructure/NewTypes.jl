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