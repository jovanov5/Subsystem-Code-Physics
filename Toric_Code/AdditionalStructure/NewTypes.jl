# New Types

struct VertexSquareLattice
    # Dimensions of lattice on the torus
    L::Int
    nbits::Int
    cell_num::Int
    # Projectors and Ranges used in TEN and TEE measurments, assuming ususal 
    subsys::Int
    subsys_size::Int
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
    subsys_size::Int
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

struct VertexHoneyLattice # Isometric to the EdgeSquareLattice up to some sheering!
    # Dimensions of lattice on the torus
    L::Int
    nbits::Int
    cell_num::Int
    # Projectors and Ranges used in TEN and TEE measurments, assuming ususal 
    subsys::Int
    subsys_size::Int
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
    subsys_size = k*L*subsys
    P_A = Diagonal(repeat(vcat(repeat([true], subsys_size), repeat([false], 3*subsys_size)), 2))
    P_B = Diagonal(repeat(vcat(repeat([false], subsys_size), repeat([true], subsys_size), repeat([false], 2*subsys_size)), 2))
    P_C = Diagonal(repeat(vcat(repeat([false], 2*subsys_size), repeat([true], subsys_size), repeat([false], subsys_size)), 2))
    S = kron([false true; true false], Diagonal(repeat([true], k*L*L)));

    # For TEE: refined qubit indices for diferent continous partitions
    r_A = 0*subsys_size+1:1*subsys_size
    r_B = 1*subsys_size+1:2*subsys_size
    r_C = 2*subsys_size+1:3*subsys_size
    r_AB = 0*subsys_size+1:2*subsys_size
    r_BC = 1*subsys_size+1:3*subsys_size
    r_ABC = 0*subsys_size+1:3*subsys_size

    return VertexSquareLattice(L, nbits, cell_num, subsys, subsys_size, P_A, P_B, P_C, S, r_A, r_B, r_C, r_AB, r_BC, r_ABC)
end

function Init_EdgeSquareLattice_SimpleAnuli(subsys::Int) # The regions A, B and C are all non-contractable anuli!
    
    # Dimensions
    L = 4*subsys
    k = 2 # Number of Qubits per unit cell
    cell_num = L*L
    nbits = k*cell_num
    
    # For TEE: refined qubit indices for diferent continous partitions
    r_A = 0*subsys_size+1:1*subsys_size
    r_B = 1*subsys_size+1:2*subsys_size
    r_C = 2*subsys_size+1:3*subsys_size
    r_AB = 0*subsys_size+1:2*subsys_size
    r_BC = 1*subsys_size+1:3*subsys_size
    r_ABC = 0*subsys_size+1:3*subsys_size

    # For TEN: define the projection onto the ususal multipartitions of the torus
    a_mask = falses(nbits)
    a_mask[r_A] .= true
    P_A = Diagonal(repeat(a_mask, 2))
    b_mask = falses(nbits)
    b_mask[r_B] .= true
    P_B = Diagonal(repeat(b_mask, 2))
    c_mask = falses(nbits)
    c_mask[r_C] .= true
    P_C = Diagonal(repeat(c_mask, 2))
    S = kron([false true; true false], Diagonal(repeat([true], k*L*L)));

    return EdgeSquareLattice(L, nbits, cell_num, subsys, subsys_size, P_A, P_B, P_C, S, r_A, r_B, r_C, r_AB, r_BC, r_ABC)
end

function Init_EdgeSquareLattice_KitaevOrg(L::Integer, d::Integer) # The regions A, B and C are all contractable disks, but they corner touch!
    
    @assert (L - d) % 2 == 0 "The lattice size L and the region size d must be such that (L - d) is even!"

    # Dimensions
    k = 2
    cell_num = L*L
    nbits = k*cell_num

    # For TEE: refined qubit indices for diferent continous partitions
    r_A_cell = [i for i in 0:(L*L-1) if mod(i, L) < (L-d)/2 && div(i, L) < (L-d)/2] # Cell index
    r_A = sort([(r_A_cell.*2).+1; (r_A_cell.*2).+2]) # Qubit index
    r_B_cell = [i for i in 0:(L*L-1) if mod(i, L) >= (L-d)/2 && mod(i, L) < (L-d) && div(i, L) < (L-d)/2]
    r_B = sort([(r_B_cell.*2).+1; (r_B_cell.*2).+2])
    r_C_cell = [i for i in 0:(L*L-1) if mod(i, L) < (L-d) && div(i, L) >= (L-d)/2 && div(i, L) < (L-d)]
    r_C = sort([(r_C_cell.*2).+1; (r_C_cell.*2).+2])
    r_AB = [r_A; r_B]
    r_BC = [r_B; r_C]
    r_ABC = [r_A; r_B; r_C]

    # Unneeded for this geometry
    subsys = d # Thickes of the region around A, B and C!
    subsys_size = nbits-length(r_ABC) # Misleading NAME, be careful! This is the number of qubits in the outside region D!

    # For TEN: define the projection onto the contractible multipartitions of the torus!
    a_mask = falses(nbits)
    a_mask[r_A] .= true
    P_A = Diagonal(repeat(a_mask, 2))
    b_mask = falses(nbits)
    b_mask[r_B] .= true
    P_B = Diagonal(repeat(b_mask, 2))
    c_mask = falses(nbits)
    c_mask[r_C] .= true
    P_C = Diagonal(repeat(c_mask, 2))
    S = kron([false true; true false], Diagonal(repeat([true], k*L*L)));

    return EdgeSquareLattice(L, nbits, cell_num, subsys, subsys_size, P_A, P_B, P_C, S, r_A, r_B, r_C, r_AB, r_BC, r_ABC)
end

function Init_EdgeSquareLattice_KitaevDoNuT(L::Integer, d::Integer) # The regions A, B and C are all contractable disks, but they corner touch!
    
    @assert (L - d) % 3 == 0 "The lattice size L and the region size d must be such that (L - d) is divisible by 3!"

    # Dimensions
    k = 2
    cell_num = L*L
    nbits = k*cell_num

    # For TEE: refined qubit indices for diferent continous partitions
    r_A_cell = [i for i in 0:(L*L-1) if mod(i, L) < (L-d) && div(i, L) < (L-d)/3] # Cell index
    r_A = sort([(r_A_cell.*2).+1; (r_A_cell.*2).+2]) # Qubit index
    r_B_cell_lower = [i for i in 0:(L*L-1) if mod(i, L) < (L-d)/3 && div(i, L) >= (L-d)/3 && div(i, L) < 2*(L-d)/3]
    r_B_cell_upper = [i for i in 0:(L*L-1) if mod(i, L) >= 2*(L-d)/3 && mod(i, L) < (L-d) && div(i, L) >= (L-d)/3 && div(i, L) < 2*(L-d)/3]
    r_B_cell = sort([r_B_cell_lower; r_B_cell_upper]) # Cell index
    r_B = sort([(r_B_cell.*2).+1; (r_B_cell.*2).+2]) # Qubit index
    r_C_cell = [i for i in 0:(L*L-1) if mod(i, L) < (L-d) && div(i, L) >= 2*(L-d)/3 && div(i, L) < (L-d)]
    r_C = sort([(r_C_cell.*2).+1; (r_C_cell.*2).+2])
    r_AB = [r_A; r_B]
    r_BC = [r_B; r_C]
    r_ABC = [r_A; r_B; r_C]

    # Unneeded for this geometry
    subsys = d # Thickes of the region around A, B and C!
    subsys_size = nbits-length(r_ABC) # Misleading NAME, be careful! This is the number of qubits in the outside region D!

    # For TEN: define the projection onto the contractible multipartitions of the torus!
    a_mask = falses(nbits)
    a_mask[r_A] .= true
    P_A = Diagonal(repeat(a_mask, 2))
    b_mask = falses(nbits)
    b_mask[r_B] .= true
    P_B = Diagonal(repeat(b_mask, 2))
    c_mask = falses(nbits)
    c_mask[r_C] .= true
    P_C = Diagonal(repeat(c_mask, 2))
    S = kron([false true; true false], Diagonal(repeat([true], k*L*L)));

    return EdgeSquareLattice(L, nbits, cell_num, subsys, subsys_size, P_A, P_B, P_C, S, r_A, r_B, r_C, r_AB, r_BC, r_ABC)
end

function Init_VertexHoneyLattice(subsys::Int)
    
    # Dimensions
    L = 4*subsys
    k = 2 # Number of Qubits per unit cell
    cell_num = L*L
    nbits = k*cell_num
    
    # For TEN: define the projection onto the ususal multipartitions of the torus
    subsys_size = k*L*subsys # Misleading NAME, be careful!
    P_A = Diagonal(repeat(vcat(repeat([true], subsys_size), repeat([false], 3*subsys_size)), 2))
    P_B = Diagonal(repeat(vcat(repeat([false], subsys_size), repeat([true], subsys_size), repeat([false], 2*subsys_size)), 2))
    P_C = Diagonal(repeat(vcat(repeat([false], 2*subsys_size), repeat([true], subsys_size), repeat([false], subsys_size)), 2))
    S = kron([false true; true false], Diagonal(repeat([true], k*L*L)));

    # For TEE: refined qubit indices for diferent continous partitions
    r_A = 0*subsys_size+1:1*subsys_size
    r_B = 1*subsys_size+1:2*subsys_size
    r_C = 2*subsys_size+1:3*subsys_size
    r_AB = 0*subsys_size+1:2*subsys_size
    r_BC = 1*subsys_size+1:3*subsys_size
    r_ABC = 0*subsys_size+1:3*subsys_size

    return VertexHoneyLattice(L, nbits, cell_num, subsys, subsys_size, P_A, P_B, P_C, S, r_A, r_B, r_C, r_AB, r_BC, r_ABC)
end

function To_VertexHoneyLattice(system::EdgeSquareLattice)
    return VertexHoneyLattice(system.L, system.nbits, system.cell_num, system.subsys, system.subsys_size, system.P_A, system.P_B, system.P_C, system.S, system.r_A, system.r_B, system.r_C, system.r_AB, system.r_BC, system.r_ABC)
end

