"""
HelperTools.jl

This file contains additional tools for the Toric Code simulation, including stabilizer visualization and other useful functions.

Author: [Your Name]
Date: [Current Date]
"""

# module HelperTools :: Need to learn how to create a module in Julia!

function edge_picker(i::Int, j::Int, direction::Int)
    """
    :param i: int, the x coordinate of the vertex
    :param j: int, the y coordinate of the vertex
    :param direction: int, the direction of the edge, 0 for horizontal, 1 for vertical
    """
    return 2*(i+j*L) + direction + 1 # Julia arrays are 1-indexed.
    
end

function visualise_the_stabiliser(stabiliser::PauliOperator, system::EdgeSquareLattice)
    stab = stab_to_gf2(stabiliser)
    L = system.L
    p = plot()
    nbits = system.nbits
    # plot!([(0, 8), (0, 0), (8, 0)], color=:black, legend= false, xticks= 0:1:L, yticks= 0:1:L)
    # plot!([(0, 8), (8, 8), (8, 0)], color=:orange, legend= false, xticks= 0:1:L, yticks= 0:1:L)

    for i in 0:L-1
        for j in 0:L-1
            if stab[edge_picker(i, j, 1)] == 1
                plot!([(i, j), (i+1, j)], color=:blue, legend= false, xticks= 0:1:L, yticks= 0:1:L)
                plot!([(i+1/2, j+1/2), (i+1/2, j-1/2)], color=:blue, linestyle= :dash, legend= false, xticks= 0:1:L, yticks= 0:1:L)
            end
            if stab[edge_picker(i, j, 1)+nbits] == 1
                plot!([(i, j), (i+1, j)], color=:red, legend= false, xticks= 0:1:L, yticks= 0:1:L)
            end
            if stab[edge_picker(i, j, 0)] == 1
                plot!([(i, j), (i, j+1)], color=:blue, legend= false, xticks= 0:1:L, yticks= 0:1:L)
                plot!([(i+1/2, j+1/2), (i-1/2, j+1/2)], color=:blue, linestyle= :dash, legend= false, xticks= 0:1:L, yticks= 0:1:L)
            end
            if stab[edge_picker(i, j, 0)+nbits] == 1
                plot!([(i, j), (i, j+1)], color=:red, legend= false, xticks= 0:1:L, yticks= 0:1:L)
            end
        end
    end
    xlims!(-1, L+1)
    ylims!(-1, L+1)
    display(p)
    # println(sum(stab))
end

function z_polarised_state(system)
    """this returns the z-polarised state for the given system."""
    nbits = system.nbits
    state = Stabilizer(zeros(UInt8, nbits), # Phases
        zeros(Bool, nbits, nbits), # X Tableau (as matrix of Bools)
        Matrix(LinearAlgebra.I, nbits, nbits)); # Z Tableau
    state = MixedDestabilizer(state)
    return state
end

function z_polarised_state_temp(nbits::Int)
    """this returns the z-polarised state for the given system."""
    state = Stabilizer(zeros(UInt8, nbits), # Phases
        zeros(Bool, nbits, nbits), # X Tableau (as matrix of Bools)
        Matrix(LinearAlgebra.I, nbits, nbits)); # Z Tableau
    state = MixedDestabilizer(state)
    return state
end

function x_polarised_state(system)
    """this returns the x-polarised state for the given system."""
    nbits = system.nbits
    state = Stabilizer(ones(UInt8, nbits), # Phases
        Matrix(LinearAlgebra.I, nbits, nbits), # X Tableau (as matrix of Bools)
        zeros(Bool, nbits, nbits)); # Z Tableau
    state = MixedDestabilizer(state)
    return state
end

function y_polarised_state(system)
    """this returns the y-polarised state for the given system."""
    nbits = system.nbits
    state = Stabilizer(ones(UInt8, nbits), # Phases
        Matrix(LinearAlgebra.I, nbits, nbits), # X Tableau (as matrix of Bools)
        Matrix(LinearAlgebra.I, nbits, nbits)); # Z Tableau
    state = MixedDestabilizer(state)
    return state
end

function maximally_mixed_state(system)
    """this returns the maximally mixed state for the given system."""
    # initial state is the maximally mixed state
    # since the package does not let you create an empty tableau, I defined the identity as the stabilizer, which is the same statement
    nbits = system.nbits
    state = Stabilizer(zeros(UInt8, 1), # Phases
        zeros(Bool, 1, nbits), # X Tableau (as matrix of Bools)
        zeros(Bool, 1, nbits)); # Z Tableau
    state = MixedDestabilizer(state)
    return state
end

function maximally_mixed_state(nbits::Int)
    """this returns the maximally mixed state for the given system."""
    # initial state is the maximally mixed state
    # since the package does not let you create an empty tableau, I defined the identity as the stabilizer, which is the same statement
    state = Stabilizer(zeros(UInt8, 1), # Phases
        zeros(Bool, 1, nbits), # X Tableau (as matrix of Bools)
        zeros(Bool, 1, nbits)); # Z Tableau
    state = MixedDestabilizer(state)
    return state
end

function get_random_pure_state(system)
    return MixedDestabilizer(random_stabilizer(system.nbits))
end

function get_subdiv_array(system, n_subdiv::Int)
    """this returns the array of subdivisions for the given system."""
    L = system.L
    subdiv_array = []
    for i_sub in 1:(n_subdiv - 1)
        push!(subdiv_array, round(Int, L * i_sub / n_subdiv))
    end
    subdiv_array = Integer.(subdiv_array)
    return subdiv_array
end

function get_subdiv_array(L::Int, n_subdiv::Int)
    """this returns the array of subdivisions for the system of a given size L."""
    subdiv_array = []
    for i_sub in 1:(n_subdiv - 1)
        push!(subdiv_array, round(Int, L * i_sub / n_subdiv))
    end
    subdiv_array = Integer.(subdiv_array)
    return subdiv_array
end

function get_subdiv_array(system, n_subdiv::String)
    """this returns the array of subdivisions for the system of a given size L."""
    @assert n_subdiv == "all" "The only valid option for n_subdiv is 'all'."
    L = system.L
    subdiv_array = Array{Int}(0:L)
    return subdiv_array
end

function get_subdiv_array(L::Int, n_subdiv::String)
    """this returns the array of subdivisions for the system of a given size L."""
    @assert n_subdiv == "all" "The only valid option for n_subdiv is 'all'."
    subdiv_array = Array{Int}(0:L)
    return subdiv_array
end

function get_t_mmt_arr_simple(t_final::Int, n_t::Int)
    """this returns the array of time steps for the simulation.
    The time steps are exponentially spaced."""
    t_mmt = [round(Int, x) for x in exp2.(range(0, log2(t_final), length=n_t))]    
    return t_mmt
end

function get_t_mmt_arr_lin(t_final::Int, d_t::Int)
    """this returns the array of time steps for the simulation.
    The time steps are linearly spaced."""
    t_mmt = Array{Int}(d_t:d_t:t_final)
    return t_mmt
end

function get_t_mmt_arr_refined(t_final::Int, n_t::Int)
    """this returns the array of time steps for the simulation.
    The time steps are exponentially spaced,
    exept that we have added a lin span between last two exponential steps."""
    t_mmt_simple = get_t_mmt_arr_simple(t_final, n_t)    
    t_mmt_final = t_mmt_simple[1:end-1]
    last_steps = t_mmt_final[end]
    while t_mmt_final[end] < t_final
        push!(t_mmt_final, t_mmt_final[end] + last_steps)
    end
    t_mmt_final[end] = t_final
    return t_mmt_final
end

# Helper routines for the Zassenhauser algorithm

function zassenhausen_map_f(g::PauliOperator)
    return g ⊗ g ⊗ g ⊗ g
end

function zassenhausen_map_h(g::PauliOperator)
    e = g * g # this is the identity operator, cool trick 
    return g ⊗ g ⊗ e ⊗ e
end

function zassenhausen_state(deformator, probe_state, nbits::Int)
    """Given two sungroup of the Pauli group, given as invalid stabilisers, we return the "Zassenhauser" state.
    """
    state = maximally_mixed_state(4*nbits)

    for g in deformator
        project!(state, zassenhausen_map_f(g), keep_result= false)
    end

    for g in probe_state
        project!(state, zassenhausen_map_h(g), keep_result= false)
    end

    return state
end

function zassenhausen_state(charge_representative, deformator, probe_state, nbits::Int)
    """Given two sungroup of the Pauli group, given as invalid stabilisers, we return the "Zassenhauser" state.
    """
    state = maximally_mixed_state(4*nbits)

    project!(state, zassenhausen_map_f(charge_representative), keep_result= false)
    
    for g in deformator
        project!(state, zassenhausen_map_f(g), keep_result= false)
    end

    for g in probe_state
        project!(state, zassenhausen_map_h(g), keep_result= false)
    end

    return state
end

function zassenhausen_alg(deformator, probe_state, nbits::Int)
    """Given two sungroup of the Pauli group, given as invalid stabilisers generators, we return the dims we need.
    """
    state = zassenhausen_state(deformator, probe_state, nbits)
    N = nbits
    dim_intersection = 2*N - entanglement_entropy(state, (2*N+1):4N, Val(:clip))
    dim_sum = 4*N - dim_intersection - entanglement_entropy(state, 1:4N, Val(:clip))
    return dim_intersection, dim_sum
end

function zassenhausen_alg(charge_representative, deformator, probe_state, nbits::Int)
    """Given two sungroup of the Pauli group, given as invalid stabilisers generators, we return the dims we need.
    """
    state = zassenhausen_state(charge_representative, deformator, probe_state, nbits)
    N = nbits
    dim_intersection = 2*N - entanglement_entropy(state, (2*N+1):4N, Val(:clip))
    dim_sum = 4*N - dim_intersection - entanglement_entropy(state, 1:4N, Val(:clip))
    return dim_intersection, dim_sum
end


function taxi_metrix_torus(x1, y1, x2, y2, Lx, Ly)
    dx = abs(x1 - x2)
    dy = abs(y1 - y2)
    dx = min(dx, Lx - dx)
    dy = min(dy, Ly - dy)
    return dx + dy
end

function taxi_metrix_torus(x1::Tuple{Int, Int}, x2::Tuple{Int, Int}, Lx, Ly)
    return taxi_metrix_torus(x1[1], x1[2], x2[1], x2[2], Lx, Ly)
end

function euclidean_metrix_torus(x1, y1, x2, y2, Lx, Ly)
    dx = abs(x1 - x2)
    dy = abs(y1 - y2)
    dx = min(dx, Lx - dx)
    dy = min(dy, Ly - dy)
    return sqrt(dx^2 + dy^2)
end

function euclidean_metrix_torus(x1::Tuple{Int, Int}, x2::Tuple{Int, Int}, Lx, Ly)
    return euclidean_metrix_torus(x1[1], x1[2], x2[1], x2[2], Lx, Ly)
end

function massage_the_zess_corr(corr_mat, system, n_t)
    """Does the time average and the distance average of the correlation matrix in the taxi cab metric on a torus.
    """
    L = system.L
    
    all_r = []
    for i = 1:L
        for j = 1:L
            x_i = i - 1
            y_i = j - 1
            push!(all_r, taxi_metrix_torus((0, 0), (x_i, y_i), L, L))
        end
    end

    r_max, _ = findmax(all_r)
    plot_x = [0:r_max]

    norm_y = zeros(Float64, r_max + 1)
    for i = 1:L
        for j = 1:L
            x_i = i - 1
            y_i = j - 1
            norm_y[taxi_metrix_torus((0, 0), (x_i, y_i), L, L) + 1] += 1
        end
    end

    corr_mat_time_avr = zeros(Float64, L, L)
    for t_index in 1:n_t
        corr_mat_time_avr += corr_mat[t_index, :, :]
    end
    
    plot_y = zeros(Float64, r_max + 1)
    for i = 1:L
        for j = 1:L
            x_i = i - 1
            y_i = j - 1
            plot_y[taxi_metrix_torus((0, 0), (x_i, y_i), L, L) + 1] += corr_mat_time_avr[i, j]/n_t
        end
    end

    return plot_x, plot_y ./ norm_y
end
