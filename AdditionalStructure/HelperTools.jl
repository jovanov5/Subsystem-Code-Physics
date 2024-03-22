"""
HelperTools.jl

This file contains additional tools for the Toric Code simulation, including stabilizer visualization and other useful functions.

Author: [Your Name]
Date: [Current Date]
"""

# module HelperTools :: Need to learn how to create a module in Julia!

function visualise_the_stabiliser(stabiliser::PauliOperator, system::EdgeSquareLattice)
    stab = stab_to_gf2(stabiliser)
    L = system.L
    p = plot()
    nbits = system.nbits
    # plot!([(0, 8), (0, 0), (8, 0)], color=:black, legend= false, xticks= 0:1:L, yticks= 0:1:L)
    # plot!([(0, 8), (8, 8), (8, 0)], color=:orange, legend= false, xticks= 0:1:L, yticks= 0:1:L)

    for i in 0:L-1
        for j in 0:L-1
            if stab[edge_picker(i+j*L, 1)] == 1
                plot!([(i, j), (i+1, j)], color=:blue, legend= false, xticks= 0:1:L, yticks= 0:1:L)
                plot!([(i+1/2, j+1/2), (i+1/2, j-1/2)], color=:blue, linestyle= :dash, legend= false, xticks= 0:1:L, yticks= 0:1:L)
            end
            if stab[edge_picker(i+j*L, 1)+nbits] == 1
                plot!([(i, j), (i+1, j)], color=:red, legend= false, xticks= 0:1:L, yticks= 0:1:L)
            end
            if stab[edge_picker(i+j*L, 0)] == 1
                plot!([(i, j), (i, j+1)], color=:blue, legend= false, xticks= 0:1:L, yticks= 0:1:L)
                plot!([(i+1/2, j+1/2), (i-1/2, j+1/2)], color=:blue, linestyle= :dash, legend= false, xticks= 0:1:L, yticks= 0:1:L)
            end
            if stab[edge_picker(i+j*L, 0)+nbits] == 1
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
