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

# end
