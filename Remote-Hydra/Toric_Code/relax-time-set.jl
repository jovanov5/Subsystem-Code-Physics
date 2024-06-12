using Base.Filesystem  # for creaating directories etc
using Random, Distributions  # random numbers
using HDF5  # hdf5 files
using QuantumClifford  # this is the stabilizer simulation package
using Plots # for plotting
# using Formatting # string formatting
using LinearAlgebra # some useful matrices etc.
using Missings # for missing values
using JSON3 # for reading the JSON file
include("../../AdditionalStructure/NewTypes.jl")
include("../../AdditionalStructure/BitStringOps.jl")
include("../../AdditionalStructure/Measurements.jl")
include("../../AdditionalStructure/MICModels.jl")
include("../../AdditionalStructure/ModelRuns.jl")
include("../../AdditionalStructure/HelperTools.jl")
include("../../AdditionalStructure/Iterators.jl");

function example_run(L::Integer, d::Integer, p_f_arr::Array{Float64}, p_b_arr::Array{Float64}, t_mmt::Array{Int}, subdiv_array::Array{Int})
    """ This is the main function on terminal call.
        The simulation is of Toric Code and we try and measure the Boson and Fermion with some probabilities.
        The main measure is the TEE with a ceratin geometry (Kiatev Donut). 
    """

    system = Init_EdgeSquareLattice_KitaevDoNuT(L, d);
    dirpath = @__DIR__
    filename = @__FILE__
    n_t = length(t_mmt)
    # t_final = t_mmt[n_t]
    n_pf = length(p_f_arr)
    p_f_indices = 1:n_pf
    n_pb = length(p_b_arr)
    p_b_indices = 1:n_pb
    n_subdiv = length(subdiv_array)

    all_p_arr = collect(Iterators.product(p_f_arr, p_b_arr))
    all_p_indices = collect(Iterators.product(p_f_indices, p_b_indices))
    # TEE_array = fill(NaN, n_t, n_pf, n_pb)
    EE_cut_array  = fill(NaN, n_t, n_pf, n_pb, n_subdiv)

    println("filename: ", filename)
    println("L: ", L)
    println("d: ", d)
    println("p_f_arr: ", p_f_arr)
    println("p_b_arr: ", p_b_arr)
    # println("t_mmt: ", t_mmt)
    println("subdiv_array: ", subdiv_array)
    # println("TEE: ", TEE_array)
    # println("EE_cut: ", EE_cut_array)

    Threads.@threads for loop_index in 1:(n_pf * n_pb)
    # for loop_index in 1:(n_pf * n_pb)
        Indices = all_p_indices[loop_index]
        Probs = all_p_arr[loop_index]
        p_f_index, p_b_index = Indices
        p_f, p_b = Probs
        p_tc = 1 - p_f - p_b
        if p_tc < 0 # Easy fix for scanning the full parameter space triangle. ToDo: Implement a better way to scan the triangle.
            continue # Also ToDO! Fix the computational precision issue, some points at p_tc = 0 are not included!
        end
        stab_distro = Categorical([p_tc/2, p_tc/2, p_b, 0, p_f])
        # state = toric_code_GS(system) # Get the pure TC ground state as the initial state
        state = get_random_pure_state(system) # Get a random pure state.
        t_old = 0
        for t_index in 1:n_t
            t_evol = t_mmt[t_index] - t_old
            t_old = t_mmt[t_index]

            # if mod(round(Int, t_evol/t_final*100), 5) == 0
            #     print("$(round(Int, t_evol/t_final*100))%")
            # end # Progress bar!

            state = iterate_measurements_only_fast!(state, system, () -> toric_code(system, stab_distro), t_evol)
            # TEE_array[t_index, p_f_index, p_b_index] = entanglement_entropy_topo(state, system)
            EE_cut_array[t_index, p_f_index, p_b_index, :] = entanglement_entropy_cut(state, system, subdiv_array)
        end
    end

    # Save the data!
    outfname = dirpath*"/data/TEE_exp:relax_4.h5"
    # write output to hdf5
    h5open(outfname, "w") do outfile
        write(outfile, "filename", filename)
        write(outfile, "L", L)
        write(outfile, "d", d)
        write(outfile, "p_f_arr", p_f_arr)
        write(outfile, "p_b_arr", p_b_arr)
        write(outfile, "t_mmt", t_mmt)
        write(outfile, "EE_cut", EE_cut_array)
    end
    return EE_cut_array#, TEE_array
end

p_f_arr = [0.55]
p_b_arr = [0.0]

L = 24
d = 3

t_mmt = Array{Int}(1:1000)
subdiv_array = [12]

EE_cut_array = example_run(L, d, p_f_arr, p_b_arr, t_mmt, subdiv_array)