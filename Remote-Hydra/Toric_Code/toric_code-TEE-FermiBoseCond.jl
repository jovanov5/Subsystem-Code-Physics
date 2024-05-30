using Base.Filesystem  # for creaating directories etc
using Random, Distributions  # random numbers
using HDF5  # hdf5 files
using QuantumClifford  # this is the stabilizer simulation package
using Plots # for plotting
# using Formatting # string formatting
using LinearAlgebra # some useful matrices etc.
using Missings # for missing values
using JSON3 # for reading the JSON file
dep_path = "../../AdditionalStructure/Julia/"
include("$(dep_path)NewTypes.jl")
include("$(dep_path)BitStringOps.jl")
include("$(dep_path)Measurements.jl")
include("$(dep_path)MICModels.jl")
include("$(dep_path)ModelRuns.jl")
include("$(dep_path)HelperTools.jl")
include("$(dep_path)Iterators.jl");

function main(L::Integer, d::Integer, p_f_arr::Array{Float64}, p_b_arr::Array{Float64}, t_mmt::Array{Int}, subdiv_array::Array{Int}, exp_index::Integer, debug::Int)
    """ This is the main function on terminal call.
        The simulation is of Toric Code and we try and measure the Boson and Fermion with some probabilities.
        The main measure is the TEE with a ceratin geometry (Kiatev Donut). 
    """

    description = "All Points in the Phase Diagram triangle, TEE and EE cut."
    system = Init_EdgeSquareLattice_KitaevDoNuT(L, d);
    sys_type = "Init_EdgeSquareLattice_KitaevDoNuT"
    # simulation = SimulationTime(t_final, t_mmt);
    dirpath = @__DIR__
    filename = @__FILE__
    # println("The file is: ", filename)
    n_t = length(t_mmt)
    n_pf = length(p_f_arr)
    p_f_indices = 1:n_pf
    n_pb = length(p_b_arr)
    p_b_indices = 1:n_pb
    n_subdiv = length(subdiv_array)

    all_p_arr = collect(Iterators.product(p_f_arr, p_b_arr))
    all_p_indices = collect(Iterators.product(p_f_indices, p_b_indices))
    TEE_array = fill(NaN, n_t, n_pf, n_pb)
    EE_cut_array  = fill(NaN, n_t, n_pf, n_pb, n_subdiv)

    println("filename: ", filename)
    println("description: ", description)
    println("L: ", L)
    println("d: ", d)
    println("sys_type: ", sys_type)
    println("p_f_arr: ", p_f_arr)
    println("p_b_arr: ", p_b_arr)
    # println("t_mmt: ", t_mmt)
    println("subdiv_array: ", subdiv_array)
    println("exp_index: ", exp_index)
    # println("TEE: ", TEE_array)
    # println("EE_cut: ", EE_cut_array)

    if debug == 1

        println("Debug mode is on!")

    end

    # Threads.@threads for loop_index in 1:(n_pf * n_pb)
    for loop_index in 1:(n_pf * n_pb)
        Indices = all_p_indices[loop_index]
        Probs = all_p_arr[loop_index]
        p_f_index, p_b_index = Indices
        p_f, p_b = Probs
        p_tc = 1 - p_f - p_b
        if p_tc < 0 # Easy fix for scanning the full parameter space triangle. ToDo: Implement a better way to scan the triangle.
            continue # Also ToDO! Fix the computational precision issue, some points at p_tc = 0 are not included!
        end
        stab_distro = Categorical([p_tc/2, p_tc/2, p_b, 0, p_f])
        state = toric_code_GS(system) # Get the pure TC ground state as the initial state
        t_old = 0
        for t_index in 1:n_t
            t_evol = t_mmt[t_index] - t_old
            t_old = t_mmt[t_index]
            state = iterate_measurements_only_fast!(state, system, () -> toric_code(system, stab_distro), t_evol)
            TEE_array[t_index, p_f_index, p_b_index] = entanglement_entropy_topo(state, system)
            EE_cut_array[t_index, p_f_index, p_b_index, :] = entanglement_entropy_cut(state, system, subdiv_array)
        end
    end

    save_data_prefix = "out/"
    if debug == 1
        
        # Debug: Plot EE vs cut for the first p_f and p_b and last time
        p = plot(subdiv_array, EE_cut_array[end, 1, 1, :], xlabel="Cut", ylabel="EE", marker=:circle)
        savefig(p, dirpath*"/data/debug-out/test_plot_1.pdf")

        # Debug: Plot TEE vs p_b for the first p_f and last time
        p = scatter(p_b_arr, TEE_array[end, 1, :], xlabel="p_b", ylabel="TEE")
        savefig(p, dirpath*"/data/debug-out/test_plot_2.pdf")

        save_data_prefix = "debug-out/"
    end

    # Save the data!
    outfname = dirpath*"/data/$(save_data_prefix)TEE_exp:$(exp_index).h5"
    # write output to hdf5
    h5open(outfname, "w") do outfile
        write(outfile, "filename", filename)
        write(outfile, "description", description)
        write(outfile, "L", L)
        write(outfile, "d", d)
        write(outfile, "sys_type", sys_type)
        write(outfile, "p_f_arr", p_f_arr)
        write(outfile, "p_b_arr", p_b_arr)
        write(outfile, "t_mmt", t_mmt)
        write(outfile, "subdiv_array", subdiv_array)
        write(outfile, "exp_index", exp_index)
        write(outfile, "TEE", TEE_array)
        write(outfile, "EE_cut", EE_cut_array)
    end
    return "done"
end

# Parse the arguments
arg_file = ARGS[1]
arg_string = read(arg_file, String)
arg_dict = JSON3.read(arg_string)

L = Int(arg_dict["L"])
d = Int(arg_dict["d"])
p_f_arr = Array{Float64}(arg_dict["p_f_arr"])
p_b_arr = Array{Float64}(arg_dict["p_b_arr"])
t_mmt = Array{Int}(arg_dict["t_mmt"])
subdiv_array = Array{Int}(arg_dict["subdiv_array"])
exp_index = Int(arg_dict["exp_index"])

if length(ARGS) > 1
    debug = parse(Int, ARGS[2])
else
    debug = 0
end

main(L, d, p_f_arr, p_b_arr, t_mmt, subdiv_array, exp_index, debug);
