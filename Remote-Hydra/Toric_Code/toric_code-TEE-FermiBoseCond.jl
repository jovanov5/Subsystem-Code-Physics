using Base.Filesystem  # for creaating directories etc
using Random, Distributions  # random numbers
using HDF5  # hdf5 files
using QuantumClifford  # this is the stabilizer simulation package
using Plots # for plotting
# using Formatting # string formatting
using LinearAlgebra # some useful matrices etc.
using Missings # for missing values
include("../../AdditionalStructure/NewTypes.jl")
include("../../AdditionalStructure/BitStringOps.jl")
include("../../AdditionalStructure/Measurements.jl")
include("../../AdditionalStructure/MICModels.jl")
include("../../AdditionalStructure/ModelRuns.jl")
include("../../AdditionalStructure/HelperTools.jl")
include("../../AdditionalStructure/Iterators.jl");

function main(L::Integer, d::Integer, p_f_max::Float64, n_pf::Integer, p_b_max::Float64, n_pb::Integer, t_final::Integer, n_t::Integer, n_subdiv::Integer, exp_index::Integer)
    """ This is the main function on terminal call.
        The simulation is of Toric Code and we try and measure the Boson and Fermion with some probabilities.
        The main measure is the TEE with a ceratin geometry (Kiatev Donut). 
    """

    system = Init_EdgeSquareLattice_KitaevDoNuT(L, d);
    sys_type = "Init_EdgeSquareLattice_KitaevDoNuT"
    t_mmt = get_t_mmt_arr_refined(t_final, n_t);
    n_t_new = length(t_mmt)
    simulation = SimulationTime(t_final, t_mmt);
    dirpath = @__DIR__
    filename = @__FILE__
    # println("The file is: ", filename)
    p_f_arr = [x for x in range(0, p_f_max, length=n_pf)];
    p_f_indices = 1:n_pf;
    p_b_arr = [x for x in range(0, p_b_max, length=n_pb)];
    p_b_indices = 1:n_pb;
    all_p_arr = collect(Iterators.product(p_f_arr, p_b_arr))
    all_p_indices = collect(Iterators.product(p_f_indices, p_b_indices))
    subdiv_array = get_subdiv_array(system, n_subdiv)
    TEE_array = fill(NaN, n_t_new, n_pf, n_pb)
    EE_cut_array  = fill(NaN, n_t_new, n_pf, n_pb, n_subdiv - 1)

    println("filename: ", filename)
    println("L: ", L)
    println("d: ", d)
    println("sys_type: ", sys_type)
    println("p_f_max: ", p_f_max)
    println("n_pf: ", n_pf)
    println("p_f_arr: ", p_f_arr)
    println("p_b_max: ", p_b_max)
    println("n_pb: ", n_pb)
    println("p_b_arr: ", p_b_arr)
    println("t_final: ", t_final)
    println("n_t: ", n_t)
    println("n_t_new: ", n_t_new)
    println("t_mmt: ", t_mmt)
    println("n_subdiv: ", n_subdiv)
    println("exp_index: ", exp_index)
    println("subdiv_array: ", subdiv_array)
    # println("TEE: ", TEE_array)
    # println("EE_cut: ", EE_cut_array)

    # Threads.@threads for loop_index in 1:(n_pf * n_pb)
    for loop_index in 1:(n_pf * n_pb)
        Indices = all_p_indices[loop_index]
        Probs = all_p_arr[loop_index]
        p_f_index, p_b_index = Indices
        p_f, p_b = Probs
        p_tc = 1 - p_f - p_b
        if p_tc < 0 # Easy fix for scanning the full parameter space triangle. ToDo: Implement a better way to scan the triangle.
            continue
        end
        stab_distro = Categorical([p_tc/2, p_tc/2, p_b, 0, p_f])
        state = toric_code_GS(system) # Get the pure TC ground state as the initial state
        t_old = 0
        for t_index in 1:n_t_new
            t_evol = t_mmt[t_index] - t_old
            t_old = t_mmt[t_index]
            state = iterate_measurements_only_fast!(state, system, () -> toric_code(system, stab_distro), t_evol)
            TEE_array[t_index, p_f_index, p_b_index] = entanglement_entropy_topo(state, system)
            EE_cut_array[t_index, p_f_index, p_b_index, :] = entanglement_entropy_cut(state, system, n_subdiv)
        end
    end

    debug = false
    save_data_prefix = ""
    if debug
        # Debug: Plot EE vs cut for the first p_f and p_b and last time
        p = plot(subdiv_array, EE_cut_array[end, 2, 2, :], xlabel="Cut", ylabel="EE", marker=:circle)
        savefig(p, dirpath*"/test_plot_1.pdf")

        # Debug: Plot TEE vs p_b for the last p_f and last time
        p = scatter(p_b_arr, TEE_array[end, end, :], xlabel="p_b", ylabel="TEE")
        savefig(p, dirpath*"/test_plot_2.pdf")

        save_data_prefix = "test_"
    end

    # Save the data!
    outfname = dirpath*"/data/$(save_data_prefix)TEE_exp:$(exp_index).h5"
    # write output to hdf5
    h5open(outfname, "w") do outfile
        write(outfile, "filename", filename)
        write(outfile, "L", L)
        write(outfile, "d", d)
        write(outfile, "sys_type", sys_type)
        write(outfile, "p_f_max", p_f_max)
        write(outfile, "n_pf", n_pf)
        write(outfile, "p_f_arr", p_f_arr)
        write(outfile, "p_b_max", p_b_max)
        write(outfile, "n_pb", n_pb)
        write(outfile, "p_b_arr", p_b_arr)
        write(outfile, "t_final", t_final)
        write(outfile, "n_t", n_t)
        write(outfile, "n_t_new", n_t_new)
        write(outfile, "t_mmt", t_mmt)
        write(outfile, "n_subdiv", n_subdiv)
        write(outfile, "subdiv_array", subdiv_array)
        write(outfile, "exp_index", exp_index)
        write(outfile, "TEE", TEE_array)
        write(outfile, "EE_cut", EE_cut_array)
    end
    return TEE_array
end

L = parse(Int, ARGS[1])
d = parse(Int, ARGS[2])
p_f = parse(Float64, ARGS[3])
n_pf = parse(Int, ARGS[4])
p_b = parse(Float64, ARGS[5])
n_pb = parse(Int, ARGS[6])
t_final = parse(Float64, ARGS[7])
t_final = Integer(t_final)
n_t = parse(Int, ARGS[8])
n_subdiv = parse(Int, ARGS[9])
exp_index = parse(Int, ARGS[10])
main(L, d, p_f, n_pf, p_b, n_pb, t_final, n_t, n_subdiv, exp_index);
