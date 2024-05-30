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

function main(L::Integer, d::Integer, p_f::Float64, p_b::Float64, p_fe::Float64, t_mmt::Array{Int}, exp_index::Integer, debug::Int)
    """ This is the main function on terminal call.
        The simulation is of Toric Code and we try and measure the Boson and Fermion with some probabilities.
        The main measure is the Zassenhausen correlator of boson-boson and fermion-fermion.
        d is not used because TEE is not being measured, included for compatibility.
    """

    description = "Single Points in the Phase Diagram, TEE and EE cut."
    system = Init_EdgeSquareLattice_KitaevDoNuT(L, d);
    sys_type = "Init_EdgeSquareLattice_KitaevDoNuT";
    dirpath = @__DIR__;
    filename = @__FILE__;
    n_t = length(t_mmt);
    subdiv_array = get_subdiv_array(L, "all");
    n_subdiv = length(subdiv_array);

    EE_cut_array  = zeros(Int, n_t, n_subdiv);
    TEE_array = zeros(Int, n_t);
    
    println("filename: ", filename);
    println("description: ", description);
    println("L: ", L);
    println("d: ", d);
    println("sys_type: ", sys_type);
    println("p_f: ", p_f);
    println("p_b: ", p_b);
    println("p_fe: ", p_fe)
    println("n_t: ", n_t);
    # println("t_mmt: ", t_mmt);
    println("subdiv_array: ", subdiv_array);
    println("exp_index: ", exp_index);

    if debug == 1

        # Debug: Plot EE vs cut for the last time
        p = plot(subdiv_array, EE_cut_array[end, :], xlabel="Cut", ylabel="EE", marker=:circle)
        savefig(p, dirpath*"/data/debug-out/test_plot_1.pdf")

        save_data_prefix = "debug-out/"
    end

    p_tc = 1 - p_f - p_b - p_fe; # Defines the measuremnt dynamics!
    stab_distro = Categorical([p_tc/2, p_tc/2, p_b, 0, p_f, p_fe]);
    state = toric_code_GS(system); # Get the pure TC ground state as the initial state

    t_old = 0;
    for t_index = 1:n_t
        t_evol = t_mmt[t_index] - t_old;
        t_old = t_mmt[t_index];
        state = iterate_measurements_only_fast!(state, system, () -> toric_code(system, stab_distro), t_evol)
        TEE_array[t_index] = entanglement_entropy_topo(state, system)
        EE_cut_array[t_index, :] = entanglement_entropy_cut(state, system, subdiv_array)
    end

    save_data_prefix = "out/"
    if debug == 1

        # Debug: Plot EE vs cut for the last time
        p = plot(subdiv_array, EE_cut_array[end, :], xlabel="Cut", ylabel="EE", marker=:circle)
        savefig(p, dirpath*"/data/debug-out/test_plot_1.pdf")

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
        write(outfile, "p_f", p_f)
        write(outfile, "p_b", p_b)
        write(outfile, "p_fe", p_fe)
        write(outfile, "t_mmt", t_mmt)
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
p_f = Float64(arg_dict["p_f"])
p_b = Float64(arg_dict["p_b"])
p_fe = Float64(arg_dict["p_fe"]) # Stabiliser for the fermionic end-points!
t_mmt = Array{Int}(arg_dict["t_mmt"])
exp_index = Int(arg_dict["exp_index"])

if length(ARGS) > 1
    debug = parse(Int, ARGS[2])
else
    debug = 0
end

main(L, d, p_f, p_b, p_fe, t_mmt, exp_index, debug);
