using Base.Filesystem  # for creaating directories etc
using JSON3 # for reading the JSON file
using ITensors
using LinearAlgebra
using SparseArrays
using IterTools
using Plots

dep_path = "../../AdditionalStructure/Julia/"
include("$(dep_path)SpinHalfModels.jl");

# Define the file path
dir_path = @__DIR__
jul_path = Sys.BINDIR

# Simulation Params to pack and save!
N_arr = [200]
cutoff_arr = [1e-10]
maxdim_arr = [64]
dt_sim_arr = [0.05]
N_A_arr = [2]
N_C_arr = [2]
Js_arr = [[1, 0, 0]]
hs_arr = [[0, 1, 0]]
all_param_array = IterTools.product(N_arr, cutoff_arr, maxdim_arr, dt_sim_arr, N_A_arr, N_C_arr, Js_arr, hs_arr)
sample_times = [LinRange(0.05, 10, 200)] # Not varied
array_of_d = [0:1:20] # Not varied

# Generate the arg files
for (exp_index, (N, cutoff, maxdim, dt_sim, N_A, N_C, Js, hs)) in enumerate(all_param_array)
    arg_dict = Dict(
        "N" => N,
        "cutoff" => cutoff,
        "maxdim" => maxdim,
        "dt_sim" => dt_sim,
        "sample_times" => sample_times,
        "N_A" => N_A,
        "N_C" => N_C,
        "array_of_d" => array_of_d,
        "Js" => Js,
        "hs" => hs,
        "exp_index" => exp_index,
        "debug" => 0
    )
    open("$(dir_path)/data/args/argfile_exp:$(exp_index).json", "w") do io
        JSON3.pretty(io, arg_dict)
    end
end

