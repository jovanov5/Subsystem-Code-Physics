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

# Define the file path
dir_path = @__DIR__
jul_path = Sys.BINDIR

# Define the number of iterations for the for loop
num_iterations = 2
L_arr = [15, 18, 21, 24]
d = 3

p_f_max = 1
n_pf = 21
p_f_arr = [x for x in range(0, p_f_max, length=n_pf)]

p_b_max = 1
n_pb = 21
p_b_arr = [x for x in range(0, p_b_max, length=n_pb)]

t_final = Int(2e4) # 1e5
n_t = 5 # 6
t_mmt = get_t_mmt_arr_refined(t_final, n_t)

n_subdiv = "all"

number_of_repetitions = 1
number_of_experiments = number_of_repetitions*length(L_arr)

# Generate the arg files
for exp_index in 1:number_of_experiments
    L = L_arr[mod(exp_index - 1, length(L_arr))+1]
    subdiv_array = get_subdiv_array(L, n_subdiv)
    arg_dict = Dict(
        "L" => L,
        "d" => d,
        "p_f_arr" => p_f_arr,
        "p_b_arr" => p_b_arr,
        "t_mmt" => t_mmt,
        "subdiv_array" => subdiv_array,
        "exp_index" => exp_index
    )
    open("$(dir_path)/data/args/argfile_exp:$(exp_index).json", "w") do io
        JSON3.pretty(io, arg_dict)
    end
end

