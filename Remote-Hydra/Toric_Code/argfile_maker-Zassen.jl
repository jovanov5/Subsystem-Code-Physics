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

# Define the file path
dir_path = @__DIR__
jul_path = Sys.BINDIR

# Define the number of iterations for the for loop
# L_arr = [30, 42, 54]
# L_arr = [24, 30]
L_arr = [31]
d = 1

# p_all_arr = [(0.55, 0.0), (0.55, 0.1), (0.2, 0.21), (0.0, 0.15)] # Defined as list of doubles (p_f, p_b)
p_all_arr = [(0.75, 0.00), (0.75, 0.10)]

# t_final = Int(2e4) # 1e5
# n_t = 5 # 6
t_mmt = Array{Int}(20:20:40).+50

n_subdiv = "all"

number_of_repetitions = 200
number_of_experiments = number_of_repetitions*length(L_arr)*length(p_all_arr)
offset = 0

# Generate the arg files
for exp_index in (1 + offset):(number_of_experiments + offset)
    usable_exp_index = exp_index - 1
    usable_L_index = mod(div(usable_exp_index, length(p_all_arr)), length(L_arr))
    usable_p_index = mod(usable_exp_index, length(p_all_arr))
    L_index = usable_L_index + 1
    p_index = usable_p_index + 1

    L = L_arr[L_index]
    p_f, p_b = p_all_arr[p_index]

    subdiv_array = get_subdiv_array(L, n_subdiv)
    arg_dict = Dict(
        "L" => L,
        "d" => d,
        "p_f" => p_f,
        "p_b" => p_b,
        "t_mmt" => t_mmt,
        "exp_index" => exp_index
    )
    open("$(dir_path)/data_z/args/argfile_exp:$(exp_index).json", "w") do io
        JSON3.pretty(io, arg_dict)
    end
end

