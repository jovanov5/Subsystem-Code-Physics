# Define the file path
dir_path = @__DIR__
file_path = dir_path*"/queue_script.sh"
jul_path = Sys.BINDIR

# Open the file in write mode
file = open(file_path, "w")

# Define the number of iterations for the for loop
num_iterations = 2
L_arr = [15, 18, 21, 24]
d = 3
p_f_max = 1
n_pf = 21
p_b_max = 1
n_pb = 21
t_final = 1e5
n_t = 6
n_subdiv = 12
number_of_repetitions = 100
number_of_experiments = number_of_repetitions*length(L_arr)

# Generate the lines for the .sh file
for exp_index in 1:number_of_experiments
    L = L_arr[mod(exp_index, length(L_arr))+1]
    line = "addqueue -c \"1 day\" -m 4 $(jul_path)/julia $(dir_path)/toric_code-TEE-FermiBoseCond.jl $(L) $(d) $(p_f_max) $(n_pf) $(p_b_max) $(n_pb) $(t_final) $(n_t) $(n_subdiv) $(exp_index)"
    println(file, line)
end

# Close the file
close(file)