# Define the file path
dir_path = @__DIR__
file_path = dir_path*"/test_script.sh"

# Open the file in write mode
file = open(file_path, "w")

# Define the number of iterations for the for loop
num_iterations = 2
L = 17
d = 2
p_f_max = 0.0
n_pf = 1
p_b_max = 0.2
n_pb = 2
t_final = 1e5
n_t = 10
exp_index = 1

# Generate the lines for the .sh file
for i in 1:num_iterations
    line = "addqueue -c \"2 min\" -m 2 julia $(dir_path)/toric_code-TEE-FermiBoseCond.jl $(L) $(d) $(p_f_max) $(n_pf) $(p_b_max) $(n_pb) $(t_final) $(n_t) $(exp_index)"
    println(file, line)
end

# Close the file
close(file)