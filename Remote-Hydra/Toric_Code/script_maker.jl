# Define the file path
dir_path = @__DIR__
file_path = dir_path*"/test_script.sh"
jul_path = Sys.BINDIR

# Open the file in write mode
file = open(file_path, "w")

# Define the number of iterations for the for loop
num_iterations = 2
L_arr = [15, 18, 21, 24]
d = 3
p_f_max = 0.75
n_pf = 21
p_b_max = 0.0
n_pb = 1
t_final = 1e5
n_t = 20
number_of_repetitions = 100
number_of_experiments = number_of_repetitions*length(L_arr)

# Generate the lines for the .sh file
for i in 1:number_of_experiments
    L = L_arr[mod(i, length(L_arr))+1]
    line = "addqueue -c \"30 min\" -m 4 $(jul_path)/julia $(dir_path)/toric_code-TEE-FermiBoseCond.jl $(L) $(d) $(p_f_max) $(n_pf) $(p_b_max) $(n_pb) $(t_final) $(n_t) $(i)"
    println(file, line)
end

# Close the file
close(file)