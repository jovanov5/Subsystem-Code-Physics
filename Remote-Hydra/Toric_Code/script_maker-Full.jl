using Glob

# Define the file path
dir_path = @__DIR__
file_path = dir_path*"/queue_script-Full.sh"
jul_path = Sys.BINDIR

# Open the file in write mode
file = open(file_path, "w")

# Get all JSONs
Files = glob("*.json", "$(dir_path)/data/args")

memory = 6 # how many gigs
time = 3 # how many days
comment = "full-PD sweep - 1 days"
cores = 1 # -n 4x4 -s : need to figure it out!!!!

# Generate the lines for the .sh file
for arg_file in Files
    line = "addqueue -c \"$(comment)\" -m $(memory) $(jul_path)/julia $(dir_path)/toric_code-Full.jl \"$(arg_file)\""
    println(file, line)
end

# Close the file
close(file)