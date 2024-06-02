using Glob

# Define the file path
dir_path = @__DIR__
file_path = dir_path*"/queue_script-FermiBoseCond.sh"
jul_path = Sys.BINDIR

# Open the file in write mode
file = open(file_path, "w")

# Get all JSONs
Files = glob("*.json", "$(dir_path)/data/args")

memory = 4 # how many gigs
time = 1 # how many days

# Generate the lines for the .sh file
for arg_file in Files
    line = "addqueue -c \"$(time) day\" -m $(memory) $(jul_path)/julia $(dir_path)/toric_code-TEE-FermiBoseCond.jl \"$(arg_file)\""
    println(file, line)
end

# Close the file
close(file)