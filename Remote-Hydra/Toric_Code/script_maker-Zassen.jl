using Glob

# Define the file path
dir_path = @__DIR__
file_path = dir_path*"/queue_script-Zassenhausen.sh"
jul_path = Sys.BINDIR

# Open the file in write mode
file = open(file_path, "w")

# Get all JSONs
Files = glob("*.json", "$(dir_path)/data/args")

# Generate the lines for the .sh file
for arg_file in Files
    line = "addqueue -c \"1 day\" -m 4 $(jul_path)/julia $(dir_path)/toric_code-Zassenhausen.jl \"$(arg_file)\""
    println(file, line)
end

# Close the file
close(file)