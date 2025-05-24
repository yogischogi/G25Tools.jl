# Compare two lists of DNA samples with G25 coordinates.

using G25Tools

# Convert main file of ancient samples into simplified Davidski format.
# The main file can be found at:
# https://www.exploreyourdna.com/ancient-samples.aspx
ancient_samples = extractG25("ancient_dna - main.csv", "G25.txt")

# Adjust file names to your needs.
# You need to create a targetfiles with a list of samples in G25 format.
# See the file G25.txt as an example.
sourcefile = "G25.txt"
targetfile = "targetsamples.txt"
resultsdirectory = "temp"

# Load samples with G25 coordinates from files.
sourceG25 = readG25(sourcefile)
targetG25 = readG25(targetfile)

# Calculate G25 distances for each entry in targetG25.
results = distances(sourceG25, targetG25)

# Write distances for each entry to file.
for r in results
    writedistances(resultsdirectory, r.targetname, r.sourcesamples)
end



