# Compare two lists of DNA samples with G25 coordinates.

using G25Tools

# Adjust file names to your needs.
# You need to create a targetfiles with a list of samples in G25 format.
# See the file G25.txt as an example. This is the G25 file we created
# in example 01.
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



