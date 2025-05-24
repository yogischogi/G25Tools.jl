# Make a PCA plot of different populations.

using G25

# Load samples with G25 coordinates from file.
ancient_samples = readG25("G25.txt")

# Select populations from ancient samples.
celtic = Population("Celtic", pick(ancient_samples, populations = ["Hallstatt", "LaTene"]))
germanic = Population("Germanic", pick(ancient_samples, populations = ["Germanic"]))

# Display PCA plot.
f = plotpca([celtic, germanic])
f

# If you like to compare your own G25 coordinates to ancient
# populations save them to a file in G25 format and load them
# from the file, for example:
# myG25 = Population("myG25", readG25("myG25.txt"))
# f = plotpca([celtic, germanic, myG25], smile = true)
# f

