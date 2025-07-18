# Compare a sample to different populations using a PCA plot.
# We start with the same code as in example No. 3.

using G25

# Insert your own G25 coordinates here.
myG25 = ["Ragnar",0.13332656233600002,0.13488497177199998,0.0690654488256,0.055520039036,0.042185878582,0.021319462551439995,0.005257583968000024,0.007148006804000101,0.003945877867999992,-0.004060018999999859,-0.006072372758060001,0.0034969998680000475,-0.0064742415880000015,-0.003566617439999642,0.018543931918000134,0.0065606472260000825,-0.006926631224000224,0.0018001823780000037,0.0036236913339999996,0.004065752831000004,0.004839795716000295,0.0019055643019999513,0.0022728944959999817,0.012700345765999974,0.00008928743999980782]

# Load samples with G25 coordinates from file.
ancient_samples = readG25("G25.txt")

# Create an empty DataFrame that we will fill with populations.
ancient_populations = similar(ancient_samples, 0)

# Create a few ancient populations by filtering the samples
# according to country and take the average.
countries = ["England", "Denmark", "Sweden", "Norway", "Germany", "Italy", "Spain"]
for country in countries
    country_samples = picksamples(ancient_samples, populations = [country])
    population = average(country_samples, name = country)
    push!(ancient_populations, population)
end

# Calculate the population distances to the target sample.
closest_populations = distances(ancient_populations, myG25)
print(closest_populations)

# Add sample to the populations and plot it.
push!(ancient_populations, myG25)
f = plotpca(ancient_populations, smile = true)
f


