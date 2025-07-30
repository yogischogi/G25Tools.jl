# Display a map of ancient samples.

using CSV, DataFrames, G25Tools
using CairoMakie, GeoMakie

# Define some haplogroups.
const U106 = "R1b1a1b1a1a1"
const P312 = "R1b1a1b1a1a2"
const S1194 = "R1b1a1b1a1a3"

# Load main file because it contains the detailed information we need.
ancient_samples = DataFrame(CSV.File("ancient_dna - main.csv"))

# Select all samples belonging to a specific haplogroup.
Yhaplo = P312
Yhaplo_samples = subset(ancient_samples, "Y-dna final" => ByRow(group -> !ismissing(group) ? occursin(Yhaplo, group) : false) )

# Simplify the table so that it is easier to work with.
# This is not needed but I want to demonstrate G25Tools here.
samples = extractG25(Yhaplo_samples, cols = ["Lat.", "Long."])

# Show samples on map.
fig = Figure()
ga = GeoAxis(
    fig[1, 1];
    dest = "+proj=wintri", 
    limits = (-20, 40, 35, 70)
)
# Draw map.
img = rotr90(GeoMakie.earth())
image!(ga, -180..180, -90..90, img; interpolate = false) 
lines!(ga, GeoMakie.coastlines())

scatter!(ga, samples[!, "Long."], samples[!, "Lat."]; color = :red)
fig


