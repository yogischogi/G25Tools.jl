# Create an animation of ancient DNA samples.
# Save animation as "R1b_migration.mp4".

using CSV, DataFrames
using CairoMakie, GeoMakie
using G25Tools  # To install: add https://github.com/yogischogi/G25Tools.jl

# Define some haplogroups.
const R1b = "R1b"
const L151 =  "R1b1a1b1a1a"
const U106 =  "R1b1a1b1a1a1"
const P312 =  "R1b1a1b1a1a2"
const S1194 = "R1b1a1b1a1a3"
const Z2103 = "R1b1a1b1b"
const NeoFarmers = "G2"

# Main file for ancient DNA samples downloaded from exploreyourdna.com.
# Filename changes sometimes.
const ancient_dna_file = "ancient_dna - Sheet1.csv"
ancient_samples = DataFrame(CSV.File(ancient_dna_file))

# Filter samples for haplogroup.
Yhaplo = R1b
Yhaplo_samples = subset(ancient_samples, "Isogg final" => ByRow(group -> !ismissing(group) ? occursin(Yhaplo, group) : false) )
R1bsamples = extractG25(Yhaplo_samples, cols = ["Lat.", "Long."])

# Create a figure that should be animated.
fig = Figure()
ga = GeoAxis(
    fig[1, 1];
    dest = "+proj=wintri", # the CRS in which you want to plot
    limits = (-20, 40, 35, 70), # Europe
    #limits = (-20, 120, 20, 70), # Europe and Asia
    subtitle = "Data provided by nomad"
)

# Add earth map.
img = rotr90(GeoMakie.earth())
image!(ga, -180..180, -90..90, img; interpolate = false) 
lines!(ga, GeoMakie.coastlines())

# Parameters for animation.
start_time = -10000
finish_time = 1500
iterator = start_time:100:finish_time

# Create an animation by modifying the figure parameters frame by frame.
record(fig, "R1b_migration.mp4", iterator; framerate = 2) do step
    println(step)
    migration_time = step
    samples = picksamples(R1bsamples, timeperiod = [-Inf, migration_time])
    ga.title = "R1b samples, year: " * string(step)
    scatter!(ga, samples[!, "Long."], samples[!, "Lat."]; color = :red)
end

