# Find ancestral matches to a given set of samples.

using CSV, DataFrames
using CairoMakie, GeoMakie
using G25Tools

# Define some haplogroups.
const L151 =  "R1b1a1b1a1a"
const U106 =  "R1b1a1b1a1a1"
const P312 =  "R1b1a1b1a1a2"
const S1194 = "R1b1a1b1a1a3"

# Display samples on a map.
function plotsamples(samples)
    fig = Figure()
    ga = GeoAxis(
        fig[1, 1]; # any cell of the figure's layout
        dest = "+proj=wintri", # the CRS in which you want to plot
        limits = (-20, 40, 35, 70)
    )
    # Draw map.
    img = rotr90(GeoMakie.earth())
    image!(ga, -180..180, -90..90, img; interpolate = false) 
    lines!(ga, GeoMakie.coastlines())

    scatter!(ga, samples[!, "Long."], samples[!, "Lat."]; color = :red)
    return fig
end

# Get samples in simplified G25 format.
ancient_samples = DataFrame(CSV.File("ancient_dna - main.csv"))
samples = extractG25(ancient_samples,  cols = ["Lat.", "Long."])

Yhaplo = L151
Yhaplo_samples = subset(ancient_samples, "Y-dna final" => ByRow(group -> !ismissing(group) ? occursin(Yhaplo, group) : false) )
L151samples = extractG25(Yhaplo_samples, cols = ["Lat.", "Long."])

# Define time periods.
sourcetime =  [-3500, -2910]
targettime =  [-2900, -2700]

# Display target samples.
target = picksamples(L151samples, timeperiod = targettime)
plotsamples(target)

# Display ancestral samples.
top = topmatches(samples, L151samples;
                 sourcetimeperiod = sourcetime, targettimeperiod = targettime,
                 maxdistance = 0.08, targetfraction = 0.68)
plotsamples(top)


