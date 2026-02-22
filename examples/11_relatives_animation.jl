# Create an animation of your historical relatives depending on genetic distance.
# The animation is saved to a file named "historical_relatives.mp4".

using CSV, DataFrames
using CairoMakie, GeoMakie
using G25Tools

"""
Return the largest column index where column[index] <= distance.
"""
function index(column, distance)
    result = 0
    for (i, val) in enumerate(column)
        if val <= distance
            result = i
        else
            break
        end
    end
    return result
end

# Insert your own G25 here. They can be obtained from
# https://g25requests.app/
# Do not forget to put your name in quotation marks.
myG25 = ["Ragnar",0.13332656233600002,0.13488497177199998,0.0690654488256,0.055520039036,0.042185878582,0.021319462551439995,0.005257583968000024,0.007148006804000101,0.003945877867999992,-0.004060018999999859,-0.006072372758060001,0.0034969998680000475,-0.0064742415880000015,-0.003566617439999642,0.018543931918000134,0.0065606472260000825,-0.006926631224000224,0.0018001823780000037,0.0036236913339999996,0.004065752831000004,0.004839795716000295,0.0019055643019999513,0.0022728944959999817,0.012700345765999974,0.00008928743999980782]

# Ancient DNA to play with. G25 coordinates by Davidski.
#myG25 = ["China_Xinjiang_Abusanteer_IA_oWestEurasian.AG.SG:C4131.AG.SG__BC_667__",0.100164,0.083273,-0.019233,0.050065,-0.046778,0.027052,0.0094,0.003461,-0.036201,-0.042279,0.006496,-0.003447,0.000149,-0.012111,0.02158,0.006895,-0.025425,0.007601,0.003268,-0.011005,-0.007487,-0.00507,0.004314,-0.005061,0.002754]
#myG25 = ["Tanzania_Zanzibar_Euro_o:I0588__AD_800__",0.106994,0.152329,0.049026,-0.00323,0.047393,-0.008088,-0.00329,0.009461,0.008999,0.037176,-0.010555,0.007943,-0.016501,-0.018854,0.011536,0.000663,0.000391,0.002787,-0.006913,-0.003001,-0.010981,-0.002597,0.005176,-0.006025,0.00491]


# Main file for ancient DNA samples downloaded from exploreyourdna.com.
# Filename changes sometimes.
const ancient_dna_file = "ancient_dna - Sheet1.csv"
ancient_samples = DataFrame(CSV.File(ancient_dna_file))

# Simplify table of samples.
ancient_G25 = extractG25(ancient_samples; cols = ["Lat.", "Long."])

# Calculate distances for ancient relatives.
relatives = distances(ancient_G25, myG25)
result_table = leftjoin(relatives, ancient_G25, on = :Name)
sort!(result_table, :Distance)

# Create a figure that should be animated.
fig = Figure()
ga = GeoAxis(
    fig[1, 1];
    dest = "+proj=wintri", # the CRS in which you want to plot
    #limits = (-20, 40, 35, 70), # Europe
    limits = (-20, 120, 20, 70), # Europe and Asia
    #limits = (-180, 179, -60, 80), # Earth
    subtitle = "Data provided by nomad"
)

# Add earth map.
img = rotr90(GeoMakie.earth())
image!(ga, -180..180, -90..90, img; interpolate = false) 
lines!(ga, GeoMakie.coastlines())

# Parameters for animation.
start_distance = 0.01
finish_distance = 0.1
step_distance = 0.001
iterator = start_distance:step_distance:finish_distance

# Create an animation by modifying the figure parameters frame by frame.
record(fig, "historical_relatives.mp4", iterator; framerate = 2) do step
    println(step)
    min = index(result_table.Distance, step - step_distance)
    if min == 0
        min = 1
    end
    max = index(result_table.Distance, step)
    samples = result_table[min:max, :]
    ga.title = "Relatives, distance: " * string(step)
    scatter!(ga, samples[!, "Long."], samples[!, "Lat."]; color = :red)
end


