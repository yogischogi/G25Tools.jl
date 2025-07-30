# Find closest relatives belonging to a certain haplogroup.

using CSV, DataFrames, G25Tools

# Insert your own G25 coordinates here.
myG25 = ["Ragnar",0.13332656233600002,0.13488497177199998,0.0690654488256,0.055520039036,0.042185878582,0.021319462551439995,0.005257583968000024,0.007148006804000101,0.003945877867999992,-0.004060018999999859,-0.006072372758060001,0.0034969998680000475,-0.0064742415880000015,-0.003566617439999642,0.018543931918000134,0.0065606472260000825,-0.006926631224000224,0.0018001823780000037,0.0036236913339999996,0.004065752831000004,0.004839795716000295,0.0019055643019999513,0.0022728944959999817,0.012700345765999974,0.00008928743999980782]

# Insert your own haplogroups here.
Yhaplo = "R1b1a1b1a1a3"
mtDNAhaplo = "H1c1"

# Load main file because it contains the detailed information we need.
ancient_samples = DataFrame(CSV.File("ancient_dna - main.csv"))

# Take a subset of all samples that fit the Yhaplo group.
# The paternal haplogroups are listed in column "Y-dna final" 
# The maternal ones are in "mtdna final"
Yhaplo_samples = subset(ancient_samples, "Y-dna final" => ByRow(group -> !ismissing(group) ? occursin(Yhaplo, group) : false) )

# Extract G25 coordinates and convert the DataFrame of samples into a simpler format.
g25_Yhaplo = extractG25(Yhaplo_samples)

# Calculate the distances to the target sample and print them to the screen.
Y_distances = distances(g25_Yhaplo, myG25)
println(first(Y_distances, 10))


# The same for our maternal haplogroup.
mtDNA_samples = subset(ancient_samples, "mtdna final" => ByRow(group -> !ismissing(group) ? occursin(mtDNAhaplo, group) : false) )
g25_mtDNA = extractG25(mtDNA_samples)
mtDNA_distances = distances(g25_mtDNA, myG25)
println(first(mtDNA_distances, 10))












