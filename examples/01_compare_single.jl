# Compare the G25 coordinates of a single individual
# to a list of ancient samples.

using G25Tools

# Insert your own G25 here. They can be obtained from
# https://g25requests.app/
# For first experiments you may want to try simulated G25 coordinates.
# https://www.exploreyourdna.com/simulated-g25.aspx
# Do not forget to put your name in quotation marks.
myG25 = ["Ragnar",0.13332656233600002,0.13488497177199998,0.0690654488256,0.055520039036,0.042185878582,0.021319462551439995,0.005257583968000024,0.007148006804000101,0.003945877867999992,-0.004060018999999859,-0.006072372758060001,0.0034969998680000475,-0.0064742415880000015,-0.003566617439999642,0.018543931918000134,0.0065606472260000825,-0.006926631224000224,0.0018001823780000037,0.0036236913339999996,0.004065752831000004,0.004839795716000295,0.0019055643019999513,0.0022728944959999817,0.012700345765999974,0.00008928743999980782]

# Get the list of ancient samples from
# https://www.exploreyourdna.com/ancient-samples.aspx
# and store them in your working directory.
# You should now have a file named "ancient_dna - main.csv".
ancient_samples = extractG25("ancient_dna - main.csv")

# Look for your closest ancient relatives by calculating G25 distances.
relatives = distances(ancient_samples, myG25)

# Print top 10 ancient relatives to screen.
print(first(relatives, 10))

# Save all results to a file.
CSV.write("my_ancient_relatives.txt", relatives)


