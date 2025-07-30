# Find your closest ancient relatives in different time periods.

using G25Tools

# Insert your own G25 coordinates here.
myG25 = ["Ragnar",0.13332656233600002,0.13488497177199998,0.0690654488256,0.055520039036,0.042185878582,0.021319462551439995,0.005257583968000024,0.007148006804000101,0.003945877867999992,-0.004060018999999859,-0.006072372758060001,0.0034969998680000475,-0.0064742415880000015,-0.003566617439999642,0.018543931918000134,0.0065606472260000825,-0.006926631224000224,0.0018001823780000037,0.0036236913339999996,0.004065752831000004,0.004839795716000295,0.0019055643019999513,0.0022728944959999817,0.012700345765999974,0.00008928743999980782]

# Read ancient samples in G25 format (see example No. 02).
sourceG25 = readG25("G25.txt")

# Define parameters for a set of time periods.
# Here we divide the last 5000 years by intervals of 500 years.
start_year = -3000
interval_length = 500
final_year = 2000
intervals = (final_year - start_year) / interval_length

# Return the interval number according to the sample age.
interval(sample_age) = floor(Integer, (sample_age - start_year) / interval_length) + 1

# Create an array that holds a DataFrame of samples for each time period.
relatives_in_time = [similar(sourceG25, 0) for _ in 1:intervals]

# Fill intervals with ancient samples.
for sample in eachrow(sourceG25)
    year = getyear(sample.Name)
    if !ismissing(year)
        i = interval(year)
        if i > 0 && i <= intervals
            push!(relatives_in_time[i], sample)
        end
    end
end

# Calculate genetic distances for each table and print them to screen.
for (i, table) in enumerate(relatives_in_time)
    dist = distances(table, myG25)
    start = start_year + (i - 1) * interval_length
    final = start + interval_length
    println()
    println("Closest matches from year: $start to: $final")
    println(first(dist, 10))
end





