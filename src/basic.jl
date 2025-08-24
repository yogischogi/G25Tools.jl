"""
    function extractG25(from::AbstractDataFrame; cols = String[])

Extract G25 data from a DataFrame that is in the form
of the main file of ancient samples.
The main file can be found at:
https://www.exploreyourdna.com/ancient-samples.aspx

Return a DataFrame in a simplified format that represents samples
as Names and their G25 coordinates.

The parameter cols can be used to add extra columns from the main file,
for example it may be useful to add GPS coordinates. Just add the
column titles as a String array like cols = ["Lat.", "Long."].
"""
function extractG25(from::AbstractDataFrame; cols = String[])
    # Include only lines which contain G25 coordinates
    source = subset(from, "g25" => ByRow(g25 -> !ismissing(g25)))

    # Convert g25 column to DataFrame.
    g25table = DataFrame(Name = String[] , PC1 = Float64[], PC2 = Float64[],
        PC3 = Float64[], PC4 = Float64[], PC5 = Float64[], PC6 = Float64[],
        PC7 = Float64[], PC8 = Float64[], PC9 = Float64[], PC10 = Float64[],
        PC11 = Float64[], PC12 = Float64[], PC13 = Float64[], PC14 = Float64[],
        PC15 = Float64[], PC16 = Float64[], PC17 = Float64[], PC18 = Float64[],
        PC19 = Float64[], PC20 = Float64[], PC21 = Float64[], PC22 = Float64[],
        PC23 = Float64[], PC24 = Float64[], PC25 = Float64[]
    )
    date = source[!, "mean ad/bc"]
    for (i, line) in enumerate(source.g25)
        datestring = ""
        if !ismissing(date[i])
            try
                year = date[i]
                if typeof(year) != Float64
                    year = parse(Float64, year)
                end
                d = round(Integer, year)
                if d < 0
                    d = -d
                    datestring = "__BC_$(d)__"
                else
                    datestring = "__AD_$(d)__"
                end
            catch e
                # Happens a lot when the date is #n/a or ? or ...
                # println("Problem extracting date: $e")
            end
        end
        csvrow = split(line, ",")
        csvrow[1] *= datestring
        row = []
        push!(row, csvrow[1])
        for i in 2:26
            push!(row, parse(Float64, csvrow[i]))
        end
        push!(g25table, row)
    end

    # Add extra columns.
    for title in cols
        g25table[!, title] = source[!, title]
    end

    return g25table
end

"""
    extractG25(filename::AbstractString; cols = String[])

Extract G25 data from the main file of ancient samples.
The main file can be found at:
https://www.exploreyourdna.com/ancient-samples.aspx

Return a DataFrame in a simplified format that represents samples
as Names and their G25 coordinates.

The parameter cols can be used to add extra columns from the main file,
for example it may be useful to add GPS coordinates. Just add the
column titles as a String array.
"""
function extractG25(filename::AbstractString; cols = String[])
    fulldata = DataFrame(CSV.File(filename))
    return extractG25(fulldata; cols)
end

"""
    extractG25(filename::AbstractString, to::AbstractString; cols = String[])

Extract G25 data from the main file and save it locally to disk.
The main file can be found at:
https://www.exploreyourdna.com/ancient-samples.aspx

The G25 data is a simplified version of the main file which ir rather
large to work with. Also some population data is only available in the
simplified format.

The parameter cols can be used to add extra columns from the main file,
for example it may be useful to add GPS coordinates. Just add the
column titles as a String array.
"""
function extractG25(filename::AbstractString, to::AbstractString; cols = String[])
    data = extractG25(filename; cols)
    CSV.write(to, data)
end

"""
    readG25(filename)

Read a DataFrame containing a list of samples with G25 coordinates.
Each row must begin with the sample name followed by their G25 coordinates.
Example: Name,PC1,PC2,...PC25
"""
function readG25(filename)
    colnames = ["Name", "PC1", "PC2", "PC3", "PC4", "PC5",
        "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12",
        "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19",
        "PC20", "PC21", "PC22", "PC23", "PC24", "PC25"
    ]
    samples = DataFrame(CSV.File(filename))
    if names(samples)[1:26] != colnames[1:26]
        error(
        """
        Invalid file format!
        Please make sure that the first line of the file begins with
        Name,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,PC11,PC12,PC13,PC14,PC15,PC16,PC17,PC18,PC19,PC20,PC21,PC22,PC23,PC24,PC25
        """)
    end
    return samples
end

"""
    write(directory, filename, distances)

Write a DataFrame in G25 format to a specific directory.
This is a convenience function that is most useful when
writing multile files to a single directory.
"""
function writedistances(directory, filename, sourcedistances::DataFrame)
    filepath = joinpath(directory, filename)
    CSV.write(filepath, sourcedistances)
end

"""
    function picksamples(samples::AbstractDataFrame; populations = String[], timeperiod = [-Inf, Inf])

Pick a table of samples belonging to certain populations and/or a time persiod.
Example: picksamples(samples, populations=["Germany", "England"], timeperiod=[0, 500])

For more advanced filtering options use DataFrames' subset() or filter() methods.
"""
function picksamples(samples::AbstractDataFrame; populations = String[], timeperiod = [-Inf, Inf])
    result = similar(samples, 0)

    # Check if sample belongs to one of the selected populations.
    for sample in eachrow(samples)
        ok = true

        if populations != []
            for p in populations
                if occursin(p, sample.Name)
                    ok = true
                    break
                end
                ok = false
            end
        end

        # Check time period.
        if ok == true
            year = getyear(sample.Name)
            if ismissing(year) || year < timeperiod[1] || year > timeperiod[2]
                ok = false
            end
        end

        if ok == true
            push!(result, sample)
        end
    end
    return result
end

"""
    getyear(samplename)

Gets the year from a sample name in Davidskis format.
Davidski encodes the year in the sample name with a pattern
like __BC_50__ or __AD_450__.

Returns the year or missing if the year is not encoded.
"""
function getyear(samplename)
    # Look for year patterns.
    r = r"__(BC|AD)_\d{1,6}__"
    m = match(r, samplename)

    if isnothing(m)
        return missing
    end
    digits = m.match[6: end-2]

    # Convert digits to time.
    year = parse(Int64, digits)

    if m[1] == "BC"
        year = -year
    end
    return year
end

"""
    distance(source::Vector, target::Vector)

Return the vector distance between two samples with G25 coordinates.
Both source and target must be in the form:
["Name", PC1,...PC25]
where Name is a String and PC1 to PC25 are the G25 coordinates
of the sample. Additional fields are allowed because only [2:26] are
used to calculate the distance.
"""
function distance(source::Vector, target::Vector)
    s = source[2:26]
    t = target[2:26]
    return norm(t - s)
end

"""
    distance(source::DataFrameRow, target::DataFrameRow)

Return the vector distance between two samples with G25 coordinates.
Both source and target must be in the form:
["Name", PC1,...PC25]
where Name is a String and PC1 to PC25 are the G25 coordinates
of the sample. Additional fields are allowed because only [2:26] are
used to calculate the distance.
"""
function distance(source::DataFrameRow, target::DataFrameRow)
    s = collect(source[1:26])
    t = collect(target[1:26])
    return distance(s, t)
end

"""
    distances(source::DataFrame, target::Vector)

Calculate the distances between samples in a DataFrame and
a single target sample.
"""
function distances(source::DataFrame, target::Vector)
    # Calculate vector distance for each sample in source.
    dist = Float64[]
    for row in eachrow(source)
        s = collect(row[1:26])
        push!(dist, distance(s, target))
    end

    # Result is a DataFrame containing samples and distances.
    result = DataFrame(Name = source.Name, Distance = dist)
    sort!(result, :Distance)
    return result
end

"""
    distances(source::DataFrame, target::DataFrameRow)

Calculate the G25 vector distances from a list of samples
(source) to a single target.

Return a DataFrame containing the distances to the target.
Each row consists of the source name and the distance to the
target. The table is ordered, smallest distances come first.
"""
function distances(source::DataFrame, target::DataFrameRow)
    t = collect(target[1:26])
    return distances(source, t)
end

"""
    distances(source::DataFrame, target::DataFrame)

Calculate the G25 vector distances from a list of samples
(source) to multiple target samples.

The source and the target contain a DataFrame in G25 format.

Return a vector of Distances where each entry contains the
Distances of a single target sample to all source samples.
"""
function distances(source::DataFrame, target::DataFrame)
    result = Distances[]
    for t in eachrow(target)
        push!(result, Distances(t.Name, distances(source, t)))
    end
    return result
end

"""
    average(samples::DataFrame, name = "Average")

Calculate the mean average of a DataFrame of samples in G25 format.

Return the average as a vector of G25 coordinates where the first
entry is the given name, usually the name of the given population.
"""
function average(samples::DataFrame; name = "Average")
    if isempty(samples)
        return nothing
    end
    result = []
    push!(result, name)
    for col in 2:26
        push!(result, mean(samples[!, col]))
    end
    return result
end

"""
    median(samples::DataFrame; name = "Median")

Calculate the median average of a DataFrame of samples in G25 format.

Return the median average as a vector of G25 coordinates where the first
entry is the given name, usually the name of the given population.

This method may be better suited than the mean aveerage for small
sample sizes that include outliers.
"""
function medianavg(samples::DataFrame; name = "Median")
    if isempty(samples)
        return nothing
    end
    result = []
    push!(result, name)
    for col in 2:26
        push!(result, median(samples[!, col]))
    end
    return result
end

"""
    topmatches(sourcesamples, targetsamples;
               sourcetimeperiod = [-Inf, Inf], targettimeperiod = [-Inf, Inf],
               maxdistance = 0.1, targetfraction = 0.68)

Find sourcesamples that had a big impact on the targetsamples.
The intention is to find the target's most important ancestral
population.

The method looks for matches within the given G25 distance (maxdistance)
to the target samples. It selects those matches that match the specified
number of targetsamples (targetfraction) or more.
"""
function topmatches(sourcesamples, targetsamples;
                    sourcetimeperiod = [-Inf, Inf], targettimeperiod = [-Inf, Inf],
                    maxdistance = 0.1, targetfraction = 0.66)
    # Filter sourcesamples for time period.
    source = similar(sourcesamples, 0)
    for s in eachrow(sourcesamples)
        year = getyear(s.Name)
        if !ismissing(year) && year >= sourcetimeperiod[1] && year <= sourcetimeperiod[2]
            push!(source, s)
        end
    end

    # Filter targetsamples for time period.
    target = similar(targetsamples, 0)
    for s in eachrow(targetsamples)
        year = getyear(s.Name)
        if !ismissing(year) && year >= targettimeperiod[1] && year <= targettimeperiod[2]
            push!(target, s)
        end
    end

    # Find all matches.
    matches = distances(source, target)

    # Consider only matches that are close to a member of the target population.
    canditates = DataFrame(Name = String[], Distance = Float64[])
    for m in matches
        for s in eachrow(m.sourcesamples)
            if s.Distance <= maxdistance
                push!(canditates, s)
            else
                # We can break here because the samples are ordered by distance.
                break
            end
        end
    end

    # Count the number of matches for each sourcesample.
    matchcounts = Dict{String, Integer}()
    for m in eachrow(canditates)
        if haskey(matchcounts, m.Name)
            matchcounts[m.Name] +=  1
        else
            matchcounts[m.Name] = 1
        end
    end

    # Sort sample names according to highest count.
    topmatches = DataFrame(Name = String[], Count = Integer[])
    mincount = targetfraction * nrow(target)
    for name in keys(matchcounts)
        if matchcounts[name] >= mincount
            push!(topmatches, (name, matchcounts[name]))
        end
    end
    sort!(topmatches, :Count, rev = true)

    # Add columns from the source table.
    topmatches = leftjoin(topmatches, source, on = :Name, makeunique = true)
    
    return topmatches
end

"""
    clusters(samples::DataFrame; distance = 0.1, neighbors = 1)

Find clusters in a set of samples that satisfy the given conditions.
Within each cluster the nearest neighbor of a sample must be within
distance and each sample must have a minumum number of neighbors that
is given by neighbors.

Return an array of Populations where each Population represents one cluster.
If no cluster could be found an empty array is returned.

XXX Currently each sample may appear in only one cluster, which
can lead to strange results. Also samples at the fringe of a cluster
are sometimes cut off if they do not have enough neighbors.
Performance is rather poor for large sample sizes. About 6 seconds
for 500 samples on my rather old computer.
24 seconds for 600 samples.
"""
function clusters(samples::DataFrame; distance = 0.05, neighbors = 1)
    result = Population[]
    n = nrow(samples)
    if n == 0 
        return result
    end

    # Calculate distances between all samples.
    dists = zeros(n, n)
    for col = 1:n
        for row = col+1:n
            dists[col, row] = G25Tools.distance(samples[col, :], samples[row, :])
            dists[row, col] = dists[col, row]
        end
    end

    # Get indices of all samples that are within a cluster.
    idxs = _select(dists, distance, neighbors)

    # Calculate clusters, each element may appear in only one cluster.
    while length(idxs) > 1
        cluster = Integer[]
        push!(cluster, idxs[1])
        push!(cluster, _cluster(idxs[1], idxs[2:end], dists, distance)...)
        if length(cluster) > 0
            # Add the cluster to the list of populations.            
            s = similar(samples, 0)
            for i in cluster
                push!(s, samples[i, 1:end])
            end
            push!(result, Population("Cluster", s))
            # Remove elements of the first cluster from the list of indices.
            idxs = setdiff(idxs, cluster)
        end
    end

    return result
end

"""
    _select(samples::DataFrame, distance, neighbors)

Select those samples that have at least the given number of
neighbors within the given distance.
The distancematrix consists of sample indices and their distances
to each other.

Return an array sample indices that satisfy the given conditions.
"""
function _select(distancematrix, distance, neighbors)
    result = Integer[]
    (n, m ) = size(distancematrix)
    for col = 1:m
        hits = 0
        for row = 1:n
            if distancematrix[row, col] <= distance
                hits += 1
            end
        end
        if hits > neighbors
            push!(result, col)
        end
    end
    return result
end

"""
    _cluster(index::Integer, indices::Integer[], distancematrix, distance)

Return all indices that cluster with the index index.    
"""
function _cluster(index::Integer, indices::Array{Integer}, distancematrix, distance)
    result = Integer[]
    n = length(indices)
    if n < 1
        return result
    end

    # Find all elements that are close to index.
    neighbors = Integer[]
    nonneighbors = Integer[]
    for i in indices
        if distancematrix[index, i] <= distance
            push!(neighbors, i)
        else
            push!(nonneighbors, i)
        end
    end
    result = union(result, neighbors)

    # Call _cluster to find all the elements that are close to the privious neighbors.
    for i in neighbors
        next_neighbors = _cluster(i, nonneighbors, distancematrix, distance)
        result = union(result, neighbors)
    end

    return result
end


