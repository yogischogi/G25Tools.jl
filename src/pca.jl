"""
    Population

    Represents a population. Each population consists of a name,
    for example "Viking", and a list of samples in G25 format.
"""
struct Population
    name::String
    samples::DataFrame
end

"""
    plotpca(samples::DataFrame)

Plot all samples in a given DataFrame. The DataFrame must
be in the G25 format: Name followed by PC1,...PC25.
If smile == true the last marker is a smiley.

This method generates a label for each sample.
"""
function plotpca(samples::DataFrame; title = "PCA Plot", smile = false)
    f = Figure()
    ax = Axis(f[1, 1],
        title = title,
        xlabel = "PC1",
        ylabel = "PC2",
    )
    m = _getmarkers(nrow(samples), smile = smile)
    scatter!(ax, samples.PC1, samples.PC2, label = samples.Name,
        color = m.colors, marker = m.markers, markersize = 15)
    marker_elements = [MarkerElement(markersize = 15, marker = m.markers[i], color = m.colors[i]) for i = 1:nrow(samples)]
    Legend(f[1, 2], marker_elements, samples.Name)
    return f
end

"""
    plotpca(populations::Array{Population}; title = "PCA Plot", smile = false)

Plot all populations. The names of the populations are used in the plot legend.
If smile == true the last marker is a smiley.

This method generates a label for each population.
"""
function plotpca(populations::Array{Population}; title = "PCA Plot", smile = false)
    f = Figure()
    ax = Axis(f[1, 1],
        title = title,
        xlabel = "PC1",
        ylabel = "PC2",
    )
    m = _getmarkers(length(populations), smile = smile)
    for (i, p) in enumerate(populations)
        scatter!(ax, p.samples.PC1, p.samples.PC2, label = p.name,
            marker = m.markers[i], color = m.colors[i])
    end
    marker_elements = [MarkerElement(marker = m.markers[i], color = m.colors[i]) for i = 1:length(populations)]
    names = [populations[i].name for i = 1:length(populations)]
    Legend(f[1, 2], marker_elements, names)
    return f
end

"""
    _getmarkers(n::Integer; smile = false)

Create a set of n markers for a plot. If smile == true
the last marker is a smiley.

Return a named tuple (markers, colors) that contains an array
of markers and an array of colors for the plot.
"""
function _getmarkers(n::Integer; smile = false)
    markers = [:circle, :rect, :diamond, :hexagon, :cross,
        :xcross, :utriangle, :dtriangle, :ltriangle, :rtriangle,
        :pentagon, :star4, :star5, :star6, :star8,
        :vline, :hline]
    smiley = '😄'

    colors = [:grey0, :green1, :blue, :chocolate, :cyan3,
        :red, :fuchsia, :orange, :darkred]

    # Calculate markers and colors.
    plotmarkers = []
    plotcolors = []
    for i in 1:n
        j = i % length(markers) + 1
        push!(plotmarkers, markers[j])
        k = i % length(colors) + 1
        push!(plotcolors, colors[k])
    end

    # Replace last element with a smiley.
    if smile == true
        pop!(plotmarkers)
        push!(plotmarkers, smiley)        
    end
    return (markers = plotmarkers, colors = plotcolors)
end

