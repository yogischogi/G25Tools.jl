module G25Tools

using CairoMakie, CSV, DataFrames, LinearAlgebra, Statistics

export readG25, extractG25, writedistances
export average, Distances, distance, distances, medianavg
export getyear, picksamples, topmatches
export clusters
export plotpca, Population

include("structs.jl")
include("basic.jl")
include("pca.jl")


end # module G25Tools
