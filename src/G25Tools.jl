module G25Tools

using CairoMakie, CSV, DataFrames, LinearAlgebra, Statistics

export readG25, extractG25
export average, Distances, distance, distances, getyear, pick
export plotpca, Population

include("basic.jl")
include("pca.jl")


end # module G25Tools
