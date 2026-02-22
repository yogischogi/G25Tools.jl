push!(LOAD_PATH, "../src/")
using Documenter, G25Tools

makedocs(
    sitename = "G25Tools.jl", 
    remotes = nothing
)

