push!(LOAD_PATH, "../src/")
using Documenter, G25Tools

makedocs(
    sitename = "G25Tools.jl",
    modules = [G25Tools]
#    remotes = nothing
)
deploydocs(;
    repo="github.com/yogischogi/G25Tools.jl.git",
    devbranch = "main"
)


