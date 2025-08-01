# Tests for G25Tools.jl
#
# The tests rely on a file "ancient_dna - main.csv" which is
# not included in the distribution.
# Thus I put it in a parent directory and use the parent directory
# for the tests.

cd("../..")
directory = pwd()
print("Directory used for testing: $directory")

include("$directory/G25Tools.jl/examples/01_compare_single.jl")
println("\nFinished 01_compare_single.jl")

include("$directory/G25Tools.jl/examples/02_compare_multi.jl")
println("\nFinished 02_compare_multi.jl")

include("$directory/G25Tools.jl/examples/03_compare_to_populations.jl")
println("\nFinished 03_compare_to_populations.jl")

include("$directory/G25Tools.jl/examples/04_relatives_through_time.jl")
println("Finished 04_relatives_through_time.jl")

include("$directory/G25Tools.jl/examples/05_filter_haplogroups.jl")
println("Finished 05_filter_haplogroups.jl")

include("$directory/G25Tools.jl/examples/06_plot_PCA.jl")
println("Finished 06_plot_PCA.jl")

include("$directory/G25Tools.jl/examples/07_plot_celtic_germanic.jl")
println("Finished 07_plot_celtic_germanic.jl")

include("$directory/G25Tools.jl/examples/08_map_samples.jl")
println("Finished 08_map_samples.jl")

include("$directory/G25Tools.jl/examples/09_ancestral_population.jl")
println("Finished 09_ancestral_population.jl")


