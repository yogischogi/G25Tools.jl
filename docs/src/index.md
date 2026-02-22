# G25Tools
```@meta
CurrentModule = G25Tools
```
## Introduction

G25Tools is a library for the [Julia](https://julialang.org/)
programming language to work with G25 coordinates.

G25 coordinates were introduced by the blogger
[Davidski](https://eurogenes.blogspot.com/2019/07/getting-most-out-of-global25_12.html)
to compare genetic samples and populations.

## First steps

### Download ancient DNA data

1. Download the
   [spreadsheet of ancient samples](https://docs.google.com/spreadsheets/d/1EJZTb34IwUNudSQ50uqvBPKDNsL9DDC90lpN586LcTY/)
   in CSV format.
   As far as I know the samples were assembled by the user `nomad` from the
   [GenArchivist](https://www.genarchivist.net/) forum.

   Alternatively you can download some
   [sample sets from Davidski](https://eurogenes.blogspot.com/2019/07/getting-most-out-of-global25_12.html).
   However these collections are not up-to-date.

2. I also highly recommend to get your G25 coordinates from
   [G25 Requests](https://g25requests.app/). 
   If you already know your K36 coordinates you can also begin
   with simulated G25 coordinates from
   [ExploreYourDNA](https://www.exploreyourdna.com/simulated-g25.aspx).
   Note that this step is not necessary to get the examples running.

## Explore the [examples](https://github.com/yogischogi/G25Tools.jl/tree/main/examples)

1. [`compare_single.jl`](https://github.com/yogischogi/G25Tools.jl/blob/main/examples/01_compare_single.jl)\
   Compares your G25 coordinates to a list of ancient samples.

2. [`compare_multi.jl`](https://github.com/yogischogi/G25Tools.jl/blob/main/examples/02_compare_multi.jl)\
   Compares a list of target samples to a list of ancient DNA source samples.

3. [`compare_to_populations.jl`](https://github.com/yogischogi/G25Tools.jl/blob/main/examples/03_compare_to_populations.jl)\
   Compares a single sample to a list of selected populations.

4. [`relatives_through_time.jl`](https://github.com/yogischogi/G25Tools.jl/blob/main/examples/04_relatives_through_time.jl)\
   Finds your closest ancient relatives for different time periods.

5. [`filter_haplogroups.jl`](https://github.com/yogischogi/G25Tools.jl/blob/main/examples/05_filter_haplogroups.jl)\
   Find your closest ancient relatives who belong to certain haplogroups.

6. [`plot_PCA.jl`](https://github.com/yogischogi/G25Tools.jl/blob/main/examples/06_plot_PCA.jl)\
   Create a simple PCA plot that compares your sample to populations from different countries.

7. [`plot_celtic_germanic.jl`](https://github.com/yogischogi/G25Tools.jl/blob/main/examples/07_plot_celtic_germanic.jl)\
   Another simple PCA plot that shows Celtic and Germanic populations.

8. [`map_samples.jl`](https://github.com/yogischogi/G25Tools.jl/blob/main/examples/08_map_samples.jl)\
   Displays ancient samples on a map.

9. [`ancestral_population.jl`](https://github.com/yogischogi/G25Tools.jl/blob/main/examples/09_ancestral_population.jl)\
   Find ancient samples within a specified time frame that are genetically close
   to a set of target samples and show them on a map.

10. [`migration_animation.jl`](https://github.com/yogischogi/G25Tools.jl/blob/main/examples/10_migration_animation.jl)\
    Create an animation that shows R1b migrations.

11. [`relatives_animation.jl`](https://github.com/yogischogi/G25Tools.jl/blob/main/examples/11_relatives_animation.jl)\
    Create an animation of your historical relatives depending on genetic distance.


## Functions
```@autodocs
Modules = [G25Tools]
```

## Index
```@index
```

