# G25Tools

G25 Tools is a library for the Julia programming language
to work with G25 coordinates.

G25 coordinates were introduced by the blogger
[Davidski](https://eurogenes.blogspot.com/2019/07/getting-most-out-of-global25_12.html)
to compare genetic samples and populations. They have become
increasingly popular in recent years and are also used by
several companies.

For users there are several online calculators available at
[ExploreYourDNA](https://www.exploreyourdna.com/calculators.aspx)
but compared to a programming language their use is rather restricted.

G25Tools is a library that gives you more freedom to play
around with DNA results and allows for more in-depth exploration.
You can:

1. Compare genetic samples and populations.
2. Filter existing samples to create new populations.
3. Generate simple PCA plots (PCA = Principal Component Analysis).

## Getting started

### Installation

1. Install [Julia](https://julialang.org/downloads/) on
   your computer.
2. At the Julia prompt type `]` to switch to the package
   manager. Type `add https://github.com/yogischogi/G25Tools.jl`
   to install the package.

### First steps

1. Download the [main file of ancient samples](https://www.exploreyourdna.com/ancient-samples.aspx)
   from Explore your DNA.
2. I also highly recommend to get your G25 coordinates from
   [G25 requests](https://g25requests.app/). 
   If you already know your K36 coordinates you can also begin
   with simulated G25 coordinates from
   [Explore your DNA](https://www.exploreyourdna.com/simulated-g25.aspx).
   Note that this step is not necessary to get the examples running.
3. Start working with the [examples](https://github.com/yogischogi/G25Tools.jl/tree/master/src/examples).




