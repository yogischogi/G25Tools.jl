"""
    Distances

    Holds the distances from multiple source samples to a single
    target sample. Each source sample consists of a name and it's
    G25 distance to the target sample.
"""
struct Distances
    targetname::String
    sourcesamples::DataFrame
end

"""
    Population

    Represents a population. Each population consists of a name,
    for example "Viking", and a list of samples in G25 format.
"""
struct Population
    name::String
    samples::DataFrame
end


