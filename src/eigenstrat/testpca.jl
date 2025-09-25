# Testing functions.
# Most of these functioons need an installed AADR database.


include("ancestrymap.jl")



function test_compare()
    file1 = "temp.geno"
    file2 = "temp2.geno"
    file3 = "temp3.geno"

    a = read(file1)
    b = read(file2)
    c = read(file3)

    count = 0
    for (i, _) in enumerate(a)
        if a[i] != c[i]
            as = string(a[i], base = 2)
            bs = string(b[i], base = 2)
            cs = string(c[i], base = 2)
            println("Difference at position: $i, a = $as, c = $cs")
            println(as)
            println(bs)
            println(cs)
            println(a[i: i + 10])
            println(c[i: i + 10])
            println()
            count += 1
        end
    end
    println("Number of differences: $count")
end


function test_read_geno()
    # 213 MB, 1/10 database: 405s read time, 4.4 GB.
    #genofile  = "../database/v62.0_1240k_public.geno"
    #indfile   = "../database/v62.0_1240k_public.ind"
    #snpfile   = "../database/v62.0_1240k_public.snp"

    #genofile  = "../NPRdata.geno"
    #indfile   = "../NPRdata.ind"
    #snpfile   = "../NPRdata.snp"

    genofile  = "temp.geno"
    indfile   = "temp.ind"
    snpfile   = "temp.snp"

    #individuals = readlines(indfile)
    #nind = length(individuals)
    nind = countlines(indfile)
    println("Number of individuals: $nind")

    #snps = readlines(snpfile)
    #nsnp = length(snps)
    nsnp = countlines(snpfile)
    println("Number of SNPs: $nsnp")

    # Choose individuals
    inds = collect(1:nind)
    println("Choosing $inds individuals.")
    geno = read_packedancestrymap(genofile, nsnp, nind, inds)
end

function test_populations()
    inds = read_eigenstrat_inds("../database/v62.0_1240k_public.ind")
    pops = populations(inds)
    n = 0
    nmin = 100
    mainpops = String[]
    for p in eachrow(pops)
        if p.Count >= nmin
            println(p.Name)
            push!(mainpops, p.Name)
            n += 1
        end
    end
    println("$n populations with at least $nmin individuals.")
    writedlm("temp_populations.txt", mainpops)
end

"""
A test example.
"""
function test_reduce_database()
    individuals = readlines("temp2.ind")
    nind = length(individuals)
    println("Number of individuals: $nind")

    snps = readlines("temp2.snp")
    nsnp = length(snps)
    println("Number of SNPs: $nsnp")

    # Choose individuals
    indies = collect(1:35)
    println("Choosing $indies individuals.")
    geno = read_packedancestrymap("temp2.geno", nsnp, nind, indies)

    # Write selected individuals to new data base.
    open("temp3.snp", create = true, write = true) do file
        for line in snps
            write(file, line)
            write(file, "\n")
        end
    end

    open("temp3.ind", create = true, write = true) do file
        for i in indies
            write(file, individuals[i])
            write(file, "\n")
        end
    end

    ihash = hash_ids("temp3.ind")
    shash = hash_ids("temp3.snp")

    write_packedancestrymap("temp3.geno", geno; ind_hash = ihash, snp_hash = shash)
end

function test_expand_reduce()
    # Add an individual to database temp and save it as temp2.
    add_individual("temp", "temp2", "dirk-autosomal.txt", "Dirk"; gender="M")
    # Remove individual and save it as temp3.
    test_reduce_database()
end


# Working with full AADR database.
function example_2()
    individuals = readlines("../database/v62.0_HO_public.ind")
    nind = length(individuals)
    println("Number of individuals: $nind")

    snps = readlines("../database/v62.0_HO_public.snp")
    nsnp = length(snps)
    println("Number of SNPs: $nsnp")

    # Choose individuals
    indies = collect(10:100:nind)
    println("Choosing $indies individuals.")
    geno = read_packedancestrymap("../database/v62.0_HO_public.geno", nsnp, nind, indies)

    # Write selected individuals to new data base.
    open("temp.snp", create = true, write = true) do file
        for line in snps
            write(file, line)
            write(file, "\n")
        end
    end

    open("temp.ind", create = true, write = true) do file
        for i in indies
            write(file, individuals[i])
            write(file, "\n")
        end
    end

    ihash = hash_ids("temp.ind")
    shash = hash_ids("temp.snp")

    write_packedancestrymap("temp.geno", geno; ind_hash = ihash, snp_hash = shash)
end




function test_extract_hgdp()
    genofile  = "../database/v62.0_1240k_public.geno"
    indfile   = "../database/v62.0_1240k_public.ind"
    snpfile   = "../database/v62.0_1240k_public.snp"

    hgdpgenofile = "hgdp.geno"
    hgdpindfile  = "hgdp.ind"
    hgdpsnpfile  = "hgdp.snp"

    # Extract modern populations from the Human Genome Diversity Project, HGDP.

    # Indices of HGDP individuals.
    hgdp_vec = Int64[]
    # Filter entries that start with HGDP.
    inds = read_eigenstrat_ind(indfile)
    hgdp_inds = similar(inds, 0)
    for (i, row) in enumerate(eachrow(inds))
        if startswith(row.ID, "HGDP")
            push!(hgdp_vec, i)
            push!(hgdp_inds, row)
        end
    end

    # Write hgdp database.
    write_eigenstrat_ind(hgdpindfile, hgdp_inds)
    cp(snpfile, hgdpsnpfile; force = true)

    nsnp = countlines(snpfile)
    nind = countlines(indfile)
    genodata = read_packedancestrymap(genofile, nsnp, nind, hgdp_vec)

    ihash = hash_ids(hgdpindfile)
    shash = hash_ids(hgdpsnpfile)
    write_packedancestrymap(hgdpgenofile, genodata; ind_hash = ihash, snp_hash = shash)
end

function test_pca()
    genofile = "hgdp.geno"
    indfile  = "hgdp.ind"
    snpfile  = "hgdp.snp"
    popsfile = "hgdp_populations.txt"

    # Load genotypes for test puprose only.
    nsnp = countlines(snpfile)
    nind = countlines(indfile)
    genodata = read_packedancestrymap(genofile, nsnp, nind, collect(1:nind))

    # Extract populations.
    inds = read_eigenstrat_ind(indfile)
    pops = populations(inds)
    writedlm(popsfile, pops.Name)
    
    # Add individual.
    add_individual("hgdp", "hgdp_dirk", "testdata/dirk-family-finder.txt", "Dirk"; gender = "M")

    run(`./smartpca -p par.SMARTPCAhgdp`)
end

function test_distances()
    # Read 25 eigenvectors from smartpca output.
    colnames = ["NAME", "PC1", "PC2", "PC3", "PC4", "PC5",
        "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12",
        "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19",
        "PC20", "PC21", "PC22", "PC23", "PC24", "PC25"
    ]
    samples = CSV.read("temp_Eigenvectors_hgdp_dirk", DataFrame; header = colnames,
        comment = "#", delim = ' ', ignorerepeated = true)

    dists = distances(samples, samples[end, :])

    # Joins samples to populations.
    inds = read_eigenstrat_ind("hgdp.ind")
    rename!(inds, [:Name, :Gender, :Status])

    result = leftjoin(dists, inds, on = :Name)
    sorted = sort(result, :Distance, rev = false)
    CSV.write("temp_xxx_distances.txt", sorted)

end


