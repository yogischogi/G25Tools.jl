# Read and write packedancestrymap files that were defined
# by David Reich's laboratory.
# https://reich.hms.harvard.edu/software/InputFileFormats
#
# The format is used for the AADR database.
# https://dataverse.harvard.edu/dataverse/reich_lab 

using CSV, DataFrames, DelimitedFiles

"""
    _hashit(sequence::String)

Calculate the hashsum of a String.

This is basically the same method as in the software from David Reich's laboratory:
https://github.com/DReichLab/EIG

However, the original C version uses 32 bit integer values and integer
overflows which may result in undefined behavior on some machines. 
This method uses 64 bit integers and is well defined.
"""
function _hashit(sequence::String)
    hash::Int64 = 0
    for c in sequence
        hash *= 23
        hash = hash % (2^32) 
        hash += Int64(c)
    end
    return hash
end

"""
    _hashsum(sequences::Array{String})

Calculate the hashsum of an array of Strings using the _hashit() method.

Again this is basically the same method as in the software from David Reich's laboratory:
https://github.com/DReichLab/EIG
"""
function _hashsum(sequences::Array{String})
    hash::Int64 = 0
    for s in sequences
        thash = _hashit(s)
        hash *= 17
        hash = hash % (2^32) 
        hash = xor(hash, thash)
    end
    return hash
end

"""
    hash_ids(filename::String)

Create the hashsum of .snp and .ind files.
SNP and Individual files contain an ID in the first column.
This method uses those IDs to calculate a hashsum.
The sums are needed to store genotype data in packed format.
"""
function hash_ids(filename::String)
    data = readdlm(filename, String)
    hash = _hashsum(data[:, 1])
    return hash
end

"""
    _bitpair(byte::UInt8, pos::Int64)

Extract 2 bits from a byte. The bit positions starts with 0.
Allowed are values 0, 1, 2, 3.
"""
function _bitpair(byte::UInt8, pos::Int64)
    if pos == 0
        (byte & 0b11000000) >> 6
    elseif pos == 1
        (byte & 0b00110000) >> 4
    elseif pos == 2
        (byte & 0b00001100) >> 2
    elseif pos == 3
        (byte & 0b00000011)
    end
end

"""
    _alleles(individual::Int64, bytes::Array{UInt8})

Return the number of variant alleles for the specified individual.
individual: Index of the individual.
bytes: Row of bytes that encode an SNP for all individuals.
"""
function _alleles(individual::Int64, bytes::Array{UInt8})
    # Fine byte that encodes the specified SNP.
    pos = Int64(ceil(individual / 4))
    byte = bytes[pos]

    # Extract bitpair that contains individual's data.
    # bitpair index starts with 0.
    bitpair_no = (individual - 1) % 4
    bits = _bitpair(byte, bitpair_no)
    return UInt8(bits)
end

"""
    read_packedancestrymap(genofile::String, nsnp::Int64, nind::Int64, indvec::Array{Int64})

Read a genofile in PackedAncestryMap format. The file must be unzipped.

genofile: filename
nsnp: number of SNPs, listed in .snp file.
nind: number of individuals, listed in .ind file.
indvec: Indices of individuals that should be read from the file.

XXX Check for comment lines in .snp and .ind files.

File description: File header starts with GENO or TGENO (transposed GENO).
So far files in the AADR archive seem to be GENO. So this method
does not support the transposed TGENO format.

The text format contains one line per genotype:
SNP_ID  Sample_ID   Number_of_variant_alleles

The packed format:
Each SNP entry has 2 bits: 0, 1, 2, 3=missing, that denote the number
of variant alleles (https://reich.hms.harvard.edu/software/InputFileFormats).
"""
function read_packedancestrymap(genofile::String, nsnp::Int64, nind::Int64,
                                indvec::Array{Int64})
    result = zeros(UInt8, nsnp, length(indvec))

    # 1 SNP value for 4 individuals is encoded as 1 byte.
    bytes_per_line = Int64(ceil(nind / 4))
    # Header size must be at least 48 bytes.
    if bytes_per_line < 48
        bytes_per_line = 48
    end

    # Read file line by line and extract SNP alleles for all individuals.
    open(genofile) do io
        buffer = Array{UInt8}(undef, bytes_per_line)
        # Skip first line because it is a header.
        read!(io, buffer)
        for snp = 1:nsnp
            read!(io, buffer)
            # Extract alleles for all individuals.
            for i in 1:length(indvec)
                result[snp, i] = _alleles(indvec[i], buffer)
            end
        end
    end
    return result
end

"""
    write_packedancestrymap(genofile::String, genomatrix::Matrix{UInt8})

Write a genotype matrix to file in PackedAncestryMap format.

ind_hash: Hashsum of .ind file.
snp_hash: Hashsum of .snp file.
"""
function write_packedancestrymap(genofile::String, genomatrix::Matrix{UInt8};
    ind_hash::Int64 = 0, snp_hash::Int64 = 0)

    (rows, cols) = size(genomatrix)

    # 4 columns are encoded as 1 byte.
    bytes_per_line = Int64(ceil(cols / 4))
    # Header size must be at least 48 bytes.
    minlen = 48
    fillbytes = zeros(UInt8, 0)
    if bytes_per_line < minlen
        fillbytes = zeros(UInt8, minlen - bytes_per_line)
        bytes_per_line = minlen
    end

    # Create header.
    ihash = string(ind_hash, base = 16)
    shash = string(snp_hash, base = 16)
    buffer = zeros(UInt8, bytes_per_line)
    io = IOBuffer(buffer, write = true)
    header = "GENO $cols $rows $ihash $shash"
    # This should never happen, but the file would still be
    # useable if the hash sums are not checked.
    if length(header) > bytes_per_line
        header = header[1: bytes_per_line]
    end
    write(io, header)

    # Write to file.
    open(genofile, create = true, write = true) do file
        write(file, buffer)
        for row in 1:rows
            byte::UInt8 = 0
            for col in 1:cols
                value = genomatrix[row, col]
                # Encode in byte.
                bitpair_pos = (col - 1) % 4
                shift = 6 - 2 * bitpair_pos
                byte = byte | (value << shift)

                # Write byte to file.
                if  bitpair_pos == 3 || col == cols
                    write(file, byte)
                    byte = 0
                end
            end
            write(file, fillbytes)
        end
        flush(file)    
    end
end

"""
    individuals(indfile::AbstractString, populations::Array{String})

Return of vector of indices that contains the individuals that
belong to populations.

This is useful if you want to reduce the large AADR database to a
smaller set of populations and run smartpca on it.
"""
function individuals(indfile::AbstractString, populations::Array{String})
    indvec = Int64[]
    inds = read_eigenstrat_ind(indfile)
    pops = Set(populatioons)
    rows = nrow(inds)
    for i in rows
        if in(row.Status, pops)
            push!(indvec, i)
        end
    end
    return indvec
end

# Functions to import data from different vendors to plink
#
# Plink only supports the 23andMe file format:
#   # Comment
#   rsid   chromosome  position    genotype
#   rs3094315   1   742429  AG
#
# So it is best to convert to 23andMe format.
# Also only the SNPs that are covered by the AADR database
# should be exported. 
# A list of all SNPs is in the Plink .bim file or in the .snp
# file of the AADR (Packed Ancestry Map format).

"""
    read_snp_file(filename)

Read file with autosomal results from FTDNA Family Finder, My Heritage
or LivingDNA. Should also work with 24andMe files but not tested.
Return a DataFrame containing the columns:
rsid  chromosome  position  genotype
"""
function read_snp_file(filename)
    snptable = CSV.read(filename, DataFrame; header = ["rsid", "chromosome", "position", "genotype"],
        types = [String, String, String, String], comment = "#")
    # Drop header line (MyHeritage, 23andMe).
    if snptable.rsid[1] == "RSID" || snptable.rsid[1] == "rsid"
        delete!(snptable, 1)
    end
    return snptable
end

"""
    write_23andMe(filename, snptable)

Write a table of SNPs in 23andMe file format.
"""
function write_23andMe(filename, snptable)
    CSV.write(filename, snptable; delim = "\t", header = ["#rsid", "chromosome", "position", "genotype"])
end

"""
    prune_snps(snpfile::AbstractString, individual_snps)

Keeps only SNPs that are listed in the snpfile and in the
individual SNPs. This method keeps the SNPs in data base order.
This method is not related to plink's prune function which takes
care of linkage disequilibrium.

Use this method before adding a new individual to the AADR database.

snpfile: File in David Reich's .bim (plink) or .snp (Ancestrymap) format

Format of the Ancestrymap .snp file, separated by one or more spaces:
    rs3094315     1        0.020130          752566 G A
Format of the Plink .bim file, tab-separated:
1	rs3094315	0.02013	752566	G	A

XXX Some SNPs are triallelic. They should be removed.
"""
function prune_snps(snpfile::AbstractString, individual_snps::DataFrame)
    # Put all individual_snps into a Dictionary.
    indiv_snps = Dict{String, DataFrameRow}()
    for row in eachrow(individual_snps)
        # Use rsid als key.
        indiv_snps[row[1]] = row
    end

    # Read SNPS that are already in data base.
    local data_base_snps::Vector{String}
    if endswith(snpfile, ".bim")
        data = CSV.read(snpfile, DataFrame; header = ["chromosome", "rsid", "cM", "position", "allele_1", "allele_2"])
        data_base_snps = data.rsid
    elseif endswith(snpfile, ".snp")
        data = CSV.read(snpfile, DataFrame;
            header = ["rsid", "chromosome", "cM", "position", "allele_1", "allele_2"],
            delim = ' ', ignorerepeated = true)
        data_base_snps = data.rsid
    else
        throw("prune_snps: Wrong file format! File must end in .bim or .snp.")
    end

    # Keep individual's SNPs that are also listed in the data base.
    result = DataFrame()
    for snp in data_base_snps
        if haskey(indiv_snps, snp)
            push!(result, indiv_snps[snp])
        end
    end
    return result
end

"""
    combine_snps(snp_sets::Array{DataFrame})

Combines all SNP entries from an array of DataFrames.
SNP are not sorted in any particular order.
"""
function combine_snps(snp_sets::Array{DataFrame})
    # Put all SNPs from first DataFrame into a dictionary.
    snps = Dict{String, DataFrameRow}()
    for row in eachrow(snp_sets[1])
        # Use rsid als key.
        snps[row[1]] = row
    end

    # Add SNPs from other data sets.
    for s in snp_sets[2:end]
        for row in eachrow(s)
            if !haskey(snps, row[1])
                snps[row[1]] = row
            end
        end
    end

    # Move SNPs from Dictionary to DataFrame.
    result = DataFrame()
    for value in values(snps)
        push!(result, value)
    end
    return result
end


function read_eigenstrat_snps(snpfile::AbstractString)
    local eigenstrat_snps::DataFrame
    if endswith(snpfile, ".bim")
        eigenstrat_snps = CSV.read(snpfile, DataFrame; header = ["chromosome", "rsid", "cM", "position", "allele_1", "allele_2"])
    elseif endswith(snpfile, ".snp")
        eigenstrat_snps = CSV.read(snpfile, DataFrame;
            header = ["rsid", "chromosome", "cM", "position", "allele_1", "allele_2"],
            delim = ' ', ignorerepeated = true)
    else
        throw("read_eigenstrat: Wrong file format! File must end in .bim or .snp.")
    end
    return eigenstrat_snps
end

"""
Read individuals from Eigenstrat .ind file.
Gender: M (male), F (Female), U (unknown).
Status: Case, Control or population label.
"""
function read_eigenstrat_ind(indfile::AbstractString)
    inds = CSV.read(indfile, DataFrame; header = ["ID", "Gender", "Status"],
        delim = ' ', ignorerepeated = true)
    return inds
end

function write_eigenstrat_ind(filename::String, inds::DataFrame)
    CSV.write(filename, inds; writeheader = false, delim = ' ')
end


"""
    populations(eigenstrat_inds::DataFrame)

Return a DataFrame of populations listed in the ind. file.
The DataFrame contains two columns [Name, Count], the population
name and the number of individuals belonging to that population.
"""
function populations(eigenstrat_inds::DataFrame)
    println("Number of individuals: $(nrow(eigenstrat_inds))")
    popcount = Dict{String, Int64}()
    for p in eigenstrat_inds.Status
        if haskey(popcount, p)
            popcount[p] += 1
        else
            popcount[p] = 1
        end
    end

    statistics = DataFrame(Name = String[], Count = Int64[])
    for (key, value) in popcount
        push!(statistics, (key, value))
    end
    sorted = sort(statistics, :Count, rev = true)

    return sorted
end

"""
Convert a list of SNPs so that in can be added to the Geno matrix.
The result is a vector where ich SNP is represnted by one value:
0: 0 copies of reference allele.
1: 1 copy of reference allele.
2: 2 copies of reference allele.
3: Missing data.
"""
function snps_to_geno(eigenstrat_snps::DataFrame, individual_snps::DataFrame)
    result = Array{UInt8}(undef, length(eigenstrat_snps))

    # Put all individual_snps into a Dictionary.
    indiv_snps = Dict{String, DataFrameRow}()
    for row in eachrow(individual_snps)
        # Use rsid als key.
        indiv_snps[row[1]] = row
    end

    # Create an entry for each eigenstrat_snp.
    for (i, snp) in enumerate(eigenstrat_snps)
        if haskey(indiv_snps, snp.rsid)
            ind = indiv_snps(snp.rsid)
            # Count reference allele.
            count = 0
            if ind.genotype[1] == snp.allele1
                count += 1
            end
            if ind.genotype[2] == snp.allele1
                count += 2
            end
            result[i] = count
        else
            # SNP is missing.
            result[i] = 3
        end
    end
    return result
end

"""
add_individual(eigenstrat_inds::DataFram, genomatrix::Matrix{UInt8},
    ID::String, gender::String, status::String, snps::Array{UInt8})

Add an individual to the eigenstrat database.
eigenstrat_inds: Individuals from the Eigenstrat .ind file.
genomatrix: Genotdata that was read from the .geno file using read_packedancestrygeno().
id: An ID for the new individual.
gender: M, F or U (Male, Female, Unknown).
snps: Number of ancestral allele for each SNP marker in the genomatrix.
    The length of the vector must be identical to the number of rows of the genomatrix.
    Use snps_to_geno() to prepare a set of SNP data.
"""
function add_individual(eigenstrat_inds::DataFrame, genomatrix::Matrix{UInt8},
    id::String, gender::String, status::String, snps::Array{UInt8})

    row = [ID, gender, status]
    push!(result_inds, row)
    result_matrix = hcat(genomatrix, snps)
    return (result_inds, result_matrix)
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
function example_reduce_database()
    individuals = readlines("../NPRdata.ind")
    nind = length(individuals)
    println("Number of individuals: $nind")

    snps = readlines("../NPRdata.snp")
    nsnp = length(snps)
    println("Number of SNPs: $nsnp")

    # Choose individuals
    indies = collect(1:20:nind)
    println("Choosing $indies individuals.")
    geno = read_packedancestrymap("../NPRdata.geno", nsnp, nind, indies)

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

























