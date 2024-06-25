module TS

using DataFrames
using Mmap
using PyCall
using Random
using Serialization
using Statistics
using xyBnG.XY
using xyBnG.xyTypes

#=
"""
    tsbits(tsf::AbstractString)

Extract the haplotypes from a tskit file. Also return the positions, the
reference alleles and the allele frequencies.
"""
function tsbits(tsf::AbstractString)
    tskit = pyimport("tskit")
    ts = tskit.load(tsf)
    gt = Int8.(ts.genotype_matrix())
    nlc = size(gt, 1)
    vld = zeros(Bool, nlc)
    frq = zeros(nlc)
    Threads.@threads for i in 1:nlc
        vld[i] = maximum(gt[i, :]) == 1
        frq[i] = mean(gt[i, :])
    end
    alt = Int8[]
    for v in ts.variants()
        push!(alt, Int8(v.alleles[2][1]))
    end
    (gt[vld, :],
     Int32.(ts.tables.sites.position[vld]),
     Char.(ts.tables.sites.ancestral_state[vld]),
     Char.(alt[vld]),
     frq[vld])
end
=#

"""
    toxy(dir::AbstractString)

Merge the ts files in `dir` to a xy file, and serialized linkage map DataFrame.
The `xy` file is ID majored, in line with the `macs` results. This is to ease 
file appending and use less memory.

This is not necessary most of the time. I can sample the results directly from
the tskit files.
"""
function toxy(dir::AbstractString)
    tskit = pyimport("tskit")
    name = readline("$dir/desc.txt")
    nchr = sum(.!isnothing.(match.(r"ts$", readdir(dir))))
    @info "Merge the $nchr ts files of pop $name in $dir to a xy file:"
    nhp = 0
    lmp = DataFrame(chr = Int8[], pos = Int32[], ref = Char[], alt = Char[],
                    frq = Float64[])
    open("$dir/$name.xy", "w") do io
        hdr = XY.header(major = 1)
        write(io, Ref(hdr), [0, 0])
        for chr in 1:nchr
            ts = tskit.load("$dir/$chr.ts")
            gt = ts.genotype_matrix()
            nlc, nhp = size(gt)
            pos = Int32.(ts.tables.sites.position)
            frq = mean(gt, dims = 2)
            vld = Vector{Bool}(undef, nlc)
            alt = zeros(Int8, nlc)
            for v in ts.variants()
                ilc = v.site.id + 1
                vld[ilc] = length(v.alleles) == 2
                alt[ilc] = v.alleles[2][1]
                write(io, Int8.(gt[ilc, :]))
            end
            append!(lmp, DataFrame(chr = Int8(chr),
                                   pos = pos[vld],
                                   ref = ts.tables.sites.ancestral_state[vld],
                                   alt = alt[vld],
                                   frq = frq[vld]))
            print(" $chr")
        end
    end
    println()
    XY.dim!("$dir/$name.xy", nhp, size(lmp, 1))
    serialize("$dir/$name.lmp", lmp)
end

"""
    tsbits(tsf::AbstractString, hps::Vector{Int}, maf::Float64)

Extract the haplotypes from a tskit file. Also return the positions, the
reference alleles and the allele frequencies.
"""
function tsbits(tsf::AbstractString, hps::Vector{Int}, maf::Float64)
    tskit = pyimport("tskit")
    ts = tskit.load(tsf)
    gt = ts.genotype_matrix(samples = hps)
    nlc = size(gt, 1)
    vld = zeros(Bool, nlc)
    frq = zeros(nlc)
    Threads.@threads for i in 1:nlc
        vld[i] = maximum(gt[i, :]) == 1
        frq[i] = mean(gt[i, :])
    end
    alt = Int8[]
    for v in ts.variants()
        push!(alt, Int8(v.alleles[2][1]))
    end
    vld = vld .&& frq .> maf .&& frq .< 1. - maf
    (Bool.(gt[vld, :]),
     Int32.(ts.tables.sites.position[vld]),
     Char.(ts.tables.sites.ancestral_state[vld]),
     Char.(alt[vld]),
     frq[vld]
    )
end

"""
    function sample2xy(dir::AbstractString, dst::AbstractString,
        nid::Int, nlc::Int...; maf = 0.)

Sample `nid` individuals from the haplotypes in `dir` and save the sample `xy`
and `lmp` files into directory `dst`.

## Todo
- Parallelize TS reading later.
"""
function sample2xy(dir::AbstractString, dst::AbstractString,
    nid::Int, nlc::Int...; maf = 0.)
    pop, tid = begin
        info = readlines("$dir/desc.txt")
        info[1], parse(Int, info[2])
    end
    0 < nid ≤ tid && 0. ≤ maf ≤ 0.4 && all(nlc .> 0) ||
        throw(ArgumentError("Invalid argument(set)"))
    
    @info "  - Sample $nid ID with $nlc loci and MAF $maf from $dir:"
    ids = shuffle(1:tid)[1:nid] # sample IDs
    hps = sort([2ids .- 1; 2ids])
    # parallelize below later
    nchr = sum(.!isnothing.(match.(r"ts$", readdir(dir))))
    gt, pos, ref, alt, frq = tsbits("$dir/1.ts", hps, maf)
    lmp = DataFrame(chr = Int8(1),
                    pos = pos,
                    ref = ref,
                    alt = alt,
                    frq = frq)
    print(1)
    for chr in 2:nchr
        a, p, r, alt, f = tsbits("$dir/$chr.ts", hps, maf)
        gt = vcat(gt, a)
        append!(lmp, DataFrame(chr = Int8(chr),
                               pos = p,
                               ref = r,
                               alt = alt,
                               frq = f))
        print(" $chr")
    end
    println()
    @info "  - writing haplotypes into $dst by haplotype"
    tlc = nrow(lmp)
    loci, set = [], 'A'
    for x in nlc
        ss =  sort(shuffle(1:tlc)[1:x]) # a SNP set
        push!(loci, ss)
        lmp[!, "ss$set"] = zeros(Bool, tlc)
        lmp[!, "ss$set"][ss] .= true
        set += 1
    end
    ass = sort(union(loci...)) # all snp sets
    isdir(dst) || mkpath(dst)
    open("$dst/$pop.xy", "w") do io
        hdr = XY.header()
        write(io, Ref(hdr), [length(ass), 2nid])
        hp = 0
        for x in eachcol(view(gt, ass, :))
            hp += 1
            write(io, Int8.(x))
            hp % 100 == 0 && print(" $hp")
        end
    end
    println()
    serialize("$dst/$pop.lmp", lmp[ass, :])
    lmp[ass, :]
end

end # module tskit
