module TS

using DataFrames
using Mmap
using PyCall
using Random
using Serialization
using Statistics
using xyBnG.XY
using xyBnG.xyTypes

"""
    tsbits(tsf::AbstractString)

Extract the haplotypes from a tskit file. Also return the positions, the
reference alleles and the allele frequencies.
"""
function tsbits(tsf::AbstractString)
    tskit = pyimport("tskit")
    ts = tskit.load(tsf)
    gt = ts.genotype_matrix()
    nlc = size(gt, 1)
    vld = zeros(Bool, nlc)
    frq = zeros(nlc)
    Threads.@threads for i in 1:nlc
        vld[i] = maximum(gt[i, :]) == 1
        frq[i] = mean(gt[i, :])
    end
    Bool.(gt[vld, :]), Int32.(ts.tables.sites.position[vld]), Char.(ts.tables.sites.ancestral_state[vld]), frq[vld]
end

"""
    toxy(dir::AbstractString)

Merge the ts files in `dir` to a xy file, and serialized linkage map DataFrame.
"""
function toxy(dir::AbstractString)
    @info "Merge the ts files of $pop in $dir to a xy file"
    name = readline("$dir/desc.txt")
    nchr = sum(.!isnothing.(match.(r"ts$", readdir(dir))))
    gt, pos, ref, frq = tsbits("$dir/1.ts")
    lmp = DataFrame(chr = Int8(1),
                    pos = pos,
                    ref = ref,
                    frq = frq)
    print(1)
    for chr in 2:nchr
        a, p, r, f = tsbits("$dir/$chr.ts")
        gt = vcat(gt, a)
        append!(lmp, DataFrame(chr = Int8(chr),
                               pos = p,
                               ref = r,
                               frq = f))
       print(" $chr")
    end
    println()
    @info "  - writing haplotypes to the disk"
    open("$dir/$name.xy", "w") do io
        hdr = XY.header()
        write(io, Ref(hdr), Ref(size(gt)))
        hp = 0
        for i in eachcol(gt)
            hp += 1
            hp % 100 == 0 && print(" $hp")
            write(io, Int8.(i))
        end
    end
    println()
    serialize("$dir/name.lmp", lmp)
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
    vld = vld .&& frq .> maf .&& frq .< 1. - maf
    (Bool.(gt[vld, :]),
     Int32.(ts.tables.sites.position[vld]),
     Char.(ts.tables.sites.ancestral_state[vld]),
     frq[vld]
    )
end

"""
    function sample2xy(dir::AbstractString, dst::AbstractString,
        nid::Int, nlc::Int...; maf = 0.)

Sample `nid` individuals from the haplotypes in `dir` and save the sample `xy`
and `lmp` files into directory `dst`.
"""
function sample2xy(dir::AbstractString, dst::AbstractString,
    nid::Int, nlc::Int...; maf = 0.)
    pop, tid = begin
        info = readlines("$dir/desc.txt")
        info[1], parse(Int, info[2])
    end
    0 < nid ≤ tid && 0. ≤ maf ≤ 0.4 && all(nlc .> 0) ||
        throw(ArgumentError("Invalid argument(set)"))
    
    @info "Sample $nid ID with $nlc loci and MAF $maf from $dir"
    ids = sort(shuffle(1:tid)[1:nid]) # sample IDs
    hps = sort([2ids .- 1; 2ids])
    nchr = sum(.!isnothing.(match.(r"ts$", readdir(dir))))
    gt, pos, ref, frq = tsbits("$dir/1.ts", hps, maf)
    lmp = DataFrame(chr = Int8(1),
                    pos = pos,
                    ref = ref,
                    frq = frq)
    print(1)
    for chr in 2:nchr
        a, p, r, f = tsbits("$dir/$chr.ts", hps, maf)
        gt = vcat(gt, a)
        append!(lmp, DataFrame(chr = Int8(chr),
                               pos = p,
                               ref = r,
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