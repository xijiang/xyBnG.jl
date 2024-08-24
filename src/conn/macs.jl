module MaCS

using DataFrames
using Mmap
using Serialization
using Statistics
using xyBnG.XY
using xyBnG.xyTypes

# ToDo: Modify according to the new xyBnG.XY module.
"""
    function make_macs(; tdir=pwd())
Clone `MaCS` into a temp dir under `tdir` and compile it into `tdir`.
## Note
- designed only for Linux systems.
- `g++`, and `boost-devel` must be installed.
- this function return the absolute path of newly compiled `macs` in the end.
"""
function make_macs(; tdir = pwd())
    macs = Sys.which("macs")
    isnothing(macs) || return macs

    macs = joinpath(abspath(tdir), "macs")
    @debug "Making macs"
    isfile(macs) && return macs
    isdir(tdir) || mkpath(tdir)
    wdir = mktempdir(tdir)
    run(`git clone https://github.com/gchen98/macs $wdir`)
    src = joinpath.(wdir, ["simulator.cpp", "algorithm.cpp", "datastructures.cpp"])
    target = joinpath(tdir, "macs")
    run(`g++ -o $target -O3 -Wall $src`)
    return macs
end

"""
    function read_macs(file)
---
Read genotypes and physical positions from simulation results of `macs`. Returns
- genotypes of Array{Int8, 2}
- physical positions and
- allele frequencies

When `trans` is `true`, the genotypes are transposed, i.e., to of `nLoci ×
nHap`.
"""
function read_macs(file, trans = false)
    gt, ps = Int8[], Int[]
    open(file, "r") do io
        nbp = parse(Float64, split(readline(io))[4])
        for line in eachline(io)
            (line[1:4] ≠ "SITE") && continue
            _, _, p, _, as = split(line)
            push!(ps, Int(ceil(parse(Float64, p) * nbp)))
            a = parse.(Int8, collect(as))
            append!(gt, a)
        end
    end
    nlc = length(ps)
    gt = reshape(gt, :, nlc)
    trans && (gt = gt')
    fq = trans ? mean(gt, dims = 2) : mean(gt, dims = 1)
    return gt, vec(ps), vec(fq)
end

"""
    function toxy(dir; keep = false)
Convert `MaCS` simulation results into `XY` format.  The `MaCS` files are of
`chr.chr-number` and `info.chr-number` in `dir`.  The merged file are stored in
`dir/(pop.name).xy` and `dir/(pop.name).lmp`. The genotypes are of BitArray of
`nLoci × nID`, or `header.r = 1`.

When `keep` is `true`, keep the original `MaCS` files.
"""
function toxy(dir; keep = false)
    isdir(dir) || error("$dir not exists")
    @info "  - Collect genotypes simulated by `MaCS` in $dir"
    fname = readline("$dir/desc.txt")
    fxy, fmd, fmp = joinpath.(dir, ["$fname.xy", "mid.xy", "$fname.lmp"])

    chrs = Int8[]
    for f in readdir(dir)
        occursin.(r"^chr", f) && push!(chrs, parse(Int8, split(f, '.')[2]))
    end
    sort!(chrs)           # chromosome number in integer, and in order
    lmp =
        DataFrame(chr = Int8[], pos = Int64[], ref = Char[], alt = Char[], frq = Float64[])
    hdr = XY.header(major = 1) # ID majored
    thp, tlc = 0, 0
    aa = (
        ('A', 'C'),
        ('A', 'G'),
        ('A', 'T'),
        ('C', 'A'),
        ('C', 'G'),
        ('C', 'T'),
        ('G', 'A'),
        ('G', 'C'),
        ('G', 'T'),
        ('T', 'A'),
        ('T', 'C'),
        ('T', 'G'),
    )
    open(fmd, "w") do io
        write(io, Ref(hdr), [thp, tlc])
        for c in chrs
            print(" $c")
            chr = joinpath(dir, "chr.$c")
            gt, ps, fq = read_macs(chr)
            thp = size(gt, 1)
            tlc += size(gt, 2)
            write(io, gt)
            ref, alt = Char[], Char[]
            for _ = 1:size(gt, 2)
                x = aa[rand(1:12)]
                push!(ref, x[1])
                push!(alt, x[2])
            end
            append!(lmp, DataFrame(chr = c, pos = ps, ref = ref, alt = alt, frq = fq))
        end
    end
    XY.dim!(fmd, thp, tlc)
    serialize(fmp, lmp)
    println()
    XY.tr8bit(fmd, fxy)
    if !keep
        rm.(joinpath(dir, "chr.") .* string.(chrs), force = true)
        rm.(joinpath(dir, "log.") .* string.(chrs), force = true)
        rm(fmd, force = true)
    end
end

end # module MaCS
