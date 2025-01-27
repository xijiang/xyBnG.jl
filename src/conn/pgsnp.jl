"""
    module PG
This is about an old program of mine: `pgsnp`.
"""
module PG
using DataFrames
using Serialization
using xyBnG.XY
function toxy(dir::AbstractString; keep = false, flip = true)
    isdir(dir) || error("$dir not exists")
    @info "  - Collect genotypes simulated by `pgsnp` in $dir"
    fname = readline("$dir/desc.txt")
    fxy, fmd, fmp = joinpath.(dir, ["$fname.xy", "mid.xy", "$fname.lmp"])
    chrs = Int8[]
    for f in readdir(dir)
        occursin.(r"^chr", f) && push!(chrs, parse(Int8, split(f, '.')[2]))
    end
    sort!(chrs)
    lmp =
        DataFrame(chr = Int8[], pos = Int64[], ref = Char[], alt = Char[], frq = Float64[])
    hdr = XY.header()
    fs = []
    for c in chrs
        push!(fs, open(joinpath(dir, "chr.$c"), "r"))
    end
    nid, tlc, fx = 0, 0, []
    for f in fs
        nid = parse(Int, readline(f))
        ilc = parse(Int, readline(f))
        tlc += ilc
        push!(fx, rand(Bool, ilc))
    end
    open(fxy, "w") do io
        write(io, Ref(hdr), [tlc, 2nid])
        for ihp = 1:2nid
            for (i, f) in enumerate(fs)
                line = readline(f)
                hap = Int8.(collect(line) .- '0')
                flip && (hap[fx[i]] .= 1 .- hap[fx[i]])
                write(io, hap)
            end
        end
    end
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
    c = 1
    for f in fs
        while !eof(f)
            line = readline(f)
            ps, fq = parse.(Int, split(line))
            ref, alt = aa[rand(1:12)]
            push!(lmp, (chrs[c], ps, ref, alt, fq / (2nid)))
        end
        c += 1
        close(f)
    end
    keep || rm.("$dir/chr." .* string.(chrs), force = true)
    serialize(fmp, lmp)
end
end
