module Founder

using DataFrames
using Distributions
using Mmap
using Serialization
using xyBnG.Conn.TS
using xyBnG.XY
import xyBnG.xyTypes: Cattle, Species, Trait

"""
    sim_base(pop::Cattle, dir::AbstractString)

Simulate a cattle population of name `pop.name`, and `pop.nid` ID in `dir`.

Note: This is a simulation with the coancestor/backward simulator msprime.
"""
function sim_base(pop::Cattle, dir::AbstractString)
    '~' ∈ dir && error("Character '~' is forbidden in dir")
    isdir(dir) || mkpath(dir)
    @info "Simulating a cattle population of name $(pop.name), and $(pop.nid) ID in $dir
      Note: This is a simulation with the coancestor/backward simulator msprime."
    Threads.@threads for chr in 1:29
        print(" $chr")
        cmd = pipeline(`stdpopsim BosTau -c $chr -o $dir/$chr.ts
            -d HolsteinFriesian_1M13 Holstein_Friesian:$(pop.nid)`, stderr=devnull)
        run(cmd)
    end
    open("$dir/desc.txt", "w") do io
        println(io, pop.name)
        println(io, pop.nid)
    end
end

"""
    sim_base(pop::Species, dir::AbstractString)

Default function to catch unsupported species.
"""
function sim_base(pop::Species, dir::AbstractString)
    @info "This species is not supported yet." # or as below
    # https://popsim-consortium.github.io/stdpopsim-docs/stable/api.html#sec-api-generic-models
end

"""
    norm_v(Q::AbstractMatrix, v; ϵ = 1e-6, μ=0., σ=1.)

Normalize QTL effect `v`, such that the `Q'v` has mean `μ` and standard deviation
`σ`.
"""
function norm_v(Q::AbstractMatrix, v; ϵ = 1e-6, μ=0., σ=1.)
    nqtl= size(Q, 1)
    gv = Q'v
    m, s = mean(gv), std(gv)
    while abs(m - μ) > ϵ || abs(s - σ) > ϵ
        v .-= (m - μ)/nqtl
        v .*= σ / s
        gv = Q'v
        m, s = mean(gv), std(gv)
    end
end

"""
    uniq(ixy::T, oxy::T) where T <: AbstractString

Make SNP of a (founder) population unique.
"""
function uniq(ixy::T, oxy::T) where T <: AbstractString
    hdr, (nlc, nhp) = XY.header(ixy), XY.dim(ixy)
    hdr.u == 0 || error("Not a SNP file")
    type = XY._type(hdr.type)
    if nhp > typemax(Int32) - 2
        error("Too many ID in the file")
    elseif nhp > typemax(Int16) - 2
        type = Int32
    else nhp > typemax(Int8) - 2
        type = Int16
    end
    hdr.type, hdr.u = XY._type(type), 1
    gt = Mmap.mmap(ixy, Matrix{Int8}, (nlc, nhp), 24)
    v = zeros(type, nlc) # for type safety
    open(oxy, "w") do io
        write(io, Ref(hdr), [nlc, nhp])
        for i in 1:nhp
            copyto!(v, gt[:, i])
            v .+= 2(i - 1)
            write(io, v)
        end
    end
end

function sample_founder(bdir::AbstractString, tdir::AbstractString,
    nid::Int, nchp::Int, nref::Int, trts::Trait...)
    nqtl = [t.nQTL for t in trts]
    name = readline("$bdir/desc.txt")
    lmp = TS.sample2xy(bdir, tdir, nid, nchp, nref, nqtl...)
    dic = Dict("ssA" => "chip", "ssB" => "dark")
    c = 'C'
    for t in trts
        dic["ss$c"] = t.name
        c += 1
    end
    rename!(lmp, dic)
    ped = DataFrame(id = Int32.(1:nid),
                    sire = Int32(0),
                    dam = Int32(0),
                    sex = rand(Int8.(0:1), nid),
                    grt = Int16(1),
                    )
    snps = mmap("$tdir/$name.xy", Matrix{Int8}, (nrow(lmp), 2nid), 24)
    for t in trts
        qtl = lmp[!, t.name]
        Q = snps[qtl, 1:2:end] + snps[qtl, 2:2:end]
        D = Q .== 1
        a, d = rand(t.da, t.nQTL), rand(t.dd, t.nQTL)
        norm_v(Q, a)
        norm_v(D, d, σ = sqrt(t.vd))
        lmp[!, "$(t.name)_a"] .= 0.
        lmp[!, "$(t.name)_d"] .= 0.
        lmp[!, "$(t.name)_a"][qtl] .= a
        lmp[!, "$(t.name)_d"][qtl] .= d
        ped[!, "tbv_$(t.name)"] = Q'a
        ped[!, "ebv_$(t.name)"] .= -1e8
        ped[!, "gt_$(t.name)"] = ped[!, "tbv_$(t.name)"] + D'd
    end
    serialize("$tdir/$name.ped", ped)
    serialize("$tdir/$name.lmp", lmp)
    ped
end

end # module Founder