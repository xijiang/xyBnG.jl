"""
    norm_v(Q::AbstractMatrix, v; ϵ = 1e-6, μ=0., σ=1.)

Normalize QTL effect `v`, such that the `Q'v` has mean `μ` and standard deviation
`σ`.
"""
function norm_v(Q::AbstractMatrix, v; ϵ = 1e-6, μ = 0.0, σ = 1.0)
    nqtl = size(Q, 1)
    gv = Q'v
    m, s = mean(gv), std(gv)
    while abs(m - μ) > ϵ || abs(s - σ) > ϵ
        v .-= (m - μ) / nqtl
        v .*= σ / s
        gv = Q'v
        m, s = mean(gv), std(gv)
    end
end

"""
    uniq(ixy::T, oxy::T) where T <: AbstractString

Make SNP of a (founder) population unique.
"""
function uniq(ixy::T, oxy::T) where {T<:AbstractString}
    hdr, (nlc, nhp) = XY.header(ixy), XY.dim(ixy)
    hdr.u == 0 || error("Not a SNP file")
    type = XY._type(hdr.type)
    if nhp > typemax(Int32) - 2
        error("Too many ID in the file")
    elseif nhp > typemax(Int16) - 2
        type = Int32
    else
        nhp > typemax(Int8) - 2
        type = Int16
    end
    hdr.type, hdr.u = XY._type(type), 1
    gt = Mmap.mmap(ixy, Matrix{Int8}, (nlc, nhp), 24)
    v = zeros(type, nlc) # for type safety
    open(oxy, "w") do io
        write(io, Ref(hdr), [nlc, nhp])
        for i = 1:nhp
            copyto!(v, gt[:, i])
            v .+= 2(i - 1)
            write(io, v)
        end
    end
end

"""
    sample_ts(bdir::AbstractString, tdir::AbstractString, nid::Int, nchp::Int, nref::Int, trts::Trait...)
Sample loci for a cattle population founder. The founder population is stored in
`tdir` with `nid` ID, `nchp` chip SNP, and `nref` reference loci. The traits
`trts` are used to sample the QTLs. The founder population is stored in `tdir`.
"""
function sample_ts(
    bdir::AbstractString,
    tdir::AbstractString,
    nid::Int,
    nchp::Int,
    nref::Int,
    trts::Trait...,
)
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
    ped = DataFrame(
        id = Int32.(1:nid),
        sire = Int32(0),
        dam = Int32(0),
        sex = rand(Int8.(0:1), nid),
        grt = Int16(1),
    )
    snps = mmap("$tdir/$name.xy", Matrix{Int8}, (size(lmp, 1), 2nid), 24)
    for t in trts
        qtl = lmp[!, t.name]
        Q = snps[qtl, 1:2:end] + snps[qtl, 2:2:end]
        D = Q .== 1
        a, d = rand(t.da, t.nQTL), rand(t.dd, t.nQTL)
        norm_v(Q, a)
        norm_v(D, d, σ = sqrt(t.vd))
        lmp[!, "$(t.name)_a"] .= 0.0
        lmp[!, "$(t.name)_d"] .= 0.0
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

"""
    add_bool!(lmp::DataFrame, name::AbstractString, pool::AbstractVector{Int}, nspl::Int)
Sample `nspl` loci from `pool` and add them to `lmp` with name `name`. The loci
are added to ``Set`` `out`. Return the loci sampled.
"""
function add_bool!(
    lmp::DataFrame,
    name::AbstractString,
    pool::AbstractVector{Int},
    nspl::Int,
)
    loci = sort(shuffle(pool)[1:nspl])
    lmp[!, name] .= false
    lmp[!, name][loci] .= true
    loci
end

"""
    sample_hps(fxy::T, fmp::T, tdir::T, nid::Int,) where {T<:AbstractString}
Sample `2nid` haplotypes from file `fxy`, and `fmp`. Returns the haplotype
matrix and their linkage map.
"""
function sample_hps(fxy::T, fmp::T, nid::Int) where {T<:AbstractString}
    isfile(fxy) && isfile(fmp) || error("Files $fxy or $fmp not found")
    hdr = XY.header(fxy)
    hdr.u ∈ (0, 1) && hdr.major == 0 || error("$fxy is not about haplotypes")
    _, nhp = XY.dim(fxy)
    aid = nhp ÷ 2
    nid ≤ aid || error("No enough ID ($aid) to sample")
    gt = XY.mapit(fxy)
    @info "  - Sampling ID"
    hps = begin
        tmp = shuffle(1:aid)[1:nid]
        sort([2tmp .- 1; 2tmp])
    end
    sgt = Int8.(gt[:, hps]) # as fxy stores a BitArray
    lmp = deserialize(fmp)
    # Randomly change the allele names
    ref, alt = lmp.ref, lmp.alt
    for ilc ∈ axes(sgt, 1)
        if rand(Bool)
            sgt[ilc, :] .= 1 .- sgt[ilc, :]
            ref[ilc], alt[ilc] = alt[ilc], ref[ilc]
        end
    end
    lmp.frq = vec(mean(sgt, dims = 2))
    sgt, lmp
end

"""
    sample_loci(
    lmp::DataFrame,
    pool::Vector{Int},
    nlc::Int,
    name::AbstractString,
    pot::Vector{Int},
    )
Sample `nlc` loci from `pool` and add them to `lmp` with name `name`. The `loci`
sampled are added to `pot`, and returned.
"""
function sample_loci(
    lmp::DataFrame,
    pool::Vector{Int},
    nlc::Int,
    name::AbstractString,
    pot::Vector{Int},
)
    length(pool) ≥ nlc || error("Number of SNP required is larger than available")
    lmp[!, name] .= false
    loci = sort(shuffle(pool)[1:nlc])
    lmp[!, name][loci] .= true
    append!(pot, loci)
    return loci
end

"""
    sample_xy(
        fxy::AbstractString,
        fmp::AbstractString,
        tdir::AbstractString,
        nid::Int,
        maf::Float64,
        nchp::Int,
        nref::Int,
        trts::Trait...;
        replace::Bool = false,
    )
Sample a founder from a SNP file `fxy` and a linkage map `fmp` into `tdir`. The
founder has `nid` ID, `maf` (exclusive) minor allele frequency, `nchp` chip SNP,
`nref` reference loci, and `trts` traits.

Note: ID are sampled first. It is better to read the genotype sequentially on a
spinning disk. Hence the SNP file needs to be loci majored.
"""
function sample_xy(
    fxy::AbstractString,
    fmp::AbstractString,
    tdir::AbstractString,
    nid::Int,
    maf::Float64,
    nchp::Int,
    nref::Int,
    trts::Trait...;
    replace::Bool = false,
)
    0.0 ≤ maf ≤ 0.5 || error("Minor allele frequency out of range")
    sgt, lmp = sample_hps(fxy, fmp, nid)
    alc = 1:size(lmp, 1) # all loci
    pool = maf .< lmp.frq .< 1.0 - maf # loci to sample
    pot = Int[]          # selected loci
    @info "  - Sampling SNP"
    sample_loci(lmp, alc[pool], nchp, "chip", pot)
    if !replace 
        pool = pool .&& .!lmp.chip
    end
    sample_loci(lmp, alc[pool], nref, "dark", pot)
    if !replace
        pool = pool .&& .!lmp.dark
    end

    ped = DataFrame(
        id = Int32.(1:nid),
        sire = Int32(0),
        dam = Int32(0),
        sex = shuffle([zeros(Int8, nid ÷ 2); ones(Int8, nid - nid ÷ 2)]),
        grt = Int16(0),
    )

    for t in trts
        qtl = sample_loci(lmp, alc[pool], t.nQTL, t.name, pot)
        if !replace
            pool = pool .&& .!lmp[!, t.name]
        end
        Q = sgt[qtl, 1:2:end] + sgt[qtl, 2:2:end]
        D = Q .== 1
        a, d = rand(t.da, t.nQTL), rand(t.dd, t.nQTL)
        norm_v(Q, a)
        norm_v(D, d, σ = sqrt(t.vd))
        lmp[!, "$(t.name)_a"] .= 0.0
        lmp[!, "$(t.name)_d"] .= 0.0
        lmp[!, "$(t.name)_a"][qtl] .= a
        lmp[!, "$(t.name)_d"][qtl] .= d
        ped[!, "tbv_$(t.name)"] = Q'a
        ped[!, "ebv_$(t.name)"] .= -1e8
        ped[!, "gt_$(t.name)"] = ped[!, "tbv_$(t.name)"] + D'd
    end
    olc = sort(unique(pot))
    open("$tdir/founder.xy", "w") do io
        hdr = XY.header(major = 1)
        write(io, Ref(hdr), [length(olc), 2nid])
        write(io, sgt[olc, :])
    end
    serialize("$tdir/founder.ped", ped)
    serialize("$tdir/founder.lmp", lmp[olc, :])
end
