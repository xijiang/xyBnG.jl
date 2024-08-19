module Founder

using DataFrames
using Distributions
using LinearAlgebra
using Mmap
using Random
using Serialization
using xyBnG.Conn.TS
using xyBnG.XY
import xyBnG.xyTypes: Cattle, Species, Trait

"""
    ts_base(pop::Cattle, dir::AbstractString)

Simulate a cattle population of name `pop.name`, and `pop.nid` ID in `dir`.

Note:
- This is a simulation with the coancestor/backward simulator msprime.
- Needs to have `tskit`, `msprime`, `scipy` and `stdpopsim` installed.
"""
function ts_base(pop::Cattle, dir::AbstractString)
    '~' ∈ dir && error("Character '~' is forbidden in dir")
    isdir(dir) || mkpath(dir)
    @info "Simulating a cattle population of name $(pop.name), and $(pop.nid) ID in $dir
      Note: This is a simulation with the coancestor/backward simulator msprime."
    Threads.@threads for chr = 1:29
        print(" $chr")
        cmd = pipeline(
            `stdpopsim BosTau -c $chr -o $dir/$chr.ts
 -d HolsteinFriesian_1M13 Holstein_Friesian:$(pop.nid)`,
            stderr = devnull,
        )
        run(cmd)
    end
    open("$dir/desc.txt", "w") do io
        println(io, pop.name)
        println(io, pop.nid)
    end
end

"""
    ts_base(pop::Species, dir::AbstractString)

Default function to catch unsupported species.
"""
function ts_base(pop::Species, dir::AbstractString)
    @info "This species is not supported yet." # or as below
    # https://popsim-consortium.github.io/stdpopsim-docs/stable/api.html#sec-api-generic-models
end

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
    snps = mmap("$tdir/$name.xy", Matrix{Int8}, (nrow(lmp), 2nid), 24)
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
    macs_base(nid::Int, dir::AbstractString)
Simulate a populatoin with `MaCS` into `dir`. The parameters are adapted from
``https://academic.oup.com/g3journal/article/2/4/425/6026056?login=true#supplementary-data``.
This simulation used the `Ne = 100` one.

Note, the command `macs` needs to be in a searchable path with no space
character in it. Files like `chr.1`, `log.1` are generated.
"""
function macs_base(nid::Int, dir::AbstractString)
    isdir(dir) || mkpath(dir)
    Ne = 100
    μ = 2.5e-8 * (4Ne)
    # [Chromosome length in bp](https://www.ncbi.nlm.nih.gov/projects/r_gencoll/ftp_service/nph-gc-ftp-service.cgi/?HistoryId=MCID_642d5d40ceff2e2c64293c60&QueryKey=1&ReleaseType=RefSeq&FileType=GENOME_FASTA&Flat=true)
    chr = [
        158534110,
        136231102,
        121005158,
        120000601,
        120089316,
        117806340,
        110682743,
        113319770,
        105454467,
        103308737,
        106982474,
        87216183,
        83472345,
        82403003,
        85007780,
        81013979,
        73167244,
        65820629,
        63449741,
        71974595,
        69862954,
        60773035,
        52498615,
        62317253,
        42350435,
        51992305,
        45612108,
        45940150,
        51098607,
    ]
    r = 1e-8 * 4Ne
    # Scaled time
    tm =
        [
            10,
            25,
            50,
            100,
            200,
            300,
            400,
            500,
            600,
            700,
            800,
            900,
            1000,
            2000,
            3000,
            4000,
            5000,
            6000,
            7000,
            8000,
            9000,
            10000,
            20000,
            40000,
            60000,
            80000,
            100000,
            200000,
            400000,
            600000,
            800000,
        ] / (4Ne)
    # Scaled population size
    ps =
        [
            175,
            200,
            350,
            500,
            700,
            820,
            850,
            900,
            1000,
            1100,
            1275,
            1300,
            1200,
            2000,
            2500,
            3000,
            3200,
            3500,
            3800,
            4000,
            4200,
            4500,
            5456,
            7367,
            9278,
            11190,
            13101,
            22658,
            41772,
            60886,
            80000,
        ] / Ne
    eN = ""
    for i = 1:length(tm)
        eN *= " -eN $(tm[i]) $(ps[i])"
    end
    macs = Sys.which("macs")
    (isnothing(macs) || any(isspace.(collect(macs)))) && error("Command `macs` error")
    @info "  - Simulating a cattle population with MaCS into $dir"
    Threads.@threads for i = 1:length(chr)
        print(" $i")
        cmd = "$macs $(2nid) $(chr[i]) -t $μ -r $r" * eN
        cmd = Cmd(convert(Vector{String}, split(cmd)))
        run(pipeline(cmd, stdout = "$dir/chr.$i", stderr = "$dir/log.$i"))
    end
    println()
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
    sample_xy(fxy::AbstractString, fmp::AbstractString, tdir::AbstractString, nid::Int, maf::Float64, nchp::Int, nref::Int, trts::Trait...)
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
    trts::Trait...,
)
    isfile(fxy) && isfile(fmp) || error("Files $fxy or $fmp not found")
    isdir(tdir) || mkpath(tdir)
    msnp = maximum((nchp, nref))
    for t in trts
        msnp < t.nQTL && (msnp = t.nQTL)
    end
    gt = XY.mapit(fxy)
    nlc, nhp = size(gt)
    aid = nhp ÷ 2 # available ID
    nid ≤ aid || error("Number of ID to sampled ($nid) is larger than available ($aid)")
    lmp = deserialize(fmp)
    nrow(lmp) == nlc || error("Linkage map and SNP file do not match")
    0.0 ≤ maf ≤ 0.5 || error("Minor allele frequency out of range")

    @info "  - Sampling ID"
    hps = begin
        tmp = shuffle(1:aid)[1:nid]
        sort([2tmp .- 1; 2tmp])
    end
    sgt = Int8.(gt[:, hps])
    lmp.frq = vec(mean(sgt, dims = 2)) # mean(gt, dims = 2) is a matrix !

    @info "  - Sampling SNP"
    vld = (lmp.frq .> maf) .&& (lmp.frq .< 1.0 - maf)
    msnp ≤ sum(vld) || error("Number of SNP required is larger than available")
    pool, out = (1:nlc)[vld], Set{Int}()
    chip = add_bool!(lmp, "chip", pool, nchp)
    out = out ∪ chip
    dark = add_bool!(lmp, "dark", pool, nref)
    out = out ∪ dark
    for t in trts
        qtl = add_bool!(lmp, t.name, pool, t.nQTL)
        out = out ∪ qtl
    end
    out = sort(collect(out))
    @info "  - Saving the sample to $tdir, with allele name randomly named"
    lmp = copy(view(lmp, out, :))
    ogt = view(sgt, out, :)
    frq, ref, alt = lmp.frq, lmp.ref, lmp.alt
    open("$tdir/founder.xy", "w") do io
        hdr = XY.header()
        write(io, Ref(hdr), [length(out), 2nid])
        for i ∈ 1:size(ogt, 1)
            if rand(Bool)
                write(io, ogt[i, :])
            else
                write(io, 1 .- ogt[i, :])
                frq[i] = 1.0 .- frq[i]
                ref[i], alt[i] = alt[i], ref[i]
            end
        end
    end
    print("    Genotypes")
    ped = DataFrame(
        id = Int32.(1:nid),
        sire = Int32(0),
        dam = Int32(0),
        sex = shuffle([zeros(Int8, nid ÷ 2); ones(Int8, nid - nid ÷ 2)]),
        grt = Int16(0),
    )

    for t in trts
        qtl = lmp[!, t.name]
        Q = ogt[qtl, 1:2:end] + ogt[qtl, 2:2:end]
        D = Q .== 1
        a = rand(t.da, t.nQTL)
        d = rand(t.dd, t.nQTL)
        norm_v(Q, a)
        norm_v(D, d, σ = sqrt(t.vd))
        lmp[!, "$(t.name)_a"] .= 0.0
        lmp[!, "$(t.name)_d"] .= 0.0
        lmp[!, "$(t.name)_a"][qtl] = a
        lmp[!, "$(t.name)_d"][qtl] = d
        ped[!, "tbv_$(t.name)"] = Q'a
        ped[!, "ebv_$(t.name)"] .= -1e8
        ped[!, "gt_$(t.name)"] = ped[!, "tbv_$(t.name)"] + D'd
    end
    serialize("$tdir/founder.ped", ped)
    serialize("$tdir/founder.lmp", lmp)
    println(", pedigree, and linkage map")
end

end # module Founder
