"""
    Sum
The ``Sum`` module provides functions to calculate summary statistics of the
simulations. It has majorly two parts. One is to summarize each repeat of a
simulation. The other contains a few notebooks for the presentation of the
summaries.
"""
module Sum

using DataFrames
using LinearAlgebra
using Serialization
using Statistics
using xyBnG.XY
using xyBnG.RS
import xyBnG.xyTypes: Trait

include("sum/savepar.jl")

"""
    tibd(mat::AbstractMatrix)
Calculate the inbreeding coefficient of each individual in with haplotype matrix
`mat`.
"""
function tibd(mat::AbstractMatrix)
    nlc, nhp = size(mat)
    F = zeros(nhp ÷ 2)
    Threads.@threads for i = 1:nhp÷2
        F[i] = sum(mat[:, 2i-1] .== mat[:, 2i]) / nlc
    end
    F
end

"""
    FFCV(mat::AbstractMatrix, grt::AbstractVector, eff::AbstractVector{Float64})
This function returns the allele frequencies of haplotype matrix `mat`,
population selection `Ceiling` and `Floor`, and genic variance by generation as
indicated by `grt` and QTL additive effects `eff`.
"""
function FFCV(mat::AbstractMatrix, grt::AbstractVector, eff::AbstractVector{Float64})
    nlc = size(mat, 1)
    ug = unique(grt) # unique generation numbers
    ng = length(ug)
    frq = zeros(Int16, nlc, ng)
    flr, clg, vgn = zeros(ng), zeros(ng), zeros(ng)
    a2 = 2eff .* eff  # => 2a²
    pp, nn = eff .> 0, eff .< 0
    aa, bb = 2sum(eff[nn]), 2sum(eff[pp])
    for i = 1:ng
        chp = grt .== ug[i]  # current haplotypes
        sub = view(mat, :, chp)
        copyto!(view(frq, :, i), sum(sub, dims = 2))
        nhp = sum(chp)
        p = frq[:, i] / nhp
        q = 1 .- p
        vgn[i] = sum(p .* q .* a2) # genic variance = ∑2pqa²
        # raising of the floor
        t = nn .&& frq[:, i] .== 0 # if negative and fixed at 0
        flr[i] = aa - 2sum(eff[t])
        t = pp .&& frq[:, i] .== nhp
        flr[i] += 2sum(eff[t])
        # lowering of the ceiling
        t = pp .&& frq[:, i] .== 0
        clg[i] = bb - 2sum(eff[t])
        t = nn .&& frq[:, i] .== nhp
        clg[i] += 2sum(eff[t])
    end
    frq, flr, clg, vgn
end

"""
    snphet(q::AbstractVector{Float64})
Calculate the heterozygosity of a diallelic locus with allele frequency `q`.
"""
function snphet(q::AbstractVector{Float64})
    H = q .* q + (1 .- q) .* (1 .- q) # homozygosity
    1 - mean(H) # heterozygosity
end

function dsnphet(q::AbstractVector{Float64})
    H = q .* q + (1 .- q) .* (1 .- q) # homozygosity
    1 .- H # heterozygosity
end

"""
    xysum(ped::DataFrame, xy::AbstractString, lmp::DataFrame, trait::Trait)
Summarize the simulation results in `xy` and `ped` files.

It is supposed that filename `xy` is of the pattern `repeat-scheme.xy`. Repeat
number and scheme are to be stored in the result DataFrame.

## Columns in the output DataFrame
- `repeat`: repeat number
- `scheme`: scheme of the simulation
- `grt`: generation number
- `nid`: number of individuals in the generation
- `nsire`: number of sires
- `ndam`: number of dams
- `mtbv`: mean TBV
- `vtbv`: genetic variance
- `mebv`: mean EBV
- `rtbv`: TBV-EBV correlation
- `fibd`: mean inbreeding coefficient by IBD
- `fped`: mean inbreeding coefficient by pedigree
- `xqtl`: number of QTL fixed
- `xchip`: number of chip loci fixed
- `xref`: number of reference loci fixed
- `xfq`: number of favorable QTL fixed
- `xuq`: number of unfavorable QTL fixed
- `genicv`: genic variance
- `floor`: floor of the population, i.e., the lowest mean TBV
- `ceiling`: ceiling of the population
- `covdq`: covariance between q₀ and qₜ
- `covdq2`: covariance between q₀ corrected and qₜ
- `fhet`: inbreeding by homozygosity
- `fdrift`: inbreeding by drift
"""
function xysum(ped::DataFrame, xy::AbstractString, lmp::DataFrame, trait::Trait)
    haps = XY.mapit(xy)
    ped.iF = tibd(view(haps, lmp.chip, :))
    frq, flr, clg, vgn =
        FFCV(isodd.(haps), repeat(ped.grt, inner = 2), lmp[!, trait.name*"_a"])
    A = RS.nrm(ped)
    ped.aF = diag(A) .- 1
    trt = trait.name
    rpt, scheme = split(split(basename(xy), '.')[1], '-')
    rpt = parse(Int, rpt)
    ss = combine(
        groupby(ped, :grt),
        :id => length => :nid,
        :sire => length ∘ unique => :nsire,
        :dam => length ∘ unique => :ndam,
        "tbv_" * trt => mean => :mtbv,
        "tbv_" * trt => var => :vtbv, # genetic variance
        "ebv_" * trt => mean => :mebv,
        ["tbv_" * trt, "ebv_" * trt] => cor => :rtbv,
        :iF => mean => :fibd,
        :aF => mean => :fped,
    )
    insertcols!(ss, 1, :repeat => rpt, :scheme => scheme)
    ng = size(frq, 2)
    xqtl = zeros(Int, ng)  # number of QTL fixed
    xchip = zeros(Int, ng) # number of chip loci fixed
    xref = zeros(Int, ng)  # number of reference loci fixed
    xfq = zeros(Int, ng)   # number of favorable QTL fixed
    xuq = zeros(Int, ng)   # number of unfavorable QTL fixed
    q = zeros(size(frq))
    for i = 1:ng
        chip = view(frq, lmp.chip, i)
        xchip[i] = sum(chip .== 0) + sum(chip .== 2ss.nid[i])
        dark = view(frq, lmp.dark, i)
        xref[i] = sum(dark .== 0) + sum(dark .== 2ss.nid[i])
        qtl = view(frq, lmp[!, trait.name], i)
        xqtl[i] = sum(qtl .== 0) + sum(qtl .== 2ss.nid[i])
        # favorable QTL fixed
        tt = lmp[!, trait.name] .&& lmp[!, trait.name*"_a"] .> 0
        tt = tt .&& frq[:, i] .== 2ss.nid[i]
        xfq[i] = sum(tt)
        tt = lmp[!, trait.name] .&& lmp[!, trait.name*"_a"] .< 0
        tt = tt .&& frq[:, i] .== 0
        xfq[i] += sum(tt)
        tt = lmp[!, trait.name] .&& lmp[!, trait.name*"_a"] .< 0
        tt = tt .&& frq[:, i] .== 2ss.nid[i]
        xuq[i] = sum(tt)
        tt = lmp[!, trait.name] .&& lmp[!, trait.name*"_a"] .> 0
        tt = tt .&& frq[:, i] .== 0
        xuq[i] += sum(tt)
        q[:, i] = frq[:, i] / 2ss.nid[i]
    end
    insertcols!(ss, :xqtl => xqtl, :xchip => xchip, :xref => xref, :xfq => xfq, :xuq => xuq)
    insertcols!(ss, :genicv => vgn, :floor => flr, :ceiling => clg)
    covdq = zeros(ng)      # covariance between q₀ and qᵢ
    covdq2 = zeros(ng)     # covariance between q₀ corrected and qᵢ
    fhet = zeros(ng)       # Inbreeding by heterozygosity
    fhet2 = zeros(ng)      # Inbreeding by heterozygosity
    fdrift = zeros(ng)     # Inbreeding by drift
    fdrift2 = zeros(ng)    # Inbreeding by drift

    q0 = q[:, 1]
    H0, H0d = snphet(q0), dsnphet(q0)
    loci = lmp.dark .&& 0 .< frq[:, 1] .< 2ss.nid[1]
    factor = sqrt.(q0[loci] .* (1 .- q0[loci]))
    qa = q0[loci] ./ factor
    qb = (q0[loci] .- 0.5) ./ factor
    for i = 1:ng
        δq = (q[loci, i] .- q0[loci]) ./ factor
        covdq[i] = cov(δq, qa)
        covdq2[i] = cov(δq, qb)
        Ht, Htd = snphet(q[:, i]), dsnphet(q[:, i])
        fhet[i] = (H0 - Ht) / H0
        fhet2[i] = mean((H0d - Htd) ./ H0d)
        fdrift[i] = 2var(q[:, i] - q0) / H0
        fdrift2[i] = mean(2(q[:, i] - q0) .^ 2 ./ H0d)
    end
    insertcols!(
        ss,
        :covdq => covdq,
        :covdq2 => covdq2,
        :fhet => fhet,
        :fhet2 => fhet2,
        :fdrift => fdrift,
        :fdrift2 => fdrift2,
    )
    ss
end

"""
    xysum(ped::AbstractString, xy::AbstractString, lmp::AbstractString, trait::Trait)
Summarize the simulation results in file `ped`, `xy` and `lmp` files.
"""
function xysum(ped::AbstractString, xy::AbstractString, lmp::AbstractString, trait::Trait)
    pd = deserialize(ped)
    lp = deserialize(lmp)
    xysum(pd, xy, lp, trait)
end

"""
    xysum(ped::AbstractString, xy::AbstractString, lmp::DataFrame, trait::Trait, ssg::Int)
Summarize the simulation results in `xy` and `ped` files. `lmp` is in memory already.
"""
function xysum(ped::AbstractString, xy::AbstractString, lmp::DataFrame, trait::Trait)
    pd = deserialize(ped)
    xysum(pd, xy, lmp, trait)
end

"""
    savesum(file::AbstractString, df::DataFrame)
Save the summary DataFrame `df` to file `file`. If the file exists, append the
DataFrame to the existing one.
"""
function savesum(file::AbstractString, df::DataFrame)
    if isfile(file)
        tsm = deserialize(file)
        append!(tsm, df)
        serialize(file, tsm)
    else
        serialize(file, df)
    end
end

"""
    cormat(ma::AbstractMatrix, mb::AbstractMatrix)
Calculate the correlation between the off-diagonal elements of two matrices `ma`
and `mb`.
"""
function cormat(ma::AbstractMatrix, mb::AbstractMatrix)
    size(ma) == size(mb) || error("Matrix size not match")
    isequal(size(ma)...) || error("Not square matrices")
    va, vb = Float64[], Float64[]
    for i ∈ 1:size(ma, 1)-1
        append!(va, ma[i+1:end, i])
        append!(vb, mb[i+1:end, i])
    end
    cor(va, vb)
end

"""
    cormat(dir::AbstractString)
Calculate the correlation between the off-diagonal elements of the `A`, `G`, and
`IBD` matrices. G and IBD are calculated with the reference loci indicated in
the linkage map file.
"""
function cormat(dir::AbstractString)
    isfile("$dir/summary.ser") || error("No summary file found in $dir")
    ss = deserialize("$dir/summary.ser")
    cag, cai, cgi, npd = Float64[], Float64[], Float64[], ndigits(ss.repeat[end])
    for (irpt, cskm) ∈ eachrow(unique(select(ss, :repeat, :scheme)))
        @info "Processing repeat $irpt and scheme $cskm"
        tag = lpad(irpt, npd, '0')

        # with linkage map
        lmp = deserialize("$dir/$tag-founder.lmp")
        ref = lmp.dark
        frq = lmp.frq[ref]

        # with pedigree
        ped = deserialize("$dir/$tag-$cskm.ped")
        grt = repeat(ped.grt, inner = 2)
        ugt = unique(grt)
        ngt = length(ugt)
        A = RS.nrm(ped)

        # with haplotypes
        hps = XY.mapit("$dir/$tag-$cskm.xy")
        for igt ∈ 1:ngt
            usnp = view(hps, ref, grt .== ugt[igt])
            D = RS.irm(usnp)
            gt = isodd.(usnp[:, 1:2:end]) + isodd.(usnp[:, 2:2:end])
            G = RS.grm(gt, p = frq)
            id = ped.grt .== ugt[igt]
            Aid = view(A, id, id)
            push!(cag, cormat(Aid, G))
            push!(cai, cormat(Aid, D))
            push!(cgi, cormat(G, D))
        end
    end
    ss.cag, ss.cai, ss.cgi = cag, cai, cgi
    serialize("$dir/summary.ser", ss)
end

"""
    resum(dir::AbstractString)
This is an amendment to the `xysum` function in directory `dir`, where the data
are kept. This is needed when the formulae of some indices are changed.
"""
function resum(dir::AbstractString, trait::Trait)
    isfile("$dir/summary.ser") || error("No summary file found in $dir")
    ss = deserialize("$dir/summary.ser")
    mv("$dir/summary.ser", "$dir/summary.bak", force = true)
    for (irpt, cskm) ∈ eachrow(unique(select(ss, :repeat, :scheme)))
        @info "Processing repeat $irpt and scheme $cskm"
        tag = lpad(irpt, ndigits(ss.repeat[end]), '0')
        ped = deserialize("$dir/$tag-$cskm.ped")
        lmp = deserialize("$dir/$tag-founder.lmp")
        xy = "$dir/$tag-$cskm.xy"
        ss = xysum(ped, xy, lmp, trait)
        savesum("$dir/summary.ser", ss)
    end
end

"""
    fgrm(dir::AbstractString)
Update `summary.ser` with `F_grm`, which are calculated from the `xy` files in
directory `dir`. GRM is calculated with the chip loci.
"""
function fgrm(dir::AbstractString)
    isfile("$dir/summary.ser") || error("No summary file found in $dir")
    ss = deserialize("$dir/summary.ser")
    F, npd = Float64[], ndigits(ss.repeat[end])
    for (irpt, cskm) ∈ eachrow(unique(select(ss, :repeat, :scheme)))
        @info "Processing repeat $irpt and scheme $cskm"
        tag = lpad(irpt, npd, '0')

        # with linkage map
        lmp = deserialize("$dir/$tag-founder.lmp")
        chp = lmp.chip
        frq = lmp.frq[chp]

        # with pedigree
        ped = deserialize("$dir/$tag-$cskm.ped")
        grt = repeat(ped.grt, inner = 2)
        ugt = unique(grt)
        ngt = length(ugt)
        tf = zeros(ngt)

        # with haplotypes
        hps = isodd.(XY.mapit("$dir/$tag-$cskm.xy"))
        for igt ∈ 1:ngt
            snps = view(hps, chp, grt .== ugt[igt])
            gt = snps[:, 1:2:end] .+ snps[:, 2:2:end]
            G = RS.grm(gt, p = frq)
            tf[igt] = mean(diag(G) .- 1)
        end
        append!(F, tf)
    end
    ss.fgrm = F
    serialize("$dir/summary.ser", ss)
end

end # module Sum
