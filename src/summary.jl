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
This function returns the allele frequencies of `mat`, population selection
`Ceiling` and `Floor`, and genic variance by generation as indicated by `grt`
and QTL additive effects `eff`.
"""
function FFCV(mat::AbstractMatrix, grt::AbstractVector, eff::AbstractVector{Float64})
    nlc = size(mat, 1)
    ug = unique(grt) # unique generation numbers
    ng = length(ug)
    frq = zeros(Int16, nlc, ng)
    flr, clg, vgn = zeros(ng), zeros(ng), zeros(ng)
    a2 = 2eff .* eff
    pp, nn = eff .> 0, eff .< 0
    aa, bb = sum(eff[nn]), sum(eff[pp])
    for i = 1:ng
        sub = view(mat, :, grt .== ug[i])
        frq[:, i] = sum(sub, dims = 2)
        nhp = sum(grt .== ug[i])
        p = frq[:, i] / nhp
        q = 1 .- p
        vgn[i] = sum(p .* q .* a2) # genic variance
        t = nn .&& frq[:, i] .== nhp
        flr[i] = aa - sum(eff[t])
        t = pp .&& frq[:, i] .== 0
        clg[i] = bb - sum(eff[t])
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

"""
    xysum(ped::DataFrame, xy::AbstractString, lmp::DataFrame, trait::Trait, ssg::Int)
Summarize the simulation results in `xy` and `ped` files. `ssg` is the selection
starting generation.

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

"""
function xysum(ped::DataFrame, xy::AbstractString, lmp::DataFrame, trait::Trait, ssg::Int)
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
    xfq = zeros(Int, ng)    # number of favorable QTL fixed
    xuq = zeros(Int, ng)    # number of unfavorable QTL fixed
    fhet = zeros(ng)        # Inbreeding by heterozygosity
    fdrift = zeros(ng)      # Inbreeding by drift
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
    covdq = zeros(ng)       # covariance between q₀ and qᵢ
    covdq2 = zeros(ng)      # covariance between q₀ corrected and qᵢ

    q0 = q[:, ssg]
    H0 = snphet(q0)
    for i = ssg+1:ng
        loci = frq[:, i] .≠ 2ss.nid[i] .&& frq[:, i] .≠ 0
        factor = sqrt.(q0[loci] .* (1 .- q0[loci]))
        δq = (q[loci, i] .- q0[loci]) ./ factor
        q₀ = q0[loci] ./ factor
        covdq[i] = cov(δq, q₀)
        q₀ = (q0[loci] .- 0.5) ./ factor
        covdq2[i] = cov(δq, q₀)
        Ht = snphet(q[:, i])
        fhet[i] = (H0 - Ht) / H0
        fdrift[i] = 2var(q[:, i] - q0) / H0
    end
    insertcols!(ss, :covdq => covdq, :covdq2 => covdq2, :fhet => fhet, :fdrift => fdrift)
    ss
end

"""
    xysum(ped::AbstractString, xy::AbstractString, lmp::AbstractString, trait::Trait, ssg::Int)
Summarize the simulation results in file `ped`, `xy` and `lmp` files.
"""
function xysum(
    ped::AbstractString,
    xy::AbstractString,
    lmp::AbstractString,
    trait::Trait,
    ssg::Int,
)
    pd = deserialize(ped)
    lp = deserialize(lmp)
    xysum(pd, xy, lp, trait, ssg)
end

"""
    xysum(ped::AbstractString, xy::AbstractString, lmp::DataFrame, trait::Trait, ssg::Int)
Summarize the simulation results in `xy` and `ped` files. `lmp` is in memory already.
"""
function xysum(
    ped::AbstractString,
    xy::AbstractString,
    lmp::DataFrame,
    trait::Trait,
    ssg::Int,
)
    pd = deserialize(ped)
    xysum(pd, xy, lmp, trait, ssg)
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
    corMat(fxy::AbstractString, fpd::AbstractString, fmp::AbstractString)
Calculate the correlation between the off-diagonal elements of the `A`, `G`, and
`IBD` matrices. G and IBD are calculated with the chip loci incidated in the
`DataFrame` stored in the `fmp` file.
"""
function corMat(fxy::AbstractString, fpd::AbstractString, fmp::AbstractString)
    ped = deserialize(fpd)
    nid = size(ped, 1)
    lmp = deserialize(fmp)
    @info "  - Calculating A"
    A = RS.nrm(ped)
    mat = XY.mapit(fxy)
    gt = isodd.(mat[:, 1:2:end]) + isodd.(mat[:, 2:2:end])
    mat = nothing
    @info "  - Calculating G"
    G = RS.grm(gt)
    @info "  - Calculating IBD matrix"
    M = RS.irm(fxy, lmp.chip, 1:nid)
    @info "  - Calculating correlations across all generations"
    x = y = z = x2 = y2 = z2 = xy = xz = yz = 0.0
    for i = 1:nid
        for j = 1:i-1
            x += A[i, j]
            y += G[i, j]
            z += M[i, j]
            x2 += A[i, j]^2
            y2 += G[i, j]^2
            z2 += M[i, j]^2
            xy += A[i, j] * G[i, j]
            xz += A[i, j] * M[i, j]
            yz += G[i, j] * M[i, j]
        end
    end
    n = nid * (nid - 1) / 2
    c1, c2, c3 = (xy - x * y / n) / sqrt((x2 - x^2 / n) * (y2 - y^2 / n)),
    (xz - x * z / n) / sqrt((x2 - x^2 / n) * (z2 - z^2 / n)),
    (yz - y * z / n) / sqrt((y2 - y^2 / n) * (z2 - z^2 / n))
    @info "  - Calculating correlations of the generations"
    id = ped.id[ped.grt.==ped.grt[end]]
    x = y = z = x2 = y2 = z2 = xy = xz = yz = 0.0
    for i in id
        for j in id
            j == i && break
            x += A[i, j]
            y += G[i, j]
            z += M[i, j]
            x2 += A[i, j]^2
            y2 += G[i, j]^2
            z2 += M[i, j]^2
            xy += A[i, j] * G[i, j]
            xz += A[i, j] * M[i, j]
            yz += G[i, j] * M[i, j]
        end
    end
    n = length(id) * (length(id) - 1) / 2
    c4, c5, c6 = (xy - x * y / n) / sqrt((x2 - x^2 / n) * (y2 - y^2 / n)),
    (xz - x * z / n) / sqrt((x2 - x^2 / n) * (z2 - z^2 / n)),
    (yz - y * z / n) / sqrt((y2 - y^2 / n) * (z2 - z^2 / n))
    return [c1, c2, c3, c4, c5, c6]
end

end # module Sum
