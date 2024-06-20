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
    Threads.@threads for i in 1:nhp ÷ 2
        F[i] = sum(mat[:, 2i-1] .== mat[:, 2i]) / nlc
    end
    F
end

"""
    gfreq(mat::AbstractMatrix, grt::AbstractVector)
Calculate allele frequencies of `mat` by generation indicator `grt`.
"""
function gfreq(mat::AbstractMatrix, grt::AbstractVector)
    nlc = size(mat, 1)
    ug = unique(grt)
    ng = length(ug)
    frq = zeros(Int16, nlc, ng)
    for i in 1:ng
        sub = view(mat, :, grt .== ug[i])
        frq[:, i] = sum(sub, dims=2)
    end
    frq
end

"""
    xysum(ped::DataFrame, xy::AbstractString, lmp::DataFrame, dsg::Int, trait::Trait)
Summarize the simulation results in `xy` and `ped` files. Arg `dsg` is the
generation number when Directional selection started.

It is supposed that filename `xy` is of the pattern `repeat-scheme.xy`. Repeat
number and scheme are to be stored in the result DataFrame.
"""
function xysum(ped::DataFrame, xy::AbstractString, lmp::DataFrame, dsg::Int, trait::Trait)
    haps = XY.mapit(xy)
    ped.iF = tibd(view(haps, lmp.chip, :))
    frq = gfreq(isodd.(haps), repeat(ped.grt, inner=2))
    A = RS.nrm(ped)
    ped.aF = diag(A) .- 1
    trt = trait.name
    ss = combine(groupby(ped, :grt),
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
    ng = size(frq, 2)
    nfxqtl = zeros(Int16, ng)
    nfxchp = zeros(Int16, ng)
    nfxref = zeros(Int16, ng)
    for i in 1:ng
        qfq = view(frq, lmp[!, trait.name], i)
        nfxqtl[i] = sum(qfq .== 0) + sum(qfq .== 2ss.nid[i])
        nfxchp[i] = sum(view(frq, lmp.chip, i) .== 0) + sum(view(frq, lmp.chip, i) .== 2ss.nid[i])
        nfxref[i] = sum(view(frq, lmp.dark, i) .== 0) + sum(view(frq, lmp.dark, i) .== 2ss.nid[i])
    end
    ss.nfxchp = nfxchp
    ss.nfxref = nfxref
    ss.nfxqtl = nfxqtl
    ss
end

end # module Sum
