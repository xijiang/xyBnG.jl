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
