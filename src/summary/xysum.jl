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
- `covdq1`: covariance between q₀ and qₜ, dark, method 1
- `covdq2`: covariance between q₀ corrected and qₜ, dark method 2
- `covdq3`: covariance between q₀ corrected and qₜ, chip method 2
- `fhet1`: inbreeding by homozygosity, dark, method 1
- `fhet2`: inbreeding by homozygosity, dark, method 2
- `fhet3`: inbreeding by homozygosity, chip, method 2
- `fdrift1`: inbreeding by drift, dark, method 1
- `fdrift2`: inbreeding by drift, dark, method 2
- `fdrift3`: inbreeding by drift, chip, method 2
"""
function xysum(ped::DataFrame, xy::AbstractString, lmp::DataFrame, trait::Trait)
    haps = XY.mapit(xy)
    #ped.iF = tibd(view(haps, lmp.chip, :))
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
        #:iF => mean => :fibd,
        #:aF => mean => :fped,
    )
    # Mean relationship using A matrix
    grt = sort(unique(ped.grt))
    aF = Float64[]
    for i in grt
        rg = ped.grt .== i
        push!(aF, mean(A[rg, rg] - I)/2)
    end
    ss.aF = aF
    # Mean relationship using IBD
    hg = repeat(ped.grt, inner = 2)
    iF = Float64[]
    for i in grt
        rg = hg .== i
        T = irm(view(haps, lmp.chip, rg))
        n = sum(rg)
        tif = (sum(T) - sum(diag(T))) / n / (n - 1)
        push!(iF, tif)
    end
    ss.iF = iF
    insertcols!(ss, 1, :repeat => rpt, :scheme => scheme)
    # Mean IBD relationship

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
    covdq1 = zeros(ng)     # covariance between q₀ and qᵢ
    covdq2 = zeros(ng)     # covariance between q₀ corrected and qᵢ, dark
    covdq3 = zeros(ng)     # covariance between q₀ corrected and qᵢ, chip
    fhet1 = zeros(ng)      # Inbreeding by homozygosity
    fhet2 = zeros(ng)      # Inbreeding by homozygosity, dark SNP
    fhet3 = zeros(ng)      # Inbreeding by homozygosity, chip SNP
    fdrift1 = zeros(ng)    # Inbreeding by drift
    fdrift2 = zeros(ng)    # Inbreeding by drift, dark SNP
    fdrift3 = zeros(ng)    # Inbreeding by drift, chip SNP

    q₀ = q[:, 1]

    dark = lmp.dark .&& 0 .< frq[:, 1] .< 2ss.nid[1]
    chip = lmp.chip .&& 0 .< frq[:, 1] .< 2ss.nid[1]
    H₀, H₀v = snphet(q₀[dark]), v_snphet(q₀)
    d_dark = sqrt.(q₀[dark] .* (1 .- q₀[dark]))
    d_chip = sqrt.(q₀[chip] .* (1 .- q₀[chip]))
    q₁ = q₀[dark] ./ d_dark
    q₂₁ = (q₀[dark] .- 0.5) ./ d_dark
    q₂₂ = (q₀[chip] .- 0.5) ./ d_chip
    q₃ = (q₀[chip] .- 0.5) ./ d_chip
    for i = 1:ng
        δq = (q[dark, i] .- q₀[dark]) ./ d_dark
        covdq1[i] = cov(δq, q₁)
        covdq2[i] = cov(δq, q₂₁)
        δq = (q[chip, i] .- q₀[chip]) ./ d_chip
        covdq3[i] = cov(δq, q₂₂)
        Ht, Htᵥ = snphet(q[dark, i]), v_snphet(q[:, i])
        fhet1[i] = (H₀ - Ht) / H₀
        fhet2[i] = mean((H₀v[dark] - Htᵥ[dark]) ./ H₀v[dark])
        fhet3[i] = mean((H₀v[chip] - Htᵥ[chip]) ./ H₀v[chip])
        fdrift1[i] = 2var(q[dark, i] - q₀[dark]) / H₀
        fdrift2[i] = mean(2(q[dark, i] - q₀[dark]) .^ 2 ./ H₀v[dark])
        fdrift3[i] = mean(2(q[chip, i] - q₀[chip]) .^ 2 ./ H₀v[chip])
    end
    insertcols!(
        ss,
        :covdq1 => covdq1,
        :covdq2 => covdq2,
        :covdq3 => covdq3,
        :fhet1 => fhet1,
        :fhet2 => fhet2,
        :fhet3 => fhet3,
        :fdrift1 => fdrift1,
        :fdrift2 => fdrift2,
        :fdrift3 => fdrift3,
    )
    ss
end
