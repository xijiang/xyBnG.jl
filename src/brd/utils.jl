"""
    sumMap(lmp; cM = 1e6)

Summary of a linkage map, which includes 3 columns, :chr, :pos(in bp), and :frq.
A DataFrame of 4-column for each chromosome is returned:
- numbering (1, 2, ...)
- length in Morgen
- number of loci
- beginning number of the first locus
"""
function sumMap(lmp; cM=1e6)
    lms = DataFrame(chr=Int8[],
          len=Float64[], # as λ
          nlc=Int[],
          bgn=Int[])
    bgn = 1
    for grp in groupby(lmp, :chr)
        chr = first(grp).chr
        len = (last(grp).pos - first(grp).pos) / (100cM)
        nlc = size(grp, 1)
        push!(lms, (chr, len, nlc, bgn))
        bgn += nlc
    end
    return lms
end

"""
    crossovers(lms)

Give a linkage map summary `lms`, which can be from `sumMap`, this function
return a vector of crossover points along the whole genome.

DataFrame `lms` has four columns, chromosome number, its length in Morgen,
number of loci it contains, and the beginning number of its first locus in the
whole genome.

The first number is `1`, or `2`, the starting alternative haplotype. The vector
is then start from locus `1`, and loci that cross-over happens.  A
``cross-over`` may also happen at the first locus of a chromsome. Otherwise, the
segment continues with the previous haplotype.  This eases the segment copy.  At
the beginning of a ``cross-over`` happens on `rand([1:2])`.  The number of
cross-over occurred on a chromosome follows a Poisson distribution.

The cross-over location is then sampled from a uniform distribution. It can can
also be U-shaped about the centromere, which might be implemented later.

The chromosome, or linkage group, should be in the order of the genotype file.
Or, unpredictable errors can happen.
"""
function crossovers(lms)
    pts = [rand(1:2)]           # crossover points
    for (_, λ, nlc, bgn) in eachrow(lms)
        bgn > 1 && rand(false:true) && push!(pts, bgn)
        nc = rand(Poisson(λ))
        append!(pts, rand(1:nlc, nc) .+ (bgn - 1))
    end
    push!(pts, last(lms).nlc - 1 + last(lms.bgn))
    pts
end

"""
    gamete(prt, hap, lms)

Generate a gamete `hap`, a **vector** view of a child's haplotype,
from `prt`, a view of a parents genotypes,
according crossovers generated from a linkage map summary, `lms`.
"""
function gamete(prt, hap, lms)
    cvs = crossovers(lms)
    h, l = cvs[1], 1            # starting
    for cv in cvs[2:end]
        copyto!(view(hap, l:cv), view(prt, l:cv, h))
        l = cv + 1
        h = 3 - h
    end
end

"""
    pair!(sirs, dams, noff::Int, ped::DataFrame, mate::Symbol)

Pair the `sirs` and `dams` to produce `noff` offspring.  The `mate` can be
`:random` or `:balanced`.  The new families will be appended to `ped`.
"""
function pair!(sirs, dams, noff::Int, ped::DataFrame, mate::Symbol)
    nid = ped.id[end]
    ids = Int32.(nid+1:nid+noff)
    sex = rand(Int8[0, 1], noff)
    nhsb = Int(ceil(noff / length(sirs)))
    nfsb = Int(ceil(noff / length(dams)))
    pa, ma = begin
        if mate == :random
            repeat(sort(sirs), outer = nhsb)[1:noff], repeat(shuffle(dams), outer = nfsb)[1:noff]
        else
            repeat(sort(sirs), inner = nhsb)[1:noff], repeat(shuffle(dams), inner = nfsb)[1:noff]
        end
        #sortslices([pa ma], dims=1, by=x -> (x[1], x[2]))
    end
    df = DataFrame(id = ids, pa = pa, ma = ma, sex = sex, grt = ped.grt[end] + 1)
    append!(ped, df, cols = :subset)
end

"""
    function drop(pg::AbstractArray, og::AbstractArray, pm, lms)
Drop haplotypes `pg` of parents into `og`, their offspring genotypes.
Parents of each offspring are defined in `pm`, which are rows of ``pa ma``.
Linkage map summary `lms` is from `summap`.

!!! ``Caution``:
- Merged data matrix from `MaCS` is `n-ID × n-loci`. Treat it with care.
- It is supposed all the genes to drop from are in `pg`.
- This function will be removed in the future.
"""
function drop(pg::AbstractArray, og::AbstractArray, pm, lms)
    nf = size(pm)[1]
    Threads.@threads for id in 1:nf
        ip = pm[id, 1]
        pa = view(pg, :, 2ip-1:2ip)
        zi = vec(view(og, :, 2id - 1))
        gamete(pa, zi, lms)
    end
    Threads.@threads for id in 1:nf
        im = pm[id, 2]
        ma = view(pg, :, 2im-1:2im)
        zi = vec(view(og, :, 2id))
        gamete(ma, zi, lms)
    end
end
