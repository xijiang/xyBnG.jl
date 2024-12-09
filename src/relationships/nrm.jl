"""
    nrm(ped::DataFrame; m = 30_000)

Given a pedigree `ped`, this function returns a full numerical relationship
matrix, `A`. This function is better for small pedigrees, and for demonstration
only. The maximal matrix size is thus limited to 30k. One can try to set `m` to
a bigger value if your RAM is enough.
"""
function nrm(ped::DataFrame; m = 30_000)
    N = size(ped, 1)
    N > m && error("Pedigree size ($N > $m) too big")
    A = zeros(N, N) + I(N)
    for (id, sire, dam) in eachrow(select(ped, :id, :sire, :dam))
        sire * dam ≠ 0 && (A[id, id] += 0.5A[sire, dam])
        for jd = 1:id-1
            sire ≠ 0 && (A[id, jd] = 0.5A[jd, sire])
            dam ≠ 0 && (A[id, jd] += 0.5A[jd, dam])
            A[jd, id] = A[id, jd]
        end
    end
    A
end

# below are some recursive functions for kinship calculation
"""
    function kinship(ped, i, j)
---
This function is handy if just to calculate relationship of a few (pairs of) ID.
The first 2 columns of ped must be `pa` and `ma`.
It can also speed up by adding `Thread.@threads` before your pair loop.
"""
function kinship(ped, i, j)
    (i == 0 || j == 0) && return 0
    ipa, ima = ped[i, :]          # used both for below and the last
    i == j && (return 1 + 0.5kinship(ped, ipa, ima))
    if i < j
        jpa, jma = ped[j, :]
        return 0.5(kinship(ped, i, jpa) + kinship(ped, i, jma))
    end
    return 0.5(kinship(ped, j, ipa) + kinship(ped, j, ima))
end

"""
    function kinship(ped, i::Int, j::Int, dic::Dict{Tuple{Int, Int}, Float64})
Recursive kinship calculation with kinship of ID pair `(i, j)`
stored in dictionary `dic`.  The first 2 columns of ped must be `pa` and `ma`.
The memory usage may be bigger than Meuwissen and Luo 1992, or Quaas 1995.
The speed is however standable.
The recursive algorithm is also easy to understand.
"""
function kinship(ped, i::Int, j::Int, dic::Dict{Tuple{Int,Int},Float64})
    (i == 0 || j == 0) && return 0
    ip, im = ped[i, :]
    if i == j
        haskey(dic, (i, i)) || (dic[(i, i)] = 1 + 0.5kinship(ped, ip, im, dic))
        return dic[(i, i)]
    end
    if i < j
        jp, jm = ped[j, :]
        haskey(dic, (i, jp)) || (dic[(i, jp)] = kinship(ped, i, jp, dic))
        haskey(dic, (i, jm)) || (dic[(i, jm)] = kinship(ped, i, jm, dic))
        return 0.5(dic[(i, jp)] + dic[(i, jm)])
    end
    return 0.5(kinship(ped, j, ip, dic) + kinship(ped, j, im, dic))
end

function fkin()
    kin = Dict((Int32(0), Int32(0)) => 0.0)

    function ck(ped, i::T, j::T) where {T<:Integer} # calculate kinship
        i == 0 || j == 0 && return 0.0
        m, n = i ≤ j ? (i, j) : (j, i)
        haskey(kin, (m, n)) && return kin[(m, n)]
        ip, im = ped[i, :]
        if i == j
            kin[(i, j)] = 1 + 0.5ck(ped, ip, im)
        elseif i < j
            jp, jm = ped[j, :]
            kin[(i, j)] = 0.5(ck(ped, i, jp) + ck(ped, i, jm))
        else
            kin[(m, n)] = 0.5(ck(ped, j, ip) + ck(ped, j, im))
        end
        return kin[(m, n)]
    end
end

ckinship = fkin()
