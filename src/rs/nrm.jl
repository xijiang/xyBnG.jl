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
        for jd in 1:id-1
            sire ≠ 0 && (A[id, jd]  = 0.5A[jd, sire])
            dam  ≠ 0 && (A[id, jd] += 0.5A[jd, dam])
            A[jd, id] = A[id, jd]
        end
    end
    A
end
