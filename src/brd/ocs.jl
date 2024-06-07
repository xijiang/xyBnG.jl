"""
    TM1997(ped, A, K)

`K` = is the restriction on average relationships (``Meuwissen, 1997, J Anim
Sci``). Let the constraint on the inbreeding in generation `t`` be:

`k(t) = k(t-1) + DF*(1-k(t-1))`

where `k(t-1)` = constraint on inbreeding in generation `t-1` and `k(0)=0`. Then
the constraint on relationship is `K=2k(t)`.

Based on Theo's function below with some simplifications.

DataFrame `ped` must have columns `idx` and `sex` and `A` must be a square
matrix of the same size as the number of rows in `ped`.

The column `idx` can be EBV of a trait, e.g., ``select(ped, "ebv_milk => "idx",
"sex")``. Or, it can be an index combining EBV of several traits. 

The function returns a vector of coefficients for the contribution of each ID.
"""
function TM1997(ped, A, K)
    ["idx", "sex"] ⊆ names(ped) || error("Not a proper pedigree")
    issymmetric(A) || error("A is not symmetric")
    nid, itr = nrow(ped), 0
    size(A, 1) == nid || error("A and ped do not match")
    id = 1:nid
    
    while true
        itr += 1
        u = ped.idx[id]
        Q = [ped.sex[id] .== 1 ped.sex[id] .== 0]
        Ai   = inv(A[id, id])
        QAQ  = Q'Ai * Q
        QAQi = inv(QAQ)
        f1 = u' * (Ai - Ai * Q * QAQi * Q'Ai) * u # numerator
        f2 = 4K - sum(QAQi) # denominator
        f2 ≤ 0 && return 0.5Ai * Q * QAQi * ones(2) # can't meet the constraint
        f1 = max(f1, 1e-10)
        λ₀ = sqrt(f1 / f2)
        λ  = QAQi * (Q'Ai * u .- λ₀)
        c  = Ai * (u - Q * λ) / (2λ₀)
        ix = findall(c .> 0) # indices of ID of next round
        if length(ix) == length(id)
            rc = zeros(nid)
            rc[id] = c
            return rc
        end
        id = id[ix]
    end
end

function TM2024()
end
