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
    rc = zeros(nid)

    while true
        itr += 1
        u = ped.idx[id]
        Q = [ped.sex[id] .== 1 ped.sex[id] .== 0]
        Ai = inv(A[id, id])
        QAQ = Q'Ai * Q
        QAQi = inv(QAQ)
        f1 = u' * (Ai - Ai * Q * QAQi * Q'Ai) * u # numerator
        f2 = 4K - sum(QAQi) # denominator
        if f2 ≤ 0
            rc[id] = 0.5Ai * Q * QAQi * ones(2) # can't meet the constraint
            @warn "Cannot achieve constraint $K"
            break
        end
        f1 = max(f1, 1e-10)
        λ₀ = sqrt(f1 / f2)
        λ = QAQi * (Q'Ai * u .- λ₀)
        c = Ai * (u - Q * λ) / (2λ₀)
        ix = findall(c .> 0) # indices of ID of next round
        if length(ix) == length(id)
            rc[id] = c
            break
        end
        id = id[ix]
    end
    return rc
end

"""
    TM2024(ped, A, K)
This is a wrapper for Theo's program with equal contribution.
"""
function TM2024(ped, A, K)
    ncd = size(A, 1) # number of candidates
    nsire = sum(ped.sex)
    ndam = ncd - nsire
    @debug nsire, ndam, ncd
    DOSop(ped.idx, A, zeros(ncd), 1.0, K, [nsire, ndam], ped.sex .+ 1) / 2.0
end

#=
"""
    konstraint(dF::Float64, k₀::Float64, igrt::Int; ong = false)
Calculate the constraint about inbreeding increase in the next generation. Note,
when `G` relationship matrix is used, there is a consideration about frequency
already. Hence `2dF` is used as `K`.
"""
function konstraint(dF::Float64, k₀::Float64, igrt::Int; ong = false)
    if ong
        2dF
    else
        2(1 - (1 - k₀) * (1 - dF)^(igrt + 1))
    end
end
=#

"""
    OC_fixsex2(EBV, A, K, c2)
This function is to optimize contribution on the male side. It is modified such
that this function has similar arguments as `TM1997` and `TM2024`. 

- `c2` are fixed contribution of females (Σc₂ = 1).
- A is relationship matrix of all animals
    - A11 = relationship of candidates
    - A22 = relationship of selected females
"""
function OC_fixsex2(EBV, A, K, c2)
    c2 = 0.5 * c2 / sum(c2)

    n = size(EBV, 1)
    n1 = n
    n2 = size(c2, 1)
    ind = collect(1:n)
    # println(" Ncandidates           = ", n)
    # println(" Nfixed_contributions  =", n2)
    # println(" contraint relationship= ", K)
    A12c2 = A[1:n1, n1+1:n1+n2] * c2
    c = zeros(n)
    ierr = 0
    for it = 1:1000
        u = EBV[ind]
        AI = inv(A[ind, ind])
        sumAI = sum(AI)
        AIA12 = AI * A[ind, n1+1:n1+n2]
        A22_11 = A[n1+1:n1+n2, n1+1:n1+n2] - A[n1+1:n1+n2, ind] * AIA12
        denominat = 4 * (K - c2' * A22_11 * c2 - ((sum(AIA12 * c2) + 0.5)^2) / sumAI)
        numerat = u' * AI * u - sum(u' * AI)^2 / sumAI
        if (denominat <= 0.0)
            println(
                "Cannot achieve constraint ",
                K,
                " Minimisation of relationships ",
                size(ind),
            )
            ierr = 0
            c = 0.5 * sum(AI, dims = 2) / sumAI  #minimise
        else
            numerat = max(numerat, 1.e-10)
            # println(" denominat ",denominat)
            # println(" numerator ",numerat)
            lamb02 = numerat / denominat
            lamb0 = sqrt(lamb02)
            println(" lamb0 ", lamb0)
            lamb = (sum(AI * (u - 2 * lamb0 * A12c2[ind])) - lamb0) / sumAI
            c = AI * (u - 2 * lamb0 * A12c2[ind] .- lamb) / (2 * lamb0)
        end
        ind2 = findall(c .>= 0.0)
        println(" iter ", it, " still in solution ", length(ind2), " (old= ", length(ind))
        if (length(ind2) == length(ind))
            #            println(" solution found; n=",length(ind))
            break
        end
        ind = ind[ind2]
    end

    if (ierr == 0)
        cc = zeros(n)
        cc[ind] = c
        println("  cAc-2cA12c2= ", c' * A[ind, ind] * c .- 2 * cc' * A12c2, " K=  ", K)
        return 2 * cc
    end
end
