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
