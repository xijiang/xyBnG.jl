"""
    snphet(q::AbstractVector{Float64})
Calculate the mean heterozygosity of diallelic loci with allele frequency vector
`q`.
"""
function snphet(q::AbstractVector{Float64})
    H = q .* q + (1 .- q) .* (1 .- q) # homozygosity
    1 - mean(H) # heterozygosity
end

"""
    v_snphet(q::AbstractVector{Float64})
Calculate the locus specific heterozygosity of a diallelic locus with allele
frequency vector `q`.
"""
function v_snphet(q::AbstractVector{Float64})
    H = q .* q + (1 .- q) .* (1 .- q) # homozygosity
    1 .- H # heterozygosity
end
