"""
    tibd(mat::AbstractMatrix)
Calculate the inbreeding coefficient of each individual in with haplotype matrix
`mat`.
"""
function tibd(mat::AbstractMatrix)
    nlc, nhp = size(mat)
    F = zeros(nhp ÷ 2)
    Threads.@threads for i = 1:nhp÷2
        F[i] = sum(mat[:, 2i-1] .== mat[:, 2i]) / nlc
    end
    F
end
