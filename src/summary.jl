"""
    Sum
The ``Sum`` module provides functions to calculate summary statistics of the
simulations. It has majorly two parts. One is to summarize each repeat of a
simulation. The other contains a few notebooks for the presentation of the
summaries.
"""
module Sum

using xyBnG.XY

"""
    inbreeding(xy::AbstractString, loci::T) where T <: AbstractVector{Int}
Calculate the inbreeding coefficient of each individual in with the coded
genotype file `xy` at `loci` loci.
"""
function inbreeding(xy::AbstractString, loci::T) where {T<:AbstractVector{Int}}
    nlc, nhp = XY.dim(xy)
    mat = XY.mapit(xy)
    F = zeros(nhp ÷ 2)
    id = 1
    for i in 1:2:nhp
        F[id] = sum(mat[:, i] .== mat[:, i + 1]) / nlc
        id += 1
    end
    F
end

end # module Sum
