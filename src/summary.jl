module Sum

using xyBnG.XY

"""
    inbreeding(xy::AbstractString, loci::T) where T <: AbstractVector{Int}
Calculate the inbreeding coefficient of each individual in with the coded
genotype file `xy` at `loci` loci.
"""
function inbreeding(xy::AbstractString, loci::T) where {T<:AbstractVector{Int}}
end

end # module Sum