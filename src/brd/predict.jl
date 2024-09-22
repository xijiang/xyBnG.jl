"""
    Predict!(ID::AbstractVector{T}, ped::DataFrame, trts::Trait...) where T <: Integer

Give random EBV of trait `trts` about ID specified with `ID` in DataFrame `ped`.
The rest of the values, and the sex not specific values are set as missing.
"""
function Predict!(ID::AbstractVector{T}, ped::DataFrame, trts::Trait...) where {T<:Integer}
    1 ≤ minimum(ID) ≤ maximum(ID) ≤ size(ped, 1) || throw(ArgumentError("ID out of range"))
    for trt in trts
        ebv = "ebv_" * trt.name
        hasproperty(ped, ebv) || (ped[!, ebv] .= missing)
        copyto!(view(ped, ID, ebv), rand(length(ID)))
    end
    ped
end

"""
    Predict!(ID::AbstractVector{T}, ped::DataFrame,
        fixed::Vector{AbstractString}, giv::Matrix{Float64},
        trts::Trait...) where T <: Integer

Predict the trait(s) `trts` EBV of `ID` in DataFrame `ped` with the fixed
effects `fixed`, which is/are also column(s) in `ped`, the inverse relationship
matrix `giv`.
"""
function Predict!(
    ID::AbstractVector{T},
    ped::DataFrame,
    fixed::Vector{S},
    giv::Matrix{Float64},
    trts::Trait...,
) where {T<:Integer,S<:AbstractString}
    1 ≤ minimum(ID) ≤ maximum(ID) ≤ size(ped, 1) || throw(ArgumentError("ID out of range"))
    X = incidence_matrix(select(ped, fixed))
    for trt in trts
        ebv, ft = "ebv_" * trt.name, "ft_" * trt.name
        _, sol = animalModel(ped[!, ft], giv, trt.h2, X)
        copyto!(view(ped, ID, ebv), sol[ID])
    end
end

#=
"""
    Predict!(ped::DataFrame, trts::Trait...; all = false)

Give random values to the :ebv of trait `trts` in DataFrame `ped`. If `all` is
true, or the ebv column the trait doesn't exist in `ped`, then all ID will be 
assigned a random value. Otherwise, only those in the last generation will be
assigned random values.
"""
function Predict!(ped::DataFrame, trts::Trait...; all = false)
    for trt in trts
        #ex = Regex("\\Q$(trt.name)\\E\$") # to find all about a trait
        tgt = all ? ped : groupby(ped, :grt)[end]
        tgt[!, "ebv_" * trt.name] = rand(nrow(tgt))
    end
    ped
end

function Predict!(ped::DataFrame, cf, giv::Matrix{Float64}, trts::Trait...;
    all = false, h2 = :known)
end

function Predict!(ID::AbstractVector{Int}, ped::DataFrame, cf, giv::Matrix{Float64},
    trts::Trait...; all = false, h2 = :known)
end
=#
