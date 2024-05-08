"""
    function predict!(ID::AbstractVector{Int}, ped::DataFrame, trts::Trait...)

Give random EBV of trait `trts` about ID specified with `ID` in DataFrame `ped`.
The rest of the values, and the sex not specific values are set as missing.
"""
function predict!(ID::AbstractVector{Int}, ped::DataFrame, trts::Trait...)
    1 ≤ minimum(ID) ≤ maximum(ID) ≤ nrow(ped) || throw(ArgumentError("ID out of range"))
    for trt in trts
        ebv = "ebv_" * trt.name
        hasproperty(ped, ebv) || (ped[!, ebv] .= missing)
        rand!(view(ped, ID, ebv))
    end
    ped
end

#=
"""
    predict!(ped::DataFrame, trts::Trait...; all = false)

Give random values to the :ebv of trait `trts` in DataFrame `ped`. If `all` is
true, or the ebv column the trait doesn't exist in `ped`, then all ID will be 
assigned a random value. Otherwise, only those in the last generation will be
assigned random values.
"""
function predict!(ped::DataFrame, trts::Trait...; all = false)
    for trt in trts
        #ex = Regex("\\Q$(trt.name)\\E\$") # to find all about a trait
        tgt = all ? ped : groupby(ped, :grt)[end]
        tgt[!, "ebv_" * trt.name] = rand(nrow(tgt))
    end
    ped
end

function predict!(ped::DataFrame, cf, giv::Matrix{Float64}, trts::Trait...;
    all = false, h2 = :known)
end

function predict!(ID::AbstractVector{Int}, ped::DataFrame, cf, giv::Matrix{Float64},
    trts::Trait...; all = false, h2 = :known)
end
=#