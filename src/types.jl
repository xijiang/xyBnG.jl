"""
    module xyTypes

Types for `xyBnG`.

## Provides
- `struct Trait`
- `function Trait`
- `abstract type Species`
- `struct Cattle`
- `function Cattle`
"""
module xyTypes

using Distributions
export Trait

"""
    struct Trait
        name::AbstractString
        type::DataType
        sex::Int     # 0, or 1, sex limited, 2 for both sexes.
        age::Float64 # age that trait is measured
        h2::Float64  # narrow-sense heritability
        nQTL::Int
        da::ContinuousUnivariateDistribution  # e.g., Normal()
        vd::Float64
        dd::ContinuousUnivariateDistribution
        rev::Bool    # default true, higher the better
    end

Note: The total additive genetic variance is normalized to mean 0, and variance
1.0. Hence vd should not be greater than 1./h2 - 1.0.
"""
struct Trait
    name::AbstractString
    type::DataType
    sex::Int     # 0, or 1, sex limited, 2 for both sexes.
    age::Float64 # age that trait is measured
    h2::Float64
    nQTL::Int
    da::ContinuousUnivariateDistribution
    vd::Float64
    dd::ContinuousUnivariateDistribution
    rev::Bool    # default true, higher the better
    function Trait(name, type, sex, age, h2, nQTL, da, vd, dd, rev)
        isempty(setdiff(name, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_")) ||
            error("Only alphabet and underscore are allowed in a name")
        type ∈ [Int, Float64] || error("Type must be Int or Float64")
        sex ∈ 0:2 || error("Sex can only be in 0:2")
        age > 0 || error("Age must be positive")
        0. < h2 ≤ 1. || error("h² must be in (0, 1]")
        0 < nQTL || error("nQTL must be in (0, ∞)")
        ul = 1. / h2 - 1.
        0. ≤ vd < ul || error("Var dominance must be in [0., (1 - h²) / h² = $ul)")
        new(name, type, sex, age, h2, nQTL, da, vd, dd, rev)
    end
end

"""
    function Trait(name, h2, nQTL; type = Float64, sex = 2, age = 1.,
            da = Normal(), vd = 0., dd = Normal())
        Trait(name, sex, age, h2, nQTL, da, vd, dd)
    end
"""
function Trait(name, h2, nQTL; type = Float64, sex = 2, age = 1.,
         da = Normal(), vd = 0., dd = Normal(), rev = true)
    Trait(name, type, sex, age, h2, nQTL, da, vd, dd, rev)
end

"""
    function name(trt::Trait)

Return the name of the trait.
"""
function name(trt::Trait)
    trt.name
end

abstract type Species end

struct Cattle <: Species
    name::AbstractString
    nid::Int
    function Cattle(name::AbstractString, nid::Int)
        nid ≤ 0 && error("nid must be positive")
        new(name, nid)
    end
end

function Cattle(nid::Int)
    Cattle("BosTau", nid)
end

struct Pig <: Species
    name::AbstractString
    nid::Int
    function Pig(name::AbstractString, nid::Int)
        nid ≤ 0 && error("nid must be positive")
        new(name, nid)
    end
end

function Pig(nid::Int)
    Pig("SusScr", nid)
end

mutable struct Plan
    npa::Int
    nma::Int
    noff::Int
    mate::Symbol
    function Plan(npa::Int, nma::Int, noff::Int; mate = :hierarchical)
        npa > 0 || error("npa must be positive")
        nma > 0 || error("nma must be positive")
        noff > 0 || error("noff must be positive")
        mate ∈ [:random, :hierarchical, :factorial] || 
            error("mate must be in [:random, :hierarchical, :factorial]")
        new(npa, nma, noff, mate)
    end
end

# ToDo: Add a generic species

end # module xyTypes
