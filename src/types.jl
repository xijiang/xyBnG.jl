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
        h2::Float64  # narrow-sense heritability
        nQTL::Int
        da::ContinuousUnivariateDistribution  # e.g., Normal()
        vd::Float64
        dd::ContinuousUnivariateDistribution
    end

Note: The total additive genetic variance is normalized to mean 0, and variance
1.0. Hence vd should not be greater than 1./h2 - 1.0.
"""
struct Trait
    name::AbstractString
    h2::Float64
    nQTL::Int
    da::ContinuousUnivariateDistribution
    vd::Float64
    dd::ContinuousUnivariateDistribution
    function Trait(name, h2, nQTL, da, vd, dd)
        isempty(setdiff(name, "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_")) ||
            error("Only alphabet and underscore are allowed in a name")
        0. < h2 ≤ 1. || error("h² must be in (0, 1]")
        0 < nQTL || error("nQTL must be in (0, ∞)")
        ul = 1. / h2 - 1.
        0. ≤ vd < ul || error("Var dominance must be in [0., (1 - h²) / h² = $ul)")
        new(name, h2, nQTL, da, vd, dd)
    end
end

"""
    function Trait(name, h2, nQTL; da = Normal(), vd = 0., dd = Normal())
        Trait(name, h2, nQTL, da, vd, dd)
    end
"""
function Trait(name, h2, nQTL; da = Normal(), vd = 0., dd = Normal())
    Trait(name, h2, nQTL, da, vd, dd)
end

end # module xyTypes
