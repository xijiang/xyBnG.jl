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
const _nvldtype = 13

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
        0.0 < h2 ≤ 1.0 || error("h² must be in (0, 1]")
        0 < nQTL || error("nQTL must be in (0, ∞)")
        ul = 1.0 / h2 - 1.0
        0.0 ≤ vd < ul || error("Var dominance must be in [0., (1 - h²) / h² = $ul)")
        new(name, type, sex, age, h2, nQTL, da, vd, dd, rev)
    end
end

function Base.show(io::IO, trt::Trait)
    println(io, "             Name: $(trt.name)")
    println(io, "             Type: $(trt.type)")
    println(io, "      Express age: $(trt.age)")
    if trt.sex == 0
        println(io, "     Expresses in: females")
    elseif trt.sex == 1
        println(io, "       Express in: males")
    else
        println(io, "       Express in: both sexes")
    end
    println(io, "      Express age: $(trt.age)")
    println(io, "     Heritability: $(trt.h2)")
    println(io, "          No. QTL: $(trt.nQTL)")
    println(io, " TBV distribution: $(trt.da)")
    println(io, "    Var dominance: $(trt.vd)")
    println(io, "   Dominance dist: $(trt.dd)")
    print(io, "Higher the better: $(trt.rev)")
end

"""
    function Trait(name, h2, nQTL; type = Float64, sex = 2, age = 1.,
            da = Normal(), vd = 0., dd = Normal())
        Trait(name, sex, age, h2, nQTL, da, vd, dd)
    end
"""
function Trait(
    name,
    h2,
    nQTL;
    type = Float64,
    sex = 2,
    age = 1.0,
    da = Normal(),
    vd = 0.0,
    dd = Normal(),
    rev = true,
)
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

function Base.show(io::IO, c::Cattle)
    println(io, "     Name: $(c.name)")
    print(io, "Pop. size: $(c.nid)")
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

function Base.show(io::IO, p::Plan)
    println(io, "    No. sires: $(p.npa)")
    println(io, "     No. dams: $(p.nma)")
    println(io, "No. offspring: $(p.noff)")
    print(io, "       Mating: $(p.mate)")
end

# ToDo: Add a generic species

"""
    _type(x::Union{Int8, DataType})

A literal bidirectional mapping between `Int8` and `DataType`.
This is used to decide the element type of the matrix,
or the type tag in the header of an `xy` file.

No range test, as this is an internal function. The valid types are
`Int8`, `Int16`, `Int32`, `Int64`, `Int128`, `UInt8`, `UInt16`, `UInt32`,
`UInt64`, `UInt128`, `Float32`, and `Float64`.
"""
function _type(x::Union{Int8,DataType})
    vt = (
        Int8,
        Int16,
        Int32,
        Int64,
        Int128,
        UInt8,
        UInt16,
        UInt32,
        UInt64,
        UInt128,
        Float32,
        Float64,
        Bool,
    )
    if x isa Int8
        vt[x]
    else
        Int8(findfirst(t -> t == x, vt))
    end
end

"""
    mutable struct header

The header struct of an `xy` file. It consists of 8 `Int8` fields.
- x::Int8
- y::Int8
- v::Int8      # x, y, v should be 'x', 'y', ' '. they are the magic chars
- flus::Int8   # FLUS matrix type: full, lower triangle, upper triangle, symmetric(lower)
- major::Int8  # 0 for loci majored, 1 for ID majored, or else
- type::Int8   # element type of the matrix, determined by function _type
- r::Int8      # 1 for BitArray. 0 for others
- u::Int8      # r and u are reserved
It is by default as (x = 'x', y = 'y', v = ' ', flus = 'F', major = 0, type = 1, r = 0, u = 0).
Or one can specify the fields by keyword arguments. For example,
```julia
header(x = 'x', y = 'y', v = ' ', flus = 'F', major = 0, type = Int8, r = 0, u = 0)
```

## Specification of `u`
- 0 for SNP coding
- 1 for IBD coding
- 2 for genotype coding
- 3+ else
"""
mutable struct header
    x::Int8
    y::Int8
    v::Int8      # x, y, v should be 'x', 'y', ' '. they are the magic chars
    flus::Int8   # FLUS matrix type: full, lower triangle, upper triangle, symmetric
    major::Int8  # 0 for loci majored, 1 for ID majored, or else
    type::Int8   # element type of the matrix, determined by function _type
    r::Int8      # 1 for BitArray. 0 for others
    u::Int8      # 0 for SNP coding, 1 for IBD coding, 2 for genotype coding, 3+ else
    function header(;
        flus = 'F',  # full, lower triangle, upper triangle, symmetric
        major = 0,
        type = 1,
        r = 0,
        u = 0,
    )
        flus ∉ "FLUS" ||
            major ∉ 0:2 ||
            type ∉ 1:_nvldtype ||
            r ∉ 0:1 ||
            u < 0 && error("Invalid header")
        new('x', 'y', ' ', flus, major, type, r, u)
    end
end

function Base.show(io::IO, h::header)
    if h.r == 1
        print(io, "A BitArray matrix")
    else
        print(io, "A normal matrix")
    end
    if h.u == 0
        println(io, " coding SNP alleles.")
    elseif h.u == 1
        println(io, " coding IBD alleles.")
    elseif h.u == 2
        println(io, " coding genotypes.")
    else
        println(io, " with unknown coding.")
    end
    println(io, "        -------:-------")
    m = join(Char.((h.x, h.y, h.v)))
    println(io, "The magic chars: '$m'.")
    if h.flus == Int8('F')
        println(io, "         Matrix: full matrix.")
    elseif h.flus == Int8('L')
        println(io, "         Matrix: lower triangle of a square matrix.")
    elseif h.flus == Int8('U')
        println(io, "         Matrix: upper triangle of a square matrix.")
    elseif h.flus == Int8('S')
        println(io, "         Matrix: symmetric matrix.")
    else
        println(io, "         Matrix: unknown.")
    end
    if h.major == 0
        println(io, "          Major: loci majored.")
    elseif h.major == 1
        println(io, "          Major: ID majored.")
    else
        println(io, "          Major: $(h.major) -- unknown.")
    end
    println(io, "   Element type: $(_type(h.type))")
end

end # module xyTypes
