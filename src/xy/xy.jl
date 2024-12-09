"""
    module XY
Functions around `xy` file format.

## Provides
- `function _type(x::Union{Int8, DataType})`
- `function isEqual(h1::T, h2::T) where T <: header`
- `function init!(xy::AbstractString, hdr::header, nrows::Int64, ncols::Int64)`
- `function header(xy::AbstractString)`
- `function dim(xy::AbstractString)`
- `function dim!(xy::AbstractString, nrows::Int64, ncols::Int64)`
- `function sub(ixy::AbstractString, rows::T, cols::T, oxy::AbstractString) where T <: AbstractVector{Int64}`
- `function sub(ixy::AbstractString, rows::T, oxy::AbstractString) where T <: AbstractVector{Int64}`
- `function mat(m::AbstractMatrix, oxy::AbstractString; flus = 'F')`
- `function mat(xy::AbstractString)`
- `function append!(ixy::T, oxy::T) where T <: AbstractString`
- `function merge(ixy::T, jxy::T, oxy::T; horizontal=true) where T <: AbstractString`
- `function code(ixy::T, oxy::T) where T <: AbstractString`
- function transpose!(ixy::T, oxy::T) where T <: AbstractString
- `function snp2gt(ixy::T, oxy::T) where T <: AbstractString`
"""
module XY

using xyBnG.Util
using Mmap
using Serialization
import xyBnG.xyTypes: header, _type, _nvldtype

include("header.jl")
include("dim.jl")
include("init-xy.jl")
include("extract-chr.jl")
include("map-it.jl")
include("recognize.jl")
include("transpose-n-code-byte-2-bit.jl")
include("transpose.jl")
include("uniq-code.jl")
include("snp2gt.jl")
include("merge.jl")
include("sub.jl")
include("xy2matrix.jl")
include("append.jl")

end # module XY
