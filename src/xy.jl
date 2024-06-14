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

const _nvldtype = 12

"""
    _type(x::Union{Int8, DataType})

A literal bidirectional mapping between `Int8` and `DataType`.
This is used to decide the element type of the matrix,
or the type tag in the header of an `xy` file.

No range test, as this is an internal function. The valid types are
`Int8`, `Int16`, `Int32`, `Int64`, `Int128`, `UInt8`, `UInt16`, `UInt32`,
`UInt64`, `UInt128`, `Float32`, and `Float64`.
"""
function _type(x::Union{Int8, DataType})
    vt = ( Int8,  Int16,  Int32,  Int64,  Int128,
          UInt8, UInt16, UInt32, UInt64, UInt128,
          Float32, Float64)
    x isa Int8 ? vt[x] : Int8(findfirst(t -> t == x, vt))
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
- r::Int8
- u::Int8      # r and u are reserved
It is by default as (x = 'x', y = 'y', v = ' ', flus = 'F', major = 0, type = 1, r = 0, u = 0).
Or one can specify the fields by keyword arguments. For example,
```julia
header(x = 'x', y = 'y', v = ' ', flus = 'F', major = 0, type = Int8, r = 0, u = 0)
```
"""
mutable struct header
    x::Int8
    y::Int8
    v::Int8      # x, y, v should be 'x', 'y', ' '. they are the magic chars
    flus::Int8   # FLUS matrix type: full, lower triangle, upper triangle, symmetric
    major::Int8  # 0 for loci majored, 1 for ID majored, or else
    type::Int8   # element type of the matrix, determined by function _type
    r::Int8      # r is reserved
    u::Int8      # 0 for SNP coding, 1 for IBD coding, 2 for genotype coding, 3+ else
    function header(;
        flus  = 'F',  # full, lower triangle, upper triangle, symmetric
        major = 0,
        type  = 1,
        u     = 0,
        )
        flus ∉ "FLUS" || major ∉ 0:1 || type ∉ 1:_nvldtype || u < 0 && 
            error("Invalid header")
        new('x', 'y', ' ', flus, major, type, 0, u)
    end
end

"""
    isEqual(h1::T, h2::T) where T <: header

Check if two headers are equal.
"""
function isEqual(h1::T, h2::T) where T <: header
    h1.x == h2.x && h1.y == h2.y && h1.v == h2.v && h1.flus == h2.flus &&
        h1.major == h2.major && h1.type == h2.type && h1.r == h2.r && h1.u == h2.u
end

"""
    init!(xy::AbstractString, hdr::header, nrows::Int64, ncols::Int64)

Initialize an `xy` file with header `hdr` and dimensions `nrows` and `ncols`.
The matrix part is filled with zeros.
"""
function init!(xy::AbstractString, hdr::header, nrows::Int64, ncols::Int64)
    isfile(xy) && rm(xy, force=true)
    block = 24
    nbyte = sizeof(_type(hdr.type))
    if hdr.flus == Int8('F')
        block += nrows * ncols * nbyte
    else
        nrows == ncols || error("nrows must equal to ncols for full matrix")
        block += nrows * (nrows + 1) ÷ 2 * nbyte
    end

    if Sys.islinux() || Sys.isbsd() || Sys.isapple() || Sys.isfreebsd() || Sys.isnetbsd() || Sys.isopenbsd() || Sys.isunix()
        run(pipeline(`head -c $block /dev/zero`, xy))
    elseif Sys.iswindows()
        run(`fsutil file createnew $xy $block`)
    else
        error("Unsupported OS")
    end

    open(xy, "r+") do io
        write(io, Ref(hdr))
        write(io, [nrows, ncols])
    end
end

"""
    header(xy::AbstractString)

Read header from file `xy` and return a `header` object.
Also check if the file is a valid `xy` file. If not, return `nothing`.
"""
function header(xy::AbstractString)
    # check if the file is a valid xy file
    sz = filesize(xy)
    sz ≤ 24 && return nothing
    hdr = header()
    read!(xy, Ref(hdr))
    (hdr.x == Int8('x') && hdr.y == Int8('y') && hdr.v == Int8(' ') && 0 < hdr.type ≤ _nvldtype) || return nothing
    nrows, ncols = dim(xy)
    nbyte = sizeof(_type(hdr.type))
    if hdr.flus == Int8('F')
        sz == 24 + nrows * ncols * nbyte || return nothing
    else
        sz == 24 + nrows * (nrows + 1) ÷ 2 * nbyte || return nothing
    end
    hdr
end

function header!(xy::AbstractString, hdr::header)
    open(xy, "r+") do io
        write(io, Ref(hdr))
    end
end

"""
    dim(xy::AbstractString)

Read dimensions from file `xy` and return a tuple `(nrows, ncols)`.
"""
function dim(xy::AbstractString)
    filesize(xy) ≤ 24 && error("file size too small")
    dm = Mmap.mmap(xy, Vector{Int64}, 2, 8)
    dm[1], dm[2]
end

"""
    dim!(xy::AbstractString, nrows::Int64, ncols::Int64)

Write dimensions `nrows` and `ncols` to file `xy`.
"""
function dim!(xy::AbstractString, nrows::Int64, ncols::Int64)
    open(xy, "r+") do io
        seek(io, 8)
        write(io, [nrows, ncols])
    end
end

"""
    sub(ixy::AbstractString, rows::AbstractVector, cols::AbstractVector)

Extract a submatrix from `ixy` and return a matrix in memory.
"""
function sub(ixy::AbstractString, rows::AbstractVector, cols::AbstractVector)
    hdr = header(ixy)
    length(rows) * length(cols) * sizeof(_type(hdr.type)) ≤ Util.memavail() || error("Submatrix too large to fit in memory.")
    isnothing(hdr) && error("$ixy is not a valid xy file.")
    hdr.flus == Int8('F') || error("Needs to be a full matrix.")
    m, n = dim(ixy)
    (lr, hr) = eltype(rows) == Bool ? (1, length(rows)) : extrema(rows)
    (lc, hc) = eltype(cols) == Bool ? (1, length(cols)) : extrema(cols)
    (lr < 1 || hr > m || lc < 1 || hc > n) && error("row or column index out of range")
    mat = Mmap.mmap(ixy, Matrix{_type(hdr.type)}, (m, n), 24)
    copy(mat[rows, cols])
end

"""
    sub(ixy::AbstractString, rows::T, cols::T, oxy::AbstractString) where T <: AbstractVector{Int64}

Extract a submatrix from `ixy` and write sub to `oxy`.

- *Note*:
  - The target file `oxy` will be overwritten if it exists.
  - Rows and columns can be not sorted, but should be in range.
  - This function deals with full matrix.
"""
function sub(ixy::AbstractString, rows::T, cols::T, oxy::AbstractString) where T <: AbstractVector{Int64}
    isfile(oxy) && rm(oxy, force=true)
    hdr = header(ixy)
    isnothing(hdr) && error("$ixy is not a valid xy file.")
    hdr.flus == Int8('F') || error("Needs to be a full matrix.")
    m, n = dim(ixy)
    (lr, hr) = eltype(rows) == Bool ? (1, length(rows)) : extrema(rows)
    (lc, hc) = eltype(cols) == Bool ? (1, length(cols)) : extrema(cols)
    (lr < 1 || hr > m || lc < 1 || hc > n) && error("row or column index out of range")
    open(oxy, "w") do io
        write(io, Ref(hdr))
        write(io, [length(rows), length(cols)])
        mat = Mmap.mmap(ixy, Matrix{_type(hdr.type)}, (m, n), 24)
        write(io, mat[rows, cols])
    end
end

"""
    sub(ixy::S, rows::T, oxy::S) where {S <: AbstractString, T <: AbstractVector}

Extract a triangle submatrix from `ixy` of a triangle matrix and write sub to `oxy`.

- *Note*:
  - The target file `oxy` will be overwritten if it exists.
  - Rows must be sorted.
  - This function deals with triangle matrix.
  - Bear in mind that matrix is stored in column major order.
"""
function sub(ixy::S, rows::T, oxy::S) where {S <: AbstractString, T <: AbstractVector}
    isfile(oxy) && rm(oxy, force=true)
    hdr = header(ixy)
    isnothing(hdr) && error("$ixy is not a valid xy file.")
    hdr.flus == Int8('F') && error("Needs to be a triangle matrix.")
    (m, _), (lr, hr), len = dim(ixy), extrema(rows), length(rows)
    (lr < 1 || hr > m) && error("Index out of range")
    issorted(rows) || error("Index must be sorted")
    nbyte = sizeof(_type(hdr.type))
    open(ixy, "r") do ii
        tvc = zeros(_type(hdr.type), m) # temp vector for column
        open(oxy, "w") do io
            write(io, Ref(hdr))
            write(io, [len, len])
            for i in 1:len
                if hdr.flus == Int8('U')
                    seek(ii, 24 + rows[i] * (rows[i] - 1) ÷ 2 * nbyte)
                    read!(ii, tvc[1:rows[i]])
                    write(io, tvc[rows[1:i]])
                else
                    seek(ii, 24 + (rows[i] + 1 + m) * (m - rows[i]) ÷ 2 * nbyte)
                    read!(ii, tvc[rows[i]:end])
                    write(io, tvc[rows[i:end]])
                end
            end
        end
    end
end

"""
    mat(m::AbstractMatrix, oxy::AbstractString; flus = 'F')

Write matrix `m` to file `oxy` in `xy` format. When doing half matrix, this
function will only check if matrix is square. It will not check if the matrix is
symmetric.
"""
function mat(m::AbstractMatrix, oxy::AbstractString; flus = 'F')
    x, y = size(m)
    type = _type(eltype(m))
    isnothing(type) && error("Element type not supported")
    x ≠ y && flus ≠ 'F' && error("Matrix must be square for lower or upper triangle")
    hdr = header(flus = flus, type = type)
    open(oxy, "w") do io
        if flus == 'F'
            write(io, Ref(hdr))
            write(io, [x, y])
            write(io, m)
        else
            write(io, Ref(hdr))
            write(io, [x, y])
            if flus == 'U'
                for i in 1:x
                    write(io, m[1:i, i])
                end
            else
                for i in 1:x
                    write(io, m[i:end, i])
                end
            end
        end
    end
end

"""
    mat(xy::AbstractString)

Read matrix from file `xy` and return a matrix in memory.
"""
function mat(xy::AbstractString)
    hdr = header(xy)
    isnothing(hdr) && error("$xy is not a valid xy file.")
    (x, y) = dim(xy)
    m = zeros(_type(hdr.type), x, y)

    open(xy, "r") do io
        seek(io, 24)
        if hdr.flus == Int8('F')
            read!(io, m)
        else
            if hdr.flus == Int8('U')
                for i in 1:x
                    read!(io, view(m, 1:i, i))
                end
            else
                for i in 1:x
                    read!(io, view(m, i:x, i))
                end
                if hdr.flus == Int8('S')
                    for i in 1:x-1
                        copy!(view(m, i, i+1:x), view(m, i+1:x, i))
                    end
                end
            end
        end
    end
    m
end

"""
    append!(ixy::T, oxy::T) where T <: AbstractString

Append matrix `oxy` to `ixy`. The two files must have the same header. Their
number of rows must be the same. This applies only to full matrix.
"""
function append!(ixy::T, oxy::T) where T <: AbstractString
    hi, ho = header(ixy), header(oxy)
    isEqual(hi, ho) || error("Header of $ixy and $oxy are different")
    hi.flus == Int8('F') || error("Needs to be a full matrix")
    (mi, ni), (mo, no) = dim(ixy), dim(oxy)
    mi == mo || error("Number of rows of $ixy and $oxy are different")
    dim!(ixy, mi, ni + no)
    open(ixy, "a") do io
        write(io, Mmap.mmap(oxy, Matrix{_type(hi.type)}, (mo, no), 24))
    end
end

"""
    append!(fxy::AbstractString, gt::Matrix)

Append matrix `gt` to `fxy`.
"""
function append!(fxy::AbstractString, gt::Matrix)
    hdr = header(fxy)
    isnothing(hdr) && error("$fxy is not a valid xy file.")
    _type(hdr.type) == eltype(gt) || error("Element type of $fxy and gt are different")
    hdr.flus == Int8('F') || error("$fxy is not a full matrix")
    (x, y) = dim(fxy)
    x == size(gt, 1) || error("Number of rows of $fxy and gt are different")
    open(fxy, "a") do io
        write(io, gt)
    end
    dim!(fxy, x, y + size(gt, 2))
end

"""
    merge(ixy::T, jxy::T, oxy::T) where T <: AbstractString

Merge two `xy` files `ixy` and `jxy` and write to `oxy`.
- The two files must have the same header.
- The merge dimension must be same.
- This applies only to full matrix.
"""
function merge(ixy::T, jxy::T, oxy::T; horizontal=true) where T <: AbstractString
    (oxy == ixy || oxy == jxy) && error("Target file should be different from the source files.")
    hi, hj = header(ixy), header(jxy)
    isEqual(hi, hj) || error("Header of $ixy and $jxy are different")
    hi.flus == Int8('F') || error("Needs to be a full matrix")
    (mi, ni), (mj, nj) = dim(ixy), dim(jxy)
    type = _type(hi.type)

    if horizontal
        mi == mj || error("Number of rows of $ixy and $jxy are different")
        open(oxy, "w") do io
            write(io, Ref(hi))
            write(io, [mi, ni + nj])
            write(io, Mmap.mmap(ixy, Matrix{type}, (mi, ni), 24))
            write(io, Mmap.mmap(jxy, Matrix{type}, (mj, nj), 24))
        end
    else
        ni == nj || error("Number of columns of $ixy and $jxy are different")
        open(oxy, "w") do io
            write(io, Ref(hi))
            write(io, [mi + mj, ni])
            a = Mmap.mmap(ixy, Matrix{type}, (mi, ni), 24)
            b = Mmap.mmap(jxy, Matrix{type}, (mj, nj), 24)
            for i in 1:ni
                write(io, a[:, i])
                write(io, b[:, i])
            end
        end
    end
end

"""
    code(ixy::T, oxy::T) where T <: AbstractString

Code the SNP alleles unique values in `ixy` and write to `oxy`. The last bit
store the SNP allele types. The matrix element may increase size. Use
`isodd.(mat)` to retrieve the SNP allele types.
"""
function code(ixy::T, oxy::T) where T <: AbstractString
    hdr = header(ixy)
    isnothing(hdr) && error("$ixy is not a valid xy file.")
    (x, y) = dim(ixy)
    hdr.flus == Int8('F') && hdr.type == 1 && y % 2 == 0 || error("Needs to be a haplotype matrix")
    hdr.type = Int8.(ceil(log2(y) / 8))

    open(oxy, "w") do io
        write(io, Ref(hdr))
        write(io, [x, y])
        imt = Mmap.mmap(ixy, Matrix{UInt8}, (x, y), 24)
        inc = 2
        for i in 1:y
            write(io, _type(hdr.type).(imt[:, i] .+ inc))
            inc += 2
        end
    end
end

"""
    transpose!(ixy::T, oxy::T) where T <: AbstractString

Transpose matrix `ixy` and write to `oxy`.
For half matrix, `U` is transposed to `L`. `L` and `S` are transposed to `U`.
"""
function transpose!(ixy::T, oxy::T) where T <: AbstractString
    hdr = header(ixy)
    isnothing(hdr) && error("$ixy is not a valid xy file.")
    (x, y) = dim(ixy)
    type = _type(hdr.type)
    open(oxy, "w") do io
        if hdr.flus == Int8('F')
            write(io, Ref(hdr))
            write(io, [y, x])
            write(io, transpose(Mmap.mmap(ixy, Matrix{type}, (x, y), 24)))
        else
            @info "Under construction"
            vec = Mmap.mmap(ixy, Vector{type}, x * (x + 1) ÷ 2, 24)
            hdr.flus = (hdr.flus == Int8('U') ? Int8('L') : Int8('U'))
            write(io, Ref(hdr))
            write(io, [y, x])
            #if hdr.flus == Int8('U')
            #    for i in 1:x
            #        write(io, transpose(view(Mmap.mmap(ixy, Matrix{_type(hdr.type)}, (x, y), 24), 1:i, i)))
            #    end
            #else
            #    for i in 1:x
            #        write(io, transpose(view(Mmap.mmap(ixy, Matrix{_type(hdr.type)}, (x, y), 24), i:end, i)))
            #    end
            #end
        end
    end
end

"""
    snp2gt(ixy::T, oxy::T) where T <: AbstractString

Convert SNP alleles to genotypes and write to `oxy`.
"""
function snp2gt(ixy::T, oxy::T) where T <: AbstractString
    hdr = header(ixy)
    isnothing(hdr) && error("$ixy is not a valid xy file.")
    (x, y) = dim(ixy)
    hdr.flus == Int8('F') && hdr.type == 1 && y % 2 == 0 || error("Needs to be a haplotype matrix")
    open(oxy, "w") do io
        write(io, Ref(hdr))
        write(io, [x, y ÷ 2])
        mat = Mmap.mmap(ixy, Matrix{UInt8}, (x, y), 24)
        write(io, mat[:, 1:2:end] + mat[:, 2:2:end])
    end
end

# some helper functions

"""
    isfull(hdr::header)
Check if the matrix is a full matrix, i.e., not a triangle one.
"""
function isfull(hdr::header)
    hdr.flus == Int8('F')
end

"""
    locMajor(hdr::header)

Check if the matrix is loci majored.
"""
function locMajor(hdr::header)
    hdr.major == 0
end

"""
    isbyte(hdr::header)

Check if the element type of the matrix is `Int8`.
"""
function isbyte(hdr::header)
    _type(hdr.type) == Int8
end

"""
    ishap(hdr::header, dm::Tuple{Int, Int})

Check if the matrix is a haplotype matrix.
"""
function ishap(hdr::header, dm::Tuple{Int, Int})
    _, n = dm
    locMajor(hdr) && n % 2 == 0
end

"""
    issnp(fxy::AbstractString)

Check if the matrix is a SNP matrix.
"""
function issnp(fxy::AbstractString)
    hdr, dm = header(fxy), dim(fxy)
    isfull(hdr) && isbyte(hdr) && ishap(hdr, dm) && begin
        gt = Mmap.mmap(fxy, Matrix{Int8}, dm, 24)
        Set(gt) == Set([0, 1])
    end
end

"""
    function mapit(fxy::AbstractString)
Map the entire matrix in `fxy` and return the Matrix.
"""
function mapit(fxy::AbstractString)
    hdr = header(fxy)
    hdr.flus == Int8('F') || error("Only a full matrix is supported")
    tp = _type(hdr.type)
    m, n = dim(fxy)
    Mmap.mmap(fxy, Matrix{tp}, (m, n), 24)
end
end # module XY
