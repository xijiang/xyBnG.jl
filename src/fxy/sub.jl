"""
    sub(ixy::AbstractString, rows::AbstractVector, cols::AbstractVector)
Extract a submatrix from `ixy` and return a matrix in memory.
"""
function sub(ixy::AbstractString, rows::AbstractVector, cols::AbstractVector)
    hdr = header(ixy)
    length(rows) * length(cols) * sizeof(_type(hdr.type)) ≤ Util.memavail() ||
        error("Submatrix too large to fit in memory.")
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
function sub(
    ixy::AbstractString,
    rows::T,
    cols::T,
    oxy::AbstractString,
) where {T<:AbstractVector{Int64}}
    isfile(oxy) && rm(oxy, force = true)
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
function sub(ixy::S, rows::T, oxy::S) where {S<:AbstractString,T<:AbstractVector}
    isfile(oxy) && rm(oxy, force = true)
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
            for i = 1:len
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
