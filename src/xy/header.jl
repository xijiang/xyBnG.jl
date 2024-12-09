"""
    isEqual(h1::T, h2::T) where T <: header
Check if two headers are equal.
"""
function isEqual(h1::T, h2::T) where {T<:header}
    h1.x == h2.x &&
        h1.y == h2.y &&
        h1.v == h2.v &&
        h1.flus == h2.flus &&
        h1.major == h2.major &&
        h1.type == h2.type &&
        h1.r == h2.r &&
        h1.u == h2.u
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
    read!(xy, hdr)
    (
        hdr.x == Int8('x') &&
        hdr.y == Int8('y') &&
        hdr.v == Int8(' ') &&
        0 < hdr.type ≤ _nvldtype
    ) || return nothing
    nrows, ncols = dim(xy)
    nbyte = sizeof(_type(hdr.type))
    if hdr.flus == Int8('F')
        sz == 24 + nrows * ncols * nbyte || nothing
    else
        sz == 24 + nrows * (nrows + 1) ÷ 2 * nbyte || return nothing
    end
    hdr
end

"""
    header!(xy::AbstractString, hdr::header)
Update header of file `xy`
"""
function header!(xy::AbstractString, hdr::header)
    open(xy, "r+") do io
        write(io, hdr)
    end
end
