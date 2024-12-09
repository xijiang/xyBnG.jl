"""
    merge(ixy::T, jxy::T, oxy::T) where T <: AbstractString
Merge two `xy` files `ixy` and `jxy` and write to `oxy`.
- The two files must have the same header.
- The merge dimension must be same.
- This applies only to full matrix.
"""
function merge(ixy::T, jxy::T, oxy::T; horizontal = true) where {T<:AbstractString}
    (oxy == ixy || oxy == jxy) &&
        error("Target file should be different from the source files.")
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
            for i = 1:ni
                write(io, a[:, i])
                write(io, b[:, i])
            end
        end
    end
end
