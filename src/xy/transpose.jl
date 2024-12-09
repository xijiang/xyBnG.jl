"""
    transpose!(ixy::T, oxy::T) where T <: AbstractString
Transpose matrix `ixy` and write to `oxy`.
For half matrix, `U` is transposed to `L`. `L` and `S` are transposed to `U`.
"""
function transpose!(ixy::T, oxy::T) where {T<:AbstractString}
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
