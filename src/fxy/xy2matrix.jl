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
                for i = 1:x
                    write(io, m[1:i, i])
                end
            else
                for i = 1:x
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
                for i = 1:x
                    read!(io, view(m, 1:i, i))
                end
            else
                for i = 1:x
                    read!(io, view(m, i:x, i))
                end
                if hdr.flus == Int8('S')
                    for i = 1:x-1
                        copy!(view(m, i, i+1:x), view(m, i+1:x, i))
                    end
                end
            end
        end
    end
    m
end
