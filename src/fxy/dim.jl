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
Update dimensions to `nrows` and `ncols` to file `xy`.
"""
function dim!(xy::AbstractString, nrows::Int64, ncols::Int64)
    open(xy, "r+") do io
        seek(io, 8)
        write(io, [nrows, ncols])
    end
end
