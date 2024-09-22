"""
    function mapit(fxy::AbstractString)
Map the entire matrix in `fxy` and return the Matrix.
"""
function mapit(fxy::AbstractString)
    hdr = header(fxy)
    hdr.flus == Int8('F') || error("Only a full matrix is supported")
    m, n = dim(fxy)
    if hdr.r == 1
        Mmap.mmap(fxy, BitArray, (m, n), 24)
    else
        Mmap.mmap(fxy, Matrix{_type(hdr.type)}, (m, n), 24)
    end
end
