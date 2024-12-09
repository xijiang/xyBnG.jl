"""
    isfull(hdr::header)
Check if the matrix is a full matrix, i.e., not a triangle one.
"""
function isfull(hdr::header)
    hdr.flus == Int8('F')
end

function isfull(hdr::Nothing)
    @error "Header is nothing"
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
function ishap(hdr::header, dm::Tuple{Int,Int})
    _, n = dm
    locMajor(hdr) && n % 2 == 0
end

"""
    issnp(fxy::AbstractString)
Check if the matrix is a SNP matrix.
"""
function issnp(fxy::AbstractString)
    hdr, dm = header(fxy), dim(fxy)
    isfull(hdr) &&
        isbyte(hdr) &&
        ishap(hdr, dm) &&
        begin
            gt = Mmap.mmap(fxy, Matrix{Int8}, dm, 24)
            Set(gt) == Set([0, 1])
        end
end
