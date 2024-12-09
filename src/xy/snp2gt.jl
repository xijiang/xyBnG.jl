"""
    snp2gt(ixy::T, oxy::T) where T <: AbstractString
Convert SNP alleles to genotypes and write to `oxy`.
"""
function snp2gt(ixy::T, oxy::T) where {T<:AbstractString}
    hdr = header(ixy)
    isnothing(hdr) && error("$ixy is not a valid xy file.")
    (x, y) = dim(ixy)
    hdr.flus == Int8('F') && hdr.type == 1 && y % 2 == 0 ||
        error("Needs to be a haplotype matrix")
    open(oxy, "w") do io
        write(io, Ref(hdr))
        write(io, [x, y ÷ 2])
        mat = Mmap.mmap(ixy, Matrix{UInt8}, (x, y), 24)
        write(io, mat[:, 1:2:end] + mat[:, 2:2:end])
    end
end
