"""
    code(ixy::T, oxy::T) where T <: AbstractString
Code the SNP alleles unique values in `ixy` and write to `oxy`. The last bit
store the SNP allele types. The matrix element may increase size. Use
`isodd.(mat)` to retrieve the SNP allele types.
"""
function code(ixy::T, oxy::T) where {T<:AbstractString}
    hdr = header(ixy)
    isnothing(hdr) && error("$ixy is not a valid xy file.")
    (x, y) = dim(ixy)
    hdr.flus == Int8('F') && hdr.type == 1 && y % 2 == 0 ||
        error("Needs to be a haplotype matrix")
    hdr.type = Int8.(ceil(log2(y) / 8))

    open(oxy, "w") do io
        write(io, Ref(hdr))
        write(io, [x, y])
        imt = Mmap.mmap(ixy, Matrix{UInt8}, (x, y), 24)
        inc = 2
        for i = 1:y
            write(io, _type(hdr.type).(imt[:, i] .+ inc))
            inc += 2
        end
    end
end
