"""
    append!(ixy::T, oxy::T) where T <: AbstractString
Append matrix `oxy` to `ixy`. The two files must have the same header. Their
number of rows must be the same. This applies only to full matrix.
"""
function append!(ixy::T, oxy::T) where {T<:AbstractString}
    hi, ho = header(ixy), header(oxy)
    isEqual(hi, ho) || error("Header of $ixy and $oxy are different")
    hi.flus == Int8('F') || error("Needs to be a full matrix")
    (mi, ni), (mo, no) = dim(ixy), dim(oxy)
    mi == mo || error("Number of rows of $ixy and $oxy are different")
    dim!(ixy, mi, ni + no)
    open(ixy, "a") do io
        write(io, Mmap.mmap(oxy, Matrix{_type(hi.type)}, (mo, no), 24))
    end
end

"""
    append!(fxy::AbstractString, gt::Matrix)
Append matrix `gt` to `fxy`.
"""
function append!(fxy::AbstractString, gt::Matrix)
    hdr = header(fxy)
    isnothing(hdr) && error("$fxy is not a valid xy file.")
    _type(hdr.type) == eltype(gt) || error("Element type of $fxy and gt are different")
    hdr.flus == Int8('F') || error("$fxy is not a full matrix")
    (x, y) = dim(fxy)
    x == size(gt, 1) || error("Number of rows of $fxy and gt are different")
    open(fxy, "a") do io
        write(io, gt)
    end
    dim!(fxy, x, y + size(gt, 2))
end
