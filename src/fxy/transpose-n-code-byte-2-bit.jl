"""
    tr8bit(fxy::AbstractString, bxy::AbstractString; bs = 2^16)
Transpose the matrix in `fxy`, and write to `bxy` as BitArray. The problem is
that reading `fxy` is not continuous. It takes a lot of time, espcielly when
`fxy` is on a spinning disk. This function convert the matrix by block. The
block is continous on disk. Hence it is dozens of times faster.

This function is specially designed for merged data from `msprime`, or, `MaCS`,
which are `n_Haplotyes x n_Loci` matrices of Int8. I don't do the type checking
here. Users are not supposed to use this function directly.
"""
function tr8bit(fxy::AbstractString, bxy::AbstractString; bs = 2^16)
    hdr, (nhp, nlc) = header(fxy), dim(fxy)
    hdr.r, hdr.major = 1, 0 # BitArray, and transposed to be locus majored
    m = mapit(fxy)
    @info "  - Transposing $fxy to $bxy as a BitArray:"
    t = nothing
    open(bxy, "w") do io
        write(io, Ref(hdr), [nlc, nhp])
        if nlc ≤ bs
            write(io, BitArray(m'))
        else
            t = BitArray(view(m, :, 1:bs)')
            rg = collect(2bs:bs:nlc)
            nlc ∈ rg || push!(rg, nlc)
            left, cnt = bs + 1, 1
            for right in rg
                print("\r$(lpad(cnt, 12)) / $(length(rg)) blocks")
                t = vcat(t, BitArray(view(m, :, left:right)'))
                left, cnt = right + 1, cnt + 1
            end
            write(io, t)
            println()
        end
    end
end
