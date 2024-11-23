import xyBnG.Util: blksz

"""
    grm(gt::AbstractArray; p::Union{Bool, AbstractVector{Float64}} = false)

Given the genotypes of `Matrix{Int8}`, this function calculate the genomic 
relationship matrix `GRM`. If the matrix is too big, the content will be 
calculated block by block and written to a file, else it will return a 
`Matrix{Float64}`.

By default the matrix needs to be 'nlc by nid' to speed up the calculation.
Such a matrix stores the genotypes continuously, which can speed up the 
matrix multiplication very much.  Another reason is that genotypes are
usually stored continuous for each individual.  They can can be read 
continuously in the `gt` columns.

If a `δ`, e.g., `δ = 0.01`, to the diagonals, you have to do this after this
function.

Important update: the function now can handle the case when the allele
frequencies are given.  In this case, the `p` must be a vector of allele
frequencies
"""
function grm(gt::AbstractArray; p::Union{Bool,AbstractVector{Float64}} = false)
    nlc, nid = size(gt)
    length(p) ≠ nlc && (p = mean(gt, dims = 2) / 2) # allele frequencies
    #p = mean(gt, dims = 2) / 2 # allele frequencies
    d = 2(1 .- p)'p             # the denominator
    mem = memavail() * 99 ÷ 100 # not all available memory
    gmt = nid^2 * 8             # memory by G
    zmt = nid * nlc * 8         # memory by Z
    if gmt + zmt < mem          # brute force
        # @info "G and Z are stored in memory"
        Z = gt .- 2p
        G = Z'Z ./ d
        return G
    else                        # minimal memory mode
        c1 = 2gt'p
        c2 = 4p'p
        if gmt < mem            # G can still be fit
            @info "only G were stored in memory"
            G = zeros(nid, nid)
            matmul!(G, gt', gt)
            G .-= c1
            G .-= c1'
            G .+= c2
            G ./= d
            return G
        else                            # G is too large
            file = basename(tempname()) # will write the result in pwd.
            @warn "G is too big to fit in memory. It is being writting into $file.
              False will be returned. $file can be read back in to memory, if enough,
              with `QTL.MIO.readmat($file)`"
            # ToDo: check disk space here
            m = mem ÷ 8 ÷ nid
            m = blksz(nid, m) # determine number of ID to be dealed a time
            stops = collect(m:m:nid)
            stops[end] == nid || push!(stops, nid)
            start = 1
            open(file, "w") do io
                write(io, [nid, nid, 13]) # 13 for Float64
                for stop in stops
                    sg = zeros(nid, stop - start + 1)
                    matmul!(sg, gt', gt[:, start:stop])
                    sg .-= c1
                    sg .-= c1[start:stop]'
                    sg .+= c2
                    sg ./= d
                    write(io, sg)
                    start = stop + 1
                end
            end
            return file
        end
    end
end

"""
    grm(xy::AbstractArray, loci)
Calculate the genomic relationship matrix `GRM` for the genotypes `xy` at
`loci`.
"""
function grm(xy::AbstractString, loci)
    hdr = XY.header(xy)
    hap = XY.mapit(xy)
    gt =
        hdr.u == 0 ? hap[loci, 1:2:end] + hap[loci, 2:2:end] :
        isodd.(hap[loci, 1:2:end]) + isodd.(hap[loci, 2:2:end])
    grm(gt)
end

"""
    grm(xy::AbstractArray, loci, frq)
Calculate the genomic relationship matrix `GRM` for the genotypes `xy` at `loci`
with allele frequencies `frq` of size(gt, 1).
"""
function grm(xy::AbstractString, loci, frq)
    hdr = XY.header(xy)
    hap = XY.mapit(xy)
    gt =
        hdr.u == 0 ? hap[loci, 1:2:end] + hap[loci, 2:2:end] :
        isodd.(hap[loci, 1:2:end]) + isodd.(hap[loci, 2:2:end])
    grm(gt, p = frq[loci])
end

"""
    grm(xy::AbstractString, loci, ids, frq)
Calculate the genomic relationship matrix `GRM` for the genotypes `xy` at `loci`
with allele frequencies `frq` of size(gt, 1) for the individuals `ids`.
"""
function grm(xy::AbstractString, loci, ids, frq)
    hdr = XY.header(xy)
    hap = XY.mapit(xy)
    hs = sort(2ids .- 1; 2ids)
    gt =
        hdr.u == 0 ? hap[loci, hs] + hap[loci, hs] :
        isodd.(hap[loci, hs]) + isodd.(hap[loci, hs])
    grm(gt, p = frq[loci])
end
