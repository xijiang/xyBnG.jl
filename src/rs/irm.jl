"""
    mibd(a::T, b::T, c::T, d::T) where T<:AbstractVector

Given the 4 haplotypes of a pair of individuals, calculate the mean IBD. This
function ignores boundery check.
"""
function mibd(a::T, b::T, c::T, d::T) where T<:AbstractVector
    r  = sum(a .== c)
    r += sum(a .== d)
    r += sum(b .== c)
    r += sum(b .== d)
    r /= 2length(a)
end

"""
    irm(xy::AbstractString, loc, ids)

Uses the IBD info stored in the `xy` file to generate a mean IBD matrix.

- `xy` is the filename of `XY` format of haplotypes, loci majored,
- `loc` is a Bool vector, specifies which loci to be used,
- `ids` specifies the ID whose gamete matrix is to be generated.
"""
function irm(xy::AbstractString,
    loc::AbstractVector{Bool},
    ids::AbstractVector{Int},
    )
    hdr, (nlc, nhp), nid = XY.header(xy), XY.dim(xy), length(ids)
    hdr.u == 1 || error("Not a unique SNP file")
    1 ≤ minimum(ids) ≤ maximum(ids) ≤ nhp ÷ 2 || error("ID number out of range")
    length(loc) ≠ nlc && error("Loci number not match")
    type = XY._type(hdr.type)
    gt = mmap(xy, Matrix{type}, (nlc, nhp), 24)
    IBD = zeros(nid, nid)
    Threads.@threads for (i, j) in [(i, j) for i in 1:nid for j in 1:i]
        IBD[i, j] = mibd(view(gt, loc, 2i-1),
                         view(gt, loc, 2i),
                         view(gt, loc, 2j-1),
                         view(gt, loc, 2j))
        IBD[j, i] = IBD[i, j]
    end
    IBD
end

"""
    irm(xy::AbstractString, loc, ids, jds)
Calculate mean IBD relationships between `ids` and `jds` using the IBD
information stored in the `xy` file. It them average the 4 cells to calculate
the relationships between ID based on IBD information. In reality, these IBD
information can be obtained from dense genotypes very accurately.
"""
function irm(xy::AbstractString, # uniquely coded genotypes
    loc::AbstractVector{Bool}, # specify which loci to be used
    ids::AbstractVector{Int}, # specify which IDs to be used
    jds::AbstractVector{Int}, # specify which IDs to be used
    )
    issorted(ids) || sort!(ids)
    issorted(jds) || sort!(jds)
    mat, mid, nid = XY.mapit(xy), length(ids), length(jds)
    (length(loc) ≠ size(mat, 1) || 
     2ids[end] > size(mat, 2) ||
     2jds[end] > size(mat, 2)) && error("Not number of loci or IDs")
    ihp = sort([2ids .- 1; 2ids])
    jhp = sort([2jds .- 1; 2jds])
    igt = copy(mat[loc, ihp])
    jgt = copy(mat[loc, jhp])
    IBD = zeros(mid, nid)
    Threads.@threads for (i, j) in [(i, j) for i in 1:mid for j in 1:nid]
        IBD[i, j] = mibd(view(igt, :, 2i-1), view(igt, :, 2i), view(jgt, :, 2j-1), view(jgt, :, 2j))
    end
    IBD
end

"""
    xirm(G, fxy, loci, mid, nid)
When a new generation is generated with `reproduce!`, this function update IRM
`G` by just calculating the relationships between the (`mid`) old and new
`(nid-mid)` ID, and those among the new ID, using the bool vector `loci`, which
specifies which loci are used in file `fxy`. It then expand `G` with the new
relationships.
"""
function xirm(G::Matrix{Float64}, fxy::AbstractString,
              loci::AbstractVector{Bool},
              mid::Int, nid::Int)
    M = zeros(nid, nid)
    ra, rb = 1:mid, mid+1:nid

    copyto!(view(M, ra, ra), G)
    T = irm(fxy, loci, ra, rb)
    copyto!(view(M, ra, rb), T)
    copyto!(view(M, rb, ra), T')
    T = irm(fxy, loci, rb)
    copyto!(view(M, rb, rb), T)
    M
end
