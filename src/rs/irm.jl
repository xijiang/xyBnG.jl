"""
    mibd(a::T, b::T, c::T, d::T) where T<:AbstractVector

Given the 4 haplotypes of a pair of individuals, calculate the mean IBD. This
function ignores boundery check.
"""
function mibd(a::T, b::T, c::T, d::T) where {T<:AbstractVector}
    r = sum(a .== c)
    r += sum(a .== d)
    r += sum(b .== c)
    r += sum(b .== d)
    r /= 2length(a)
end

"""
    irm(hps::AbstractMatrix)
Calculate the mean IBD matrix for a locus majored haplotype matrix.
"""
function irm(hps::AbstractMatrix)
    iseven(size(hps, 2)) || error("Not a haplotype matrix")
    nid = size(hps, 2) ÷ 2
    IBD = zeros(nid, nid)
    Threads.@threads for (i, j) in [(i, j) for i = 1:nid for j = 1:i]
        IBD[i, j] = mibd(hps[:, 2i - 1], hps[:, 2i], hps[:, 2j - 1], hps[:, 2j])
        IBD[j, i] = IBD[i, j]
    end
    IBD
end

"""
    irm(xy::AbstractString, loc::Vector{Bool}, id::AbstractVector{Int})

Uses the IBD info stored in the `xy` file to generate a mean IBD matrix.

- `xy` is the filename of `XY` format of haplotypes, loci majored,
- `loc` is a Bool vector, specifies which loci to be used,
- `id` specifies the ID whose gamete matrix is to be generated.
"""
function irm(xy::AbstractString, loc::Vector{Bool}, id::AbstractVector{Int})
    hdr, (nlc, nhp), nid = XY.header(xy), XY.dim(xy), length(id)

    hdr.u == 1 || error("Not a uniquely coded SNP file")
    1 ≤ id[1] ≤ id[end] ≤ nhp ÷ 2 || error("ID number out of range")
    length(loc) ≠ nlc && error("Loci number not match")

    type = XY._type(hdr.type)
    gt = mmap(xy, Matrix{type}, (nlc, nhp), 24)
    idx = sort([collect(2id .- 1); collect(2id)])
    irm(view(gt, loc, idx))
end

"""
    irm(xy::AbstractString loc::Vector{Bool}, id::UnitRange{Int}, jd::UnitRange{Int})

Calculate mean IBD relationships between `id` and `jd` using the IBD information
stored in the `xy` file. It them average the 4 cells to calculate the
relationships between ID based on IBD information. In reality, these IBD
information can be obtained from dense genotypes very accurately.
"""
function irm(
    xy::AbstractString, # uniquely coded genotypes
    loc::Vector{Bool}, # specify which loci to be used
    id::UnitRange{Int}, # specify which IDs to be used
    jd::UnitRange{Int}, # specify which IDs to be used
)
    hdr, (nlc, nhp) = XY.header(xy), XY.dim(xy)
    hdr.u == 1 || error("Not a uniquely coded SNP file")
    1 ≤ id[1] ≤ id[end] ≤ nhp ÷ 2 || error("ID number out of range")
    1 ≤ jd[1] ≤ jd[end] ≤ nhp ÷ 2 || error("ID number out of range")
    length(loc) == nlc || error("Loci number not match")
    gt = mmap(xy, Matrix{XY._type(hdr.type)}, (nlc, nhp), 24)

    IBD = zeros(length(id), length(jd))
    Threads.@threads for (i, j) in [(i, j) for i in eachindex(id) for j in eachindex(jd)]
        IBD[i, j] = mibd(
            view(gt, :, 2id[i] - 1),
            view(gt, :, 2id[i]),
            view(gt, :, 2jd[j] - 1),
            view(gt, :, 2jd[j]),
        )
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
function xirm(
    G::Matrix{Float64},
    fxy::AbstractString,
    loci::AbstractVector{Bool},
    mid::Int,
    nid::Int,
)
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
