"""
    fileIRM(
    file::AbstractString,
    xy::AbstractString,
    chp::AbstractVector{Bool},
    id::AbstractVector{Int64};
    ε = 0.0
)
Read IRM from `file` if it exists and of the correct size, otherwise calculate
IRM with the genotypes in `xy`, chip information `chp`, and individual IDs `id`.
New IRM is saved to `file` for future use. Return the IRM.
"""
function fileIRM(
    file::AbstractString,
    xy::AbstractString,
    chp::AbstractVector{Bool},
    id::AbstractVector{Int64};
    ε = 0.0,
)
    n = length(id)
    if isfile(file) && filesize(file) == sizeof(Float64) * n * n
        G = zeros(Float64, n, n)
        read!(file, G)
        return G
    end
    G = irm(xy, chp, id) + ε * I
    write(file, G)
    return G
end
