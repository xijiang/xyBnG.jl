"""
    function phenotype!(xy::AbstractString, ped::DataFrame, lmp::DataFrame,
        trts::Trait...) # method 01

Phenotype trait `trts` with QTL information in `lmp` and SNP information in
`xy`. The phenotype columns about the traits will be added/updated to `ped`. The
phenotype column name for a trait `trt` is `ft_\$(trt.name)`.
"""
function phenotype!(xy::AbstractString, ped::DataFrame, lmp::DataFrame,
                    trts::Trait...)
    @info "under construction"
    hdr, (nlc, nhp) = XY.header(xy), XY.dim(xy)
    #hps = Mmap.mmap(xy, Matrix{XY._type(hdr.type)}, (nlc, nhp), 24)
    ped
end

"""
    function phenotype!(xy::AbstractString, ped::AbstractString,
        lmp::AbstractString, trts::Trait...)

This method reads `ped` and `lmp` from files and call `phenotype!` method 01.
"""
function phenotype!(xy::AbstractString, ped::AbstractString,
                    lmp::AbstractString, trts::Trait...)
    ped = deserialize(ped)
    lmp = deserialize(lmp)
    phenotype!(xy, ped, trts...)
end

"""
    function phenotype!(ID::T, xy::AbstractString, ped::DataFrame, lmp::DataFrame,
        trts::Trait...) where T <: AbstractVector{Int} # method 02

As method phenotype! 01, but only phenotype the ID in `ID`.
"""
function phenotype!(ID::T, xy::AbstractString, ped::DataFrame, lmp::DataFrame,
                    trts::Trait...) where T <: AbstractVector{Int}
    @info "under construction"
    hdr, (nlc, nhp) = XY.header(xy), XY.dim(xy)
    #hps = Mmap.mmap(xy, Matrix{XY._type(hdr.type)}, (nlc, nhp), 24)
    ped
end