#=
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
=#
"""
    phenotype!(ID::AbstractVector{T}, ped::DataFrame, trts::Trait...) where T <: Integer

As method phenotype! 01, but only phenotype the ID in `ID`.
"""
function phenotype!(ID::AbstractVector{T}, ped::DataFrame, trts::Trait...) where T <: Integer
    @debug "Generate phenotypes of trait $(join(name.(trts), ", "))"
    for trt in trts
        if trt.type == Int 
            @info "Non continuous trait not implemented"
            continue
        end
        ft = "ft_$(trt.name)"
        gt = "gt_$(trt.name)"
        if !hasproperty(ped, ft)
            tmp = Vector{Union{Missing, trt.type}}(missing, nrow(ped))
            ped[!, ft] = tmp
        end
        sde = sqrt(1.0 / trt.h2 - 1.0 - trt.vd)
        id = zeros(Bool, nrow(ped))
        id[ID] .= true
        if trt.sex == 2 # all ID have phenotypes
            tgt = view(ped, id, ft)
        else
            id = zeros(Bool, nrow(ped))
            id[ID] .= true
            id = id .&& ped.sex .== trt.sex
            tgt = view(ped, id, ft)
        end
        copyto!(tgt, rand(Normal(0, sde), sum(id)) + view(ped, id, gt))
    end
    ped
end