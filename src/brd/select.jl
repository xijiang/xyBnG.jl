"""
    function mate(pas::Vector{T}, mas::Vector{T}, plan::Plan, ped::DataFrame) where T <: Integer

Mate `pas`, `mas` according to `plan` to produce `noff` offspring in DataFrame
`ped`.
"""
function mate(pas::AbstractVector{T}, mas::AbstractVector{T}, plan::Plan) where T <: Integer
    if plan.mate == :random
        sire = rand(pas, plan.noff)
        dam  = rand(mas, plan.noff)
    else
        st = Int(ceil(plan.noff / plan.npa))
        dt = Int(ceil(plan.noff / plan.nma))
        if plan.mate == :hierarchical
            sire = repeat(shuffle(pas), inner = st)[1:plan.noff]
            dam  = repeat(mas, inner = dt)[1:plan.noff]
        elseif plan.mate == :factorial
            sire = repeat(shuffle(pas), outer = st)[1:plan.noff]
            dam  = repeat(mas, outer = dt)[1:plan.noff]
        else
            error("mate must be in [:random, :hierarchical, :factorial]")
        end
    end
    pm = sortslices([sire dam], dims=1, by=x -> (x[1], x[2]))
    DataFrame(
        id = 1:plan.noff,
        sire = pm[:, 1],
        dam = pm[:, 2],
        sex = rand(Int8.(0:1), plan.noff),
        grt = 1,
    )
end

"""
    Select(ID::AbstractVector{T}, plan::Plan, ped::DataFrame, trt::Trait) where T <: Integer

Selection `ID` on their EBV of trait `trt` in DataFrame `ped` according to plan
`plan`.
"""
function Select(ID::AbstractVector{T}, plan::Plan, ped::DataFrame, trt::Trait) where T <: Integer
    @debug "Selection on $(trt.name)"
    df = select(view(ped, ID, :), :id, :sex, r"ebv_")
    sort!(df, "ebv_" * trt.name, rev = trt.rev)
    gps = groupby(df, :sex)
    mas = gps[1].id[1:plan.nma]
    pas = gps[2].id[1:plan.npa]
    ng  = mate(pas, mas, plan)
    ng.id .+= nrow(ped)
    ng.grt .+= maximum(ped.grt[ID])
    ng
end

"""
    Select(ID::AbstractVector{T}, plan::Plan, ped::DataFrame, dic::Dict; rev = true) where T <: Integer

Selection `ID` on their weighted EBV of traits in dictionary `dic`.
"""
function Select(ID::AbstractVector{T}, plan::Plan, ped::DataFrame, dic::Dict; rev = true) where T <: Integer
    @debug "Selection on $(join(keys(dic), ", "))"
    df = select(view(ped, ID, :), :id, :sex, r"ebv_")
    index = zeros(nrow(df))
    for trt in keys(dic)
        index += df[!, "ebv_" * trt] * dic[trt]
    end
    df.index = index
    sort!(df, :index, rev = rev)
    gps = groupby(df, :sex)
    mas = gps[1].id[1:plan.nma]
    pas = gps[2].id[1:plan.npa]
    ng = mate(pas, mas, plan)
    ng.id .+= nrow(ped)
    ng.grt .+= ped.grt[end]
    ng
end

function Select(ID::AbstractVector{T}, plan::Plan, ped::DataFrame, c::AbstractVector{Float64}) where T <: Integer
    @debug "Select `ID` according to their contribution `c`"
    pm = begin
        df = DataFrame(id = ID, sex = ped.sex[ID], c = c,
                        n = Int.(round.(c * plan.noff)))
        df = view(df, df.c .> 0, :)
        groupby(df, :sex)
    end
    if plan.mate == :random
        pa = sample(pm[2].id, Weights(pm[2].c), plan.noff)
        ma = sample(pm[1].id, Weights(pm[1].c), plan.noff)
    else
        pa, ma = Int[], Int[]
        # Note that the number of offspring may not be exactly `noff`
        for (id, _, _, n) in eachrow(pm[2])
            append!(pa, fill(id, n))
        end
        for (id, _, _, n) in eachrow(shuffle(pm[1]))
            append!(ma, fill(id, n))
        end
        pa = sort(pa[1:plan.noff])
        ma = ma[1:plan.noff] # already for hierarchical
        if plan.mate == :factorial
            ma = shuffle(ma)
        end
    end
    nid = nrow(ped)
    DataFrame(
        id = nid + 1:nid + plan.noff,
        sire = pa,
        dam = ma,
        sex = rand(Int8.(0:1), plan.noff),
        grt = maximum(ped.grt[ID]) + 1,
    )
end

#=
"""
    Select!(ped::DataFrame, ppo::Tuple{Int, Int, Int}, trt::Trait; mate = :random, rev = true)

Selection `npa = ppo[1]`, `nma = ppo[2]` on EBV of trait `trt` in DataFrame
`ped` to reproduce `noff = ppo[3]` offspring. The `mate` can be `:random` or
`:balanced`. The new generation will be appended to `ped`. If `rev` is true,
which is by default, select the parents of the highest EBV. Otherwise, select
the parents of the lowest EBV.
"""
function Select!(ped::DataFrame, ppo::Tuple{Int, Int, Int}, trt::Trait; mate = :random, rev=true)
    npa, nma, noff = ppo
    tgt = groupby(ped, :grt)[end]
    ebv = "ebv_" * trt.name
    ♀s, ♂s = groupby(select(tgt, :id, :sex, ebv), :sex)
    sirs = sort(♂s, ebv, rev = rev).id[1:npa]
    dams = sort(♀s, ebv, rev = rev).id[1:nma]
    pair!(sirs, dams, noff, ped, mate)
end

"""
    Select!(ped::DataFrame, ppo::Tuple{Int, Int, Int}, dic::Dict{String, Float64}; mate = :random, rev = true)

Selection on dictionary `Dict{String, Float64}` `dic`, which defines the traits
to select on, and their weight. If the key is in the pedigree, then that column
is included in the weighted sum. If the key is not in the pedigree, then
`ebv_key` is included in the weighted sum. Otherwise, an error occurrs.
"""
function Select!(ped::DataFrame, ppo::Tuple{Int, Int, Int},
    dic::Dict{String, Float64}; mate = :random, rev = true)
    npa, nma, noff = ppo
    tgt = groupby(ped, :grt)[end]
    vals = zeros(nrow(tgt))
    for c in keys(dic)
        if c ∈ names(tgt)
            vals += tgt[!, c] * dic[c]
        elseif "ebv_" * c ∈ names(tgt)
            vals += tgt[!, "ebv_" * c] * dic[c]
        else
            error("$c or ebv_$c not found in `ped`")
        end
    end
    tdf = DataFrame(id = tgt.id, sex = tgt.sex, vals = vals)
    ♀s, ♂s = groupby(tdf, :sex)
    sirs = sort(♂s, :vals, rev = rev).id[1:npa]
    dams = sort(♀s, :vals, rev = rev).id[1:nma]
    pair!(sirs, dams, noff, ped, mate)
end

"""
    function Select!(ped::DataFrame, ppo::Tuple{Int, Int, Int}, dic::Dict{Symbol, Float64}; mate = :random)

Convert `Symbol` keys to `String` keys and call `Select!(ped, ppo, dic; mate = mate)`.
"""
function Select!(ped::DataFrame, ppo::Tuple{Int, Int, Int},
    dic::Dict{Symbol, Float64}; mate = :random)
    dd = Dict(zip(string.(keys(dic)), values(dic)))
    Select!(ped, ppo, dd; mate = mate)
end

# under construction
function Select!(ped::DataFrame, ppo::Tuple{Int, Int, Int},
    trt::Trait, on::Union{AbstractString, Symbol}; mate = :random)
    @debug "Selection on $on of $(trt.name)"
end

function Select!(ped::DataFrame, ppo::Tuple{Int, Int, Int},
    c::AbstractVector{Float64}; mate = :random)
    @debug "Selection last generation of `ped` with contribution $c"
end
=#