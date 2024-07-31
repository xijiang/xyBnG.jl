function reproduce!(
    ng::DataFrame,
    ped::DataFrame,
    xy::AbstractString,
    lmp::DataFrame,
    trts::Trait...,
)
    @debug "Reproduce the $(nrow(ng)) offspring"
    lms = sumMap(lmp)
    hdr, (nlc, nhp) = XY.header(xy), XY.dim(xy)
    etp = XY._type(hdr.type)
    pg = mmap(xy, Matrix{etp}, (nlc, nhp), 24)
    og = zeros(etp, nlc, 2nrow(ng))
    drop(pg, og, [ng.sire ng.dam], lms)
    XY.append!(xy, og)
    for trt in trts
        qtl = lmp[!, trt.name]
        a = lmp[!, trt.name*"_a"][qtl]
        d = lmp[!, trt.name*"_d"][qtl]
        qh = view(og, qtl, :)
        qg =
            hdr.u == 0 ? qh[:, 1:2:end] + qh[:, 2:2:end] :
            isodd.(qh[:, 1:2:end]) + isodd.(qh[:, 2:2:end])
        dg = qg .== 1
        ng[!, "tbv_"*trt.name] = qg'a
        ng[!, "gt_"*trt.name] = ng[!, "tbv_"*trt.name] + dg'd
    end
    append!(ped, ng, cols = :subset)
end
