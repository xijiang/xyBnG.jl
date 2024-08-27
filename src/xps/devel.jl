function devel()
    data = "../pgsnp"
    rst = "$data/rst"
    species, trait = Cattle(5_000), Trait("growth", 0.25, 10_000)
    nchp, nref, nrng, nsel = 25_000, 5_000, 5, 20
    plan, fixed = Plan(25, 50, 200), ["grt"]
    dF, nrpt = 0.011, 50

    isdir(rst) && rm(rst, force = true, recursive = true)
    mkpath(rst)
    CULLS = (gblup, ablup, iblup)
    OCSS = (aaocs, iiocs, ggocs, agocs, igocs)
    scenario = (
        Data = data,
        Result = rst,
        Species = species,
        Trait = trait,
        Nchp = nchp,
        Nref = nref,
        Nrng = nrng,
        Nsel = nsel,
        Plan = plan,
        Fixed = fixed,
        ΔF = dF,
        Nrpt = nrpt,
        Schemes = union(CULLS, OCSS),
    )
    savepar(scenario, "$rst/scenario.par")

    for irpt = 1:nrpt
        tag = lpad(irpt, ndigits(nrpt), '0')
        fxy, fmp, maf =
            "$data/$irpt/$(species.name).xy", "$data/$irpt/$(species.name).lmp", 0.0
        @info "==========> Repeat: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"
        lmp, F0 = initPop(fxy, fmp, rst, plan, maf, nchp, nref, nrng, trait, tag)
        for scheme in CULLS
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(rst, foo, bar, lmp, nsel, trait, fixed, plan)
            summary = Sum.xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
            Sum.savesum("$rst/summary.ser", summary)
        end
        for scheme in OCSS
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(rst, foo, bar, lmp, nsel, trait, fixed, plan, dF, F0)
            summary = Sum.xysum("$rst/$bar.ped", "$rst/$bar.xy", lmp, trait)
            Sum.savesum("$rst/summary.ser", summary)
        end
    end
    open("$test/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end
