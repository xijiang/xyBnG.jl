function dosblup(;
    data = "rst",
    baseDir = "tskit",
    testDir = "cattle",
    species = Cattle(5_000),
    trait = Trait("milk", 0.25, 10000; sex = 0),
    nchp = 50_000,
    nref = 10_000,
    nrng = 5,
    nsel = 20,
    plan = Plan(25, 50, 200),
    fixed = ["grt"],
    dF = 0.011,
    nrpt = 1,
    keep = true,
)

    # Scenario recording
    base, test = "$data/$baseDir", "$data/$testDir"
    isdir("$test") || mkpath("$test")
    scenario = (
        Data = data,
        BaseDir = baseDir,
        TestDir = testDir,
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
    )
    savepar(scenario, "$test/scenario.par")
    isfile("$test/summary.ser") && rm("$test/summary.ser", force = true)

    # Prepare a base population
    if !isfile("$base/desc.txt")
        ts_base(species, base)
    else
        desc = readlines("$base/desc.txt")
        parse(Int, desc[2]) < species.nid ||
            desc[1] ≠ species.name && ts_base(cattle, baseDir)
    end
    sname = desc[1]
    fxy, fmp, maf = "$base/$sname.xy", "$base/$sname.lmp", 0.0
    #schemes = (aaocs, iiocs, iidos, ggocs, agocs, igocs, gblup, ablup, iblup)
    schemes = (aaocs, iiocs, ggocs, agocs, igocs, gblup, ablup, iblup)

    # Simulations
    for irpt = 1:nrpt
        tag = lpad(irpt, ndigits(nrpt), '0')
        println()
        @info "==========> Repeat: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"

        lmp, F0 = initPop(fxy, fmp, test, plan, maf, nchp, nref, nrng, trait, tag, true)
        for scheme in schemes
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            if occursin("blup", bar)
                scheme(test, foo, bar, lmp, nsel, trait, fixed, plan)
            else
                scheme(test, foo, bar, lmp, nsel, trait, fixed, plan, dF, F0)
            end
            summary = Sum.xysum("$test/$bar.ped", "$test/$bar.xy", lmp, trait)
            Sum.savesum("$test/summary.ser", summary)
        end
        if !keep
            for f in readdir("$test")
                occursin(Regex("^$(tag)"), f) && rm("$test/$f", force = true)
            end
        end
    end
    open("$test/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end
