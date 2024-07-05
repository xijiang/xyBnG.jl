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
        keep = true)

    # Scenario recording
    base, test = "$data/$baseDir", "$data/$testDir"
    isdir("$test") || mkpath("$test")
    scenario = (Data=data, BaseDir = baseDir, TestDir = testDir, Species = species,
                Trait = trait, Nchp = nchp, Nref = nref, Nrng = nrng, Nsel = nsel,
                Plan = plan, Fixed = fixed, ΔF = dF, Nrpt = nrpt,)
    savepar(scenario, "$test/scenario.par")
    isfile("$test/summary.ser") && rm("$test/summary.ser", force=true)

    # Prepare a base population
    if !isfile("$base/desc.txt")
        ts_base(species, base)
    else
        desc = readlines("$base/desc.txt")
        parse(Int, desc[2]) < species.nid || desc[1] ≠ species.name && ts_base(cattle, baseDir)
    end
    sname = desc[1]
    fxy, fmp, maf = "$base/$sname.xy", "$base/$sname.lmp", 0.0
    schemes = (aaocs, iiocs, iidos, ggocs, agocs, igocs, gblup, ablup, iblup)

    F0 = 0.027

    # Simulations
    for irpt in 1:nrpt
        tag = lpad(irpt, ndigits(nrpt), '0')
        println()
        @info "==========> Repeat: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"

        sample_xy(fxy, fmp, test, plan.noff, maf, nchp, nref, trait)
        mv("$test/founder.xy", "$test/snp.xy", force=true)
        uniq("$test/snp.xy", "$test/founder.xy")
        lmp = deserialize("$test/founder.lmp")

        # The starting point: random selection
        randbrd(test, "founder", "$tag-rand", lmp, nrng, trait, plan; ibd=true)
        for scheme in schemes
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            occursin("blup", bar) ?
                scheme(test, foo, bar, lmp, nsel, trait, fixed, plan) :
                scheme(test, foo, bar, lmp, nsel, trait, fixed, plan.noff, dF, F0)
            summary = Sum.xysum("$test/$bar.ped", "$test/$bar.xy", lmp, trait, nrng + 1)
            Sum.savesum("$test/summary.ser", summary)
        end
        if !keep
            for f in readdir("$test")
                occursin(Regex("^$(tag)"), f) && rm("$test/$f", force=true)
            end
        end
    end
    open("$test/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end
