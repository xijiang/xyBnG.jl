function dosblup(;
        data = "rst",
        baseDir = "base",
        testDir = "allschmes",
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
        sim_base(species, base)
    else
        desc = readlines("$base/desc.txt")
        parse(Int, desc[2]) < species.nid || desc[1] ≠ species.name && sim_base(cattle, baseDir)
    end
    sname = desc[1]
    schemes = [aaocs, iiocs, iidos, ggocs, agocs, igocs, gblup, ablup, iblup]
    lbls = ["aaocs", "iiocs", "iidos", "ggocs", "agocs", "igocs", "gblup", "ablup", "iblup"]
    F0 = 0.027

    # Simulations
    for irpt in 1:nrpt
        tag = lpad(irpt, ndigits(nrpt), '0')
        println()
        @info "==========> Repetition: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"

        sample_founder(base, test, plan.noff, nchp, nref, trait)
        for ext in ["ped", "lmp"]
            mv("$test/$sname.$ext", "$test/founder.$ext", force=true)
        end
        uniq("$test/$sname.xy", "$test/founder.xy")
        rm("$test/$sname.xy", force=true)
        lmp = deserialize("$test/founder.lmp")

        # The starting point: random selection
        randbrd(test, "founder", "$tag-rand", lmp, nrng, trait, plan; ibd=true)
        for i in 1:9
            scheme = schemes[i]
            lbl = lbls[i]
            if i < 7 # OCS
                scheme(test, "$tag-rand", "$tag-$lbl", lmp, nsel, trait, fixed, plan.noff, dF, F0)
            else
                scheme(test, "$tag-rand", "$tag-$lbl", lmp, nsel, trait, fixed, plan)
            end
            summary = Sum.xysum("$test/$tag-$lbl.ped", "$test/$tag-$lbl.xy", lmp, trait, nrng + 1)
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
