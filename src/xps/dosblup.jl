function dosblup(;
        data = "rst",
        baseDir = "base",
        testDir = "dosblup",
        species = Cattle(5_000),
        trait = Trait("milk", 0.25, 10000; sex = 0),
        nchp = 50_000,
        nref = 10_000,
        nrng = 5,
        nsel = 20,
        plan = Plan(25, 50, 200),
        fixed = ["grt"],
        dF = 0.011,
        nrpt = 1)

    # Scenario recording
    base, test = "$data/$baseDir", "$data/$testDir"
    isdir("$test") || mkpath("$test")
    scenario = (Data=data, BaseDir = baseDir, TestDir = testDir, Species = species,
                Trait = trait, Nchp = nchp, Nref = nref, Nrng = nrng, Nsel = nsel,
                Plan = plan, Fixed = fixed, ΔF = dF, Nrpt = nrpt,)
    savepar(scenario, "$test/scenario.par")

    # Prepare a base population
    if !isfile("$base/desc.txt")
        sim_base(species, base)
    else
        desc = readlines("$base/desc.txt")
        parse(Int, desc[2]) < species.nid || desc[1] ≠ species.name && sim_base(cattle, baseDir)
    end
    sname = desc[1]

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
        F0 = 0.027


        # The starting point: random selection
        randbrd(test, "founder", "$tag-rand", lmp, nrng, trait, plan; ibd=true)
        # OCS
        aaocs(test, "$tag-rand", "$tag-aaocs", lmp, nsel, trait, fixed, plan.noff, dF, F0)
        iiocs(test, "$tag-rand", "$tag-iiocs", lmp, nsel, trait, fixed, plan.noff, dF, F0)
        iidos(test, "$tag-rand", "$tag-iidos", lmp, nsel, trait, fixed, plan.noff, dF, F0)
        #tgocs(test, "$tag-rand", "$tag-tgocs", lmp, nsel, trait, fixed, plan.noff, dF, F0)
        ggocs(test, "$tag-rand", "$tag-ggocs", lmp, nsel, trait, fixed, plan.noff, dF, F0)
        agocs(test, "$tag-rand", "$tag-agocs", lmp, nsel, trait, fixed, plan.noff, dF, F0)
        igocs(test, "$tag-rand", "$tag-igocs", lmp, nsel, trait, fixed, plan.noff, dF, F0)
        # Ordinary BLUP
        gblup(test, "$tag-rand", "$tag-gblup", lmp, nsel, trait, fixed, plan)
        ablup(test, "$tag-rand", "$tag-ablup", lmp, nsel, trait, fixed, plan)
        iblup(test, "$tag-rand", "$tag-iblup", lmp, nsel, trait, fixed, plan)
    end
    open("$test/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end
