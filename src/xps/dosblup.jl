function dosblup(;
        data = "rst",
        baseDir = "base",
        testDir = "dosblup",
        species = Cattle("cattle", 5_000),
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
    savepar(Dict(:Data => data, :BaseDir => baseDir, :TestDir => testDir,
            :Species => species, :Trait => trait, :Nchp => nchp, :Nref => nref,
            :Nrng => nrng, :Nsel => nsel, :Plan => plan, :Fixed => fixed, :dF => dF, :nrpt => nrpt),
            "$test/scenario.par")
    open("$test/par.csv", "w") do io
        println(io, "Parameter,Value")
        println(io, "Species,$(species.name)")
        println(io, "Base pop size,$(species.nid)")
        println(io, "Trait,$(trait.name)")
        println(io, "No. QTL,$(trait.nQTL)")
        println(io, "Heritability,$(trait.h2)")
        println(io, "Sex,$(trait.sex)")
        println(io, "No. chip SNP,$nchp")
        println(io, "No. reference SNP,$nref")
        println(io, "No. rnd. sel. generations,$nrng")
        println(io, "No. directional sel. generations,$nsel")
        println(io, "No. of selected sires,$(plan.npa)")
        println(io, "No. of selected dams,$(plan.nma)")
        println(io, "No. of offspring,$(plan.noff)")
        println(io, "Mating type,$(plan.mate)")
        println(io, "Fixed effects,$(join(fixed, " "))")
    end

    # Prepare a base population
    if !isfile("$base/desc.txt")
        sim_base(cattle, base)
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

        @goto founder
        sample_founder(base, test, plan.noff, nchp, nref, trait)
        for ext in ["ped", "lmp"]
            mv("$test/$sname.$ext", "$test/founder.$ext", force=true)
        end
        uniq("$test/$sname.xy", "$test/founder.xy")
        rm("$test/$sname.xy", force=true)
        @label founder

        lmp = deserialize("$test/founder.lmp")
        F0 = 0.027


        # The starting point: random selection
        #randbrd(test, "founder", "$tag-rand", lmp, nrng, trait, plan)
        #aaocs(test, "$tag-rand", "$tag-aablup", lmp, nsel, trait, fixed, plan.noff, dF, F0)
        #iiocs(test, "$tag-rand", "$tag-iiblup", lmp, nsel, trait, fixed, plan.noff, dF, F0)
        #iidos(test, "$tag-rand", "$tag-iidos",  lmp, nsel, trait, fixed, plan.noff, dF, F0)
        #ggocs(test, "$tag-rand", "$tag-ggblup", lmp, nsel, trait, fixed, plan.noff, dF, F0)
        #agocs(test, "$tag-rand", "$tag-agblup", lmp, nsel, trait, fixed, plan.noff, dF, F0)
        #igocs(test, "$tag-rand", "$tag-igblup", lmp, nsel, trait, fixed, plan.noff, dF, F0)
        #gblup(test, "$tag-rand", "$tag-gblup", lmp, nsel, trait, fixed, plan)
        ablup(test, "$tag-rand", "$tag-gblup", lmp, nsel, trait, fixed, plan)
    end
end
