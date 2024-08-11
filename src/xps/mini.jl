"""
    mini()
This function is to run a very small example of the package `xyBnG`. It is
for my laptop to run in a reasonable time.
"""
function mini(;
    data = "rst",
    baseDir = "tskit",
    testDir = "mini",
    species = Cattle(5_000),
    trait = Trait("milk", 0.25, 1_000; sex = 0),
    nchp = 3_000,
    nref = 1_000,
    nrng = 5,
    nsel = 20,
    plan = Plan(5, 10, 40),
    fixed = ["grt"],
    dF = 0.011,
    nrpt = 1,
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
    #schemes = (aaocs, iiocs, ggocs, agocs, igocs, gblup, ablup, iblup)
    cullSchemes = (gblup, ablup)
    F0 = 0.027

    # Simulations
    for irpt = 1:nrpt
        tag = lpad(irpt, ndigits(nrpt), '0')
        println()
        @info "==========> Repeat: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"

        sample_xy(fxy, fmp, test, plan.noff, maf, nchp, nref, trait)
        mv("$test/founder.xy", "$test/snp.xy", force = true)
        uniq("$test/snp.xy", "$test/founder.xy")
        lmp = deserialize("$test/founder.lmp")
        # The starting point: random selection
        randbrd(test, "founder", "$tag-rand", lmp, nrng, trait, plan; ibd = true)
        for scheme in cullSchemes
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(test, foo, bar, lmp, nsel, trait, fixed, plan)
            summary = Sum.xysum("$test/$bar.ped", "$test/$bar.xy", lmp, trait)
            Sum.savesum("$test/summary.ser", summary)
        end
    end
    open("$test/scenario.par", "a") do io
        println(io, "Ended: ", time())
    end
end
