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
    keep = true,
)
    # Scenario recording
    base, test = "$data/$baseDir", "$data/$testDir"
    isdir("$test") && rm("$test", force = true, recursive = true)
    mkpath("$test")
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

    chkbase(base, species) # Prepare/verify a base population
    fxy, fmp, maf = "$base/$(species.name).xy", "$base/$(species.name).lmp", 0.0
    #schemes = (aaocs, iiocs, ggocs, agocs, igocs, gblup, ablup, iblup)
    culls = (gblup, ablup)
    # Simulations
    for irpt = 1:nrpt
        tag = lpad(irpt, ndigits(nrpt), '0')
        println()
        @info "==========> Repeat: $tag / $nrpt <=========="
        @info "  - Prepare a founder population"

        lmp, F0 = initPop(fxy, fmp, test, plan, maf, nchp, nref, nrng, trait, tag)
        for scheme in culls
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            scheme(test, foo, bar, lmp, nsel, trait, fixed, plan)
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
