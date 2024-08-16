"""
    rumigen(;
        data = "rst",
        baseDir = "tskit",
        testDir = "cattle",
        species = Cattle(5_000),
        trait = Trait("milk", 0.25, 10000),
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
Function to repeat the simulation for the rumigen project. This is to verify if
my new codes are right. Several issues need to pay attention to.
    
Previously, I forgot to mask the male phenotypes. The number of cadidates, i.e.,
`c` is non-zero, will have more males than females. This is because that EBV for
males are not as accurate as those of females. Male EBV has thus smaller
vairance. 

In the new version of `xyBnG`, I was trying to limit the number of males as set
in `plan`. So I have two `plan`s in this function. If OCS, then `plan` 2 asks
for 50 males and 50 females. The real numbers selected in the simulation are
around 30. `plan` 2 guanruntees the candidates are to be selected.

The other issue is about the GRM. Before I was using allele frequencies accross
generations. Here I only use that of the generation when directional selection
started. Hence selection schemes involving GRM will be different. To verify the
codes, only OCS without GRM are to be considered.
"""
function rumigen(;
    data = "rst",
    baseDir = "tskit",
    testDir = "cattle",
    species = Cattle(5_000),
    trait = Trait("milk", 0.25, 10000),
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
    schemes = (aaocs, iiocs, ggocs, agocs, igocs, gblup, ablup, iblup)
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
        Schemes = schemes,
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
    pln2 = Plan(50, 50, 200)

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
        for scheme in schemes
            foo, bar = "$tag-rand", tag * '-' * string(scheme)
            if occursin("blup", bar)
                scheme(test, foo, bar, lmp, nsel, trait, fixed, plan)
            else
                scheme(test, foo, bar, lmp, nsel, trait, fixed, pln2, dF, F0)
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

"""
    run_rumigen(op)
Run the `rumigen` function to repeat previous results.
"""
function run_rumigen(op; dir = "/mnt/a/store/xybng")
    if op == 1
        # Simulation started: 2024-08-09T10:55:59.562
        # Simulation finised: 2024-08-11T20:19:59.645
        rumigen(data = dir, testDir = "rumigen/01", nrpt = 100)
    elseif op == 2
        # Simulation started: 2024-08-12T09:54:56.039
        # Simulation finised: 2024-08-12T18:23:53.339
        rumigen(data = dir, testDir = "rumigen/02", nrpt = 100, nrng = 0, nsel = 10)
    elseif op == 3
        rumigen(data = dir, testDir = "rumigen/03", nrpt = 20)
        rumigen(data = dir, testDir = "rumigen/04", nrpt = 20, nrng = 0, nsel = 10)
    else
        @info "1 ≤ op ≤ 2. No operation is performed."
    end
end
