"""
    a_tale_of_selection()

A tale of selection. Many of us are using script languages to do our work. A
script also means the written text of a play, movie, or broadcast. It tells a
story. Here is a story of selection, using my package of `xyBnG`.
"""
function a_tale_of_selection()
    ########################################
    # Scenarios, or characters in the play
    baseDir, tstDir = "rst/base", "rst/tale"
    cattle = Cattle("cattle", 5_000)
    milk   = Trait("milk",   0.3, 10000; sex = 0, vd = 0.1)
    growth = Trait("growth", 0.5, 10000; sex = 2, vd = 0.2)
    nid, nchp, nref = 200, 50_000, 10_000
    nrng, nsel = 5, 10
    plan = Plan(25, 50, 200) # default :hierarchical
    fixed = ["grt"]

    ########################################
    title = md"""# Simple selections on a small cattle population"""
    tprintln(title)
    synopsis = md"""## Synpsis:

    In this play, we will simulate a catle base population.
    We will then sample a founder population from it.
    Several selection schemes will be carried out on the
    founder population."""
    tprintln(synopsis, "\n")

    ########################################
    act = md"""
    ## Act I: simulate the base population

    Simulate a cattle base population begins here."""
    tprintln(act)

    isfile("$baseDir/desc.txt") || sim_base(cattle, baseDir)
    line = readlines("$baseDir/desc.txt")
    desc = md"""
    This population is named $(line[1]), and of $(line[2]) individuals
    in directory $baseDir.

    !!! tip "Sex specific traits"
        The sex prorperty for trait milk is 0, meaning that it is only expressed in cows.

    When you see above, you know the base population is ready.
    """
    tprintln(desc, "\n")

    ########################################
    act = md"""
    ## Act II: sample a founder population

    Sample a founder population from the base population.
    This may take a while.
    """
    tprintln(act, "\n")
    isfile("$tstDir/founder.xy") || sample_founder(baseDir, tstDir, nid, nchp, nref, milk, growth)
    for ext in ["ped", "lmp", "xy"]
        cp("$tstDir/cattle.$ext", "$tstDir/founder.$ext", force=true)
    end
    desc = md"""
    A total of $nid ID, or, $(2nid) haplotypes are sampled from the base population.
    These ID have $nchp chip SNPs, $nref reference/hidden SNPs, $(milk.nQTL) QTL
    for trait milk, and $(growth.nQTL) QTL for trait growth.
    The TBV variances of the traits are standardized to 1.0. Narrow sense
    heritability of milk is $(milk.h2), and ``h^2`` of growth is $(growth.h2).
    """
    tprintln(desc, "\n")

    ########################################
    act = md"""
    ## Act III: Selection on milk of the founder population
    """
    tprintln(act, "\n")
    foo, bar = "founder", "rand"
    for ext in ["ped", "xy"]
        cp("$tstDir/$foo.$ext", "$tstDir/$bar.$ext", force = true)
    end
    lmp = deserialize("$tstDir/$foo.lmp")
    ped = deserialize("$tstDir/$bar.ped")
    xy  = "$tstDir/$bar.xy"
    for igrt in 1:nrng
        println()
        @info lpad("<--- Generation $igrt / $nrng --->", 40)
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, milk, growth)
        predict!(ids, ped, milk, growth)
        ng = Select(ids, plan, ped, milk)
        reproduce!(ng, ped, xy, lmp, milk, growth)
    end
    serialize("$tstDir/$bar.ped", ped)

    ########################################
    act = md"""
    ## Act IV: Selection on milk and growth of the founder population
    """
    tprintln("\n", act, "\n")
    bar = "both"
    for ext in ["ped", "xy"]
        cp("$tstDir/$foo.$ext", "$tstDir/$bar.$ext", force = true)
    end
    ped = deserialize("$tstDir/$bar.ped")
    xy  = "$tstDir/$bar.xy"
    wgt = Dict{String, Float64}("milk" => 0.6, "growth" => 0.4)
    for igrt in 1:nrng
        println()
        @info lpad("<--- Generation $igrt / $nrng --->", 40)
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, milk, growth)
        predict!(ids, ped, milk, growth)
        ng = Select(ids, plan, ped, wgt)
        reproduce!(ng, ped, xy, lmp, milk, growth)
    end
    serialize("$tstDir/$bar.ped", ped)

    ########################################
    act = md"""
    ## Act V: Directional selection on EBV of milk

    This selection is based on the $nrng generations of random selection
    results on milk. The ID in the last generation are as candidates.
    """
    tprintln("\n", act, "\n")
    foo, bar = "rand", "ablup"
    for ext in ["ped", "xy"]
        cp("$tstDir/$foo.$ext", "$tstDir/$bar.$ext", force = true)
    end
    ped = deserialize("$tstDir/$bar.ped")
    xy  = "$tstDir/$bar.xy"
    for igrt in 1:nsel
        println()
        @info lpad("<--- Generation $igrt / $nsel --->", 40)
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, milk, growth)
        G = nrm(ped)
        giv = inv(G)
        predict!(ids, ped, fixed, giv, milk)
        ng = Select(ids, plan, ped, milk)
        reproduce!(ng, ped, xy, lmp, milk, growth)
    end

    ########################################
    act = md"""
    ## Act VI: Optimal contribution selection on milk

    This selection is based on the random selection results also.
    This is AABLUP.
    """
    tprintln("\n", act, "\n")
    bar = "ocs"
    for ext in ["ped", "xy"]
        cp("$tstDir/$foo.$ext", "$tstDir/$bar.$ext", force = true)
    end
    ped = deserialize("$tstDir/$bar.ped")
    xy  = "$tstDir/$bar.xy"
    plan = Plan(25, 50, 200, mate = :random)
    for igrt in 1:nsel
        println()
        @info lpad("<--- Generation $igrt / $nsel --->", 40)
        ids = view(ped, ped.grt .== ped.grt[end], :id)
        phenotype!(ids, ped, milk, growth)
        giv, R = begin
            G = nrm(ped)
            R = G[ids, ids]
            inv(G), R
        end
        predict!(ids, ped, fixed, giv, milk)
        tgt = select(ped, "ebv_milk" => "idx", "sex")[ids, :]
        K = 0.09
        c = TM1997(tgt, R, K)
        #= with Theo's or
        dat = Float64.([tgt.idx  tgt.sex .== 1 tgt.sex .== 0])
        cb = fungencont(dat, R, K)
        return [c cb]  # passed
        =#
        ng = Select(ids, plan, ped, c)
        reproduce!(ng, ped, xy, lmp, milk, growth)
    end
    serialize("$tstDir/$bar.ped", ped)

    ########################################
    act = md"""
    ## Act VII: Sum test
    """
    tprintln("\n", act, "\n")
end
